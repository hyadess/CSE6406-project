from __future__ import annotations

import time
from collections import deque
from concurrent.futures import Future, ProcessPoolExecutor
from multiprocessing import get_context
from typing import Callable, Iterator, Sequence


def _format_duration(seconds: float) -> str:
    seconds = int(max(seconds, 0))
    days, seconds = divmod(seconds, 86_400)
    hours, seconds = divmod(seconds, 3_600)
    minutes, seconds = divmod(seconds, 60)
    clock = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
    return f"{days}d {clock}" if days else clock


class ProgressReporter:
    """Print a single-line completion, rate, and ETA trace for one phase."""

    def __init__(self, label: str, total: int, *, every: float = 15.0):
        self.label, self.total, self.every = label, total, every
        self.started = self.last = time.monotonic()
        self.done = 0
        if total:
            print(f"[{label}] starting {total} task(s)", flush=True)

    def advance(self, count: int = 1) -> None:
        self.done += count
        now = time.monotonic()
        if now - self.last < self.every and self.done != self.total:
            return
        self.last = now
        elapsed = now - self.started
        rate = self.done / elapsed if elapsed > 0 else 0.0
        remaining = (self.total - self.done) / rate if rate > 0 else float("inf")
        percent = 100.0 * self.done / self.total if self.total else 100.0
        print(
            f"[{self.label}] {self.done}/{self.total} ({percent:.1f}%) "
            f"elapsed {_format_duration(elapsed)} "
            f"eta {_format_duration(remaining) if rate > 0 else '--'} "
            f"rate {rate:.2f}/s",
            flush=True,
        )

    def finish(self) -> None:
        elapsed = time.monotonic() - self.started
        if self.total:
            state = "completed" if self.done == self.total else "stopped after"
            print(f"[{self.label}] {state} {self.done}/{self.total} task(s) in "
                  f"{_format_duration(elapsed)}", flush=True)


class WorkerPool:
    """Run picklable task functions across processes, yielding results in order.

    Submission is bounded so that a phase with hundreds of thousands of tasks
    never materializes every pending result at once, and completion order is
    restored to submission order so every CSV table keeps the exact row order
    that the sequential pipeline produces.
    """

    def __init__(
        self, workers: int, *,
        initializer: Callable[..., None] | None = None,
        initargs: tuple = (),
        pending_per_worker: int = 4,
    ):
        self.workers = max(1, int(workers))
        self.initializer, self.initargs = initializer, initargs
        self.pending_per_worker = max(1, int(pending_per_worker))
        self._pool: ProcessPoolExecutor | None = None

    def __enter__(self) -> "WorkerPool":
        if self.workers > 1:
            self._pool = ProcessPoolExecutor(
                max_workers=self.workers,
                mp_context=get_context("spawn"),
                initializer=self.initializer,
                initargs=self.initargs,
            )
        elif self.initializer is not None:
            # One worker means in-process execution, which still needs the
            # module-level context that the task functions read.
            self.initializer(*self.initargs)
        return self

    def __exit__(self, *exception) -> None:
        if self._pool is not None:
            self._pool.shutdown(wait=True, cancel_futures=True)
            self._pool = None

    def map_ordered(
        self, function: Callable, tasks: Sequence[tuple], *, label: str,
    ) -> Iterator:
        progress = ProgressReporter(label, len(tasks))
        try:
            if self._pool is None:
                for task in tasks:
                    yield function(*task)
                    progress.advance()
            else:
                yield from self._stream(function, iter(tasks), progress)
        finally:
            progress.finish()

    def _stream(self, function: Callable, tasks: Iterator[tuple], progress) -> Iterator:
        pool = self._pool
        assert pool is not None
        limit = self.workers * self.pending_per_worker
        pending: deque[Future] = deque()
        exhausted = False
        try:
            while True:
                while not exhausted and len(pending) < limit:
                    task = next(tasks, None)
                    if task is None:
                        exhausted = True
                        break
                    pending.append(pool.submit(function, *task))
                if not pending:
                    return
                result = pending.popleft().result()
                progress.advance()
                yield result
        except BaseException:
            for future in pending:
                future.cancel()
            raise
