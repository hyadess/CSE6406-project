import unittest

from cse6406.core.errors import PipelineError
from cse6406.parallel.executor import WorkerPool, _format_duration
from cse6406.parallel.pipeline import default_workers


def _double(value):
    return value * 2


def _fail_on(value, bad):
    if value == bad:
        raise PipelineError(f"task {value} failed")
    return value


def _needs_context(_value):
    from cse6406.parallel import tasks
    return tasks.CONTEXT is not None


def _install_marker(marker):
    from cse6406.parallel import tasks
    tasks.CONTEXT = marker


class WorkerPoolTests(unittest.TestCase):
    def test_results_keep_submission_order(self):
        tasks = [(value,) for value in range(200)]
        with WorkerPool(4, pending_per_worker=2) as pool:
            results = list(pool.map_ordered(_double, tasks, label="test"))
        self.assertEqual(results, [value * 2 for value in range(200)])

    def test_single_worker_runs_in_process(self):
        tasks = [(value,) for value in range(10)]
        with WorkerPool(1) as pool:
            self.assertIsNone(pool._pool)
            results = list(pool.map_ordered(_double, tasks, label="test"))
        self.assertEqual(results, [value * 2 for value in range(10)])

    def test_worker_failure_propagates_to_the_caller(self):
        tasks = [(value, 7) for value in range(50)]
        with WorkerPool(4) as pool:
            with self.assertRaises(PipelineError):
                list(pool.map_ordered(_fail_on, tasks, label="test"))

    def test_initializer_runs_in_every_worker(self):
        with WorkerPool(3, initializer=_install_marker, initargs=("marker",)) as pool:
            results = list(pool.map_ordered(
                _needs_context, [(value,) for value in range(30)], label="test"
            ))
        self.assertTrue(all(results))

    def test_default_worker_count_is_usable(self):
        self.assertGreaterEqual(default_workers(), 1)

    def test_durations_are_human_readable(self):
        self.assertEqual(_format_duration(0), "00:00:00")
        self.assertEqual(_format_duration(3_661), "01:01:01")
        self.assertEqual(_format_duration(90_061), "1d 01:01:01")


if __name__ == "__main__":
    unittest.main()
