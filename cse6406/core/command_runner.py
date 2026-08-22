from __future__ import annotations

import shlex
import subprocess
from pathlib import Path

from .errors import PipelineError


class CommandRunner:
    """Execute one auditable external command and retain its complete log."""

    def run(
        self,
        command: list[object],
        *,
        log: Path,
        cwd: Path | None = None,
        timeout: int = 3600,
    ) -> str:
        argv = [str(part) for part in command]
        display = shlex.join(argv)
        try:
            result = subprocess.run(
                argv,
                cwd=str(cwd) if cwd else None,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            raise PipelineError(f"Could not run {display}: {error}") from error

        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(f"$ {display}\n\n{result.stdout}")
        if result.returncode:
            raise PipelineError(
                f"Command failed with exit code {result.returncode}: {display}; see {log}"
            )
        return result.stdout
