from __future__ import annotations

import os


def aster_threads(threads: str) -> str:
    """Convert IQ-TREE's AUTO spelling to the integer ASTER requires."""
    if threads.upper() == "AUTO":
        return str(os.cpu_count() or 1)
    return threads
