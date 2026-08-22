import csv
from pathlib import Path


class CSVResultWriter:
    """Write stable, headered CSV output without a dataframe dependency."""

    def write(self, path: Path, rows: list[dict[str, object]]) -> None:
        if not rows:
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
