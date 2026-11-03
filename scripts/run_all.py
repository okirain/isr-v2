"""Run the complete ISR-V2 analysis chain.

Usage
-----
    python scripts/run_all.py

The script regenerates every number and figure used in the thesis:

* KRUSTY 28-hour validation and cooldown benchmark
* ISR-V2 heat-pipe array design scan
* fuel-temperature, failure, mass and operating-limit figures
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.isr_design.comparison import main as run_comparison  # noqa: E402
from src.validation.krusty_benchmark import main as run_validation  # noqa: E402


def main() -> dict:
    print("=" * 72)
    print("1/2  KRUSTY 28-hour validation")
    print("=" * 72)
    validation = run_validation()
    print(json.dumps(validation, indent=2))

    print("=" * 72)
    print("2/2  ISR-V2 design and comparison")
    print("=" * 72)
    comparison = run_comparison()
    print(json.dumps(comparison, indent=2))

    summary = {"validation": validation, "comparison": comparison}
    path = ROOT / "results" / "run_summary.json"
    path.write_text(json.dumps(summary, indent=2))
    print(f"\nSummary written to {path}")
    return summary


if __name__ == "__main__":
    main()
