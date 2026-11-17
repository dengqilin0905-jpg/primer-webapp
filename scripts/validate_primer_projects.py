#!/usr/bin/env python3
"""Run the fail-closed portfolio completion gate for all primer projects."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "app"))

from project_migration import DEFAULT_MIGRATION_REGISTRY  # noqa: E402
from project_validator import DEFAULT_PRIMER_SKILL_ROOT, audit_primer_portfolio  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the complete primer project portfolio")
    parser.add_argument("--migration-registry", type=Path, default=DEFAULT_MIGRATION_REGISTRY)
    parser.add_argument("--skill-root", type=Path, default=DEFAULT_PRIMER_SKILL_ROOT)
    args = parser.parse_args()
    report = audit_primer_portfolio(args.migration_registry, skill_root=args.skill_root)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
