#!/usr/bin/env python3
"""Audit historical-project migration classification without modifying projects."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "app"))

from project_migration import DEFAULT_MIGRATION_REGISTRY, audit_migration_registry  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit primer project migration registry")
    parser.add_argument("--registry", type=Path, default=DEFAULT_MIGRATION_REGISTRY)
    args = parser.parse_args()
    report = audit_migration_registry(args.registry)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
