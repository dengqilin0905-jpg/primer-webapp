#!/usr/bin/env python3
"""Run one canonical primer project through the shared contract orchestration."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "app"))

from project_orchestration import run_project_directory  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Build and validate one Primer Project Contract project")
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--runtime-inputs", type=Path)
    parser.add_argument("--order-template", type=Path)
    parser.add_argument("--sequencing-template", type=Path)
    parser.add_argument("--overwrite-outputs", action="store_true")
    args = parser.parse_args()
    result = run_project_directory(
        args.project_root,
        runtime_input_path=args.runtime_inputs,
        order_template_path=args.order_template,
        sequencing_template_path=args.sequencing_template,
        overwrite_outputs=args.overwrite_outputs,
    )
    print(json.dumps({
        "summary": result["summary"],
        "conformance": result["conformance"],
        "normalized_canonical_result_sha256": result[
            "normalized_canonical_result_sha256"
        ],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
