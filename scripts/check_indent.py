from pathlib import Path
import re
import sys


def main() -> int:
    root = Path(__file__).resolve().parents[1] / "app"
    bad = []
    for path in sorted(root.glob("*.py")):
        for ln, line in enumerate(path.read_text(encoding="utf-8", errors="ignore").splitlines(), start=1):
            m = re.match(r"^[ \t]+", line)
            if not m:
                continue
            indent = m.group(0)
            if "\t" in indent:
                bad.append((str(path), ln))
    if bad:
        print("Indent check failed: tabs found in leading indentation.")
        for p, ln in bad:
            print(f"{p}:{ln}")
        return 1
    print("Indent check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

