#!/usr/bin/env python3

from __future__ import annotations

import re
import sys
from pathlib import Path


def extract_restart_tag(restart_file: str) -> str:
    name = Path(restart_file).name
    match = re.search(r"(\d+)(?!.*\d)", name)
    if match:
        return match.group(1)

    stem = Path(restart_file).stem
    sanitized = re.sub(r"[^A-Za-z0-9]+", "_", stem).strip("_")
    return sanitized or "restart"


def sanitize_optional_tag(tag: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9._-]+", "_", tag).strip("_")
    return sanitized


def main() -> int:
    if len(sys.argv) not in (6, 7):
        raise SystemExit(
            "Usage: build_output_root.py <script_root> <restart_file> <target_T> <force_mag> [extra_tag] <output_file>"
        )

    script_root = Path(sys.argv[1]).resolve()
    restart_file = sys.argv[2]
    target_t = sys.argv[3]
    force_mag = sys.argv[4]
    if len(sys.argv) == 6:
        extra_tag = ""
        output_file = Path(sys.argv[5]).resolve()
    else:
        extra_tag = sanitize_optional_tag(sys.argv[5])
        output_file = Path(sys.argv[6]).resolve()

    restart_tag = extract_restart_tag(restart_file)
    output_root_name = f"output_restart{restart_tag}_Tstar_{target_t}_F{force_mag}"
    if extra_tag:
        output_root_name = f"{output_root_name}_{extra_tag}"
    output_root = script_root / output_root_name
    output_file.write_text(f"variable output_root string {output_root}\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
