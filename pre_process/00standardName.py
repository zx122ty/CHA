#!/usr/bin/env python3
"""
Standardize file names in a folder:
If a file's basename ends with a letter immediately followed by digits
(no "-" separator), insert "-" between them.

Examples:
  N-Blank1.mzML  →  N-Blank-1.mzML
  N-QC6.mzML     →  N-QC-6.mzML
  N-Blank10.mzML →  N-Blank-10.mzML
"""

import os
import re
import sys


def standardize_filename(filename: str) -> str:
    """
    Insert '-' between a letter and trailing digits if they are directly adjacent
    in the basename (before the extension).
    """
    name, ext = os.path.splitext(filename)
    # Match: letter followed by one or more digits at the end of basename
    new_name = re.sub(r'([a-zA-Z])(\d+)$', r'\1_\2', name)
    return new_name + ext


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python 00standardName.py <folder_path>")
        sys.exit(1)

    folder_path = sys.argv[1]

    if not os.path.isdir(folder_path):
        print(f"Error: '{folder_path}' is not a valid directory.")
        sys.exit(1)

    # Collect all files (skip subdirectories)
    entries = sorted(os.listdir(folder_path))
    files = [e for e in entries if os.path.isfile(os.path.join(folder_path, e))]

    if not files:
        print(f"No files found in '{folder_path}'.")
        return

    print(f"Folder: {folder_path}  ({len(files)} files)\n{'─' * 60}")

    renamed = 0
    for fname in files:
        new_fname = standardize_filename(fname)

        if new_fname == fname:
            print(f"  ✓  {fname}")
            continue

        old_path = os.path.join(folder_path, fname)
        new_path = os.path.join(folder_path, new_fname)

        if os.path.exists(new_path):
            print(f"  ⚠  {fname}  →  SKIP (target '{new_fname}' already exists)")
            continue

        os.rename(old_path, new_path)
        print(f"  ✗  {fname}  →  {new_fname}")
        renamed += 1

    print(f"{'─' * 60}\nRenamed {renamed} file(s).")


if __name__ == "__main__":
    main()