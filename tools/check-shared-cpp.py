#!/usr/bin/env python3
"""Detect drift between misha's and pymisha's shared C++ sources.

The two packages share their scan/format/PSSM core by manual copy, and nothing
has ever checked that the copies stay in sync. They have silently diverged more
than once - misha's 5.10.3 below-N skip guard never reached pymisha, and a
change to gtrackimportwig's arity would have broken pymisha's binding the next
time the sources were synced.

Comparison ignores comments and whitespace, so reformatting does not register as
drift. Files that are legitimately adapted for the Python binding are recorded in
the baseline with the misha-side hash they were last reconciled against; when
that side changes, the check asks a human whether the change needs porting. It
does not try to reconcile them automatically - most of the adaptations are
deliberate.

Usage:
    check-shared-cpp.py --misha <dir> --pymisha <dir> [--update-baseline]
"""

import argparse
import difflib
import hashlib
import json
import os
import re
import sys

BLOCK_COMMENT = re.compile(r"/\*.*?\*/", re.S)
LINE_COMMENT = re.compile(r"//[^\n]*")
BASELINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "shared-cpp-baseline.json")


def normalize(path):
    """Source text with comments stripped and whitespace collapsed."""
    with open(path, "r", errors="replace") as fh:
        text = fh.read()
    text = LINE_COMMENT.sub("", BLOCK_COMMENT.sub("", text))
    lines = (re.sub(r"\s+", " ", line).strip() for line in text.split("\n"))
    return "\n".join(line for line in lines if line)


def digest(text):
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def collect(misha_src, pymisha_src):
    r = {f for f in os.listdir(misha_src) if f.endswith((".cpp", ".h"))}
    p = {f for f in os.listdir(pymisha_src) if f.endswith((".cpp", ".h"))}
    out = {}
    for name in sorted(r & p):
        a = normalize(os.path.join(misha_src, name))
        b = normalize(os.path.join(pymisha_src, name))
        out[name] = {"identical": a == b, "misha": digest(a), "pymisha": digest(b),
                     "_a": a, "_b": b}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--misha", required=True)
    ap.add_argument("--pymisha", required=True)
    ap.add_argument("--update-baseline", action="store_true")
    ap.add_argument("--show-diff", action="store_true")
    args = ap.parse_args()

    current = collect(os.path.join(args.misha, "src"), os.path.join(args.pymisha, "src"))

    if args.update_baseline:
        payload = {n: {k: v for k, v in d.items() if not k.startswith("_")}
                   for n, d in current.items()}
        with open(BASELINE, "w") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
            fh.write("\n")
        ident = sum(1 for d in current.values() if d["identical"])
        print(f"baseline written: {len(current)} shared files, "
              f"{ident} identical, {len(current) - ident} adapted")
        return 0

    if not os.path.exists(BASELINE):
        print("no baseline; run with --update-baseline first", file=sys.stderr)
        return 2
    with open(BASELINE) as fh:
        base = json.load(fh)

    broke, unreviewed, added = [], [], []
    for name, cur in current.items():
        prev = base.get(name)
        if prev is None:
            added.append(name)
        elif prev["identical"] and not cur["identical"]:
            broke.append(name)
        elif not prev["identical"] and cur["misha"] != prev["misha"]:
            unreviewed.append(name)
    removed = sorted(set(base) - set(current))

    for name in broke:
        print(f"DRIFT   {name}: was identical in both packages, no longer is")
        if args.show_diff:
            d = difflib.unified_diff(current[name]["_a"].split("\n"),
                                     current[name]["_b"].split("\n"),
                                     "misha", "pymisha", lineterm="", n=2)
            print("\n".join(f"        {l}" for l in d))
    for name in unreviewed:
        print(f"REVIEW  {name}: adapted in pymisha, and misha's side changed - "
              f"check whether the change needs porting")
    for name in added:
        print(f"NEW     {name}: newly shared; add it to the baseline")
    for name in removed:
        print(f"GONE    {name}: no longer shared; drop it from the baseline")

    if broke or unreviewed or added or removed:
        print(f"\n{len(broke)} broke parity, {len(unreviewed)} need review, "
              f"{len(added)} new, {len(removed)} gone.")
        print("If the change is intentional and pymisha has been updated to match, "
              "re-run with --update-baseline and commit the result.")
        return 1

    print(f"{len(current)} shared C++ files: no unreviewed drift.")
    return 0


sys.exit(main())
