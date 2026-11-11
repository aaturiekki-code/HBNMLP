#!/usr/bin/env python3
"""Generate missing QE sweep configurations without submitting jobs."""

import argparse
import itertools
from pathlib import Path
import re
import subprocess
import sys

import numpy as np


def parse_range(values, default=None):
    if values is None:
        return [default]

    if not np.isfinite(values).all():
        raise ValueError("Range values must be finite")

    start, end, count_value = values
    count = int(count_value)

    if count < 1 or count != count_value:
        raise ValueError("Point count must be a positive integer")

    return np.linspace(start, end, count).tolist()


def remove_overrides(arguments):
    arities = {
        "--distance": 1,
        "--scale-xy": 1,
        "--slide": 2,
        "-o": 1,
        "--output": 1,
        "--name": 1,
        "--slurm-name": 1,
        "--overwrite": 0,
        "--no-overwrite": 0,
    }

    result = []
    index = 0

    while index < len(arguments):
        argument = arguments[index]
        flag = argument.split("=", 1)[0]

        if flag in arities:
            index += 1 if "=" in argument else 1 + arities[flag]
        else:
            result.append(argument)
            index += 1

    return result


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, allow_abbrev=False
    )
    parser.add_argument("--parent-dir", required=True, type=Path)
    parser.add_argument("--distance-range", nargs=3, type=float)
    parser.add_argument("--scale-range", nargs=3, type=float)
    parser.add_argument("--slide-x-range", nargs=3, type=float)
    parser.add_argument("--slide-y-range", nargs=3, type=float)
    parser.add_argument("--variant-tag")

    args, forwarded = parser.parse_known_args()
    forwarded = remove_overrides(forwarded)

    if args.variant_tag and not re.fullmatch(
        r"[A-Za-z0-9_-]+", args.variant_tag
    ):
        parser.error("Invalid variant tag")

    try:
        distances = parse_range(args.distance_range)
        scales = parse_range(args.scale_range, 1.44)
        xs = parse_range(args.slide_x_range)
        ys = parse_range(args.slide_y_range)
    except ValueError as error:
        parser.error(str(error))

    if xs == [None] and ys == [None]:
        slides = [(None, None)]
    else:
        slides = list(itertools.product(
            [0.0] if xs == [None] else xs,
            [0.0] if ys == [None] else ys,
        ))

    parent = args.parent_dir.expanduser().resolve()
    parent.mkdir(parents=True, exist_ok=True)

    generator = Path(__file__).resolve().parent / "createscf.py"
    if not generator.is_file():
        parser.error(f"Missing {generator}")

    configurations = []
    names = set()

    for distance, scale, (sx, sy) in itertools.product(
        distances, scales, slides
    ):
        parts = []
        if distance is not None:
            parts.append(f"dist_{distance:.3f}")
        parts.append(f"scale_{scale:.3f}")
        if sx is not None:
            parts.append(f"slide_{sx:.3f}_{sy:.3f}")
        if args.variant_tag:
            parts.append(f"variant_{args.variant_tag}")

        name = "_".join(parts)
        if name in names:
            parser.error(f"Rounded folder-name collision: {name}")

        names.add(name)
        configurations.append((name, distance, scale, sx, sy))

    added = skipped = failed = 0

    for name, distance, scale, sx, sy in configurations:
        folder = parent / name

        if folder.exists():
            skipped += 1
            continue

        command = [
            sys.executable, str(generator),
            "-o", str(folder),
            "--no-overwrite",
            "--scale-xy", f"{scale:.6f}",
            "--slurm-name", f"{parent.name}_{name}",
        ]

        if distance is not None:
            command += ["--distance", f"{distance:.6f}"]
        if sx is not None:
            command += ["--slide", f"{sx:.6f}", f"{sy:.6f}"]

        command += forwarded
        result = subprocess.run(
            command, capture_output=True, text=True
        )

        complete = all(
            (folder / filename).is_file()
            for filename in (
                "pw.in", "run.sh", "createscf_meta.json"
            )
        )

        if result.returncode or not complete:
            failed += 1
            print(
                f"[FAILED] {folder}\n{result.stdout}\n{result.stderr}",
                file=sys.stderr,
            )
        else:
            added += 1
            print(f"[ADDED] {folder}")

    print(f"Added={added}, existing/skipped={skipped}, failed={failed}")
    print("No jobs submitted.")

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
