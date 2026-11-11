#!/usr/bin/env python3
"""Generate QE structures with selectable calculation workflows."""

import argparse
import json
import sys

import createscf_base as base
from qe_response import write_ph


def write_workflow_runner(
    job_name,
    filename,
    cores,
    time_limit,
    mem_limit,
    workflow,
):
    text = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --time={time_limit}
#SBATCH --mem={mem_limit}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={cores}
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out

set -e
module restore qemod

echo "Starting SCF"
srun pw.x -np {cores} -inp pw.in > pw.out
"""

    if workflow != "scf-only":
        text += f"""
echo "Starting response calculation"
srun ph.x -np {cores} -inp ph.in > ph.out
"""

    if workflow == "full-response":
        text += """
echo "Starting dynmat"
srun --ntasks=1 dynmat.x -in dynmat.in > dynmat.out
"""

    text += '\necho "Workflow completed."\n'
    filename.write_text(text)
    filename.chmod(filename.stat().st_mode | 0o111)


def main():
    parser = argparse.ArgumentParser(add_help=False)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--scf-only", action="store_true")
    modes.add_argument("--electronic-only", action="store_true")
    modes.add_argument("--full-response", action="store_true")
    options, remaining = parser.parse_known_args()

    if options.scf_only:
        workflow = "scf-only"
    elif options.electronic_only:
        workflow = "electronic-only"
    else:
        # Preserve legacy behavior unless explicitly changed.
        workflow = "full-response"

    sys.argv = [sys.argv[0]] + remaining

    original_write_pw = base.write_pw
    output_folders = []

    def capture_write_pw(*args, **kwargs):
        filename = kwargs["filename"]
        output_folders.append(filename.parent)
        return original_write_pw(*args, **kwargs)

    def response_writer(prefix, outdir, filename):
        if workflow != "scf-only":
            write_ph(
                prefix,
                outdir,
                filename,
                electronic_only=(workflow == "electronic-only"),
            )

    def runner_writer(job_name, filename, cores, time_limit, mem_limit):
        write_workflow_runner(
            job_name,
            filename,
            cores,
            time_limit,
            mem_limit,
            workflow,
        )

    base.write_pw = capture_write_pw
    base.write_ph = response_writer
    base.write_bash_runner = runner_writer

    if workflow != "full-response":
        base.write_dynmat = lambda *args, **kwargs: None

    base.main()

    for folder in output_folders:
        path = folder / "createscf_meta.json"
        if not path.is_file():
            raise RuntimeError(f"Missing generated metadata: {path}")

        metadata = json.loads(path.read_text())
        metadata["workflow"] = workflow
        metadata["electronic_only"] = workflow == "electronic-only"

        path.write_text(
            json.dumps(metadata, indent=2, allow_nan=False) + "\n"
        )


if __name__ == "__main__":
    main()
