#!/usr/bin/env python3

import subprocess
import sys
import re

jobid = sys.argv[1]

try:
    # Try querying sacct first (preferred)
    result = subprocess.run(
        ["sacct", "-j", jobid, "--format=State", "--noheader"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    output = result.stdout.strip().lower()

    if not output:
        # fallback to scontrol if sacct is empty
        result = subprocess.run(
            ["scontrol", "show", "job", jobid],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        output = result.stdout

        if "JobState=" in output:
            state = re.search(r"JobState=(\w+)", output).group(1).lower()
        else:
            raise Exception("Could not parse scontrol output")
    else:
        # Use first line of sacct output
        state = output.split()[0]

    # Map SLURM states to Snakemake states
    if state in ["pending", "configuring", "completing"]:
        print("running")  # Snakemake treats this as still active
    elif state in ["running", "resizing", "suspended"]:
        print("running")
    elif state in ["cancelled", "failed", "timeout", "node_fail", "preempted"]:
        print("failed")
    elif state in ["completed"]:
        print("success")
    else:
        print("unknown")

except Exception as e:
    print("failed")
    sys.exit(1)