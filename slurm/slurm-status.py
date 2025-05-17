#!/usr/bin/env python3
import subprocess
import sys

jobid = sys.argv[1]

try:
    # sacct gives one line:  <JobID>  <State>
    out = subprocess.check_output(
        ["sacct", "-P", "-b", "-n", "-o", "State", "-j", jobid],
        stderr=subprocess.DEVNULL,
        text=True
    ).strip().split()[0]
except subprocess.CalledProcessError:
    # sacct cannot find the job – assume vanished/failed
    print("failed")
    sys.exit(0)

state = out.split("+")[0]       # strip sub-states (e.g. COMPLETED+)
state = state.upper()

if state in {"PENDING", "CONFIGURING", "RUNNING", "COMPLETING", "SUSPENDED"}:
    print("running")
elif state in {"COMPLETED"}:
    print("success")
else:                            # CANCELLED, FAILED, TIMEOUT, NODE_FAIL, etc.
    print("failed")
