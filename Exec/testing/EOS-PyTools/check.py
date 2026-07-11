#!/usr/bin/env python3
"""Gate on the pytest log written by ./run.

Pass requires a pytest summary line reporting only passes (and optionally
skips), e.g. '15 passed in 1.4s'. Any 'failed'/'error' token, a missing
pytest module, or an absent/empty log fails the case.
"""

import os
import re
import sys

LOG = "pytest_log.txt"


def fail(msg):
    print("EOS-PyTools: FAIL — %s" % msg)
    sys.exit(1)


def main():
    if not os.path.exists(LOG):
        fail("no %s (did ./run execute?)" % LOG)
    with open(LOG) as f:
        log = f.read()
    # strip ANSI colour codes (pytest may colourise depending on the env)
    log = re.sub(r"\x1b\[[0-9;]*m", "", log)
    if not log.strip():
        fail("empty pytest log")
    if "No module named pytest" in log:
        fail("pytest is not installed (see Exec/python_analysis/environment.yml)")

    # last pytest summary line, e.g. '15 passed in 1.39s' or
    # '1 failed, 14 passed in 1.42s'
    summaries = re.findall(r"^([0-9].*(?:passed|failed|error).*)$", log,
                           re.MULTILINE)
    if not summaries:
        fail("no pytest summary line found")
    summary = summaries[-1]
    if re.search(r"\b(failed|error)s?\b", summary):
        fail(summary)
    m = re.search(r"(\d+) passed", summary)
    if not m or int(m.group(1)) == 0:
        fail("no tests passed: " + summary)

    print("EOS-PyTools: PASS — %s" % summary)
    sys.exit(0)


if __name__ == "__main__":
    main()
