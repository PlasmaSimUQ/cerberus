#!/usr/bin/env python3
"""Gate for the SESAME-table one-zone self-tests (doc/eos_sesame_plan.md).

Parses EOSTAB-SELFTEST lines from run_log.txt. Every table here is real
conditioned SESAME data, so the tier-2 policy applies (the EOS-Table
precedent): the robustness checks are gated, 'fd-vs-blocks' is reported
but not gated — the cavitated-response derivative floor (plan §3.0) makes
the blocks deliberately stiffer than the value surface inside the
crossover band.
"""

import re
import sys

TABLES = [
    "data/copper_3337_s311.eostab",
    "data/deuterium_5267_s301.eostab",
    "data/diamond_7834_s301.eostab",
    "data/hydrogen_5251_s301.eostab",
    "data/ti-beta-21s_2963_s311.eostab",
    "data/ti-beta-21s_2963_trackP.eostab",
    "data/ti-beta-21s_2963_coldext.eostab",
]
GATED = {"reader", "roundtrip-e", "roundtrip-p", "identities", "hull",
         "corner"}

line_re = re.compile(r"EOSTAB-SELFTEST\[([\w-]+)\] (PASS|FAIL) (.*)")
overall_re = re.compile(r"EOSTAB-SELFTEST OVERALL (PASS|FAIL) table=(\S+)")


def main():
    try:
        log = open("run_log.txt").read().splitlines()
    except FileNotFoundError:
        print("check.py: run_log.txt not found")
        return 1

    blocks, pending = {}, {}
    for ln in log:
        m = line_re.search(ln)
        if m:
            pending[m.group(1)] = (m.group(2), m.group(3))
            continue
        m = overall_re.search(ln)
        if m:
            blocks[m.group(2)] = pending
            pending = {}

    rc = 0
    for table in TABLES:
        checks = blocks.get(table)
        if checks is None:
            print("FAIL  %s: no self-test block found" % table)
            rc = 1
            continue
        for name, (verdict, detail) in sorted(checks.items()):
            gated = name in GATED
            ok = (verdict == "PASS") or not gated
            tag = "PASS" if verdict == "PASS" else (
                "FAIL" if gated else "report")
            print("%-6s %s %-14s %s" % (tag, table, name, detail))
            if not ok:
                rc = 1
    print("check.py: %s" % ("PASS" if rc == 0 else "FAIL"))
    return rc


if __name__ == "__main__":
    sys.exit(main())
