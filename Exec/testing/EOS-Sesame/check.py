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
    # common-energy-reference set (eref295): one shared gauge across the
    # three mixture members (HANDOFF_common_energy_reference.md + addendum)
    "data/ti-beta-21s_2963_coldext_eref295.eostab",
    "data/deuterium_5267_s301_eref295.eostab",
    "data/dry-air_5031_s301_eref295.eostab",
    # eref295 solid extensions on the same gauge (make_eref295_solids.sh)
    "data/aluminum_3720_coldext_eref295.eostab",
    "data/diamond_7834_s301_eref295.eostab",
    # sub-floor T extension of the eref295 member with the highest native
    # floor (air_lowT_extension/make_air_lowT.sh): same gauge, T floor
    # 100 K -> 17.7 K, the original table kept alongside (D5-a)
    "air_lowT_extension/data/dry-air_5031_s301_eref295_Tf1p25.eostab",
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
