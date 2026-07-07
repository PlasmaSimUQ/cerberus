#!/usr/bin/env python3
"""Gate for the one-zone tabulated-EOS self-test (Stage 2, W5).

Parses the EOSTAB-SELFTEST lines from run_log.txt.

Tier 1 (data/ideal_synthetic.eostab) is the hard gate: every check must
PASS — the table is closed-form ideal gas, so any failure is a code bug.

Tier 2 (data/D_fpeos.eostab) gates only the robustness checks (reader,
round trips, hull behaviour, degenerate corner). The 'fd-vs-blocks' check
compares the offline PCHIP-conditioned derivative blocks against finite
differences of the bilinear value surface; on real conditioned data the
two may legitimately disagree beyond the tier-1 tolerance, so it is
REPORTED but not gated (plan: tier 2 is consistency/robustness only).
"""

import re
import sys

TIER1 = "data/ideal_synthetic.eostab"
TIER2 = "data/D_fpeos.eostab"

# checks gated per tier
GATED = {
    TIER1: {"reader", "roundtrip-e", "roundtrip-p", "identities",
            "fd-vs-blocks", "hull", "corner"},
    TIER2: {"reader", "roundtrip-e", "roundtrip-p", "identities",
            "hull", "corner"},
}
# committed iteration ceiling (review amendment: fixed the first time the
# harness ran; see STAGE2.md check 2)
ITERS_MAX_CEILING = 60

line_re = re.compile(r"EOSTAB-SELFTEST\[([\w-]+)\] (PASS|FAIL) (.*)")
overall_re = re.compile(r"EOSTAB-SELFTEST OVERALL (PASS|FAIL) table=(\S+)")


def main():
    try:
        log = open("run_log.txt").read().splitlines()
    except FileNotFoundError:
        print("check.py: run_log.txt not found")
        return 1

    blocks = {}  # table path -> {check name: (PASS/FAIL, detail)}
    pending = {}
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
    for table, gated in GATED.items():
        checks = blocks.get(table)
        if checks is None:
            print("check.py: FAIL — no self-test output for %s" % table)
            rc = 1
            continue
        for name in sorted(gated):
            status, detail = checks.get(name, ("MISSING", ""))
            ok = status == "PASS"
            if not ok:
                rc = 1
            print("check.py: %-28s %-12s %s  %s"
                  % (table.split("/")[-1], name, status if ok else "** " + status,
                     detail))
            # iteration ceiling on the round-trip checks
            if name.startswith("roundtrip") and ok:
                m = re.search(r"iters_max=(\d+)", detail)
                if m and int(m.group(1)) > ITERS_MAX_CEILING:
                    print("check.py: %s %s iters_max=%s exceeds ceiling %d"
                          % (table, name, m.group(1), ITERS_MAX_CEILING))
                    rc = 1
        # ungated checks: report only
        for name, (status, detail) in sorted(checks.items()):
            if name not in gated:
                print("check.py: %-28s %-12s %s (reported, not gated)"
                      % (table.split("/")[-1], name, status))

    print("check.py:", "PASS" if rc == 0 else "FAIL")
    return rc


if __name__ == "__main__":
    sys.exit(main())
