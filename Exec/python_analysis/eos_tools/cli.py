"""Offline EOS table conditioning tool for Cerberus tabulated-EOS (plan W3).

Produces canonical ``.eostab`` files from raw source data (FPEOS first),
plus QA plots. All curation happens HERE, offline and inspectable; the C++
reader (Stage 2, MFP_eos_table) parses only the canonical format (frozen
spec: see eos_tools/formats/eostab.py and Exec/testing/EOS-Table/README.md).

Commands:
    synthetic  — analytic ideal-gas table (Stage-2 tier-1 gate + fallback)
    fpeos      — ingest FPEOS H table, isotope-scale to D, regrid,
                 condition, write, QA
    qa         — QA plots for an existing .eostab
    sources    — verify/fetch the raw-data manifest (data/raw/sources.yaml)

Examples (from Exec/testing/EOS-Table/):
    python3 ../../python_analysis/eos_table_prep.py synthetic \\
        --out data/ideal_synthetic.eostab
    python3 ../../python_analysis/eos_table_prep.py fpeos \\
        --src data/raw/FPEOS/H_EOS_09-18-20.txt \\
        --out data/D_fpeos.eostab --qa qa
    python3 ../../python_analysis/eos_table_prep.py sources --verify
"""

import argparse
import datetime
import os
import subprocess

import numpy as np

from .condition import condition, inverse_maps, shift_energy
from .constants import KB, M_D
from .formats.eostab import read_eostab, write_eostab
from .formats.fpeos import fpeos_to_deuterium, read_fpeos
from .grids import regrid
from .sources import main_sources


def git_sha():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=os.path.dirname(os.path.abspath(__file__)),
        ).decode().strip()
    except Exception:
        return "unknown"


def provenance(material, source, comp, e_shift, stats):
    now = datetime.date.today().isoformat()
    return [
        ("material", material),
        ("source", source),
        ("generator", "eos_table_prep.py %s, run %s" % (git_sha(), now)),
        ("composition", comp),
        ("units", "cgs"),
        ("e_shift", "%.10e" % e_shift),
        ("conditioning", "cv_floor=%.4e cv_floored=%d monotonised=%d maxwell=%s"
         % (stats["cv_floor"], stats["cv_floored"], stats["monotonised"],
            stats["maxwell"])),
    ]


def cmd_synthetic(args):
    """Ideal-gas gamma-law table with analytic blocks (closed-form truth)."""
    g = args.gamma
    R = KB / M_D
    lrho = np.linspace(np.log10(args.rho_min), np.log10(args.rho_max), args.n_rho)
    lT = np.linspace(np.log10(args.T_min), np.log10(args.T_max), args.n_T)
    rho = 10.0 ** lrho[:, None]
    T = 10.0 ** lT[None, :]
    e = R * T / (g - 1.0) * np.ones_like(rho)
    p = rho * R * T
    blocks = {
        "p": p, "e": e,
        "dpdT": rho * R * np.ones_like(T),
        "dpdrho": R * T * np.ones_like(rho),
        "cv": R / (g - 1.0) * np.ones_like(p),
        "dedrho": np.zeros_like(p),
        "hull": np.ones_like(p),
    }
    le, lp, T_of_e, T_of_p = inverse_maps(lrho, lT, p, e)
    blocks["T_of_e"], blocks["T_of_p"] = T_of_e, T_of_p
    stats = dict(cv_floor=0.0, cv_floored=0, monotonised=0, maxwell="none-needed")
    prov = provenance("ideal-D", "synthetic gamma-law gamma=%g (this tool)" % g,
                      "A=2.014 Z=1", 0.0, stats)
    write_eostab(args.out, prov, lrho, lT, blocks, {"le": le, "lp": lp})
    if args.qa:
        from .qa import qa_plots
        qa_plots(args.qa, lrho, lT, p, e, blocks["cv"], blocks["hull"],
                 tag="_synthetic")


def cmd_fpeos(args):
    raw = read_fpeos(args.src)
    pts = fpeos_to_deuterium(raw)
    print("source points: %d  (rho %.4g..%.4g g/cc D-equivalent, T %.4g..%.4g K)"
          % (len(pts["rho"]), pts["rho"].min(), pts["rho"].max(),
             pts["T"].min(), pts["T"].max()))
    lrho, lT, p, e, hull, loo = regrid(pts, args.n_rho, args.n_T)
    e, e_shift = shift_energy(e, hull)
    pts_plot = dict(pts, e=pts["e"])  # raw (unshifted) for overlays
    dpdT, dpdrho, cv, dedrho, stats = condition(lrho, lT, p, e)
    le, lp, T_of_e, T_of_p = inverse_maps(lrho, lT, p, e)
    blocks = dict(p=p, e=e, dpdT=dpdT, dpdrho=dpdrho, cv=cv, dedrho=dedrho,
                  T_of_e=T_of_e, T_of_p=T_of_p, hull=hull)
    prov = provenance(
        "D",
        "FPEOS %s (militzer.berkeley.edu, PRE 103 013203 (2021); H table "
        "isotope-scaled to D, see data/raw/README.md)" % os.path.basename(args.src),
        "A=2.014 Z=1", e_shift, stats)
    write_eostab(args.out, prov, lrho, lT, blocks, {"le": le, "lp": lp})
    print("conditioning: %s ; hull coverage %.1f%%"
          % (stats, 100.0 * hull.mean()))
    if args.qa:
        from .qa import qa_plots
        hug = dict(rho0=pts["hug_rho0"], e0=pts["hug_e0"] + e_shift,
                   p0=pts["hug_p0"])
        if args.hug_ref and os.path.exists(args.hug_ref):
            hug["ref"] = np.loadtxt(args.hug_ref)
            hug["ref_label"] = os.path.basename(args.hug_ref)
        qa_plots(args.qa, lrho, lT, p, e, cv, hull, pts=pts_plot, hug=hug,
                 loo=loo, tag="_D_fpeos")


def cmd_qa(args):
    from .qa import qa_plots
    prov, meta, blocks = read_eostab(args.table)
    qa_plots(args.qa, meta["lrho"], meta["lT"], blocks["p"], blocks["e"],
             blocks["cv"], blocks["hull"], tag="_" + os.path.basename(args.table))
    print("provenance:", {k: prov[k] for k in
                          ("material", "source", "conditioning") if k in prov})


def cmd_splice(args):
    import datetime
    import sys

    gen = "eos_table_prep.py %s, run %s" % (
        git_sha(), datetime.date.today().isoformat())
    if args.material == "Ti":
        from .qa_splice import run_qa_ti
        from .splice_ti import splice_titanium, write_spliced_ti
        res = splice_titanium(n_rho=args.n_rho, n_T=args.n_T)
        rc, _ = run_qa_ti(res, qa_dir=args.qa)
        if args.out:
            write_spliced_ti(res, args.out, gen)
        sys.exit(rc)

    from .qa_splice import run_qa
    from .splice import splice_deuterium, write_spliced

    res = splice_deuterium(n_rho=args.n_rho, n_T=args.n_T,
                           use_coolprop=not args.no_coolprop)
    for a in res["align_stats"]:
        print("splice: align %s cE=%.4e (n=%d, |dev|/kT=%.4f)"
              % (a["pair"], a["cE"], a["n"], a["std_over_kT"]))
    rc, _ = run_qa(res, qa_dir=args.qa)
    if args.out:
        write_spliced(res, args.out, gen)
    sys.exit(rc)


def cmd_coldmodel(args):
    import sys

    from .coldmodel_qa import run_qa
    sys.exit(run_qa(args))


def main():
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("synthetic", help="analytic ideal-gas .eostab")
    s.add_argument("--out", required=True)
    s.add_argument("--qa", default=None)
    s.add_argument("--gamma", type=float, default=1.4)
    s.add_argument("--n-rho", type=int, default=64)
    s.add_argument("--n-T", type=int, default=64)
    s.add_argument("--rho-min", type=float, default=1e-4)
    s.add_argument("--rho-max", type=float, default=10.0)
    s.add_argument("--T-min", type=float, default=1e3)
    s.add_argument("--T-max", type=float, default=1e7)
    s.set_defaults(func=cmd_synthetic)

    s = sub.add_parser("fpeos", help="condition an FPEOS element table -> D")
    s.add_argument("--src", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--qa", default=None)
    s.add_argument("--n-rho", type=int, default=96)
    s.add_argument("--n-T", type=int, default=96)
    s.add_argument("--hug-ref", default=None,
                   help="published Hugoniot points file (T rho P_GPa "
                        "compression) to overlay")
    s.set_defaults(func=cmd_fpeos)

    s = sub.add_parser("qa", help="QA plots for an existing .eostab")
    s.add_argument("--table", required=True)
    s.add_argument("--qa", required=True)
    s.set_defaults(func=cmd_qa)

    s = sub.add_parser("coldmodel",
                       help="build + QA the semi-analytic cold/low-T model")
    s.add_argument("--material", default="D", choices=["D", "Ti"])
    s.add_argument("--qa", default=None, help="QA plot output directory")
    s.add_argument("--n-rho", type=int, default=192)
    s.add_argument("--n-T", type=int, default=192)
    s.add_argument("--no-coolprop", action="store_true",
                   help="skip the CoolProp fluid piece (degraded; QA only)")
    s.set_defaults(func=cmd_coldmodel)

    s = sub.add_parser("splice",
                       help="build + QA the spliced wide-range table (SS3)")
    s.add_argument("--material", default="D", choices=["D", "Ti"])
    s.add_argument("--out", default=None,
                   help=".eostab output path (also writes .gz)")
    s.add_argument("--qa", default=None, help="QA output directory")
    s.add_argument("--n-rho", type=int, default=None)
    s.add_argument("--n-T", type=int, default=None)
    s.add_argument("--no-coolprop", action="store_true")
    s.set_defaults(func=cmd_splice)

    s = sub.add_parser("sources", help="verify/fetch the raw-data manifest")
    s.add_argument("--manifest", default=None,
                   help="path to sources.yaml (default: "
                        "Exec/testing/EOS-Table/data/raw/sources.yaml)")
    s.add_argument("--verify", action="store_true",
                   help="checksum every manifest file (default action)")
    s.add_argument("--fetch", action="store_true",
                   help="download absent files of scripted sources, then verify")
    s.add_argument("--quiet", action="store_true", help="report failures only")
    s.set_defaults(func=main_sources)

    args = ap.parse_args()
    args.func(args)
