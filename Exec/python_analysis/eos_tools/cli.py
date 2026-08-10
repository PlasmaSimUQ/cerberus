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
    sesame     — extract + condition one SESAME ASCII2 material/table
                 (doc/eos_sesame_plan.md; Maxwell + monotone crossover,
                 no tension clip)
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
from .constants import AMU_G, KB, M_D
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
    if args.mass_amu is None:
        R = KB / M_D
        material, comp = "ideal-D", "A=2.014 Z=1"
    else:
        R = KB / (args.mass_amu * AMU_G)
        material = args.material
        comp = "A=%g (mean particle mass, amu)" % args.mass_amu
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
    prov = provenance(material, "synthetic gamma-law gamma=%g (this tool)" % g,
                      comp, 0.0, stats)
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


def cmd_sesame(args):
    """SESAME ASCII2 -> .eostab (Track S; doc/eos_sesame_plan.md §3)."""
    import hashlib
    import re

    from .condition import monotonise_rho, monotonise_T
    from .constants import GPA_CGS
    from .formats.sesame_ascii2 import MJKG_CGS, index, read_sesame
    from .maxwell import condition_surface
    from .scan import g1_scan, report as g1_report

    if args.list:
        mats = index(args.src)
        pat = (args.grep or "").lower()
        for mid, m in mats.items():
            name = m["name"] or "(no 101 name)"
            if pat and pat not in name.lower() and pat != str(mid):
                continue
            print("%9d  %-55s tables: %s"
                  % (mid, name[:55], " ".join(str(t) for t in sorted(m["tables"]))))
        return
    if args.mat is None:
        raise SystemExit("--mat is required (use --list to browse materials)")

    # --- stage 0: source identity (S1, ti_splice_plan_v2) ----------------
    # One hash pass, reused for provenance. Verified against any 64-hex
    # sha256 recorded in the sources.yaml manifest (the sesame-unc entry
    # keeps its hash in free text, so we scan each entry's text).
    with open(args.src, "rb") as f:
        h = hashlib.sha256()
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    src_sha = h.hexdigest()
    if not args.no_verify_src:
        from .sources import default_manifest, load_manifest
        manifest = default_manifest()
        known = {}
        if os.path.exists(manifest):
            try:
                for ent in load_manifest(manifest)["sources"]:
                    for hx in re.findall(r"\b[0-9a-f]{64}\b", str(ent)):
                        known[hx] = ent.get("id", "?")
            except Exception as ex:
                print("WARNING: could not read manifest %s (%s); source "
                      "hash not verified" % (manifest, ex))
        if src_sha in known:
            print("source sha256 verified against sources.yaml (id: %s)"
                  % known[src_sha])
        elif known:
            raise SystemExit(
                "source hash GATE FAILED: %s\n  sha256 %s\n  matches no "
                "manifest entry in %s\n  (pass --no-verify-src to use an "
                "unmanifested library)" % (args.src, src_sha, manifest))
        else:
            print("WARNING: no sha256 records found in %s; source hash "
                  "not verified" % manifest)

    # --- stage 1: ingest (default table: 311 if present, else 301) -------
    table = args.table
    if table is None:
        try:
            raw = read_sesame(args.src, args.mat, 311)
            table = 311
        except ValueError:
            raw = read_sesame(args.src, args.mat, 301)
            table = 301
        print("table auto-selected: %d" % table)
    else:
        raw = read_sesame(args.src, args.mat, table)
        if table in (303, 304, 305):
            print("WARNING: table %d is a partial EOS component (2-T feed); "
                  "as a single-table closure it is physically incomplete" % table)
    name = re.split(r"[\s(]", raw["comment101"].split("material.")[-1].strip()
                    or "mat")[0].strip().lower() or "mat"

    # zero rows/columns cannot live on log axes; --T-min/--rho-min trim more
    kr = raw["rho"] > max(args.rho_min or 0.0, 0.0)
    kt = raw["T"] > max(args.T_min or 0.0, 0.0)
    rho_s, T_s = raw["rho"][kr], raw["T"][kt]
    p_s = raw["p"][np.ix_(kr, kt)] * GPA_CGS
    e_s = raw["e"][np.ix_(kr, kt)] * MJKG_CGS
    print("source: mat %d table %d  %dx%d kept (of %dx%d)  rho %.3g..%.3g "
          "g/cc  T %.6g..%.3g K  Helmholtz=%s"
          % (args.mat, table, len(rho_s), len(T_s), len(raw["rho"]),
             len(raw["T"]), rho_s[0], rho_s[-1], T_s[0], T_s[-1],
             raw["has_helmholtz"]))

    # --- stage 2: G1 pre-scan (routing) ----------------------------------
    g1_report(g1_scan(p_s), "raw")

    # --- stages 3-4: Maxwell placement + monotone crossover --------------
    p_c, e_c, band, mx = condition_surface(rho_s, T_s, p_s, e_s)
    print("crossover: %(n_maxwell)d maxwell / %(n_flat)d flat / "
          "%(n_ramp)d ramp isotherms; %(n_band_cells)d band cells" % mx)

    # --- stage 6: monotone-in-T value surfaces ---------------------------
    # cells the enforcement moves significantly are no longer source
    # physics -> they join the band (hull 0), like the crossover cells
    def _sig(a, b, rel=1e-6):
        return np.abs(a - b) > rel * np.maximum(np.abs(b), 1e-300)

    p_pre, e_pre = p_c.copy(), e_c.copy()
    p_c, nTp, rTp = monotonise_T(p_c)
    e_c, nTe, rTe = monotonise_T(e_c)
    n_monoT_out = int(np.sum((_sig(p_c, p_pre) | _sig(e_c, e_pre)) & ~band))
    band |= _sig(p_c, p_pre) | _sig(e_c, e_pre)
    if n_monoT_out:
        print("monoT moved %d cells outside the crossover band -> banded"
              % n_monoT_out)

    # --- stage 7: regrid onto log-uniform axes (source-native span, O3) --
    pts = {"rho": np.repeat(rho_s, len(T_s)), "T": np.tile(T_s, len(rho_s)),
           "p": p_c.ravel(), "e": e_c.ravel()}
    if args.lrho or args.lT:
        # Track-P mode: PRESCRIBED log10 axes (e.g. the shared mixture
        # T-bracket, requirement H3). regrid_onto refuses to extrapolate:
        # target cells outside the source span are nearest-filled and
        # marked hull 0 — the constant-in-T extension below the SESAME
        # T floor comes from exactly this semantics.
        from .grids import regrid_onto
        lrho = (np.linspace(args.lrho[0], args.lrho[1], int(args.lrho[2]))
                if args.lrho else
                np.linspace(np.log10(rho_s[0]), np.log10(rho_s[-1]),
                            args.n_rho or len(rho_s)))
        lT = (np.linspace(args.lT[0], args.lT[1], int(args.lT[2]))
              if args.lT else
              np.linspace(np.log10(T_s[0]), np.log10(T_s[-1]),
                          args.n_T or len(T_s)))
        p, e, hull, loo = regrid_onto(pts, lrho, lT)
        n_ext = int((hull < 0.5).sum())
        print("prescribed grid: lrho [%g, %g] x%d, lT [%g, %g] x%d; "
              "%d cells outside the source span (nearest-filled, hull 0)"
              % (lrho[0], lrho[-1], len(lrho), lT[0], lT[-1], len(lT), n_ext))
    else:
        n_rho = args.n_rho or len(rho_s)
        n_T = args.n_T or len(T_s)
        lrho, lT, p, e, hull, loo = regrid(pts, n_rho, n_T)
    if len(loo):
        print("regrid LOO: p rel median %.2e max %.2e; e span-normed "
              "median %.2e max %.2e" % (np.median(loo[:, 0]), loo[:, 0].max(),
                                        np.median(loo[:, 1]), loo[:, 1].max()))
    # band mask -> target grid (nearest source cell); band cells are hull 0
    ii = np.abs(np.log10(rho_s)[None, :] - lrho[:, None]).argmin(axis=1)
    jj = np.abs(np.log10(T_s)[None, :] - lT[:, None]).argmin(axis=1)
    band_g = band[np.ix_(ii, jj)]
    hull = np.where(band_g, 0.0, hull)

    # --- stage 7b (opt-in): solid-model cold extension (T3/T4, plan v2) --
    ce = None
    if args.cold_extend:
        from .coldext import cold_extend_stage
        p, e, hull, band_g, ce = cold_extend_stage(
            args.src, args.mat, raw, lrho, lT, p, e, hull, band_g,
            f_melt=args.melt_frac, align_tol=args.align_tol)

    # post-regrid enforcement (PCHIP can ripple): rho first, then T —
    # cummax along T preserves rho-monotonicity elementwise
    p_pre, e_pre = p.copy(), e.copy()
    p, nRp, rRp = monotonise_rho(p)
    p, nTp2, _ = monotonise_T(p)
    e, nTe2, _ = monotonise_T(e)
    hull = np.where(_sig(p, p_pre) | _sig(e, e_pre), 0.0, hull)

    # --- stage 9a: G1 acceptance on the emitted surface ------------------
    out_scan = g1_scan(p)
    g1_report(out_scan, "emitted")
    if not (out_scan["ok_rho"] and out_scan["ok_T"] and p.min() > 0.0):
        raise SystemExit("G1 gate FAILED on the emitted surface")

    # --- stage 8: gauge (common e-ref + shift), blocks, response, maps ---
    # Common-energy-reference recipe (HANDOFF addendum, 2026-08-10):
    #   (i)  e -> e - e(rho_ref, T_ref)   one common physical state
    #   (ii) e -> e + S                   ONE shared S for the whole set
    # e_shift in the header is documentary — the C++ reader never applies
    # it; the gauge must be baked into the stored values here.
    e_ref_val = None
    if args.e_ref_state:
        from .condition import sample_bilinear
        rho_ref, T_ref = args.e_ref_state
        try:
            e_ref_val = sample_bilinear(lrho, lT, e, rho_ref, T_ref)
        except ValueError as ex:
            raise SystemExit("--e-ref-state GATE FAILED: %s" % ex)
        h_ref = sample_bilinear(lrho, lT, hull, rho_ref, T_ref)
        e = e - e_ref_val
        print("e-ref (common gauge): e(%g g/cc, %g K) = %.10e erg/g "
              "subtracted; ref-point hull weight %.2f%s"
              % (rho_ref, T_ref, e_ref_val, h_ref,
                 "" if h_ref > 0.99 else
                 " (WARNING: reference sits on constructed fill — the "
                 "gauge is generator-dependent there)"))
    if args.probe_shift:
        # machine-readable minima for the set driver: the inverse-map le
        # axis takes log10 of the FULL-array minimum (hull-0 fills
        # included), so S must clear min_full, not just min_hull
        print("PROBE-SHIFT min_full=%.10e min_hull=%.10e max=%.10e"
              % (e.min(), e[hull > 0.5].min(), e.max()))
        return
    if args.e_shift is not None:
        e_shift = float(args.e_shift)
        if e.min() + e_shift <= 0.0:
            # NEVER silently top up a forced shift: a per-table top-up
            # would re-break the common gauge this flag exists to enforce
            raise SystemExit(
                "forced --e-shift GATE FAILED: min(e) + S = %.6e <= 0 "
                "(full-array min %.6e; re-probe the set and raise S)"
                % (e.min() + e_shift, e.min()))
        e = e + e_shift
    else:
        e, e_shift = shift_energy(e, hull)
        if e.min() <= 0.0:  # hull-0 fills may undershoot the hull-based shift
            extra = -e.min() + 1e-6 * (e.max() - e.min())
            e += extra
            e_shift += extra
    m_ref = (raw["abar"] * AMU_G) if raw["abar"] else M_D
    dpdT, dpdrho, cv, dedrho, stats = condition(
        lrho, lT, p, e, m_ref=m_ref,
        cs_floor_cms=(args.cs_floor * 1e5 if args.cs_floor else None))
    if stats["cs_floored"]:
        print("cs floor (W-B): dpdrho >= (%.3g km/s)^2 on %d cells"
              % (args.cs_floor, stats["cs_floored"]))

    # thermal stiffness floor: no single-phase fluid is isothermally softer
    # than ideal gas, dpdrho >= kB*T/m. Cells below it are tie-line residue
    # the band mapping missed (measured on 2963: dome-top cells flattened
    # by the regrid) -> floored, joined to the band, hull 0. Cells between
    # kT/m and c_cav^2 outside the band are GENUINE near-critical
    # softening and are left untouched. (Scope: all cells, demote failures
    # — correct for SESAME sources; see condition.thermal_floor.)
    #
    # MOLECULAR fluids (air, D2): the 201 abar is per ATOM, so kT/m_atom
    # over-floors the cold molecular region by the association factor
    # (air: 28.97/14.80 ~ 2x) and would demote the ambient itself.
    # --floor-mass-amu supplies the molecular mass; an ideal molecular gas
    # then sits ON the bound (dpdrho == kT/m exactly), so demotion is
    # restricted to cells genuinely below it (5% tolerance) — marginal
    # cells are floored to the exact bound but stay in-hull.
    from .condition import thermal_floor
    m_floor = (args.floor_mass_amu * AMU_G) if args.floor_mass_amu else m_ref
    dpdrho_f, kTm, sub = thermal_floor(dpdrho, lT, m_floor)
    demote = (dpdrho < 0.95 * kTm) if args.floor_mass_amu else sub
    dpdrho = dpdrho_f
    n_sub = int(demote.sum())
    if args.floor_mass_amu:
        print("thermal floor mass: %.4g amu (molecular, cli) vs abar %.4g"
              % (args.floor_mass_amu, raw["abar"] or 0.0))
    if n_sub or sub.any():
        band_g = band_g | demote
        hull = np.where(demote, 0.0, hull)
        print("thermal floor: %d cells floored to kT/m, %d sub-thermal "
              "(banded, hull 0)" % (int(sub.sum()), n_sub))

    if args.c_cav:
        c_cav = args.c_cav * 1e5  # km/s -> cm/s
        c_src = "cli"
    elif raw["bs0"] and raw["rho0"] and raw["bs0"] > 0:
        c_cav = float(np.sqrt(raw["bs0"] * GPA_CGS / raw["rho0"]))
        c_src = "201:BS/rho0"
    else:
        # coldest-isotherm slope at the first NON-band cell at/above rho0
        # (sampling inside the band would measure the soft bridge, not the
        # condensed branch — circular)
        i0 = int(np.abs(lrho - np.log10(raw["rho0"] or 10 ** lrho[-1])).argmin())
        while i0 < len(lrho) - 1 and band_g[i0, 0]:
            i0 += 1
        c_cav = float(np.sqrt(max(dpdrho[i0, 0], 0.0)))
        c_src = "coldest-isotherm slope at rho=%.3g g/cc" % 10 ** lrho[i0]
    c2 = c_cav ** 2
    n_resp = int(np.sum(band_g & (dpdrho < c2)))
    dpdrho = np.where(band_g, np.maximum(dpdrho, c2), dpdrho)
    print("cavitated response (D-e): c_cav = %.3f km/s (%s); dpdrho floored "
          "on %d band cells" % (c_cav / 1e5, c_src, n_resp))

    le, lp, T_of_e, T_of_p = inverse_maps(lrho, lT, p, e)
    blocks = dict(p=p, e=e, dpdT=dpdT, dpdrho=dpdrho, cv=cv, dedrho=dedrho,
                  T_of_e=T_of_e, T_of_p=T_of_p, hull=hull)

    # --- write -----------------------------------------------------------
    base = "%s_%d_s%d" % (name, args.mat, table)
    out = args.out or ("data/%s.eostab" % base)
    # QA tag follows the actual output name (an --out rename like the
    # Track-P table must not overwrite another table's QA plots)
    base = os.path.splitext(os.path.basename(out))[0]
    comment = " ".join(raw["comment101"].split())[:160]
    cond = ("cv_floor=%.4e cv_floored=%d monotonised=%d "
            "maxwell=constructed:%d/flats:%d/ramp:%d band_cells=%d "
            "c_cav=%.4e resp_floored=%d thermal_floored=%d tension_clip=0 "
            "monoT_p=%d(%.1e) monoT_e=%d(%.1e) monoRho_p=%d(%.1e)"
            % (stats["cv_floor"], stats["cv_floored"], stats["monotonised"],
               mx["n_maxwell"], mx["n_flat"], mx["n_ramp"], mx["n_band_cells"],
               c_cav, n_resp, n_sub, nTp + nTp2, rTp, nTe + nTe2, rTe,
               nRp, rRp))
    if stats.get("cs_floored"):
        cond += " cs_floor=%.4e cs_floored=%d" % (stats["cs_floor"],
                                                  stats["cs_floored"])
    if args.floor_mass_amu:
        cond += (" floor_mass_amu=%g floor_floored=%d floor_demote_tol=0.95"
                 % (args.floor_mass_amu, int(sub.sum())))
    if ce:
        v0, B0, B0p = ce["fit"]["vinet"]
        cond += (" cold_extend=vinet306(rho0K=%.4f,B0=%.4e,B0p=%.3f,"
                 "rms=%.2e,n=%d,theta0=%.0fK) melt_cap=411 f=%.2f "
                 "rho_max=%g cols=%d zone1=%d blend=%d bridge=%d trunc=%d "
                 "align_e=%.6e align_stdkT=%.3f pmis=%.2e cx2_band=%d"
                 % (1.0 / v0, B0, B0p, ce["fit"]["rms"], ce["fit"]["n"],
                    ce["model"]["theta0"], ce["f_melt"], ce["rho_max"],
                    ce["cols"], ce["zone1"], ce["blend"], ce["bridge"],
                    ce["blend_trunc"], ce["align_cE"], ce["align_std_kT"],
                    ce["p_mismatch_max"], ce["cx2_new_band"]))
    prov = [
        ("material", name),
        ("source", "SESAME ASCII2 %s (sha256 %s) material %d table %d; "
         "101: %s" % (os.path.basename(args.src), src_sha[:12],
                      args.mat, table, comment)),
        ("generator", "eos_table_prep.py %s, run %s" % (
            git_sha(), datetime.date.today().isoformat())),
        ("composition", "A=%s Z=%s rho0=%s (SESAME 201)" % (
            raw["abar"], raw["zbar"], raw["rho0"])),
        ("units", "cgs"),
        ("e_shift", "%.10e" % e_shift),
        ("conditioning", cond),
    ]
    if e_ref_val is not None:
        # common-gauge record: stored e = e_raw - e_ref + e_shift, so the
        # emitted table reads e = e_shift exactly at the reference state
        prov.insert(6, ("e_ref_state",
                        "rho=%.10g T=%.10g e_ref=%.10e"
                        % (args.e_ref_state[0], args.e_ref_state[1],
                           e_ref_val)))
    write_eostab(out, prov, lrho, lT, blocks, {"le": le, "lp": lp})
    print("hull coverage %.1f%%; conditioning: %s"
          % (100.0 * (hull > 0.5).mean(), prov[-1][1]))
    if args.qa:
        from .qa import qa_plots
        qa_plots(args.qa, lrho, lT, p, e, cv, hull, tag="_" + base)


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
    s.add_argument("--mass-amu", type=float, default=None,
                   help="mean particle mass in amu (default: deuterium, "
                        "2.0136, with legacy ideal-D provenance)")
    s.add_argument("--material", default="ideal-gas",
                   help="provenance material label (used with --mass-amu)")
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

    s = sub.add_parser("sesame",
                       help="extract + condition a SESAME ASCII2 table "
                            "(doc/eos_sesame_plan.md)")
    s.add_argument("--src", required=True, help="path to the .ascii2 file")
    s.add_argument("--mat", type=int, default=None, help="SESAME material ID")
    s.add_argument("--list", action="store_true",
                   help="list materials in the file (with --grep filter) "
                        "and exit")
    s.add_argument("--grep", default=None,
                   help="with --list: case-insensitive name filter")
    s.add_argument("--table", type=int, default=None,
                   help="SESAME table ID (default: 311 if present, else 301)")
    s.add_argument("--out", default=None,
                   help="output path (default data/<name>_<mat>_s<table>.eostab)")
    s.add_argument("--qa", default=None, help="QA plot output directory")
    s.add_argument("--n-rho", type=int, default=None,
                   help="target grid points (default: source count)")
    s.add_argument("--n-T", type=int, default=None)
    s.add_argument("--rho-min", type=float, default=None,
                   help="trim source densities below this (g/cc); the rho=0 "
                        "column is always dropped")
    s.add_argument("--T-min", type=float, default=None,
                   help="trim source temperatures below this (K); the T=0 "
                        "row is always dropped")
    s.add_argument("--c-cav", type=float, default=None,
                   help="cavitated-response sound speed in km/s (default: "
                        "sqrt(BS/rho0) from the 201 table, else cold-curve "
                        "slope at rho0)")
    s.add_argument("--cs-floor", type=float, default=None,
                   help="opt-in W-B sound-speed safety net in km/s "
                        "(dpdrho >= floor^2 on ALL cells); off by default "
                        "to preserve baseline byte-identity")
    s.add_argument("--e-ref-state", type=float, nargs=2, default=None,
                   metavar=("RHO", "T"),
                   help="common-energy-reference state (g/cc, K): subtract "
                        "e(RHO,T), sampled bilinearly off the finished "
                        "surface with the C++ reader's convention, before "
                        "the positivity shift — step (i) of the mixture-set "
                        "common-gauge recipe")
    s.add_argument("--e-shift", type=float, default=None,
                   help="FORCE the positivity constant S (erg/g) instead of "
                        "the per-table automatic shift — step (ii): one "
                        "shared S across the set. Hard-fails if min(e)+S "
                        "<= 0 anywhere (no silent per-table top-up)")
    s.add_argument("--probe-shift", action="store_true",
                   help="print the post-reference e minima (PROBE-SHIFT "
                        "line: full array + hull) and exit before the "
                        "derivative/inverse-map build; used by the set "
                        "driver to choose the shared S")
    s.add_argument("--floor-mass-amu", type=float, default=None,
                   help="mass (amu) for the thermal stiffness floor kT/m "
                        "(default: 201 abar). Pass the MOLECULAR mass for "
                        "molecular fluids (air 28.97, D2 4.028) — the "
                        "atomic abar over-floors the cold molecular region "
                        "by ~2x and demotes the ambient; with this flag "
                        "only cells < 0.95*kT/m are demoted")
    s.add_argument("--no-verify-src", action="store_true",
                   help="skip the S1 sha256 gate against sources.yaml")
    s.add_argument("--cold-extend", action="store_true",
                   help="T3/T4 (ti_splice_plan_v2): replace the sub-hull "
                        "cold fill with the 306-fitted solid model, capped "
                        "at the 411 melt line, bridged to the hull edge")
    s.add_argument("--melt-frac", type=float, default=0.9,
                   help="solid-model weight ends at this fraction of "
                        "T_m(rho) (default 0.9)")
    s.add_argument("--align-tol", type=float, default=5.0,
                   help="H5 offset-constancy gate: max std(e-offset)/kT "
                        "over the solid<->SESAME overlap (default 5.0)")
    s.add_argument("--lrho", type=float, nargs=3, default=None,
                   metavar=("LO", "HI", "N"),
                   help="prescribed log10-rho axis (Track P); cells outside "
                        "the source span are nearest-filled and hull-0")
    s.add_argument("--lT", type=float, nargs=3, default=None,
                   metavar=("LO", "HI", "N"),
                   help="prescribed log10-T axis (Track P; e.g. the shared "
                        "mixture bracket '1.25 9.0 384')")
    s.set_defaults(func=cmd_sesame)

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
