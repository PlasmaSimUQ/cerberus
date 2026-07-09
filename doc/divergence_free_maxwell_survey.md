# Divergence-free-by-construction Maxwell updates for multifluid plasma — literature survey

Status: research survey (2026-07-06). Scope: whether constrained-transport-like
(divergence-free by construction) update algorithms exist for the **full Maxwell
equations** in a multifluid/five-moment plasma solver, and how they would map onto
Cerberus (cell-centred, block-structured AMReX finite-volume, CTU hydro).

This is a *landscape* document, not an engineering plan. It records the candidate
algorithms, their data-layout requirements, AMR compatibility, and the trade-off
against the hyperbolic-cleaning approach that essentially every production
multifluid code actually ships. Method: fan-out web search across five angles,
22 primary sources fetched, 240 claims extracted, adversarially verified
(2 of 103 verification votes refuted — both scope nuances, noted inline below).

---

## 0. TL;DR for Cerberus

- **Yes**, exact `div B = 0` **and** `div E = ρ_c/ε₀` by construction exist for the
  *full* Maxwell system coupled to a two-fluid/five-moment plasma — this is the
  **Balsara line** (§3). It is a solved research problem.
- The classic **constrained transport (CT)** you know from MHD (§2) only solves
  `div B = 0` for the ideal-MHD *induction* equation — no displacement current,
  no Gauss's law. It is a **building block**, not a Maxwell solver.
- The **frontier gap** for us specifically: exact schemes either want **staggered
  face/edge storage** (AMR-proven, but a data-model change from Cerberus's
  cell-centred `field` state) **or** are **co-located** (natural AMReX fit) but
  have **not been demonstrated on AMR**.
- **What the peer group ships**: hyperbolic divergence cleaning (Perfectly
  Hyperbolic Maxwell / PHM). Approximate, cheap, cell-centred, AMR-trivial. The
  Shumlak/Hakim/Gkeyll lineage uses this (§5).
- **WarpX** (itself an AMReX code) proves staggered Yee E/B live comfortably in
  the AMReX AMR machinery (§6) — de-risks the "add staggered storage" concern —
  but its constraint recipe is PIC-flavoured (Esirkepov deposition), not a
  finite-volume fluid mechanism.

Recommendation ordering in §7.

---

## 1. Framing: the two constraints and three camps

Maxwell has two involution (divergence) constraints that a scheme must respect:

- **Gauss's law for B:** `∇·B = 0`
- **Gauss's law for E:** `∇·E = ρ_c/ε₀` (with charge density — *not* just `∇·E = 0`)

If the update does not respect these, errors accumulate and, as Hakim–Loverich–
Shumlak note (citing Jiang et al.), for initial-boundary-value problems the
divergence equations are **not** redundant — ignoring them yields spurious
solutions and charge non-conservation.

Three families of remedy appear in the literature:

| Camp | What it constrains | Maturity |
|---|---|---|
| **Classic CT** (§2) | `div B = 0`, induction eq. only | very mature, AMR-proven |
| **Full-Maxwell constraint-preserving** (§3) | both, full Maxwell, plasma-coupled | mature (Balsara), AMR partly open |
| **Hyperbolic cleaning / PHM** (§5) | both, approximately | ubiquitous in production multifluid codes |

---

## 2. Classic constrained transport (div B = 0 only — building blocks)

All induction-equation-only (no displacement current, no Gauss's law for E), so
none is a complete Maxwell solver. They supply reconstruction/AMR machinery the
full-Maxwell schemes reuse.

- **Evans & Hawley (1988)** — the original staggered CT idea (face-centred B,
  edge-centred EMF); `div B = 0` via Stokes' theorem.
- **Balsara & Spicer (1999)**, *JCP* 149:270 — "flux-CT": builds edge EMFs
  directly from the upwinded Godunov fluxes so face B stays solenoidal by
  construction. States plainly that *zone-centred Godunov schemes cannot hold
  `div B = 0` without a divergence-cleaning step.*
  https://www.sciencedirect.com/science/article/abs/pii/S0021999198961538
- **Gardiner & Stone (2005)**, *JCP* 205:509 — single-step 2nd-order **unsplit
  Godunov + CTU + CT**. This is the **same CTU (corner-transport-upwind) family
  Cerberus already uses for hydro**. Corner EMFs built from 1D Riemann fluxes;
  warns PPM and CTU need **nontrivial multidimensional extensions** to marry with
  CT. https://arxiv.org/abs/astro-ph/0501557
- **Londrillo & Del Zanna (2004)**, *JCP* 195:17 — **Upwind Constrained Transport
  (UCT)**: general Godunov-CT framework, `div B = 0` to machine precision; edge
  EMF is a genuine **4-state multidimensional** function. Warns that naive
  cell-centred upwinding without CT-consistent reconstruction produces **O(1)
  numerical magnetic monopoles**.
- **Balsara (2001)**, *JCP* 174:614 — the **AMR** piece; see §4.
  https://ui.adsabs.harvard.edu/abs/2001JCoPh.174..614B/abstract

Key structural takeaway: classic CT stores **B on cell faces, E (EMF) on cell
edges, fluid vars at cell centres**. That staggering is what makes `div B = 0`
exact — and is exactly what a cell-centred code must add.

---

## 3. Full-Maxwell constraint-preserving family (the direct answer)

Balsara and collaborators generalised CT to the *complete* Maxwell system,
including displacement current and Gauss's law with charge density. These are the
"constrained-transport-for-Maxwell" answers.

### 3.1 ⭐ Balsara, Amano, Garain & Kim (2016), *JCP* 318:169 — the canonical match
"*A high-order relativistic two-fluid electrodynamic scheme with consistent
reconstruction of electromagnetic fields and a multidimensional Riemann solver
for electromagnetism.*"
https://ui.adsabs.harvard.edu/abs/2016JCoPh.318..169B/abstract

The paper that most precisely matches the goal. Two-fluid plasma coupled to the
**full set of Maxwell's equations including displacement current**:

- Reconstructs **B divergence-free within each zone** *and* **E consistent with
  Gauss's law** (`div E = ρ_c/ε₀`) within each zone — both constraints by
  construction, not by cleaning.
- Updates B via a **discrete Faraday's law** and E via a **discrete generalised
  Ampère's law**, driven by a **multidimensional Riemann solver for
  electromagnetism** that supplies edge-centred E (for the Stokes-law B update)
  and edge-centred B (for the E update).
- Handles the charge source honestly: charge is a zone-averaged FV quantity, and
  the higher moments of `div E` are constrained to match the reconstructed
  moments of charge density (shown at 2nd/3rd/4th order). **The fluid solver's
  numerical fluxes supply a self-consistent face-centred current density** for
  Ampère's law — the fluid-coupling half that PIC codes get from particle
  deposition.
- Couples to the FV fluid solver at high order via **IMEX** time integration
  (implicit on the stiff Lorentz/current source, explicit on the flux);
  **explicitly extensible to multiple ion species** — i.e. beyond two fluids to a
  general five-moment multifluid set.
- Cost premium over a plain zone-centred Godunov scheme: **< 7% CPU/step**.
- Validated to design accuracy in both the MHD limit and the EM-wave-propagation
  limit; runs a relativistic GEM reconnection problem.

**Requirements / caveats:** **face-collocated (Yee-type) primary E and B** plus
**edge-based multidimensional Riemann solves** — staggered face/edge storage, not
pure cell-centred. **The paper does not present an AMR treatment** (cites the
divergence-free-AMR machinery rather than demonstrating it).

### 3.2 ⭐ Kumar, Chandrashekar & Balsara (2024), *JSC* 101:46 — the co-located variant
"*Second order divergence constraint preserving entropy stable finite difference
schemes for ideal two-fluid plasma flow equations.*"
https://link.springer.com/article/10.1007/s10915-024-02685-0 ·
https://arxiv.org/abs/2409.16004

Arguably the **more AMReX-friendly** option:

- Governs exactly our system: ideal two-fluid (ion + electron) with a **linear
  Maxwell flux coupled only through source terms** — structurally identical to
  Cerberus's five-moment formulation (three flux blocks: ion, electron, Maxwell;
  coupling via sources).
- Preserves **both** E and B divergence constraints by construction, via a
  **multidimensional Riemann solver evaluated at cell vertices**.
- **Co-located (cell-centred, non-staggered)** while still guaranteeing exactly
  divergence-free B, entropy stability, 2nd-order accuracy. Directly refutes the
  assumption that constraint preservation *requires* a staggered Yee mesh.
- **Explicit and IMEX** time integration for the stiff source coupling.
- **Directly benchmarks against PHM cleaning and against no cleaning** — gives
  head-to-head trade-off data, not just assertion.

Relativistic follow-up: **arXiv:2503.20372**, "*Second order divergence constraint
preserving schemes for two-fluid relativistic plasma flow equations*" — same
vertex-based multidimensional-Riemann mechanism, modular w.r.t. any
consistent/stable (entropy-stable) fluid discretisation.

**Caveats:** these are **finite-difference** (Cerberus is finite-volume — usually
a modest reinterpretation since both are conservative and share the Riemann
machinery), demonstrated in **1D/2D**, and **AMR at refinement boundaries is not
addressed**.

### 3.3 Other full-Maxwell constraint-preserving schemes

- **Balsara & Käppeli (2017)**, *JCP* — "*Computational electrodynamics in
  material media… Part I: 2nd-order FVTD.*" Explicit **synthesis of FDTD/Yee
  staggering + Godunov upwinding**: face-collocated normal D and B, edge-collocated
  E and H, multidimensional (MuSIC) Riemann solvers at edges. Preserves both
  Gauss-law constraints exactly **including `div D = ρ_E` with nonzero charge**.
  Uses a **single control volume** (vs FDTD's two linked ones); authors state this
  makes it **well-suited to AMR** (promised in a sequel). Descends from an earlier
  Balsara et al. formulation "*focused on electrodynamics in plasma.*"
  https://www.sciencedirect.com/science/article/abs/pii/S0021999117305326
- **Hazra, Chandrashekar & Balsara (2019)**, *JCP* — "*Globally constraint-preserving
  FR/DG for Maxwell at all orders.*" DG up to **5th order**, mimetically preserving
  both involutions; dispersion error ~75× smaller than FDTD.
  ⚠️ **Verification flag (refuted-as-overclaim):** *this scheme as published
  handles only the charge-free case (`div D = 0`); nonzero charge density and
  nonzero conductivity are explicitly deferred to future work* (quote from §3 of
  arXiv:1809.03816). So it is **not** directly usable for a charged plasma where
  `div E = ρ_c/ε₀` — despite being cited as a "constraint-preserving Maxwell"
  result. Also: at **4th order and above it must additionally evolve some
  cell-centred modes** (hybrid layout).
  https://www.sciencedirect.com/science/article/abs/pii/S0021999119304048
- **Non-staggered central FV Maxwell on AMR** (arXiv:1511.04994) — cell-centred FV
  central scheme preserving both `div B = 0` and `div E = 4πρ` discretely via
  **vertex-based auxiliary-potential flux redistribution** (deliberately chosen
  *over* GLM cleaning), with IMEX for stiff plasma currents. **Reflection-free at
  coarse-fine interfaces** (Yee/FDTD produces spurious reflections there). As of
  that paper, **not yet coupled to a multifluid/PIC solver** (future work) — a
  candidate, not a demonstration.
- **DDFV / covolume mimetic 2D Maxwell on arbitrary non-conforming meshes**
  (*JCP* 2008, S0021999108002945) — preserves Gauss's law + `div B = 0`,
  energy-conserving, 2nd-order on non-conforming/distorted meshes (AMR-like), but
  **2D only, primal-dual mesh** (mismatch with plain cell-centred `MultiFab`s),
  no fluid coupling.

---

## 4. AMR compatibility (the AMReX-specific question)

**Staggered CT + AMR is a solved problem — with machinery:**

- **Balsara (2001)**, *JCP* 174:614 — the foundational reference: divergence-free
  **prolongation** (coarse→fine reconstruction that reduces to TVD in 1D) *and*
  **restriction** operators for face-centred fields, plus an **electric-field
  (edge-EMF) correction at coarse-fine boundaries that keeps `div B = 0` exact
  even under subcycling** (refined levels taking smaller timesteps). Validated on
  3D AMR-MHD in the RIEMANN framework.
- **CHARM** (Miniati & Martin 2011, *ApJS* 195:5) — production CT+AMR via face
  restriction/prolongation + a **"reflux-curl"** operation at refinement
  boundaries — the CT analogue of AMReX's flux refluxing.
  https://iopscience.iop.org/article/10.1088/0067-0049/195/1/5
- **BHAC** (Olivares et al. 2019, *A&A*) — coordinate-independent family of
  divergence-preserving prolongation operators for face fields (Tóth & Roe 2002 a
  special case).
- **WarpX** — staggered-Yee Maxwell + mesh refinement in an AMReX code (§6).

⚠️ **The key AMR warning (BHAC):** *cell-centred ("flux-CT") variants of CT are
**incompatible with AMR** because the face-representation information needed to
enforce the constraint across coarse-fine interfaces is lost.* The classic route
to AMR-safe exact `div B = 0` **requires staggered face/edge storage** — in
AMReX, face- and edge-centred `MultiFab`s plus divergence-preserving prolong/
restrict + an edge-EMF sync at coarse-fine boundaries.

**Open question for the co-located 2024 schemes (§3.2):** they sidestep staggering
but **neither 2024/2025 paper demonstrates AMR** — whether their vertex-based
constraint survives coarse-fine interfaces cleanly is unproven. This is the
genuine research gap for "exact constraints *and* co-located AMReX data."

---

## 5. What multifluid codes actually use: hyperbolic cleaning (PHM)

Essentially every established five-moment/two-fluid code uses **hyperbolic
divergence cleaning**, not exact CT.

- **Hakim, Loverich & Shumlak (2006)**, *JCP* — the foundational five-moment
  two-fluid solver: uses **mixed-potential** *or* **PHM** cleaning; **PHM chosen
  for all production runs**. Cell-centred LeVeque wave-propagation FV — same data
  model as Cerberus.
  https://www.sciencedirect.com/science/article/abs/pii/S0021999106001707
- **Gkeyll Moments app** — production five-/ten-moment multifluid-Maxwell → **PHM**,
  cell-centred high-resolution wave-propagation FV.
  https://gkyl.readthedocs.io/en/latest/gkyl/App/Moment/Moments.html
- **Abgrall & Kumar (2013)**, *JSC* — cell-centred MUSCL two-fluid, **PHM cleaning**
  + positivity preservation.
  https://link.springer.com/article/10.1007/s10915-013-9809-6

**How PHM works** (Munz et al. 2000): augment Maxwell with two correction
potentials ψ (magnetic) and φ (electric):

```
∂B/∂t + ∇×E + γ∇ψ = 0            (1/χ)∂φ/∂t + ∇·E = ρ_c/ε₀
ε₀μ₀ ∂E/∂t − ∇×B + χ∇φ = −μ₀J    (ε₀μ₀/γ)∂ψ/∂t + ∇·B = 0
```

Divergence *errors* become extra hyperbolic waves at speeds `±cγ, ±cχ` that advect
and damp the error away. 1D eigenvalues become `{−cγ, cγ, −cχ, cχ, −c, −c, c, c}`
(vs `{−c,−c,c,c,0,0}` for unmodified Maxwell — the two static zero-speed modes are
the un-propagating divergence error PHM targets). **Constraints hold exactly only
in the limit γ, χ → ∞**; in practice γ, χ = c or 2c → constraints held only
*approximately* (to the scheme's order of accuracy). It is a hyperbolic
generalisation of Hodge projection.

**Trade-offs — exact (CT / Balsara-Maxwell) vs PHM cleaning:**

| | Exact | PHM cleaning |
|---|---|---|
| Constraint | machine precision by construction | approximate (order-accurate) |
| Data layout | staggered face/edge (classic) — *or* co-located (2024) | **cell-centred — trivially fits AMReX** |
| AMR | solved but needs special prolong/restrict + edge-EMF sync; **cell-centred flux-CT is AMR-incompatible** | **couples to AMR straightforwardly** |
| Implementation cost | high (descriptors, reconstruction, multidim Riemann, AMR operators) | **low — two extra scalar fields + standard Riemann solver** |
| Known failure modes | — | spurious divergence can **persist in periodic geometries and at stagnation points**; adds error-wave dynamics |

⚠️ **Caution against over-claiming CT superiority:** the often-cited Balsara & Kim
(2004) "*all divergence-cleaning schemes show deficiencies*" result was **flagged
in verification as an overgeneralisation** — that study tested only **Poisson/
projection (elliptic) cleaning**, *not* hyperbolic GLM/PHM. It is **not** evidence
that exact CT beats PHM for our use case. The fair, directly-relevant comparison
is the head-to-head in **Kumar-Balsara (2024)** (§3.2).

---

## 6. WarpX confirmation (staggered Yee on AMReX AMR)

Checked directly against WarpX docs and the `BLAST-WarpX/warpx` repository, since
it is itself an AMReX code and a useful de-risking data point.

- **Staggered Yee grid — confirmed.** Docs: *"quantities are given on a staggered
  (or 'Yee') grid, where the electric field components are located between nodes
  and the magnetic field components are located in the center of the cell faces"* —
  i.e. **E edge-centred, B face-centred**. Source: `Source/FieldSolver/
  FiniteDifferenceSolver/` has `EvolveE.cpp`, `EvolveB.cpp` (curl FDTD updates),
  CKC/Lehe stencil variants, plus a PSATD spectral solver.
- **Constraint preservation — by construction (Yee/FDTD mechanism):**
  - `div B = 0` preserved automatically by the staggered curl update
    (`div(curl E) ≡ 0` discretely).
  - `div E = ρ/ε₀` preserved via **charge-conserving current deposition
    (Esirkepov)** — the discrete continuity equation is made exact, so the Yee
    Ampère update keeps Gauss's law satisfied.
  - Optional explicit cleaning also present: `EvolveF.cpp`/`EvolveG.cpp`
    (**Langdon–Marder** F/G correction fields), and a PSATD-JRhom *propagative
    divergence cleaning* option.
- **AMR-ready — confirmed.** Built on `amrex::AmrCore`; staggered Yee fields on
  **every refinement level**, split into **coarse-patch (`cp`) / fine-patch (`fp`)**
  `MultiFab`s per MR level. Coarse-fine handling uses an **additive substitution
  method** `F(a) = F(r) + I[F(s) − F(c)]` with PMLs terminating each patch.

**Caveats — why WarpX is only a partial analogue for Cerberus:**

- WarpX is a **kinetic PIC** code, not a fluid code. Its Gauss's-law preservation
  rests on **Esirkepov particle→grid deposition** — a *particle* mechanism. A
  five-moment **fluid** solver has no particles; it needs a charge-conserving
  *flux* for the current density instead (exactly what Balsara-2016 does with the
  fluid fluxes as a self-consistent face current). WarpX validates "staggered Yee
  EM on AMReX AMR works," **not** the fluid-coupling half.
- Its mesh refinement is **additive substitution** (Vay et al.), **not** a
  divergence-preserving prolongation in the Balsara-2001 / CHARM reflux-curl
  sense. Docs note it can produce *"spurious multipole fields from imperfect
  cancellation between patches,"* mitigated with transition cells. So WarpX is
  AMR-*capable* with staggered EM, but its MR does not claim machine-precision
  constraint preservation across coarse-fine interfaces the way CT-MHD codes do.

**Net:** WarpX is solid proof-of-existence that **edge/face-staggered E/B fields
live comfortably in AMReX AMR** (multi-level face/edge `MultiFab`s, per-level
solves, coarse-fine handling) — de-risking the "add staggered storage to AMReX"
concern. But WarpX's specific recipe (Yee + Esirkepov) is PIC-flavoured; the
finite-volume fluid analogue we want is the Balsara-2016 line.

Sources: WarpX PIC theory (Yee grid, charge conservation)
https://warpx.readthedocs.io/en/25.11/theory/pic.html ; WarpX AMR theory
https://warpx.readthedocs.io/en/latest/theory/amr.html ; WarpX source
https://github.com/BLAST-WarpX/warpx/tree/development/Source/FieldSolver/FiniteDifferenceSolver

---

## 7. Recommendation ordering for Cerberus

Given Cerberus is a **cell-centred AMReX five-moment multifluid FV solver with
CTU-based hydro actions**:

1. **Least disruption to the data model, exact constraints:** study
   **Kumar/Chandrashekar/Balsara (2024)** + relativistic follow-up (2503.20372) —
   **co-located, vertex multidimensional Riemann, IMEX**, explicitly our two-fluid
   system, explicitly benchmarked vs PHM. Unknowns to solve: **AMR at refinement
   boundaries** (unaddressed) and FD→FV port.
2. **Most complete, battle-tested full-Maxwell + multifluid, if staggered storage
   is acceptable:** **Balsara et al. (2016)** — face-collocated E/B, edge
   multidimensional Riemann, multi-species-extensible, self-consistent current
   from fluid fluxes. Add face/edge `MultiFab`s and reuse **Balsara (2001) +
   CHARM reflux-curl** for AMR. WarpX proves staggered-EM-on-AMReX-AMR is viable
   in-framework.
3. **Pragmatic default (what the field does):** implement **PHM cleaning** —
   cheap, cell-centred, AMR-clean, the acknowledged production choice of the
   Shumlak/Hakim/Gkeyll lineage. Give up exactness for a straightforward path:
   two extra scalar potentials + a modified flux on the existing `field` state.

**Strategic read:** exact div-free-by-construction Maxwell for multifluid plasma
is a solved research problem (Balsara line), but coupling it to **cell-centred
block-structured AMR specifically** remains the frontier — exact schemes either
want staggering (AMR-proven, data-model change) or are co-located but
AMR-undemonstrated. PHM is the low-risk default the peer group ships.

---

## 8. Staggered-grid implementation of constraint-preserving Maxwell's

This section captures the concrete pathway discussed for putting a *staggered*
(Yee-type) constraint-preserving Maxwell solver into Cerberus's AMReX data model —
i.e. option 2 of §7 (Balsara 2016 line), which is the best-tested route for AMR.
It is an engineering sketch, not a committed design.

### 8.1 Why staggered, and what "best-tested for AMR" actually means

Taking as given that swapping the `field` state from cell-centred to staggered
storage is acceptable (co-located-vs-staggered *not* a decision factor), the
question becomes: which exact scheme has the most mature AMR track record?

- The **best-*tested* AMR constraint machinery is the CT-MHD lineage** —
  **Balsara (2001)** divergence-free prolongation/restriction + edge-EMF
  coarse-fine sync, with **CHARM** (Miniati & Martin 2011) the most complete
  production implementation and **BHAC** (Olivares 2019) an independent check.
  Caveat: this lineage is proven for **`div B = 0` only** (ideal-MHD induction) —
  no displacement current, no Gauss's law for E.
- **No full-Maxwell constraint-preserving scheme has published AMR results.**
  Balsara (2016) and Balsara & Käppeli (2017) present the *uniform-grid* solver
  and cite the divergence-free-AMR machinery rather than demonstrating it.
- Therefore the **best bet** is Balsara (2016)/Balsara & Käppeli (2017) precisely
  because they are built on the **same face-staggered Stokes'-law update** as
  CT-MHD, so they should **inherit the Balsara-2001/CHARM AMR operators** for the
  B-field half largely unchanged. The genuinely new AMR work is the **E-field /
  Gauss's-law** analogue with a fluid-supplied charge density (§8.5).
- **WarpX** is the best-tested staggered-Maxwell-*on-AMReX implementation*, but its
  mesh refinement is **additive substitution** (`F(a)=F(r)+I[F(s)−F(c)]`), which
  can leave **spurious multipole fields** at coarse-fine interfaces — *not*
  machine-precision constraint preservation. It de-risks the storage question, not
  the constraint-across-refinement question.

### 8.2 The staggered data layout in AMReX terms

Classic CT stores, and Balsara-2016 reuses:

- **B on cell faces** — in AMReX, `AMREX_D_DECL` face-centred `MultiFab`s
  (`IndexType` with one `NODE` direction), i.e. `Bx` on x-faces, `By` on y-faces,
  `Bz` on z-faces. Already a first-class AMReX index type (used for fluxes).
- **E / EMF on cell edges** — edge-centred `MultiFab`s (two `NODE` directions in
  3D; nodal in 2D). Also a native AMReX index type.
- **Fluid variables (and charge density) at cell centres** — unchanged from today.

In 2D the edge quantity collapses to a **node-centred** scalar `MultiFab` (the
out-of-plane EMF `Ez`), which is the simplest first target. So the storage change
is "add face- and edge-typed `MultiFab`s to the `field` state descriptor," not a
new framework capability — AMReX already carries these index types.

### 8.3 The update (uniform grid): Stokes'-law Faraday + Ampère

- **B update (Faraday, discrete Stokes' law):** the change of magnetic flux
  through a face equals the circulation of the edge EMF around that face's four
  edges. Because it is a discrete curl, `div B` (sum of face fluxes over a cell)
  is preserved **identically** — the defining CT property.
- **E update (generalised Ampère):** symmetric — change of `E`/`D` flux driven by
  circulation of edge `B`/`H` plus the current-density source. Balsara-2016
  constrains the higher moments of `div E` to match the reconstructed moments of
  the zone charge density, so `div E = ρ_c/ε₀` holds by construction.
- **Edge quantities come from a multidimensional Riemann solver** at edges (2D) /
  vertices (their co-located variant) — not 1D face Riemann fluxes. This is the
  single biggest new numerical component; Cerberus's existing 1D Riemann solvers
  do not suffice for the edge EMF.
- **Current density:** the fluid solver's numerical fluxes supply a self-consistent
  face-centred current for Ampère's law — the FV analogue of PIC charge deposition.
- **Stiff Lorentz/current coupling:** IMEX time integration (implicit on the source,
  explicit on the flux), consistent with how Cerberus already treats stiff plasma
  source terms.

### 8.4 AMR: the reflux-curl operator and its AMReX mapping

**What CHARM is.** CHARM (Miniati & Martin 2011, *ApJS* 195:5) is a
constrained-transport MHD code on **block-structured AMR** — built on **Chombo**,
LBNL's block-structured AMR framework and a close sibling of AMReX (same
coarse-fine bookkeeping model, so its operators map almost one-to-one). Layout:
face-B, edge-E(EMF), cell-fluid; `div B = 0` to machine precision via Stokes' law.

**What the reflux-curl operator is.** It is the CT analogue of the ordinary
**refluxing** correction that AMReX already performs for conserved quantities:

- *Ordinary reflux (AMReX `FluxRegister` today):* at a coarse-fine boundary a
  coarse cell was advanced with a coarse-level face flux, while the fine grid
  computed a more accurate flux across the same interface. The register stores the
  mismatch `F_coarse − Σ F_fine` and adds it back to the adjacent **cell-centred**
  conserved quantity → conservation restored.
- *Reflux-curl (the CT version):* the managed "flux" is the **edge-EMF**, and the
  corrected quantity is **face-B**. The coarse edge-EMF disagrees with the sum of
  fine sub-edge EMFs tiling it: `ΔE_edge = E_coarse_edge − Σ E_fine_subedges`. The
  correction applied to the bounding faces is the **curl of that EMF mismatch**,
  `ΔB_face = curl(ΔE_edge)`. Because the correction is itself a discrete curl, it
  carries **zero divergence** — so it repairs the coarse-fine flux mismatch *while
  leaving `div B = 0` exactly intact*. A naive additive B correction would fix the
  values but reintroduce numerical monopoles at the boundary; the curl form is what
  avoids that.

**One line:** *reflux-curl = refluxing where the register accumulates edge-EMF
mismatches instead of face-flux mismatches, and the correction is applied as the
curl of that register onto face-B — so the coarse-fine interface stays both
conservative and divergence-free to machine precision.*

**Mapping onto AMReX.** The scaffolding largely exists:

- AMReX has `FluxRegister` (face) and `EdgeFluxRegister` / `YAFluxRegister`
  machinery plus the coarse-fine iteration used by refluxing. Reflux-curl is a
  register that **accumulates edge quantities** and applies a **curl stencil**
  rather than a straight additive face correction.
- **Balsara (2001)** supplies the other two operators: **divergence-free
  prolongation** (coarse→fine face-B interpolation that stays solenoidal; reduces
  to TVD in 1D) and **restriction** (fine→coarse face averaging that preserves the
  flux), plus the **edge-EMF sync under subcycling** (refined levels taking smaller
  timesteps) — all directly portable to AMReX face/edge `MultiFab`s.
- WarpX already demonstrates multi-level face/edge `MultiFab`s with per-level
  solves and coarse-fine handling in-framework — proof the storage and iteration
  work; only the *constraint-preserving* correction differs.

### 8.5 The remaining gap (the genuinely new work)

The B-field reflux-curl is a solved problem (CHARM / Balsara 2001). The open piece
for a *full-Maxwell* multifluid solver is the **E-field / Gauss's-law analogue**:
an edge-`B`-register reflux-curl for the Ampère update **together with** a
coarse-fine treatment of the **fluid-supplied charge density** so that
`div E = ρ_c/ε₀` is preserved across refinement boundaries — not just `div B = 0`.
No published scheme demonstrates this on AMR. Balsara & Käppeli (2017) argue their
single-control-volume FVTD formulation is well-suited to AMR (promised in a sequel)
and is the most likely template, but it remains to be shown.

### 8.6 Suggested implementation ordering

1. **Add staggered storage** to the `field` state: face-B and edge/nodal-E
   `MultiFab`s alongside (or replacing) the cell-centred field variables. Start in
   **2D**, where the EMF is a single node-centred scalar `Ez`.
2. **Uniform-grid CT update** first: discrete Faraday (B) + Ampère (E) via Stokes'
   law, with the multidimensional edge Riemann solver. Validate `div B` and
   `div E` to machine precision on a single level (EM-wave + MHD-limit tests).
3. **Fluid coupling:** wire the self-consistent face current from the fluid fluxes
   into Ampère; IMEX for the stiff Lorentz/current source.
4. **AMR for the B half:** port Balsara-2001 prolong/restrict + the reflux-curl
   register (B-field) onto AMReX `FluxRegister`/`EdgeFluxRegister`. Validate
   `div B = 0` across coarse-fine interfaces and under subcycling.
5. **AMR for the E half (research):** the Gauss's-law reflux-curl analogue with
   fluid charge density (§8.5) — the novel contribution.

Fallback at any step: **PHM cleaning** (§5) remains the cell-centred, AMR-clean,
low-risk default the peer group ships.

---

## Appendix: primary sources

Full-Maxwell / multifluid constraint-preserving:
- Balsara, Amano, Garain & Kim (2016), *JCP* 318:169 — https://ui.adsabs.harvard.edu/abs/2016JCoPh.318..169B/abstract
- Kumar, Chandrashekar & Balsara (2024), *JSC* 101:46 — https://link.springer.com/article/10.1007/s10915-024-02685-0 · https://arxiv.org/abs/2409.16004
- Relativistic two-fluid follow-up (2025) — https://arxiv.org/abs/2503.20372
- Balsara & Käppeli (2017), *JCP*, FVTD Part I — https://www.sciencedirect.com/science/article/abs/pii/S0021999117305326
- Hazra, Chandrashekar & Balsara (2019), FR/DG — https://www.sciencedirect.com/science/article/abs/pii/S0021999119304048 · https://arxiv.org/pdf/1809.03816 (charge-free only)
- Non-staggered central FV Maxwell on AMR — https://arxiv.org/pdf/1511.04994
- DDFV 2D Maxwell, non-conforming meshes (2008) — https://www.sciencedirect.com/science/article/abs/pii/S0021999108002945

Classic CT (div B = 0, foundations):
- Balsara & Spicer (1999), *JCP* 149:270 — https://www.sciencedirect.com/science/article/abs/pii/S0021999198961538
- Gardiner & Stone (2005), *JCP* 205:509 — https://arxiv.org/abs/astro-ph/0501557
- Londrillo & Del Zanna (2004), UCT, *JCP* 195:17 — https://www.researchgate.net/publication/222562319
- Balsara (2001), Divergence-Free AMR, *JCP* 174:614 — https://ui.adsabs.harvard.edu/abs/2001JCoPh.174..614B/abstract
- Balsara & Kim (2004), cleaning-vs-CT intercomparison — https://arxiv.org/pdf/astro-ph/0310728 (tested projection cleaning only)

CT + AMR:
- Miniati & Martin (2011), CHARM, *ApJS* 195:5 — https://iopscience.iop.org/article/10.1088/0067-0049/195/1/5
- Olivares et al. (2019), BHAC, *A&A* — https://www.aanda.org/articles/aa/full_html/2019/09/aa35559-19/aa35559-19.html
- WarpX relativistic reconnection w/ MR — https://arxiv.org/pdf/2408.08960

Hyperbolic cleaning / multifluid practice:
- Hakim, Loverich & Shumlak (2006), *JCP* — https://www.sciencedirect.com/science/article/abs/pii/S0021999106001707
- Abgrall & Kumar (2013), *JSC* — https://link.springer.com/article/10.1007/s10915-013-9809-6
- Gkeyll Moments app docs — https://gkyl.readthedocs.io/en/latest/gkyl/App/Moment/Moments.html
- Hakim, Maxwell eigensystem / PHM notes — https://ammar-hakim.org/sj/maxwell-eigensystem.html
- Munz et al. (2000), *JCP* 161:484 — hyperbolic divergence correction (via citations above)
