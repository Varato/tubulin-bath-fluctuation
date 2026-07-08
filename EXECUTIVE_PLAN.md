# Executive Plan — Operational Dashboard

**Project:** Multi-scale site-energy fluctuation analysis of the tubulin
tryptophan (Trp) network. See `RESEARCH_PLAN.md` for the full scientific
specification; this file tracks execution status, decisions made, and
critical parameters flowing downstream.

**Code:** `scripts/` (flat layout, shared `utils.py`)
**Data:** `/Volumes/SandiskData/tubulin-md/noise_analysis/data/`
**Results:** `results/phaseN_<name>/` (one folder per phase)
**Python:** `/Users/alan/miniforge3/envs/qbio/bin/python` (qbio env)

---

## Phase Status

| Phase | RESEARCH_PLAN § | Status | Report |
|---|---|---|---|
| 0 — Validation gate | §1.6 | ✅ done | (inline) |
| 1 — Basic statistics & variance decomp | §2 | ✅ done | `PHASE1_REPORT.md` |
| 2 — Multi-timescale dynamics | §3 | ✅ done | `PHASE2_REPORT.md` |
| 3 — Source-resolved attribution | §4 | ✅ done | `PHASE3_REPORT.md` |
| 4 — Spatial correlation | §5 | ✅ done | `PHASE4_REPORT.md` |
| 5 — Exciton dynamics model | §6 | ⬜ pending | — |
| 6 — Methodological generalization | §7 | ⬜ pending | — |

---

## Phase 0 — Validation (done)

Two gates, both pass:
- **Linearity:** max|E_total − Σ E_sources| / max|E_total| = 4×10⁻¹⁴ (slow),
  8×10⁻¹⁴ (fast). Confirms the source decomposition is numerically exact.
- **Spectral normalization (periodogram):** ∫S df / σ² = 1.0000 (fast),
  0.9999 (slow). Confirms the PSD code.
- **Finding:** Welch underestimates σ² by 6–13% (Hann window + per-segment
  demean). Carry forward to Phase 3: use periodogram as ground truth.

Files: `results/phase0_validation/{linearity,spectral_norm}.csv`

---

## Phase 1 — Basic Statistics (done)

**Task 1.1 — disorder strength.** σ_total ≈ 860–915 cm⁻¹ (mean over sites),
σ_total / J ≈ **19–20** (J = 45 cm⁻¹). Confirms the σ/J ≫ 1 strong-disorder
regime.

**Task 1.2 — variance decomposition (focal result).**
- Screening ratio R_screen = 1 − σ_total²/Σσᵢ² = **0.66 (slow), 0.56 (fast)**.
  Cross-covariances cancel 56–66% of the naive variance sum.
- Protein–water correlation = **−0.39 (slow), −0.23 (fast)** — the dielectric
  screening signature.
- **Implication for Phase 5:** the Monte Carlo sampling in Task 5.2 MUST draw
  from the full 8×8 covariance of `delta_s_total`, not from per-component σ².
  Using Σσᵢ² overestimates σ by ~70%.

**Task 1.3 — Gaussianity.** |skew| < 0.33, |excess kurtosis| < 0.69 across
sites. Formally rejects Gaussianity at N = 2×10⁵ (over-powered test) but the
deviations are small. → ACF/PSD are nearly complete descriptors; Gaussian
sampling in Phase 5 is acceptable. Heaviest tails: αW388, αW346, βW344.

Files: `results/phase1_basic_stats/` (9 files: σ matrix, variance decomp,
screening ratio, covariance matrices, Gaussianity).

---

## Phase 2 — Multi-Timescale Dynamics (done → `PHASE2_REPORT.md`)

**Task 2.1 — stitched ACF.** Crossover at dt_slow = 10 ps. The two trajectories
disagree in their overlap: σ²_slow / σ²_fast = **1.21** (fast traj misses slow
modes). Discontinuity at crossover = 0.08 (small).

**Task 2.2 — fitting (key deviation from plan).** The plan expected
bi-exponential optimality. **It is rejected**: ΔAIC = **−331**, ΔBIC = −324.
Tri-exp reveals three physical timescales:

| component | τ | A | origin |
|---|---|---|---|
| τ₁ | 0.044 ps | 0.50 | water libration |
| τ₂ | 1.70 ps | 0.33 | water rotation / H-bond |
| τ₃ | 1140 ps | 0.16 | protein conformational |

Authoritative τ_slow = **2663 ps** (slow-only fit on slow ACF); the stitched
bi-exp τ_slow = 671 ps is biased by the discontinuity.

**Task 2.3 — stitched PSD.** Periodogram normalization exact (1.0000). Crossover
at 1 cm⁻¹. Shape reliable; absolute level ambiguous by ~20% (mixed σ² refs).

**Task 2.4 — sampling bias.** 10 ps sampling destroys sub-ps timescale info.
The dt/2 trapezoidal spike replaces the true fast integral (~0.4 ps) with a
~5 ps artifact. τ_fast is fundamentally unresolvable below dt. τ_int per-site
bias 0.66–4.74 (largest for fast-dominated sites).

Files: `results/phase2_timescale/` (17 files).

---

## Phase 3 — Source-Resolved Attribution (done)

**Goal (§4):** attribute fluctuations to physical sources (protein, water,
nucleotide, ions) and quantify the dipole-orientation contribution.

**Task 3.1 — component-wise ACF & τ (validates Phase 2 τ assignments).**
Computed per-component ACFs on *both* trajectories (plan said fast-only; but
protein is slow and needs the slow traj). Single-exp + plateau fit, mean over
sites:

| component | τ_fast (ps) | offset | τ_slow (ps) | Phase 2 τ match |
|---|---|---|---|---|
| protein | 1.84 | 0.29 | **2452** | ≈ τ₃ (2663 slow-only) ✓ |
| water | **1.60** | 0.40 | 1172 | ≈ τ₂ (1.70) ✓ |
| nucleotide | 0.31 | 0.35 | 1777 | minor / intermediate |
| ions | 6.04 | 0.75 | 973 | minor / intermediate |

Bi-exp fast-band fit gives τ₁ ≈ 0.12–0.17 ps across all components — sub-ps but
**not** matching Phase 2's τ₁ = 0.044 ps. This was a model-resolution artifact
(bi-exp on the narrow fast band lacks leverage). The **tri-exp per component**
(the same model Phase 2 used, on the stitched ACF) cleanly validates all three:
water τ₁ = **0.054 ps** (≈ 0.044 ✓), τ₂ = **1.12 ps** (≈ 1.70 ✓); protein
τ₃ = **2501 ps** (≈ 2663 ✓). **Conclusion: water owns τ₁ + τ₂, protein owns τ₃.**

**Task 3.2 — component-wise PSD & band attribution.** Stitched per component at
1 cm⁻¹. Band-resolved power fractions (against uncorrelated Σ_comp):

| band | range (cm⁻¹) | protein | water | nucleotide | ions |
|---|---|---|---|---|---|
| slow | 0.02–1 | 26.8% | **54.3%** | 1.3% | 17.6% |
| mid | 1–50 | 39.2% | **52.1%** | 0.7% | 8.0% |
| fast | 50–500 | 32.8% | **63.8%** | 0.4% | 3.1% |
| ultrafast | 500–1668 | 35.3% | **64.5%** | 0.0% | 0.2% |

- Water dominates every band (52–65%); protein is #2 everywhere (27–39%);
  nucleotide negligible at all frequencies (<1.3%, 1.6–3.2 decades below water);
  ions are a real slow-band contributor (17.6%) but negligible above 50 cm⁻¹
  (<3%, 1.3+ decades below water).
- **Cross-spectral check:** ∫Σ PSD_comp / ∫PSD_total = **1.995 ≈ 2.0**. The
  component PSDs sum to *double* the measured total — the protein-water
  anti-correlation causes massive cross-spectral cancellation. This is the
  frequency-domain confirmation of the Phase 1 screening ratio (R_screen ≈
  0.5–0.66); independently consistent with σ²_total/Σσ²_comp² ≈ 0.5.

**Task 3.3 — dipole-orientation contribution.**
- Projection verified to machine precision: δs = **−μ̂·E** (corr = −1.0,
  rel err 5×10⁻¹⁷).
- Indole reorientation τ_μ ≈ **17 ns** mean (3.6–25 ns range; βW397 outlier at
  3.6 ns). Fit on slow traj; fast traj barely decays (C_μ(500 ps) = 0.88).
- **Exciton-timescale relevance:** C_μ(T_obs = 2 ps) = **0.984** → dipole is
  98.4% correlated on T_obs → effectively **static** for exciton dynamics.
- Control test: σ(fixed-μ)/σ(total) = **1.01 (fast), 0.998 (slow)** → freezing
  μ̂ changes variance by <1%. Dipole reorientation contributes ≈ **0%** to
  site-energy noise. **Field fluctuations fully dominate.**
- **Implication for Phase 5:** treat μ̂ as a fixed per-site vector (the
  time-average) in the exciton Hamiltonian. Dipole-orientation dynamics can be
  neglected — only δs = −μ̄·E(t) matters, where E(t) carries all the dynamics.

Files: `results/phase3_source_attribution/` — `comp_acf.png`, `comp_tau.csv`,
`comp_tau_biexp_fast.csv`, `comp_psd.png`, `comp_psd_bysite.png`,
`comp_band_power.csv`, `dipole_acf.png`, `dmu_reorientation.csv`,
`dipole_control.csv`.

---

## Phase 4 — Spatial Correlation (done → `PHASE4_REPORT.md`)

**Result:** spatial correlations are negligible. Mean off-diagonal |r| =
**0.019** (slow traj) — indistinguishable from zero. The independent-site
assumption holds.

**Phase 5 implication:** treat sites as independent. Use per-site σ_total
(diagonal covariance). Covariance matrix at
`results/phase4_spatial_correlation/corr_cov_slow.npz`.

---

## Phase 5 — Exciton Dynamics (pending — see critical parameters below)

**Goal (§6):** map fluctuation statistics → exciton model; separate static
slow disorder from dynamic fast dephasing; verify strong localization.

**Tasks:**
- 5.1: Fast/slow variance split (σ_slow² = A_slow·σ²_total, σ_fast² = A_fast·σ²_total).
- 5.2: Monte Carlo (5000–10000 realizations) of the 8-site exciton Hamiltonian
  with slow disorder → participation ratio (localization).
- 5.3: Markovian dephasing rate γ_φ = σ_fast² · τ_fast; compare τ_φ to τ_J = 74 fs.
- 5.4: Static+dynamic dual-mechanism model vs traditional pure-static model.

---

## Phase 6 — Methodological Generalization (pending)

**Goal (§7):** derive general criteria for when the static-disorder
approximation is valid, and MD sampling-frequency guidelines.

**Inputs from Phase 2:** the sampling-bias result (10 ps destroys sub-ps info)
and the bi-exp rejection (ΔAIC = −331) directly support the dual-parameter
criterion (R_slow, τ_slow/T_obs) proposed in §7.2.

---

## Critical Parameters for Phase 5

These are the load-bearing numbers from Phases 0–4 that Phase 5 depends on.
**(Also saved to memory for cross-conversation recall.)**

### σ_slow and the localization question
- A_slow ≈ 0.16 (tri-exp) → σ_slow ≈ √A_slow · σ_total ≈ **0.40 × σ_total**
- σ_slow / J ≈ 0.40 × 19 ≈ **8 ≫ 1** → strong exciton localization survives
  even after removing all fast dynamics. The core conclusion (§8.2 of
  RESEARCH_PLAN) is robust.
- For the Monte Carlo (Task 5.2): sample from the **full 8×8 covariance of
  `delta_s_total`** (Phase 4 result: well-conditioned, cond = 4.7, nearly
  diagonal with off-diag/diag = −0.003). The independent-site approximation
  (diagonal only) is excellent — but do NOT reconstruct σ from Σσᵢ² — the
  screening ratio is 0.56–0.66.

### τ_fast for the dephasing rate γ_φ = σ_fast² · τ_fast
- The system has **two** sub-ps timescales. The choice matters:
  - τ₁ = 0.044 ps (libration) → γ_φ too small by ~40× vs τ₂.
  - τ₂ = 1.70 ps (water rotation) → **physically relevant** for exciton
    dephasing (τ_J = 74 fs sits in this regime).
- **Recommendation: use τ₂ = 1.70 ps for γ_φ.** τ₁ is too fast to couple
  efficiently to the exciton transfer dynamics.
- σ_fast² ≈ (A₁ + A₂)·σ_total² ≈ 0.84 × σ_total².

### τ_slow for the static criterion (Phase 6)
- τ_slow = **2663 ps** (slow-only fit, authoritative).
- τ_slow / T_obs = 2663 / 2 ≈ **1300 ≫ 1** → slow disorder is quasi-static
  on the exciton timescale. Static-disorder approximation is valid for the
  slow component.

### Sampling for Monte Carlo
- Gaussianity is acceptable (Phase 1.3: |skew| < 0.33, |kurt| < 0.69).
- Watch the heavy-tailed sites (αW388, αW346, βW344) for localization-tail
  statistics.
- Fixed seed (RESEARCH_PLAN §1.6 requires it).

---

## Methodological Decisions Log

| # | Decision | Rationale |
|---|---|---|
| 1 | `qbio` conda env for all Python | has numpy 2.4, scipy 1.17, matplotlib 3.10 |
| 2 | Flat `scripts/` + per-phase `results/` folders | user preference; traceability |
| 3 | Phase 0 validation gate before any analysis | RESEARCH_PLAN §1.6; catches data/code errors early |
| 4 | Periodogram (not Welch) as PSD ground truth | Welch underestimates by 6–13% (Hann windowing) |
| 5 | ACF crossover at dt_slow = 10 ps | slow traj has no data below 10 ps; fast traj biased above |
| 6 | PSD crossover at 1 cm⁻¹ | matches the boundary between slow/fast baths |
| 7 | Report tri-exp as the physical model | ΔAIC = −331 rejects bi-exp; three timescales are physical |
| 8 | Slow-only fit for authoritative τ_slow | stitched-ACF discontinuity biases the joint fit |
| 9 | First-zero-crossing truncation for τ_int | robust to long-lag noise; consistent across ACFs |
