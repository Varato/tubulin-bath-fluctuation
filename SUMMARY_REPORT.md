# Summary Report — Tubulin Trp Network: From Noise Characterization to Exciton Dynamics

**System:** 8 tryptophan (Trp) sites in the tubulin dimer. Two MD trajectories
resolve complementary timescales: fast (10 fs, 0–2 ns) and slow (10 ps, 10–50 ns).
Exciton coupling from Craddock 2014 Hamiltonian; strongest pair Trp4–βW344
(J = −59 cm⁻¹).

---

## The Story in Three Acts

1. **The noise has three timescales** (Phase 2), each owned by a distinct
   molecular process (Phase 3): water libration, water rotation, protein
   conformational dynamics.
2. **The disorder is strong** (Phase 1): σ/J ≈ 19, suggesting Anderson
   localization. Dielectric screening cancels half the variance (Phases 1+3).
3. **But fast noise breaks the localization** (Phase 5): water fluctuations
   drive noise-assisted transport at 50% efficiency, while static disorder
   alone traps at 17%. This is a non-Markovian effect invisible to Lindblad.

---

## Act I — Three Timescales and Their Origins

### The ACF is tri-exponential (Phase 2)

The stitched ACF is decisively tri-exponential (ΔAIC = **−331** vs bi-exp).
Three physical timescales:

| mode | τ | amplitude | origin |
|---|---|---|---|
| τ₁ | **44 fs** | 0.50 | water O–H libration |
| τ₂ | **1.7 ps** | 0.33 | water rotation / H-bond |
| τ₃ | **2.7 ns** | 0.16 | protein conformational |

Authoritative τ_slow = **2663 ps** (slow-only fit); the stitched value (1140 ps)
is biased by the 10 ps crossover discontinuity.

### Source decomposition confirms the assignment (Phase 3)

Fitting the same tri-exp model to each source independently:

| Phase 2 τ | water (Phase 3) | protein (Phase 3) | verdict |
|---|---|---|---|
| τ₁ = 44 fs | **54 fs ✓** | — | water owns τ₁ |
| τ₂ = 1.7 ps | **1.1 ps ✓** | — | water owns τ₂ |
| τ₃ = 2.7 ns | — | **2.5 ns ✓** | protein owns τ₃ |

This mutual validation — temporal decomposition (Phase 2) and source decomposition
(Phase 3) converging from orthogonal directions — confirms the three timescales
are physically real, not fitting artifacts.

### Water dominates every frequency band

| band (cm⁻¹) | water | protein | nucleotide | ions |
|---|---|---|---|---|
| 0.02–1 | **54%** | 27% | 1.3% | 18% |
| 1–50 | **52%** | 39% | 0.7% | 8% |
| 50–500 | **64%** | 33% | 0.4% | 3% |
| 500–1668 | **65%** | 35% | 0% | 0.2% |

---

## Act II — Strong Disorder and Dielectric Screening

### σ/J ≈ 19: localization expected

Per-site σ_total ≈ 675–1453 cm⁻¹ (mean 915). With J = 45–59 cm⁻¹, every site
has σ/J > 11. Even removing fast dynamics (σ_slow ≈ 0.40 × σ_total):

> σ_slow / J ≈ **8 ≫ 1** → strong localization survives in the static limit.

### Dielectric screening cancels half the variance

Protein-water anti-correlation (r = **−0.39**) causes cross-covariances to cancel
**56–66%** of the naive variance sum. The component PSDs sum to **double** the
total (∫Σ PSD_comp / ∫PSD_total = **2.0**) — the frequency-domain confirmation.

Any model treating sources as independent overestimates σ by ~2×.

### Simplifications verified

- **Dipole frozen:** τ_μ ≈ 17 ns, C_μ(2 ps) = 0.984. Only E(t) carries dynamics. (Phase 3)
- **Spatial independence:** mean off-diagonal |r| = 0.019. Sites are independent. (Phase 4)
- **Near-Gaussian:** |skew| < 0.33, |kurt| < 0.69. Gaussian sampling justified. (Phase 1)

---

## Act III — Noise-Assisted Transport (Phase 5)

### Method: exact trajectory Monte Carlo

For each realization: generate per-site colored noise (three OU processes with
per-site Aₖ, τₖ from Phase 2), solve the TDSE with H(t) = J σₓ + ΔE/2 σ_z +
δs(t)·|n⟩⟨n|, average over 300–500 realizations.

This is exact for classical Gaussian noise — no Born-Markov, no perturbation
theory, no rotating-wave approximation. OU noise verified against theory:
ACF RMS = 0.028 over 10 fs – 7.5 ns.

### The result: fast noise enables transfer, static traps

Trp4 → βW344 dimer (J = −59, ΔE = 41 cm⁻¹, per-site noise):

| noise config | P₇ at 2 ps | character |
|---|---|---|
| τ₁ only (libration) | **0.51** | smooth monotonic rise, no oscillation |
| τ₂ only (rotation) | **0.49** | same |
| τ₁+τ₂ | 0.52 | same |
| full (τ₁+τ₂+τ₃) | 0.50 | same |
| τ₃ only (static) | **0.17** | weak oscillation, localized |

**Fast noise drives ~50% incoherent transfer** — the fluctuating energy gap
crosses resonance, and J transfers population during each crossing window.
Static disorder (τ₃) keeps the gap fixed far from resonance → trapped at 17%.

This is **environment-assisted quantum transport** (ENAQT): the noise is not
merely decoherence to be minimized — it is the transport mechanism.

τ₁+τ₃ ≈ full (0.51 vs 0.50): **τ₂ barely matters**. The transport is driven
almost entirely by the fastest component (water libration).

### All-pairs scan: βW344 is the exclusive partner

Using the full Craddock 2014 Hamiltonian, Trp4 transfers only to βW344
(J = −59). All other 6 Trps have |J| ≤ 6 and P(2 ps) < 0.05.

### Lindblad fails by 3×

The Markovian approximation underestimates transfer: Lindblad predicts P₇ = 0.17
vs trajectory MC = 0.49 for τ₁-only. A sigma scan confirms MC matches Lindblad
at κ ≪ 1 and diverges at κ ≫ 1 — the expected Markov breakdown. All Kubo numbers
κₖ = σₖτₖ/ℏ are ≫ 1 (5–10 for τ₁, 165 for τ₂, 10⁵ for τ₃).

---

## Data and Code Review

### Numerical consistency (all verified)

| check | result |
|---|---|
| Tri-exp amplitudes sum to 1 (all 8 sites) | ✓ |
| Band power fractions sum to 100% per band | ✓ |
| Mean τ₁, τ₂, τ₃ from per-site CSVs match reported means | ✓ |
| Mean ΔAIC = −331 (per-site: −84 to −474) | ✓ |
| σ²_slow/σ²_fast = 1.21 recomputed from raw data | ✓ |
| OU sum ACF vs theory (30 ns trajectory): RMS = 0.028 | ✓ |
| Lindblad mesolve matches analytic formula exactly | ✓ |
| MC matches Lindblad at κ ≪ 1 (sigma scan) | ✓ |
| P₄ + P₇ = 1 to machine precision (conservation) | ✓ |
| Noiseless sesolve matches sin²(Jt) | ✓ |

### Code review (no bugs found)

- AIC/BIC formulas, free-parameter counts (k=3 bi-exp, k=5 tri-exp) ✓
- Amplitude constraints (third derived, not free) ✓
- ACF/PSD normalization (Parseval = 1.0000 for periodogram) ✓
- OU exact discretization (analytic SDE solution, not Euler-Maruyama) ✓
- Qutip sesolve with time-dependent H: verified against analytic static case ✓

---

## Phase Reports

| Phase | Report |
|---|---|
| 0 — Validation | (inline in EXECUTIVE_PLAN.md) |
| 1 — Basic statistics | `PHASE1_REPORT.md` |
| 2 — Multi-timescale dynamics | `PHASE2_REPORT.md` |
| 3 — Source-resolved attribution | `PHASE3_REPORT.md` |
| 4 — Spatial correlation | `PHASE4_REPORT.md` |
| 5 — Exciton dynamics | `PHASE5_REPORT.md` |
