# Phase 5 Report: Exciton Dynamics — Noise-Assisted Transport in the Trp4–Trp7 Dimer

**System:** two-site dimer Trp4 (αW407) – Trp7 (βW397), J = −59 cm⁻¹.
σ₄ = 999 cm⁻¹, σ₇ = 1453 cm⁻¹. Tri-exp noise from Phase 2.

**Method:** stochastic Schrödinger equation — 500 noise-trajectory Monte Carlo,
component ablation (each timescale isolated). Scripts in `scripts/phase5_task1_exciton_mc.py`,
`scripts/test_ou_sum_acf.py`.

---

## The finding

Fast environmental noise **enables** exciton transfer. Static disorder
**prevents** it. The hierarchy:

| config | P₇ at 2 ps | character |
|---|---|---|
| τ₁ only (libration) | **0.499** | smooth monotonic rise, no oscillation |
| τ₁+τ₂ (all dynamic) | **0.506** | same |
| full (τ₁+τ₂+τ₃) | 0.478 | same, slightly suppressed |
| τ₂ only (rotation) | 0.412 | same, slower |
| τ₃ only (static) | **0.092** | weak damped oscillation, trapped |

Noiseless Rabi transfer reaches 0.63 at 200 fs but oscillates back.

The fast-noise curves (τ₁, τ₂, τ₁+τ₂, full) share a universal shape:
rapid monotonic rise to a plateau near 0.5, no oscillation, log-like
saturation. The static-disorder curve (τ₃) shows weak oscillation below 0.2.

---

## Method: stochastic Schrödinger Monte Carlo

### The algorithm

For each noise realization $i = 1, \dots, 500$:

**Step 1 — generate noise trajectories.** Two independent colored-noise
time series, one per site, each a sum of three Ornstein-Uhlenbeck
processes:

$$\delta s_k^{(i)}(t) = \sum_{j=1}^{3} x_{k,j}^{(i)}(t), \qquad
\langle x_{k,j}(0)\, x_{k,j}(t) \rangle = A_j \sigma_k^2\, e^{-t/\tau_j}$$

where $k \in \{4, 7\}$ indexes the site, and $(A_j, \tau_j)$ are the Phase 2
tri-exp parameters. The OU discretization is **exact** (analytic solution of
the linear SDE, not Euler-Maruyama):

$$x[n+1] = x[n]\, e^{-\Delta t/\tau} + \sqrt{A_j \sigma_k^2\,(1 - e^{-2\Delta t/\tau})}\;\mathcal{N}(0,1)$$

verified in `test_ou_sum_acf.py`: the sum's ACF matches the theoretical
$\sum A_j e^{-t/\tau_j}$ to RMS = 0.028 over 10 fs – 7.5 ns.

**Step 2 — construct the time-dependent Hamiltonian.**

$$H^{(i)}(t) = \underbrace{J\,\sigma_x}_{\text{coupling}} +
\underbrace{\delta s_4^{(i)}(t)\,|0\rangle\langle 0| + \delta s_7^{(i)}(t)\,|1\rangle\langle 1|}_{\text{classical noise}}$$

In matrix form (site basis, $|0\rangle$ = Trp4, $|1\rangle$ = Trp7):

$$H^{(i)}(t) = \begin{pmatrix} \delta s_4^{(i)}(t) & J \\ J & \delta s_7^{(i)}(t) \end{pmatrix}$$

**Step 3 — solve the Schrödinger equation** (unitary, time-dependent):

$$i\hbar\frac{d}{dt}|\psi^{(i)}(t)\rangle = H^{(i)}(t)\,|\psi^{(i)}(t)\rangle,
\qquad |\psi^{(i)}(0)\rangle = |0\rangle$$

Implemented via Qutip `sesolve` with an adaptive Runge-Kutta (RK45)
integrator. The noise arrays are linearly interpolated at the solver's
internal sub-steps (which may be smaller than $\Delta t$).

**Step 4 — extract observables** per realization, then ensemble-average:

$$\langle P_7(t) \rangle = \frac{1}{N_\text{real}} \sum_{i=1}^{N_\text{real}}
\langle \psi^{(i)}(t) | 1 \rangle \langle 1 | \psi^{(i)}(t) \rangle$$

The standard error of the mean is $\text{SEM} = \text{std}_i(P_7^{(i)}) / \sqrt{N_\text{real}}$.

### Parameters

| parameter | value | justification |
|---|---|---|
| $\Delta t$ | **2 fs** | resolves $\tau_1 = 44$ fs with 22 points/decay |
| $T_\text{max}$ | 2 ps | exciton observation window ($T_\text{obs}$) |
| $N_\text{steps}$ | 1001 | — |
| $N_\text{real}$ | 500 | SEM $\approx$ 0.01 for $P_7$ |
| $\hbar$ | 1 | natural units; energy in rad/ps |

$\Delta t = 2$ fs is the **output grid** only. The RK45 solver takes
internal sub-steps as small as needed for the ODE tolerance. The fastest
scale in $H(t)$ is $\sigma_7 / \hbar = 274$ rad/ps → phase accumulation
$0.55$ rad per output step — well within the adaptive solver's reach.

### Why this approach is exact (for our assumptions)

The MD electric-field fluctuations are **classical** $c$-number trajectories.
They are not quantum bath operators — no back-action, no entanglement with
the bath. The correct treatment of classical noise on a quantum system is:

1. For each noise realization, evolve the quantum state under the
   resulting $H(t)$ — this is exact unitary dynamics, no approximation.
2. Average over noise realizations — this gives the reduced density matrix
   $\rho(t) = \overline{|\psi^{(i)}(t)\rangle\langle\psi^{(i)}(t)|}$.

No Born-Markov approximation. No perturbation theory in $\sigma/J$. No
rotating-wave approximation. The only assumptions are:

- **Gaussianity:** the noise is Gaussian (Phase 1: $|\text{skew}| < 0.33$,
  $|\text{excess kurt}| < 0.69$ → acceptable).
- **ACF model:** the tri-exp captures the noise statistics (Phase 2:
  $\Delta\text{AIC} = -331$ vs bi-exp → three timescales are physical).
- **Site independence:** $\langle \delta s_4(t)\, \delta s_7(t')\rangle = 0$
  for $t \neq t'$ (Phase 4: mean off-diagonal $|r| = 0.019$ → negligible).

### Why simpler approaches fail here

| method | requires | our system | verdict |
|---|---|---|---|
| Lindblad (Markovian) | $\kappa = \sigma\tau/\hbar \ll 1$ | $\kappa = 5.3,\; 165,\; 1.8\times10^5$ | ✗ wrong lineshape for all three |
| Redfield (perturbative) | $\sigma \ll J$ | $\sigma/J \approx 17\text{–}25$ | ✗ not perturbative |
| Static MC (freeze all noise) | $\tau \gg T_\text{obs}$ for all components | $\tau_1 = 44$ fs $\ll T_\text{obs}$ | ✗ misses fast dynamics |
| **Trajectory MC (ours)** | Gaussian noise | verified Phase 1 | ✓ exact |

The trajectory MC is the gold standard for classical Gaussian noise at any
$\sigma/J$ and any $\kappa$. Its cost is $N_\text{real}$ times one `sesolve`
call — trivial for a 2-site system (ms per call), feasible for 8 sites.

---

## Physical explanation

### Static disorder (τ₃ only): the localization trap

Each Monte Carlo realization draws a **fixed** energy offset
Δε = δs₄ − δs₇ from a Gaussian with σ_Δ = √(A₃) · √(σ₄²+σ₇²) ≈ 705 cm⁻¹.
Since σ_Δ / J ≈ 12, most realizations have |Δε| ≫ J. The eigenstates are
nearly site-localised, and the per-realisation population transfer is:

$$P_7(t) = \frac{J^2}{J^2 + (\Delta\varepsilon/2)^2}\sin^2\!\left(\frac{\sqrt{J^2+(\Delta\varepsilon/2)^2}\;t}{\hbar}\right)$$

- **Amplitude** ≈ (2J/Δε)² ≪ 1 for most realizations → very little transfer per shot.
- **Frequency** ≈ |Δε|/2ℏ, different for each realization → ensemble dephasing
  kills the oscillation.
- **Net**: P₇ oscillates weakly and stays below ~0.1. The excitation is
  trapped on Trp4 — this is the Anderson-localisation limit.

### Fast noise (τ₁ and/or τ₂): noise-assisted incoherent hopping

Now Δε(t) fluctuates. When it crosses zero, the sites come into transient
resonance and J drives efficient transfer. The noise randomises the quantum
phase, so there is no coherent Rabi oscillation — but energy still flows via
**incoherent hopping**. The ensemble-averaged dynamics become classical
rate-equation-like:

$$P_7(t) \approx \tfrac{1}{2}(1 - e^{-kt}), \qquad k \sim \frac{2J^2 \tau_c}{\hbar^2}$$

- **Monotonic, no oscillation** — the phase is scrambled faster than the
  Rabi period (τ₁ = 44 fs < τ_Rabi = 141 fs).
- **Plateau at 0.5** — the two sites have no mean energy bias (both ε = 0),
  so incoherent hopping equilibrates at 50/50.
- **Log-like saturation** — the exponential approach to equilibrium.

This is **environment-assisted quantum transport** (ENAQT): the noise is not
merely a source of decoherence to be minimised — it actively drives the system
through resonance configurations that static disorder would forbid.

### Why τ₁ transfers more than τ₂

Both cause incoherent hopping, but τ₁ = 44 fs crosses resonance more
frequently (σ₁τ₁ ≈ 5 in Kubo units — fast enough to sample many resonance
events within τ_Rabi). τ₂ = 1.7 ps is slower (σ₂τ₂ ≈ 165 — quasi-static on
the transfer timescale), so it acts partly as additional static disorder.
The transfer hierarchy τ₁ > τ₂ > τ₃ reflects the **rate of resonance
sampling**: faster fluctuations explore configuration space more thoroughly.

### Full noise: fast wins, slow hurts

Adding τ₃ to the dynamic noise (τ₁+τ₂ → full) introduces a per-realisation
static bias that shifts the time-averaged gap away from zero. This slightly
reduces the hopping rate, lowering the plateau from 0.51 to 0.48. The
fast components still dominate the transport mechanism.

---

## Coherence dynamics

| config | |ρ₄₇| at 20 fs | |ρ₄₇| at 100 fs |
|---|---|---|
| noiseless | 0.215 | 0.397 |
| τ₃ only | 0.168 | 0.210 |
| τ₁ only | 0.126 | **0.257** (rising!) |
| full | 0.102 | 0.197 |

The τ₁-only coherence **rises** from 20 to 100 fs. This is motional narrowing:
after the initial dephasing, the fast bath averages out and the J coupling
regenerates coherence. The τ₃-only case cannot average out (it's frozen),
so coherence decays monotonically.

With σ/J ≈ 17–25, all three Kubo numbers κₖ = σₖτₖ/ℏ are ≫ 1 (5.3, 165,
180 000). A Lindblad / Markovian treatment would give the wrong lineshape
for every component. The trajectory approach is essential.

---

## Key numbers

| quantity | value |
|---|---|
| Rabi transfer time (noiseless) | 141 fs |
| P₇ plateau, fast noise | ~0.50 |
| P₇ at 2 ps, static only | 0.092 |
| P₇ at 2 ps, full noise | 0.478 |
| Coherence peak (noiseless) | 0.50 at 71 fs |
| Coherence at 100 fs, τ₁ only | 0.257 (motional narrowing) |

---

## Files

```
results/phase5_exciton_dynamics/
    mc_results.npz              ensemble observables (P4, P7, coh per config)
    ablation_coherence.png      |rho_47(t)| for 5 configs + noiseless
    ablation_population.png     P4 (dashed) + P7 (solid) for 5 configs

results/phase5_noise_check/
    ou_sum_acf.png              OU sum ACF vs theory (verification)
```
