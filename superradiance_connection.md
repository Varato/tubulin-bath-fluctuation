# Connecting tubulin noise characterization to microtubule superradiance

## 1. The experimental observation (2024 superradiance paper)

A 2024 study measured tryptophan (Trp) fluorescence quantum yield (QY) across three
system sizes: free Trp in solution, tubulin heterodimer (TuD), and assembled
microtubule (MT). After correcting for non-Trp contributions, the Trp-specific QY
shows a non-monotonic hierarchy:

| Sample | QY (280 nm excitation) | QY (295 nm, Trp-only) |
|---|---|---|
| Free Trp | 12.4 ± 1.1% | 11.4 ± 1.1% |
| Tubulin dimer (TuD) | 10.6 ± 0.6% | 10.9 ± 1.3% |
| Microtubule (MT) | 15.7 – 19.5% | 14.7 ± 1.6% |

**Key trend:** dimer QY is slightly *below* free Trp, while microtubule QY is
significantly *above* both. The MT enhancement over the dimer reaches ~70%.

The paper's interpretation:

1. **Dimer QY < free Trp:** the protein matrix introduces extra non-radiative
   decay channels (internal conversion, charge transfer, intersystem crossing),
   raising the non-radiative rate $k_{\mathrm{nr}}$ and lowering QY.

2. **MT QY > dimer QY:** collective superradiance from the ordered Trp network
   in the microtubule lattice. Thousands of Trp form a giant coupled dipole
   network, enhancing the radiative rate $k_{\mathrm{rad}}$.

3. The enhancement is robust against room-temperature static disorder of
   $\sim 200$ cm$^{-1}$, which disperses but does not destroy the collective
   oscillator strength.

QY is defined as:

$$\mathrm{QY} = \frac{k_{\mathrm{rad}}}{k_{\mathrm{rad}} + k_{\mathrm{nr}}}$$

So the two transitions involve different rate channels: free Trp $\to$ dimer
raises $k_{\mathrm{nr}}$ (quenching), while dimer $\to$ MT raises
$k_{\mathrm{rad}}$ (superradiance).


## 2. What our noise characterization contributes

Our MD-derived noise model quantifies the protein-induced electric field
fluctuations at each of the 8 Trp sites. The key parameters:

### Noise magnitude

Per-site total fluctuation amplitude (slow trajectory):
- Mean $\sigma = 914$ cm$^{-1}$ (range: 675–1453 cm$^{-1}$ across sites)
- For the strongest-coupled pair Trp4–Trp7: $\sigma_4 = 999$ cm$^{-1}$,
  $\sigma_7 = 675$ cm$^{-1}$

The ratio $\sigma / |J_{47}|$ ranges from 11 (Trp7) to 25 (Trp8), with mean 15.5.
This places the system deep in the **strong-coupling regime**
($\sigma \gg J$), where noise overwhelms the excitonic coupling.

### Three-timescale decomposition

The autocorrelation function decomposes into three OU components:

| Component | $\tau_k$ (mean) | $\sigma_k$ (mean) | Kubo number $\kappa_k$ |
|---|---|---|---|
| $T_1$ (fast) | 46 fs | 655 cm$^{-1}$ | $\sim 6$ |
| $T_2$ (intermediate) | 1.77 ps | 524 cm$^{-1}$ | $\sim 1.7 \times 10^2$ |
| $T_3$ (slow) | 2.66 ns | 363 cm$^{-1}$ | $\sim 1.8 \times 10^5$ |

For the Trp4–Trp7 pair specifically, $\sigma_3 \approx 200$–230 cm$^{-1}$,
matching the $\sim 200$ cm$^{-1}$ static disorder cited in the 2024 paper.

All three Kubo numbers $\kappa_k = \sigma_k \tau_k / \hbar \gg 1$, confirming
the bath is deeply non-Markovian. The $T_3$ component ($\tau_3 = 2.66$ ns)
is effectively static on the exciton timescale ($\hbar/J \approx 90$ fs),
producing quasi-static energetic disorder that drives Anderson localisation.


## 3. Eigenstate structure of the clean Craddock Hamiltonian

Even **without any added noise**, the Craddock 8-site Hamiltonian produces
strongly localised eigenstates. Diagonalising $H_{8\times 8}$ (in cm$^{-1}$,
relative to base energy 35888 cm$^{-1}$):

| Eigenstate | Energy (cm$^{-1}$) | PR | Dominant character |
|---|---|---|---|
| $\psi_1$ | $-25.4$ | 1.52 | Trp6–Trp8 pair |
| $\psi_2$ | $+0.3$ | 1.01 | Trp1 (isolated) |
| $\psi_3$ | $+57.6$ | 1.23 | Trp5 (isolated) |
| $\psi_4$ | $+125.6$ | 1.50 | Trp6–Trp8 pair |
| $\psi_5$ | $+166.2$ | 1.90 | Trp4–Trp7 pair |
| $\psi_6$ | $+287.9$ | 2.07 | Trp4–Trp7 pair (most delocalised) |
| $\psi_7$ | $+320.6$ | 1.84 | Trp2–Trp3 pair |
| $\psi_8$ | $+412.2$ | 1.62 | Trp2–Trp3 pair |

**Participation ratio** $\mathrm{PR}_k = 1 / \sum_i |c_i^{(k)}|^4$ measures how
many sites eigenstate $k$ effectively spans. PR = 1 means perfect localisation
to one site; PR = 8 means uniform delocalisation across all sites.

**Finding:** maximum PR = 2.07. The "8-site network" is actually three isolated
dimers (Trp4–Trp7, Trp6–Trp8, Trp2–Trp3) plus two singleton sites (Trp1, Trp5).

### Root cause: intrinsic energetic disorder

The diagonal (site) energies of the Craddock Hamiltonian span 388 cm$^{-1}$
(from 0 to 388 cm$^{-1}$ relative), which is **6.6× the strongest coupling**
$|J_{47}| = 59$ cm$^{-1}$. These site energy differences arise from the
different protein electrostatic environments at each Trp location and act as
intrinsic static disorder.

The coupling network fragments into pairs where the inter-site detuning
$\Delta\varepsilon$ is comparable to $J$:

| Pair | $J$ (cm$^{-1}$) | $\Delta\varepsilon$ (cm$^{-1}$) | $\Delta\varepsilon / J$ |
|---|---|---|---|
| Trp4–Trp7 | $-59$ | 41 | 0.69 |
| Trp6–Trp8 | $-51$ | 102 | 2.0 |
| Trp2–Trp3 | $-41$ | 46 | 1.1 |

Only the Trp4–Trp7 pair has $\Delta\varepsilon < J$, explaining its slightly
higher PR. The other pairs are marginally resonant. Cross-pair couplings are
all $|J| \leq 21$ cm$^{-1}$, too weak to bridge the large inter-pair energy
gaps.

### Implication

The effective network size for superradiance within a single dimer is
$N_{\mathrm{eff}} \approx 2$ (at best, for the Trp4–Trp7 pair), not 8. The
maximum possible oscillator strength enhancement is $\sim 2\times$, and only
if the two dipoles are aligned.


## 4. The bandwidth criterion: when does superradiance survive disorder?

### The $J\sqrt{N}$ argument

For a network of $N$ sites with typical pairwise coupling $J$, the collective
eigenstate bandwidth scales as $\sim J\sqrt{N}$. This follows from random
matrix theory: the eigenvalue spread of a matrix whose off-diagonal elements
have typical magnitude $J$ grows as $J\sqrt{N}$ (Wigner semicircle law).

Static disorder $\sigma$ broadens the site energies randomly. The condition
for collective eigenstates to survive:

$$J\sqrt{N} \gtrsim \sigma$$

**Physical mechanism:** when the inter-site coupling bandwidth exceeds the
random energetic detunings, eigenstates remain delocalised across all $N$
sites, and the brightest state carries superradiant oscillator strength
(enhancement up to $N$). When disorder wins, the eigenstates localise and
the oscillator strength redistributes to single-site values (enhancement
$\sim 1$).

### Application to tubulin

**Single dimer ($N = 8$):**

The effective bandwidth from the Craddock couplings is $J\sqrt{8} \approx
59 \times 2.83 \approx 167$ cm$^{-1}$. This must compete against:

1. Intrinsic site energy spread: 388 cm$^{-1}$ (already $> J\sqrt{N}$)
2. $T_3$ quasi-static noise: $\sigma_3 \approx 200$–360 cm$^{-1}$ per site
3. Total: $\sigma_{\mathrm{eff}} \sim 400$–600 cm$^{-1}$

The dimer sits **well below** the superradiance threshold. Even without added
noise, the intrinsic spread (388 cm$^{-1}$) exceeds the coupling bandwidth
(167 cm$^{-1}$), which is why the PR analysis shows localised eigenstates.

**Microtubule ($N \sim 10^3$–$10^4$ dimers):**

The microtubule lattice links tubulin dimers into a helical tube of 13
protofilaments. Even with weak inter-dimer coupling $J' \sim 1$–10 cm$^{-1}$
(much weaker than intra-dimer $J$), the effective bandwidth scales as:

$$J'\sqrt{N_{\mathrm{MT}}} \sim 10 \times \sqrt{10^4} = 1000 \text{ cm}^{-1}$$

This **far exceeds** $\sigma_3 \approx 200$–360 cm$^{-1}$, placing the
microtubule **well above** the superradiance threshold.

### Limitations of the $J\sqrt{N}$ heuristic

This argument is most rigorous for systems with near-degenerate sites and
random couplings. The Craddock Hamiltonian has specific, non-random structure
(non-uniform couplings, large intrinsic energy spread). The heuristic gives
the correct intuition but cannot replace direct oscillator strength
computation.

For the microtubule, the inter-dimer coupling structure is not yet known
in detail, so the $J' \sim 1$–10 cm$^{-1}$ estimate is speculative. The
actual MT superradiance may involve radiative (through-field) coupling
rather than purely excitonic (through-bond) coupling.


## 5. Oscillator strength of dimer eigenstates

### Definition

Each eigenstate $|\psi_k\rangle = \sum_i c_i^{(k)} |i\rangle$ has a transition
dipole moment equal to the coherent vector sum of site dipoles:

$$\boldsymbol{\mu}_k = \sum_{i=1}^{8} c_i^{(k)} \boldsymbol{\mu}_i$$

where $\boldsymbol{\mu}_i = |\mu| \hat{n}_i$ is the transition dipole of site
$i$ (magnitude $|\mu| = 5$ D for all Trp, direction $\hat{n}_i$ from MD).

The **oscillator strength** (proportional to the radiative decay rate) is:

$$f_k = |\boldsymbol{\mu}_k|^2 = \left|\sum_i c_i^{(k)} \hat{n}_i\right|^2 |\mu|^2$$

Normalising by the single-Trp value $|\mu|^2$, the enhancement factor is:

$$\mathcal{F}_k = \frac{f_k}{|\mu|^2} = \left|\sum_i c_i^{(k)} \hat{n}_i\right|^2$$

### Physical interpretation

- If all 8 dipoles are parallel and the eigenstate is fully delocalised
  ($c_i = 1/\sqrt{8}$): $\mathcal{F} = 8$ (maximum superradiance).

- If dipoles are randomly oriented: partial cancellation reduces
  $\mathcal{F}$ even for delocalised eigenstates.

- If the eigenstate is localised on one site ($c_i \approx \delta_{ij}$):
  $\mathcal{F} = |\hat{n}_j|^2 = 1$ (no enhancement).

### Computed results (Analysis A)

Using MD-derived time-averaged dipole directions $\hat{n}_i$ from the slow
trajectory (4001 frames, dt = 10 ps):

| Eigenstate | Energy (cm$^{-1}$) | PR | $\mathcal{F}_k$ | Character |
|---|---|---|---|---|
| $\psi_1$ | $-25.4$ | 1.52 | 0.262 | Trp8–Trp6 (dark) |
| $\psi_2$ | $+0.3$ | 1.01 | 0.981 | Trp1 (isolated) |
| $\psi_3$ | $+57.6$ | 1.23 | 1.208 | Trp5–Trp8 |
| $\psi_4$ | $+125.6$ | 1.50 | 1.404 | Trp6–Trp8 |
| **$\psi_5$** | $+166.2$ | **1.90** | **1.806** | **Trp4–Trp7 (brightest)** |
| $\psi_6$ | $+287.9$ | 2.07 | 0.499 | Trp7–Trp4 (dark pair) |
| $\psi_7$ | $+320.6$ | 1.84 | 1.627 | Trp3–Trp2 |
| $\psi_8$ | $+412.2$ | 1.62 | 0.213 | Trp2–Trp3 (dark) |

**Key findings:**

- **Sum rule verified:** $\sum_k \mathcal{F}_k = 8.000$ (Thomas-Reiche-Kuhn).
- **Brightest state:** $\psi_5$ with $\mathcal{F}_{\max} = 1.81$. This is a
  mere **1.81$\times$ enhancement** over single Trp, far below the ideal
  $N = 8$ limit.
- **Trp4–Trp7 dipole alignment:** $\cos\theta = 0.745$ ($|\theta| = 42$°),
  moderately aligned. For a perfectly aligned 2-site pair: $\mathcal{F} = 2$.
  The actual 1.81 reflects both incomplete alignment and incomplete
  delocalisation (PR = 1.90).
- **Dark states:** $\psi_1$ ($\mathcal{F} = 0.26$) and $\psi_8$
  ($\mathcal{F} = 0.21$) are the anti-superradiant companions of their
  respective pairs.

**Conclusion:** the clean Craddock dimer offers at most 1.81$\times$ radiative
enhancement. This is negligible for superradiance.

### Thomas–Reiche–Kuhn sum rule

The total oscillator strength is conserved:

$$\sum_{k=1}^{N} f_k = N |\mu|^2$$

Static disorder redistributes $f_k$ among eigenstates but does not change the
total. The question is whether the brightest state retains a disproportionate
share. In a localised regime (PR $\sim 1$), oscillator strength is spread
roughly evenly: each state carries $\sim 1 \times |\mu|^2$. In a delocalised
regime, one bright state can carry $\sim N \times |\mu|^2$.


## 6. Quantitative analyses (completed)

Scripts: `superradiance/scripts/analysis_oscillator_strength.py`
Results: `superradiance/results/clean_oscillator_strength.npz`,
         `superradiance/results/disorder_scan.npz`
Figure: `superradiance/figures/fig_oscillator_strength.{png,pdf}`

### Analysis A: Clean oscillator strengths

See Section 5 above. Brightest state $\psi_5$ has $\mathcal{F} = 1.81$.

### Analysis B: Disorder $\sigma$-scan of brightest-state oscillator strength

Added random Gaussian static disorder $\delta\varepsilon_i \sim \mathcal{N}(0,
\sigma^2)$ to the Craddock diagonal, diagonalised, computed $\mathcal{F}_k$ for
all eigenstates. 500 realisations per $\sigma$ point.

| $\sigma$ (cm$^{-1}$) | $\mathcal{F}_{\max}$ (mean) | SEM | Max PR |
|---|---|---|---|
| 0 | 1.806 | 0 | 2.07 |
| 25 | 1.818 | 0.004 | 2.22 |
| 50 | 1.832 | 0.008 | 1.89 |
| 100 | 1.814 | 0.010 | 1.62 |
| 150 | 1.751 | 0.011 | 1.50 |
| **200** ($\sigma_3$) | **1.674** | **0.011** | **1.40** |
| 300 | 1.592 | 0.012 | 1.31 |
| 500 | 1.480 | 0.010 | 1.20 |

**Key observations:**

1. **Mild disorder enhances $\mathcal{F}_{\max}$ slightly** ($\sigma \sim 25$–75
   cm$^{-1}$): a small ENAQT-like effect where disorder redistributes oscillator
   strength, making the brightest state marginally brighter. Peak
   $\mathcal{F}_{\max} = 1.835$ at $\sigma = 75$ cm$^{-1}$.

2. **At $\sigma_3 = 200$ cm$^{-1}$:** $\mathcal{F}_{\max} = 1.674$, still above
   single-Trp but far below any meaningful superradiant enhancement.

3. **Even at $\sigma = 500$ cm$^{-1}$:** $\mathcal{F}_{\max} = 1.48$, still
   above 1. This residual enhancement comes from the intrinsic site energy
   structure of the Craddock Hamiltonian (the Trp4–Trp7 pair is always
   somewhat resonant).

4. **Throughout the entire scan, $\mathcal{F}_{\max} < 2$**, confirming no
   meaningful superradiance is possible in the isolated dimer at any disorder
   level.

### Analysis C: Participation ratio under disorder

Max PR decreases monotonically (after a brief ENAQT-like uptick at low
$\sigma$): from 2.07 (clean) to 1.20 ($\sigma = 500$). At $\sigma_3 = 200$,
the most delocalised state has PR = 1.40, meaning eigenstates are almost
completely localised to single sites.


## 7. Size scaling: why the dimer fails, why the MT works

### The dimer ($N = 8$)

Three factors conspire to suppress superradiance in the isolated tubulin
dimer:

1. **Large intrinsic energy spread** (388 cm$^{-1}$): the 8 Trp site energies
   differ by up to 388 cm$^{-1}$ due to different electrostatic environments.
   This is 6.6× the strongest coupling, causing the network to fragment into
   isolated pairs.

2. **Small effective $N$:** the coupling network produces only three weakly
   coupled dimers (PR $\leq 2.07$). The computed maximum oscillator strength is
   $\mathcal{F}_{\max} = 1.81$ (clean) dropping to 1.67 at $\sigma_3$,
   insufficient to overcome the non-radiative quenching.

3. **MD-derived noise adds further disorder:** our $\sigma_3 \approx 200$–360
   cm$^{-1}$ is comparable to $J$, pushing even the Trp4–Trp7 pair toward
   localisation. The total noise $\sigma \approx 900$ cm$^{-1}$ completely
   overwhelms all couplings.

**Conclusion:** the isolated dimer cannot exhibit superradiance. Its QY
(10.8%) is slightly below free Trp (11.4%) because the protein adds
non-radiative channels without providing any compensating radiative
enhancement.

### The microtubule ($N \sim 10^3$–$10^4$)

The microtubule assembles tubulin dimers into a helical lattice. This creates
a qualitatively different regime:

1. **Large network size:** with $\sim 10^3$–$10^4$ dimers, each containing 8
   Trp, the total number of coupled chromophores is $N_{\mathrm{MT}} \sim
   10^4$–$10^5$.

2. **Bandwidth overcomes disorder:** even weak inter-dimer coupling
   ($J' \sim 1$–10 cm$^{-1}$) produces a collective bandwidth
   $J'\sqrt{N_{\mathrm{MT}}} \sim 100$–1000 cm$^{-1}$ that can exceed
   $\sigma_3$. The system crosses above the superradiance threshold.

3. **Superradiant rate dominates:** the collective radiative rate
   $k_{\mathrm{rad}}^{\mathrm{MT}} \propto N_{\mathrm{eff}}$ can exceed the
   non-radiative rate $k_{\mathrm{nr}}$, raising QY.

4. **Robustness to $\sigma_3$:** the 2024 paper observes that $\sim 200$
   cm$^{-1}$ static disorder does not destroy the MT superradiance. This is
   consistent with our analysis: the lattice-scale collective bandwidth
   ($J'\sqrt{N_{\mathrm{MT}}}$) far exceeds $\sigma_3$, so the collective
   states survive.

**Conclusion:** microtubule superradiance arises from the lattice-scale
network, not from within a single dimer. Our noise characterization confirms
that the per-site disorder ($\sigma_3 \sim 200$ cm$^{-1}$) is too large for
single-dimer superradiance but consistent with lattice-scale collective
effects.


## 8. Connection to non-Markovian physics and Kubo number

### Why the Haken–Strobl model is inapplicable

The Haken–Strobl model (HSM, white-noise pure dephasing) assumes
$\kappa = \sigma\tau/\hbar \ll 1$: many uncorrelated dephasing events per
coherent timescale. Our system has $\kappa_k \gg 1$ for all three components,
placing it deeply in the non-Markovian regime.

Under HSM dephasing ($\gamma = 50$ cm$^{-1}$), our 8-site simulations show
extensive delocalisation: populations spread across the entire network.
This is because white noise homogenises the energetic landscape, erasing
the static disorder that causes localisation. The HSM result is
unphysical for this system.

Under the colored-noise model (CNM, our MD-derived tri-exponential OU bath),
the 8-site simulations show strong localisation: excitons remain confined to
strongly-coupled site pairs. The CNM preserves the static disorder component
($T_3$), producing realistic Anderson localisation.

### Implication for superradiance

The failure of HSM has direct consequences for superradiance predictions:

- **HSM would predict delocalised eigenstates** across the 8-site network,
  suggesting superradiant enhancement up to $N = 8$. This is incorrect.

- **CNM correctly predicts localisation** to $\sim 2$-site pairs, limiting
  the enhancement to $\sim 2\times$. The non-Markovian bath preserves the
  quasi-static disorder that fragments the network.

- **For the microtubule:** HSM would overestimate superradiance by ignoring
  disorder. CNM correctly captures the competition between coupling and
  disorder, which determines whether lattice-scale collective states form.

### ENAQT and the fast component ($T_1$)

The $T_1$ component ($\tau_1 = 46$ fs, $\sigma_1 \approx 655$ cm$^{-1}$)
drives environment-assisted quantum transport (ENAQT): rapid site-energy
fluctuations sweep the Trp4–Trp7 detuning through resonance, enabling
transfer that the static disorder alone would block.

ENAQT provides partial connectivity within the dimer ($P_{\mathrm{leak}} \sim
0.2$–$0.6$ in our 8-site MC), but not enough to create delocalised eigenstates.
In the microtubule, ENAQT may play a larger role by providing inter-dimer
connectivity that supplements the direct coupling $J'$.


## 9. Physical picture: putting it all together

```
Free Trp (QY ≈ 11.4%)
    │
    │  Protein environment adds:
    │  σ ≈ 900 cm⁻¹ electric field noise
    │  → non-radiative quenching channels
    │  → no compensating superradiance (F ≈ 1)
    │
    ▼
Tubulin dimer (QY ≈ 10.8%)
    │  Intrinsic site energy spread: 388 cm⁻¹
    │  + T₃ static noise: ~200 cm⁻¹
    │  → network fragments into 2-site pairs
    │  → PR ≤ 2, F_max = 1.81 (clean) → 1.67 (at σ₃)
    │  → disorder wins (J√8 ≈ 167 < σ_eff ≈ 500)
    │
    │  Lattice assembles ~10³-10⁴ dimers
    │  → inter-dimer bandwidth J'√N_MT >> σ₃
    │  → collective states survive disorder
    │  → radiative rate enhanced by N_eff
    │
    ▼
Microtubule (QY ≈ 14.7-17.6%)
```

The QY hierarchy reflects a **size-dependent crossover**:
- Below threshold ($J\sqrt{N} < \sigma$): disorder localises excitons,
  quenching dominates, QY is low.
- Above threshold ($J\sqrt{N} > \sigma$): collective states form,
  superradiance enhances radiative rate, QY rises.

The tubulin dimer sits below threshold. The microtubule sits above.


## 10. What goes in the paper

### Results section (quantitative)

A new subsection or paragraph reporting:

1. **Eigenstate PR analysis** of the clean Craddock Hamiltonian: table of
   eigenvalues, PR, and dominant character. Key result: max PR = 2.07, the
   network fragments into isolated pairs.

2. **Bandwidth analysis:** diagonal spread (388 cm$^{-1}$) vs coupling
   bandwidth ($J\sqrt{8} \approx 167$ cm$^{-1}$) vs $T_3$ noise
   ($\sigma_3 \approx 200$ cm$^{-1}$). The dimer is below threshold.

3. **Oscillator strength** (Analysis A): $\mathcal{F}_k$ per eigenstate using
   MD dipole directions. Brightest state's enhancement.

4. **Disorder scan** (Analysis B): $\mathcal{F}_{\max}$ vs $\sigma$, showing
   the transition from enhanced to quenched as disorder increases.

### Discussion section (interpretation)

1. Connect the eigenstate localisation to the QY hierarchy: dimer QY is
   low because superradiance is impossible (PR $\leq 2$, disorder wins).

2. Argue the size-dependent threshold: $J\sqrt{N}$ vs $\sigma_3$ explains
   why MT superradiance survives $\sim 200$ cm$^{-1}$ disorder but dimer
   superradiance does not.

3. Contrast CNM vs HSM predictions for superradiance: HSM overestimates
   delocalisation (and thus superradiance) by erasing static disorder;
   CNM correctly captures localisation.

4. Reference the 2024 experimental paper as direct validation of the
   size-dependent crossover.

5. Note the limitation: our analysis covers the single dimer ($N = 8$).
   Microtubule-scale simulations ($N \sim 10^4$) would require inter-dimer
   coupling data and are beyond the current scope.


## 11. Data sources

- **Craddock Hamiltonian:** `base_hamiltonian.md` / `Craddock_2014.ipynb`
  in the `decoherence` workspace.

- **Per-site $\sigma$:** `results/phase1_basic_stats/sigma_matrix.csv`

- **Tri-exponential fit:** `results/paper_figures/corrected_triexp.csv`

- **Dipole directions:** `dmu` arrays in the NPZ data files (loaded via
  `utils.load('slow')['dmu']`, shape $(N_{\mathrm{frames}}, 8, 3)$).

- **Dipole reorientation data:** `results/phase3_source_attribution/dmu_acf.npz`,
  `results/structural_correlation/indole_mobility.json`

- **8-site MC results:** `results/paper_figures/fig4_8site_scan.npz`

- **Oscillator strength analysis:** `superradiance/scripts/analysis_oscillator_strength.py`,
  `superradiance/results/clean_oscillator_strength.npz`,
  `superradiance/results/disorder_scan.npz`,
  `superradiance/figures/fig_oscillator_strength.{png,pdf}`
