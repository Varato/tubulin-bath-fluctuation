# Detailed Research Plan: Multi\-Scale Site Energy Fluctuation Analysis \(Coding\-Agent Ready\)

## 1\. General Overview

### 1\.1 Research Objective

Quantify multi\-timescale environmental noise on the 8 tryptophan \(Trp\) sites in the tubulin network based on dual\-resolution MD simulation data\. Validate the fast/slow two\-component fluctuation model, correct the bias of the traditional pure static disorder approximation, and systematically reassess exciton localization, dephasing and long\-range coherent energy transfer feasibility in the Trp network\. Provide rigorous numerical evidence for the core conclusion that exciton long\-range coherent transport is suppressed in tubulin, and summarize general methodological criteria for static disorder approximation in protein exciton research\.

### 1\.2 Data Source \& Basic Attributes

Unified NPZ noise analysis workspace containing two complementary MD trajectory datasets with unified data schema, covering sub\-ps fast dynamics to ns slow dynamics\.

|Dataset|Sampling Interval|Total Frames|Effective Frequency Range|Core Function|
|---|---|---|---|---|
|Slow trajectory|10 ps|4001|0\.02 – 1\.67 cm⁻¹|Capture ns\-scale slow conformational fluctuations|
|Fast trajectory|10 fs \(0\.01 ps\)|200001|0\.02 – 1668 cm⁻¹|Capture sub\-ps fast environmental fluctuations|

### 1\.3 Unified NPZ Data Schema \& Unit Conversion

#### 1\.3\.1 Core Shared Keys \(Both Datasets\)

|Key|Shape|Units|Description|
|---|---|---|---|
|time\_ps|\(N,\)|ps|Trajectory frame time axis|
|dt\_ps|Scalar|ps|Simulation sampling interval|
|dmu|\(N, 8, 3\)|Unit vector|Trp transition dipole direction vector|
|E\_total / E\_protein / E\_water / E\_nucleotide / E\_ions|\(N, 8, 3\)|V/m|Total and component\-resolved electric field at indole center|
|delta\_s\_total / delta\_s\_protein / delta\_s\_water / delta\_s\_nucleotide / delta\_s\_ions|\(N, 8\)|V/m|Site energy fluctuation precursor: $-\Delta\boldsymbol{\mu}\cdot\mathbf{E}$|

#### 1\.3\.2 Exclusive Keys \(50 ns Slow Dataset Only\)

`com`, `indole_center` \(scripts must tolerate missing values in fast dataset\)

#### 1\.3\.3 Unit Conversion Constant

```Plain Text
V_TO_CM = 8.397e-7  # Convert V/m to cm⁻¹ (|Δμ| = 5 D for Trp S₀→S₁ transition)
```

### 1\.4 Trp Site Index Mapping \(0–7\)

|Index|Residue Label|Protein Chain|
|---|---|---|
|0|αW21|A|
|1|αW346|A|
|2|αW388|A|
|3|αW407|A|
|4|βW21|B|
|5|βW101|B|
|6|βW344|B|
|7|βW397|B|

### 1\.5 Fixed Physical Parameters

- Trp exciton coupling strength: $J = 45\ \text{cm}^{-1}$

- Exciton coherent time window: $T_\text{obs} = 2\ \text{ps}$

- Exciton characteristic transfer time: $\tau_J \approx 74\ \text{fs}$

### 1\.6 Coding Implementation Specifications

- All time units unified to **ps**; all energy/frequency units unified to **cm⁻¹**

- Modular code structure: shared reusable functions for ACF, PSD, data stitching, curve fitting; independent script for each analysis phase

- Output standardization: numerical results saved as NPZ/CSV; visualization results saved as PNG/PDF

- Full reproducibility: fixed random seeds for all Monte Carlo sampling processes

- Pre\-analysis validation: reproduce linearity and spectral normalization checks in README first

## 2\. Phase 1: Basic Statistics \& Variance Decomposition

### 2\.1 Research Goal

Quantify the amplitude of total and component\-resolved site energy fluctuations, verify the strong disorder regime of the system, decompose variance contribution of each environmental source, and test Gaussian statistics of fluctuations\.

### 2\.2 Task 1\.1: Fluctuation Amplitude \& Disorder Strength Calibration

#### Input

`delta_s_total, delta_s_protein, delta_s_water, delta_s_nucleotide, delta_s_ions` \(both datasets, converted to cm⁻¹ via `V_TO_CM`\)

#### Methodology

1. Demean all time series to eliminate baseline offset\.

2. Calculate standard deviation for each Trp site \(0–7\) and each component: $\sigma_X = \sqrt{\langle \Delta E_X^2 \rangle - \langle \Delta E_X \rangle^2}$

3. Compute disorder\-coupling ratio: $\sigma_\text{total}/J$

4. Count mean and standard deviation of σ values across 8 sites for statistical summary\.

#### Output

- 8×5 σ matrix \(per site \& per component\)

- System\-level average disorder strength $\overline{\sigma_\text{total}/J}$

- Validation of $\sigma_\text{total}/J \gg 1$ \(strong disorder confirmation\)

### 2\.3 Task 1\.2: Full Variance Decomposition with Cross\-Covariance

#### Input

Four component\-resolved site energy time series \(protein, water, nucleotide, ions\) in cm⁻¹

#### Methodology

Based on linear superposition principle: $\Delta E_\text{total} = \sum \Delta E_i$

Total variance formula: $\sigma_\text{total}^2 = \sum_i \sigma_i^2 + 2\sum_{i<j} \text{Cov}(i,j)$

1. Calculate 4×4 covariance matrix of four environmental components\.

2. Quantify individual variance contribution and pairwise cross\-covariance contribution to total variance\.

3. Focus on protein\-water covariance sign and magnitude \(dielectric screening effect\)\.

4. Calculate screening ratio: $1 - \sigma_\text{total}^2 / \sum\sigma_i^2$

#### Output

- Component variance fraction table

- Cross\-covariance contribution statistics

- Average dielectric screening ratio of the system

### 2\.4 Task 1\.3: Gaussianity Test

#### Input

`delta_s_total` time series of 8 Trp sites \(fast dataset, full resolution\)

#### Methodology

1. Plot normalized fluctuation histogram with Gaussian fitting curve\.

2. Calculate skewness and excess kurtosis for each site\.

3. Conduct Kolmogorov\-Smirnov normality test\.

#### Output

- Skewness \& kurtosis table for all sites

- Normality test result, judge whether Gaussian statistics hold \(support ACF/PSD as complete fluctuation descriptors\)

## 3\. Phase 2: Full\-Band Multi\-Timescale Dynamics Characterization

### 3\.1 Research Goal

Quantify timescale separation of site energy fluctuations via ACF and PSD, verify the validity of fast/slow dual\-component model, and quantitatively evaluate the systematic bias caused by low 10 ps sampling resolution\.

### 3\.2 Task 2\.1: ACF Calculation \& Dual\-Trajectory Stitching

#### Input

`delta_s_total`, `time_ps` \(fast \& slow datasets\)

#### Methodology

Normalized autocorrelation function definition: $C(t) = \frac{\langle \Delta E(0)\Delta E(t) \rangle}{\langle \Delta E^2 \rangle}$

1. Compute FFT\-based ACF for fast dataset \(0\.01 ps \~ 500 ps\) and slow dataset \(10 ps \~ 25 ns\) respectively\.

2. Stitch unified ACF: adopt fast ACF for high\-frequency short\-time region, slow ACF for low\-frequency long\-time region\.

3. Verify curve consistency in overlapping frequency band\.

#### Output

- Independent ACF arrays for fast/slow trajectories

- Unified full\-time\-scale stitched ACF \(10 fs – 25 ns\)

### 3\.3 Task 2\.2: Bi\-Exponential ACF Fitting

#### Input

Stitched full\-scale ACF per Trp site

#### Methodology

Fitting model: $C(t) = A_\text{fast} e^{-t/\tau_\text{fast}} + A_\text{slow} e^{-t/\tau_\text{slow}},\ A_\text{fast}+A_\text{slow}=1$

1. Non\-linear least squares fitting for each site\.

2. Compare tri\-exponential fitting via AIC/BIC to confirm two\-component optimality\.

3. Extract key parameters: $\tau_\text{fast}, \tau_\text{slow}, A_\text{fast}, A_\text{slow}$

4. Calculate integrated total correlation time: $\tau_\text{int} = A_\text{fast}\tau_\text{fast} + A_\text{slow}\tau_\text{slow}$

#### Output

- Per\-site dual\-component timescale \& amplitude parameter table

- System\-averaged fluctuation timescale parameters

### 3\.4 Task 2\.3: Full\-Band PSD Calculation \& Stitching

#### Input

Dual\-trajectory site energy time series, sampling interval `dt_ps`

#### Methodology

1. Compute PSD via Welch’s method for fast and slow datasets separately\.

2. Convert frequency unit: $1\ \text{ps}^{-1} = 33.356\ \text{cm}^{-1}$

3. Stitch full\-band PSD: \<1 cm⁻¹ from slow trajectory, \>1 cm⁻¹ from fast trajectory\.

4. Verify spectral normalization: $\int_0^\infty S(f)df \approx \sigma^2$

#### Output

- Full log\-scale PSD curve \(0\.02 – 1668 cm⁻¹\)

- Spectral normalization error per site

### 3\.5 Task 2\.4: Low\-Sampling Resolution Bias Quantification

#### Methodology

1. Downsample fast 10 fs dataset to 10 ps sampling interval\.

2. Calculate ACF and $\tau_\text{int}$ of downsampled data\.

3. Compare with true high\-resolution $\tau_\text{int}$, quantify overestimation bias factor\.

#### Output

Per\-site and average timescale overestimation bias of 10 ps low sampling

## 4\. Phase 3: Source\-Resolved Dynamical Attribution

### 4\.1 Research Goal

Clarify the physical origin of fast/slow fluctuations, distinguish the dynamic characteristics of different environmental components, and quantify the contribution of dipole orientation fluctuation to total site energy noise\.

### 4\.2 Task 3\.1: Component\-Wise ACF \& Timescale Fitting

#### Input

Four component\-resolved `delta_s_*` time series \(fast dataset\)

#### Methodology

1. Calculate normalized ACF for protein, water, nucleotide, ions components respectively\.

2. Fit single\-exponential decay to extract characteristic correlation time of each component\.

#### Expected Physical Result

Water dominates sub\-ps fast fluctuations; protein dominates ns\-scale slow fluctuations; nucleotide/ions contribute minor intermediate timescale noise\.

#### Output

Component\-resolved characteristic fluctuation timescale table

### 4\.3 Task 3\.2: Component\-Wise PSD Spectral Decomposition

#### Methodology

1. Compute PSD for each environmental component\.

2. Superimpose all component PSDs for comparative visualization\.

3. Quantify spectral power fraction of each component in high/low frequency bands\.

#### Output

Frequency band\-resolved component power contribution matrix

### 4\.4 Task 3\.3: Dipole Orientation Fluctuation Contribution Analysis

#### Input

`dmu`, `E_total`

#### Methodology

1. Calculate dipole orientation autocorrelation: $C_\mu(t) = \langle \hat{\boldsymbol{\mu}}(0)\cdot\hat{\boldsymbol{\mu}}(t) \rangle$, extract indole ring reorientation time\.

2. Control variable test: fix dipole direction to time\-averaged value, recalculate site energy fluctuation variance\.

3. Calculate ratio $\sigma_\text{fixed-\mu}/\sigma_\text{total}$ to quantify field vs\. dipole fluctuation contribution\.

#### Output

- Trp indole reorientation timescale

- Variance contribution fraction of dipole orientation fluctuation

## 5\. Phase 4: Spatial Correlation Analysis of Site Fluctuations

### 5\.1 Research Goal

Verify the common assumption of independent site disorder, quantify spatial correlation of Trp site energy fluctuations, and distinguish correlation characteristics of different environmental components\.

### 5\.2 Task 4\.1: 8×8 Site Cross\-Correlation Matrix

#### Input

`delta_s_total, delta_s_protein, delta_s_water`

#### Methodology

1. Calculate Pearson correlation matrix of 8 Trp sites for total, protein\-only, water\-only fluctuations\.

2. Classify and count correlation strength: intra\-A chain, intra\-B chain, inter\-chain, adjacent/distant sites\.

#### Output

- Three groups of 8×8 correlation matrices

- Statistical summary of spatial correlation distribution

### 5\.3 Task 4\.2: Spatial Cross\-Correlation Time Evolution

#### Methodology

1. Select representative site pairs \(nearest neighbor, farthest inter\-chain pair\)\.

2. Compute temporal cross\-correlation function and decay timescale\.

3. Compare with self\-ACF decay characteristics\.

#### Output

Cross\-correlation decay curves and timescale parameters of typical site pairs

## 6\. Phase 5: Exciton Dynamics Model Calibration \& Core Conclusion Verification

### 6\.1 Research Goal

Map multi\-scale fluctuation statistics to exciton model parameters, separate static slow disorder and dynamic fast dephasing effects, quantitatively compare traditional pure static model and revised dual\-mechanism model, and verify the core conclusion of exciton strong localization\.

### 6\.2 Task 5\.1: Fast/Slow Variance Splitting

#### Input

Bi\-exponential fitting parameters, total fluctuation variance

#### Methodology

$\sigma_\text{slow}^2 = A_\text{slow}\cdot\sigma_\text{total}^2,\ \sigma_\text{fast}^2 = A_\text{fast}\cdot\sigma_\text{total}^2$

Auxiliary verification: low\-pass filter fast trajectory to extract slow envelope, calculate direct $\sigma_\text{slow}$ for cross validation\.

#### Output

$\sigma_\text{slow}, \sigma_\text{fast}$ per site and $\sigma_\text{slow}/J$ ratio

### 6\.3 Task 5\.2: Static Localization Induced by Slow Disorder Only

#### Methodology

1. Monte Carlo sampling \(5000–10000 realizations\): construct static disorder Hamiltonian via $\mathcal{N}(0,\sigma_\text{slow}^2)$, support spatial correlation correction\.

2. Build 8\-site exciton Hamiltonian with fixed nearest\-neighbor coupling $J=45\ \text{cm}^{-1}$\.

3. Diagonalize Hamiltonian, calculate eigenstate participation ratio: $\text{PR} = 1/\sum_i|c_i|^4$

4. Compare PR results of slow\-disorder\-only model and original total\-disorder pure static model\.

#### Output

- Ensemble\-averaged participation ratio

- Localization degree comparison of two static models

- Verification: whether slow disorder alone is sufficient for strong exciton localization

### 6\.4 Task 5\.3: Dynamic Dephasing Rate Calculation

#### Methodology

Markovian pure dephasing rate: $\gamma_\phi \approx \sigma_\text{fast}^2 \cdot \tau_\text{fast}$

Coherence lifetime: $\tau_\phi = 1/\gamma_\phi$

Compare with exciton transfer timescale $\tau_J$ to evaluate dephasing suppression efficiency\.

#### Output

Per\-site $\gamma_\phi, \tau_\phi$ and timescale ratio $\tau_\phi/\tau_J$

### 6\.5 Task 5\.4: Static \+ Dynamic Dual\-Mechanism Model Comparison

#### Methodology

1. Simulate exciton population and coherence dynamics for two models:


    - Model A: Traditional pure static disorder model \(total σ\)

    - Model B: Revised model \(slow static disorder \+ fast dynamic dephasing\)

2. Quantify differences in coherence decay rate, long\-time localization degree and energy transfer efficiency\.

3. Evaluate the additional coherence suppression effect of fast dynamic fluctuations\.

#### Output

- Exciton population \& coherence evolution curves of two models

- Quantitative comparison table of model deviations

## 7\. Phase 6: Methodological Generalization \& Boundary Discussion

### 7\.1 Research Goal

Summarize universal methodological criteria for static disorder approximation, clarify MD sampling resolution selection principles for protein exciton quantum dynamics research\.

### 7\.2 Task 6\.1: Validity Criterion of Static Disorder Approximation

#### Key Dimensionless Parameters

- Slow variance fraction: $R_\text{slow} = \sigma_\text{slow}^2/\sigma_\text{total}^2$

- Timescale ratio: $r_\tau = \tau_\text{slow}/T_\text{obs}$

#### Methodology

Establish quantitative criterion to judge the applicability of static disorder approximation, and apply it to the tubulin Trp network system for verification\.

### 7\.3 Task 6\.2: MD Sampling Frequency Guideline

Summarize the systematic errors of 10 ps low sampling \(missing fast sub\-ps fluctuations, overestimating correlation time\), propose universal sampling principle:**MD sampling interval must be far smaller than the fastest environmental fluctuation correlation time for ps\-scale exciton coherence research**\.

## 8\. Final Core Expected Conclusions

1. Tubulin Trp network site energy fluctuations have obvious fast/slow timescale separation; water dominates sub\-ps fast dephasing, protein dominates ns static disorder\.

2. Even excluding fast dynamic fluctuations, the residual slow static disorder is still far larger than exciton coupling strength, leading to intrinsic strong exciton localization\.

3. Fast environmental fluctuations further accelerate coherence decay and suppress long\-range energy transfer, without reversing the core localization conclusion\.

4. The traditional pure static disorder approximation has quantitative bias but qualitative validity for this system; its applicable boundary is clarified via the proposed dual\-parameter criterion\.

> （注：部分内容可能由 AI 生成）
