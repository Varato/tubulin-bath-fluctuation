# OU sampling consistency with tri-exponential ACF

## The target: tri-exponential ACF

From the MD fit, each site m has normalized ACF:

```
C_m(t) = f₁ exp(-t/T₁) + f₂ exp(-t/T₂) + f₃ exp(-t/T₃),   Σfₖ = 1
```

The actual (unnormalized) ACF is σ_m² C_m(t), where σ_m is the total
fluctuation amplitude for site m.

Each component k contributes a variance of fₖ σ_m² to the total σ_m².

## Standard OU process

The SDE `dx = -x/τ dt + σ_SDE dW` has:

- Stationary variance: (1/2) σ_SDE² τ
- ACF: ⟨x(0)x(t)⟩ = (1/2) σ_SDE² τ · exp(-|t|/τ)

The amplitude (1/2)σ_SDE²τ carries dimensions of [amplitude² × time].

## Our code parameterizes by variance, not SDE amplitude

```python
def gen_ou(n, dt, var, tau, rng):
    decay = exp(-dt / tau)
    noise = sqrt(var * (1 - decay**2)) * rng.standard_normal(n)
    x[0] = sqrt(var) * rng.standard_normal()
    for i in range(1, n):
        x[i] = x[i-1] * decay + noise[i]
    return x
```

Recursion: x(t+Δt) = x(t) e^{-Δt/T} + √[var (1 - e^{-2Δt/T})] ξ

Stationary check: Var[x] = e^{-2Δt/T} Var[x] + var (1 - e^{-2Δt/T})
                 → Var[x] = var  ✓

ACF: ⟨x(0)x(t)⟩ = var · exp(-|t|/T)

**Mapping to standard OU:** var = (1/2) σ_SDE² T, i.e. σ_SDE = √(2·var/T)

Our code never needs σ_SDE; it works directly with the target variance.

## How gen_noise calls gen_ou

```python
def gen_noise(n, dt, sigma, amps, taus, rng):
    x = zeros(n)
    for A, tau in zip(amps, taus):
        x += gen_ou(n, dt, A * sigma**2, tau, rng)
    return x
```

Here A = fₖ (tri-exp weight), sigma = σ_m (total amplitude for site m).

So each OU component has:
- var = fₖ σ_m²
- ACF = fₖ σ_m² · exp(-|t|/Tₖ)

## Sum of three OU components

The total noise is δε(t) = x₁(t) + x₂(t) + x₃(t).

ACF of the sum (components are independent):

⟨δε(0) δε(t)⟩ = Σₖ fₖ σ_m² exp(-|t|/Tₖ)
              = σ_m² · C_m(t)  ✓

This matches the MD target exactly. The verification figure (OU sum ACF
overlaid on the tri-exp target) confirms this numerically.

## For the paper appendix

The current formula is:

```
x(t+Δt) = x(t) e^{-Δt/T} + √(1 - e^{-2Δt/T}) · σ√f · ξ(t)
```

where σ = σ_m and f = fₖ.

The factor σ√f is the standard deviation of the stationary distribution
of this OU component (variance = f σ_m² = fₖ σ_m²).

Suggested clarifying sentence after the formula:

"Here σ is the total site-energy fluctuation amplitude (Table 1) and f
is the fractional weight of the k-th mode, so that the stationary
variance of each OU component is fₖ σ_m² and the sum of three
independent components reproduces the target ACF σ_m² C_m(t)."
