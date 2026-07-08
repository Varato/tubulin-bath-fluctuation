#!/usr/bin/env python
"""Phase 5 / Task 5.1 — Stochastic Schrödinger MC for Trp4-Trp7 dimer.

Two-site exciton with tri-exp colored noise. Component ablation isolates
each timescale's effect on coherence decay and population transfer.

Hamiltonian:
    H(t) = J sigma_x + ds4(t) |0><0| + ds7(t) |1><1|

where ds4(t), ds7(t) are independent tri-exp OU noise trajectories.

Outputs (results/phase5_exciton_dynamics/):
    mc_results.npz              ensemble-averaged observables per config
    ablation_coherence.png      |rho_47(t)| for each noise config
    ablation_population.png     P7(t) for each noise config
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import setup_style, save_fig, phase_dir, TRP_PRETTY

# ── Constants ──
CM_TO_RADPS = 2 * np.pi * 2.998e-2     # cm^-1 -> rad/ps (hbar=1)

# ── System: Trp4 (aW407) — Trp7 (bW397) ──
J_CM    = -59.0
J       = J_CM * CM_TO_RADPS
SIG4_CM = 998.5       # slow-traj sigma_total, Trp4
SIG7_CM = 1453.2      # slow-traj sigma_total, Trp7
SIG4    = SIG4_CM * CM_TO_RADPS
SIG7    = SIG7_CM * CM_TO_RADPS

# ── Tri-exp noise model (Phase 2 system means) ──
AMPS_FULL = np.array([0.50, 0.33, 0.16])
TAUS_FULL = np.array([0.044, 1.70, 2663.0])

# ── Simulation grid ──
DT      = 0.002       # ps (2 fs)
T_MAX   = 2.0         # ps
TLIST   = np.arange(0, T_MAX + DT / 2, DT)
N_STEPS = len(TLIST)
N_REAL  = 500

# ── Operators ──
ket0 = qt.basis(2, 0)    # |Trp4>
ket1 = qt.basis(2, 1)    # |Trp7>
N4_OP = ket0 * ket0.dag()
N7_OP = ket1 * ket1.dag()
SX_OP = ket0 * ket1.dag() + ket1 * ket0.dag()
SM_OP = ket0 * ket1.dag()    # <sm> = rho_10 -> |rho_47|
H0    = J * SX_OP
PSI0  = ket0


# ── Noise generation ──
def gen_ou(n, dt, var, tau, rng):
    decay = np.exp(-dt / tau)
    noise = np.sqrt(var * (1 - decay**2)) * rng.standard_normal(n)
    x = np.empty(n)
    x[0] = np.sqrt(var) * rng.standard_normal()
    for i in range(1, n):
        x[i] = x[i - 1] * decay + noise[i]
    return x


def gen_noise(n, dt, sigma, amps, taus, rng):
    x = np.zeros(n)
    for A, tau in zip(amps, taus):
        if A > 0:
            x += gen_ou(n, dt, A * sigma**2, tau, rng)
    return x


# ── Ensemble runner ──
def run_ensemble(sig4, sig7, amps, taus, n_real, seed=42):
    """Run n_real noise trajectories, return averaged P4, P7, |rho_47|."""
    P4  = np.zeros((n_real, N_STEPS))
    P7  = np.zeros((n_real, N_STEPS))
    coh = np.zeros((n_real, N_STEPS))

    for i in range(n_real):
        rng = np.random.default_rng(seed + i)
        noise4 = gen_noise(N_STEPS, DT, sig4, amps, taus, rng)
        noise7 = gen_noise(N_STEPS, DT, sig7, amps, taus, rng)
        H_t = [H0, [N4_OP, noise4], [N7_OP, noise7]]
        res = qt.sesolve(H_t, PSI0, TLIST, e_ops=[N4_OP, N7_OP, SM_OP])
        P4[i]  = res.expect[0]
        P7[i]  = res.expect[1]
        coh[i] = np.abs(res.expect[2])

    return dict(
        P4=P4.mean(axis=0), P7=P7.mean(axis=0), coh=coh.mean(axis=0),
        P4_sem=P4.std(axis=0) / np.sqrt(n_real),
        P7_sem=P7.std(axis=0) / np.sqrt(n_real),
        coh_sem=coh.std(axis=0) / np.sqrt(n_real),
    )


def main():
    P = phase_dir(5)
    setup_style()

    print(f"Trp4-Trp7 dimer  J={J_CM} cm^-1  ({J:.2f} rad/ps)")
    print(f"Rabi transfer: {np.pi/(2*abs(J))*1e3:.0f} fs")
    print(f"sigma4={SIG4_CM}, sigma7={SIG7_CM} cm^-1")
    print(f"Grid: {N_STEPS} steps, dt={DT*1e3:.0f} fs, T={T_MAX} ps")
    print(f"Realizations: {N_REAL}\n")

    # ── Noiseless reference ──
    res0 = qt.sesolve(H0, PSI0, TLIST, e_ops=[N4_OP, N7_OP, SM_OP])
    print("Noiseless reference done.\n")

    # ── Ablation configs ──
    configs = [
        ('tau1 only',       np.array([0.50, 0., 0.])),
        ('tau2 only',       np.array([0., 0.33, 0.])),
        ('tau3 only',       np.array([0., 0., 0.16])),
        ('tau1+tau2',       np.array([0.50, 0.33, 0.])),
        ('full',            np.array([0.50, 0.33, 0.16])),
    ]

    results = {}
    for name, amps in configs:
        print(f"  running {name:15s} ...", end=' ', flush=True)
        results[name] = run_ensemble(SIG4, SIG7, amps, TAUS_FULL, N_REAL)
        sig_eff = np.sqrt(amps.sum()) * SIG4_CM
        print(f"sigma_eff={sig_eff:.0f} cm^-1")

    # ── Save ──
    save_dict = dict(tlist=TLIST,
                     P4_noiseless=np.array(res0.expect[0]),
                     P7_noiseless=np.array(res0.expect[1]),
                     coh_noiseless=np.abs(np.array(res0.expect[2])))
    for name, r in results.items():
        tag = name.replace(' ', '_').replace('+', '_')
        for key in ['P4', 'P7', 'coh', 'P4_sem', 'P7_sem', 'coh_sem']:
            save_dict[f'{tag}_{key}'] = r[key]
    np.savez_compressed(P / 'mc_results.npz', **save_dict)

    # ── Key numbers ──
    print(f"\n{'='*55}")
    print(f"{'config':<16} {'coh@20fs':>8} {'coh@100fs':>9} "
          f"{'P7@200fs':>9} {'P7@2ps':>8}")
    print(f"{'-'*55}")
    i20  = int(0.02 / DT)
    i100 = int(0.10 / DT)
    i200 = int(0.20 / DT)
    iend = -1
    for name, r in results.items():
        print(f"{name:<16} {r['coh'][i20]:8.4f} {r['coh'][i100]:9.4f} "
              f"{r['P7'][i200]:9.4f} {r['P7'][iend]:8.4f}")
    print(f"{'noiseless':<16} {np.abs(res0.expect[2])[i20]:8.4f} "
          f"{np.abs(res0.expect[2])[i100]:9.4f} "
          f"{res0.expect[1][i200]:9.4f} {res0.expect[1][iend]:8.4f}")

    # =====================================================================
    # FIGURE 1: Coherence decay
    # =====================================================================
    fig1, ax1 = plt.subplots(figsize=(9, 5.5))
    t_fs = TLIST * 1e3

    style = {
        'tau1 only':   dict(c='#e41a1c', ls='-',  lw=1.3),
        'tau2 only':   dict(c='#377eb8', ls='-',  lw=1.3),
        'tau3 only':   dict(c='#4daf4a', ls='-',  lw=1.3),
        'tau1+tau2':   dict(c='#ff7f00', ls='--', lw=1.5),
        'full':        dict(c='black',   ls='-',  lw=2.0),
    }
    ax1.plot(t_fs, np.abs(res0.expect[2]), ':', color='gray', lw=1.5,
             alpha=0.5, label='noiseless')
    for name, r in results.items():
        ax1.plot(t_fs, r['coh'], label=name, **style[name])

    ax1.set_xlabel('t (fs)')
    ax1.set_ylabel(r'$|\rho_{47}(t)|$')
    ax1.set_title('Coherence decay — component ablation')
    ax1.legend(fontsize=8, ncol=2)
    for tau in TAUS_FULL[:2]:
        ax1.axvline(tau * 1e3, color='gray', ls=':', lw=0.5, alpha=0.3)
    save_fig(fig1, P / 'ablation_coherence.png')

    # =====================================================================
    # FIGURE 2: Populations P4 (dashed) and P7 (solid) on same plot
    # =====================================================================
    # line styles: color per config, solid=P7, dashed=P4 (no ls in style dict)
    pop_style = {k: dict(c=v['c'], lw=v['lw']) for k, v in style.items()}

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    for name, r in results.items():
        c = pop_style[name]['c']
        lw = pop_style[name]['lw']
        ax2.plot(t_fs, r['P7'], '-', color=c, lw=lw, label=name)
        ax2.fill_between(t_fs, r['P7'] - r['P7_sem'], r['P7'] + r['P7_sem'],
                         alpha=0.15, color=c)
        ax2.plot(t_fs, r['P4'], '--', color=c, lw=lw)

    ax2.set_xlabel('t (fs)')
    ax2.set_ylabel('Population')
    ax2.set_title('Population dynamics — excitation starts on Trp4\n'
                  '(solid = $P_7$ Trp7, dashed = $P_4$ Trp4)')
    ax2.legend(fontsize=8, loc='center left', framealpha=0.9)
    ax2.set_ylim(-0.02, 1.02)
    for tau in TAUS_FULL[:2]:
        ax2.axvline(tau * 1e3, color='gray', ls=':', lw=0.5, alpha=0.3)
    save_fig(fig2, P / 'ablation_population.png')

    print(f"\nFigures: {P}/ablation_{{coherence,population}}.png")
    print(f"Data:    {P}/mc_results.npz")


if __name__ == '__main__':
    main()
