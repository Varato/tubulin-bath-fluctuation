#!/usr/bin/env python
"""Two diagnostic MC runs for Fig 3(b):

(1) White-noise (Markovian) limit: single OU process per site with tau -> 0
    and sigma^2*tau matched to the Lindblad dephasing rate gamma=50 cm^-1.
    The trajectory MC must converge onto the Haken-Strobl Lindblad curve --
    a validation that the MC code reduces correctly to the Markovian result.

    Matching: per-site white-noise intensity sigma^2*tau equals the Lindblad
    coherence-decay rate gamma_ps (gamma=50 cm^-1 in rad/ps), because the
    site-basis coherence |0><1| decays at rate gamma_ps in both descriptions.

(2) sigma-scaling scan: realistic tri-exp noise but sigma -> alpha*sigma_MD,
    to find the alpha at which P7(2 ps) matches the Lindblad value (0.50).
    Shows how much weaker than MD the bath must be to imitate gamma=50.

Outputs: results/paper_figures/whitenoise_mc.npz, sigma_scan.npz
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import CM_TO_RADPS, PAPER_DIR
import run_exciton_mc as M   # reuse operators, gen_ou, grid, Hamiltonian pieces

GAMMA_CM = 50.0
GAMMA_PS = GAMMA_CM * CM_TO_RADPS            # 9.42 rad/ps (coherence decay rate)
TAU_WN = 0.005                               # ps — white-noise OU tau (<< 1/J ~ 0.09 ps)
SIG_WN = np.sqrt(GAMMA_PS / TAU_WN)          # per-site sigma (rad/ps) so sigma^2 tau = gamma_ps
N_REAL_WN = 1500                             # extra realizations for a clean validation curve


def run_whitenoise():
    """Single-OU-per-site MC in the Markovian limit, matched to gamma=50."""
    J  = M.J_CM * CM_TO_RADPS
    dE = M.D_E_CM * CM_TO_RADPS
    H0 = J * M.SX + (dE / 2) * M.SZ
    amps1 = np.array([1.0, 0.0, 0.0]); taus1 = np.array([TAU_WN, 1.0, 1.0])
    r = M.run_ensemble(SIG_WN, SIG_WN, amps1, taus1, amps1, taus1, N_REAL_WN, H0)
    print(f"[white-noise] sigma={SIG_WN:.2f} rad/ps, tau={TAU_WN} ps, "
          f"sigma^2*tau={SIG_WN**2*TAU_WN:.3f} (=gamma_ps={GAMMA_PS:.3f})")
    print(f"[white-noise] P7(2ps)={r['P7'][-1]:.4f} +- {r['P7_sem'][-1]:.4f}  "
          f"(Lindblad target 0.500)")
    np.savez_compressed(PAPER_DIR / 'whitenoise_mc.npz',
                        t=M.TLIST, P7=r['P7'], P7_sem=r['P7_sem'],
                        P4=r['P4'], coh=r['coh'], sigma_wn=SIG_WN, tau_wn=TAU_WN,
                        gamma_ps=GAMMA_PS)


def run_sigma_scan():
    """Realistic tri-exp noise scaled by alpha; find alpha for P7(2ps)=0.50."""
    sig, tri = M.load_site_params()
    s4, s7 = sig[M.I4], sig[M.I7]
    a4, t4 = tri[M.I4]; a7, t7 = tri[M.I7]
    J  = M.J_CM * CM_TO_RADPS
    dE = M.D_E_CM * CM_TO_RADPS
    H0 = J * M.SX + (dE / 2) * M.SZ
    alphas = [0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.30, 0.50, 1.00]
    P7end = []
    for a in alphas:
        r = M.run_ensemble(a*s4, a*s7, a4, t4, a7, t7, 300, H0)
        P7end.append(float(r['P7'][-1]))
        print(f"[sigma-scan] alpha={a:.2f}  P7(2ps)={P7end[-1]:.3f}")
    # interpolate alpha where P7 = 0.50
    a_arr = np.array(alphas); p_arr = np.array(P7end)
    # find crossing of 0.5 on the ascending branch (small->mid alpha)
    np.savez_compressed(PAPER_DIR / 'sigma_scan.npz',
                        alpha=a_arr, P7_2ps=p_arr)
    # crude interpolation around the 0.50 level
    for i in range(len(a_arr)-1):
        if (p_arr[i]-0.5)*(p_arr[i+1]-0.5) < 0:
            a0,a1 = a_arr[i],a_arr[i+1]; p0,p1=p_arr[i],p_arr[i+1]
            a_cross = a0 + (a1-a0)*(0.5-p0)/(p1-p0)
            print(f"\n=> P7=0.50 at alpha ~ {a_cross:.3f}  "
                  f"(sigma_eff4={a_cross*s4:.0f}, sigma_eff7={a_cross*s7:.0f} cm^-1)")
            return a_cross
    print("\n(no crossing of 0.50 found in scanned range)")
    return None


def main():
    print("=== (1) white-noise Markovian-limit MC ===")
    run_whitenoise()
    print("\n=== (2) sigma-scaling scan ===")
    run_sigma_scan()


if __name__ == '__main__':
    main()
