"""Shared utilities for tubulin Trp site-energy fluctuation analysis.

Single shared module for the flat-scripts layout (RESEARCH_PLAN.md §1.6).
All time units: ps. All energy/frequency units: cm^-1.
"""
from __future__ import annotations
import os
from pathlib import Path
from dotenv import load_dotenv
import numpy as np

# numpy 2.0 renamed np.trapz -> np.trapezoid; keep a single shim everywhere
trapz = getattr(np, 'trapezoid', None) or np.trapz

# ---------------------------------------------------------------- constants
V_TO_CM = 8.397e-7          # V/m -> cm^-1 (assumes |Δμ| = 5 D for Trp S0->S1)
PS_INV_TO_CM = 33.356       # ps^-1 -> cm^-1
J_COUPLING = 45.0           # cm^-1, fixed Trp exciton nearest-neighbour coupling
T_OBS_PS = 2.0              # ps, exciton observation window

N_TRP = 8
TRP_LABELS = ['aW21', 'aW346', 'aW388', 'aW407',
              'bW21', 'bW101', 'bW344', 'bW397']
TRP_PRETTY = ['αW21', 'αW346', 'αW388', 'αW407',
              'βW21', 'βW101', 'βW344', 'βW397']
TRP_CHAIN = ['A', 'A', 'A', 'A', 'B', 'B', 'B', 'B']

COMPONENTS = ['protein', 'water', 'nucleotide', 'ions']

# ---------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parent.parent
load_dotenv(ROOT / '.env')   # machine-specific paths live in .env (gitignored)

DATA_DIR = Path(os.environ.get('TUBULIN_DATA_DIR', ''))
if not DATA_DIR.exists():
    raise FileNotFoundError(
        f"TUBULIN_DATA_DIR is unset or points to a missing path: '{DATA_DIR}'\n"
        f"  Set it in {ROOT}/.env  (e.g. TUBULIN_DATA_DIR=/path/to/noise_analysis/data)"
    )
RESULTS = ROOT / 'results'
RESULTS.mkdir(parents=True, exist_ok=True)

# canonical phase folder names (keep in sync with RESEARCH_PLAN.md §2–§7)
PHASES = {
    0: 'phase0_validation',
    1: 'phase1_basic_stats',
    2: 'phase2_timescale',
    3: 'phase3_source_attribution',
    4: 'phase4_spatial_correlation',
    5: 'phase5_exciton_dynamics',
    6: 'phase6_methodology',
}


def phase_dir(phase) -> Path:
    """Return (creating) the results folder for a phase.

    phase : int (0–6, see PHASES) or an explicit folder-name string.
    All outputs for that phase — CSV, NPZ, PNG — go in this single folder.
    """
    name = PHASES.get(phase, phase) if not isinstance(phase, str) else phase
    p = RESULTS / name
    p.mkdir(parents=True, exist_ok=True)
    return p

FILES = {
    'slow': DATA_DIR / 'noise_50ns@10ps.npz',
    'fast': DATA_DIR / 'noise_2ns@10fs.npz',
}

# ---------------------------------------------------------------- io
def load(tag: str) -> dict:
    """Load a trajectory dataset as a dict of raw arrays.

    tag : 'slow' (50 ns @ 10 ps) or 'fast' (2 ns @ 10 fs).
    Arrays are returned raw (V/m for E-fields and delta_s; convert with V_TO_CM).
    Adds convenience scalars: dt_ps, n_frames, label.
    """
    if tag not in FILES:
        raise KeyError(f"unknown tag '{tag}', choose from {list(FILES)}")
    if not FILES[tag].exists():
        raise FileNotFoundError(f"missing dataset: {FILES[tag]}")
    npz = np.load(FILES[tag])
    out = {k: npz[k] for k in npz.files}
    out['dt_ps'] = float(out['dt_ps'])
    out['n_frames'] = int(out['time_ps'].shape[0])
    out['label'] = tag
    return out


def delta_s_cm(d: dict, comp: str = 'total') -> np.ndarray:
    """Return delta_s_<comp> in cm^-1, shape (N, 8)."""
    return d[f'delta_s_{comp}'] * V_TO_CM


def write_csv(path, header: str, rows) -> Path:
    """Write a simple CSV. `rows` is an iterable of tuples; values are
    formatted with {!s} (so pass pre-formatted strings/numbers)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as fh:
        fh.write(header.rstrip('\n') + '\n')
        for r in rows:
            fh.write(','.join(str(v) for v in r) + '\n')
    return path


def acf_fft(x, dt):
    """FFT-based normalized autocorrelation function (linear / non-circular).

    x : (N,) or (N, M) array. ACF computed along axis 0.
    dt: sampling interval (same units as returned t_lag).
    Returns (t_lag, acf) where acf has the same shape as x and acf[0] = 1.
    Zero-padded to >= 2N to avoid wraparound.
    """
    x = np.asarray(x, dtype=float)
    x = x - x.mean(axis=0)
    N = x.shape[0]
    nfft = 1 << int(np.ceil(np.log2(2 * N)))
    F = np.fft.rfft(x, n=nfft, axis=0)
    acf = np.fft.irfft(F * np.conj(F), n=nfft, axis=0)[:N]
    acf /= acf[0]                           # normalize per column → C(0)=1
    t = np.arange(N) * dt
    return t, acf


def log_interp(t_new, t_orig, y_orig, axis=0):
    """Linear interpolation that works for monotonic (increasing or decreasing)
    t_orig. 1D in t but y may be multi-column along `axis`."""
    t_orig = np.asarray(t_orig)
    if t_orig[0] > t_orig[-1]:
        sl = slice(None, None, -1)
        t_orig = t_orig[sl]
        y_orig = np.take(y_orig, np.arange(len(t_orig))[sl], axis=axis)
        axis = 0
    t_new = np.asarray(t_new)
    if y_orig.ndim == 1:
        return np.interp(t_new, t_orig, y_orig)
    out = np.empty((len(t_new),) + y_orig.shape[1:])
    for j in range(y_orig.shape[1]):
        out[:, j] = np.interp(t_new, t_orig, y_orig[:, j])
    return out



# ---------------------------------------------------------------- plotting
def setup_style():
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.use('Agg')
    plt.rcParams.update({
        'figure.dpi': 130,
        'savefig.dpi': 200,
        'savefig.bbox': 'tight',
        'font.size': 10,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'axes.axisbelow': True,
        'lines.linewidth': 1.2,
        'legend.frameon': False,
    })


def save_fig(fig, path) -> Path:
    """Save a figure. `path` may be a Path (preferred, from phase_dir) or a
    bare filename (then it goes into the current phase_dir set via
    set_phase(), or RESULTS root if none)."""
    p = Path(path)
    if not p.is_absolute() and len(p.parts) == 1:
        p = (_CURRENT_PHASE or RESULTS) / p
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p)
    return p


_CURRENT_PHASE = None


def set_phase(phase) -> Path:
    """Set the current phase folder (returns it) so save_fig(fig, 'x.png')
    resolves there. Most scripts will just call phase_dir() explicitly."""
    global _CURRENT_PHASE
    _CURRENT_PHASE = phase_dir(phase)
    return _CURRENT_PHASE
