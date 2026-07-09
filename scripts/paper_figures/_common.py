"""Shared conventions for paper figures.

Naming: Trp1-Trp8 throughout (per results_outline.md). The alpha/beta
correspondence is stated once in the paper text, not repeated on figures.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# make scripts/utils.py importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils import RESULTS, V_TO_CM, J_COUPLING  # noqa: F401  (re-export)

PAPER_DIR = RESULTS / 'paper_figures'
PAPER_DIR.mkdir(parents=True, exist_ok=True)

# ── Trp naming ──
# index 0..7  ->  Trp1..Trp8 ; alpha/beta correspondence (stated once in text)
TRP_NAMES = [f'Trp{i+1}' for i in range(8)]
TRP_ALPHABETA = ['αW21', 'αW346', 'αW388', 'αW407',
                 'βW21', 'βW101', 'βW344', 'βW397']

# ── Colours (from results_outline.md) ──
SOURCE_COLORS = {
    'protein':    '#7E57C2',
    'water':      '#29B6F6',
    'nucleotide': '#FF7043',
    'ions':       '#2CA02C',   # green: high contrast on white, distinct from blue/orange/purple
}
SOURCE_LABELS = ['protein', 'water', 'nucleotide', 'ions']

_TRP_CMAP = plt.cm.tab10
TRP_COLORS = [_TRP_CMAP(i) for i in range(8)]
TRP_COLORS[6] = _TRP_CMAP(9)   # keep Trp7 visually distinct

# ── unit conversions ──
CM_TO_RADPS = 2 * np.pi * 2.998e-2   # cm^-1 -> rad/ps  (hbar = 1 natural units)
CM_TO_THZ = 2.998e-2                # cm^-1 -> THz


def setup():
    """Paper figure style: clean, compact, no chartjunk."""
    mpl.use('Agg')
    plt.rcParams.update({
        'figure.dpi': 130, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
        'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 11,
        'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 8.5,
        'axes.grid': True, 'grid.alpha': 0.25, 'axes.axisbelow': True,
        'lines.linewidth': 1.3, 'legend.frameon': False,
        'font.family': 'DejaVu Sans',
    })


def save(fig, name: str) -> Path:
    """Save a figure as both PNG (preview) and PDF (for LaTeX)."""
    stem = PAPER_DIR / name
    if stem.suffix == '.png':
        stem = stem.with_suffix('')
    p_png = stem.with_suffix('.png')
    p_pdf = stem.with_suffix('.pdf')
    fig.savefig(p_png)
    fig.savefig(p_pdf)
    plt.close(fig)
    print(f"  saved {p_png}  +  {p_pdf}")
    return p_png
