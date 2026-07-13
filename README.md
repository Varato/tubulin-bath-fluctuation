# Tubulin Trp Exciton Coloured-Noise Analysis

Analysis code for the paper:

> **Molecular Dynamics-Derived Coloured Noise Mediates Anderson Localisation and Environment-Assisted Transport of Tryptophan Excitons in Tubulin**

The scripts in this repository quantify multi-timescale site-energy fluctuations at the eight tryptophan (Trp) sites of the tubulin dimer, validate a fast/slow coloured-noise model, and assess exciton localisation and environment-assisted quantum transport on the Trp network.

## Data

Running the analysis requires **two molecular dynamics trajectories** of the tubulin dimer that are **not bundled with this repository** (~100 GB total). The data is available from the author on request. Once obtained, the two pre-processed NPZ files should be placed in a directory whose path is set as `TUBULIN_DATA_DIR` (see below):

| File | Sampling | Frames | Band |
|------|----------|--------|------|
| `noise_50ns@10ps.npz` | 10 ps | 4001 | ns-scale slow conformational fluctuations |
| `noise_2ns@10fs.npz` | 10 fs | 200001 | sub-ps fast environmental fluctuations |

The two trajectories resolve complementary frequency bands and are stitched together to reconstruct the full fluctuation spectrum. Schema and Trp site indexing are documented in `RESEARCH_PLAN.md` §1.3-1.4.

## Setup

The project uses a flat-scripts layout with a single shared module (`scripts/utils.py`). Dependencies: numpy 2.x, scipy, matplotlib, python-dotenv.

```
conda create -n qbio python=3.12 numpy scipy matplotlib python-dotenv
conda activate qbio
```

Create a `.env` file in the repo root pointing at the data directory:

```
# .env (not checked in)
TUBULIN_DATA_DIR=/path/to/noise_analysis/data
```
