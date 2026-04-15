# BBGKY for 1D Spin-Boson Model

This directory contains a BBGKY-based hierarchy propagator for a 1D spin-boson model.

## Files

- `BBGKY.py`: core hierarchy builder, filtering utilities, and time propagator.
- `run.py`: main driver script. Builds the hierarchy, propagates the density matrix, and writes populations to `Pt.txt`.
- `model/spinboson_1D.py`: 1D spin-boson model definition, including Hamiltonian, coupling operator, and default parameters.
- `plotting.py`: simple plotting script for data in `Pt.txt`.
- `Pt.txt`: latest output written by `run.py`.
- `Pt1_6.txt`, `Pt3_50.txt`, `Pt3_75.txt`: older reference/output data files.
- `results/`: auxiliary output directory.

## Model

The model in `model/spinboson_1D.py` is a two-level system coupled to one bosonic mode.

Default parameters:

- total time: `t = 30`
- time step: `dt = 0.001`
- number of electronic states: `2`
- number of bath modes: `1`
- default coupling: `sqrt(2) / 10`

The system Hamiltonian is

```text
Hs = [[0, 0.5],
      [0.5, 0]]
```

and the system coupling operator is

```text
Q = [[ 1, 0],
     [ 0,-1]]
```

## How to Run

From this directory:

```powershell
python run.py
```

The script prints progress to the terminal and writes sampled populations to `Pt.txt`.

To visualize the output:

```powershell
python plotting.py
```

## Main Runtime Settings

These are set in `run.py`:

- `TIERS`: hierarchy truncation level for the L1 hierarchy.
- `FILTER_THRESHOLD`: norm threshold for pruning small ADOs.
- `FILTER_INTERVAL`: how often filtering is applied. Set `0` to disable filtering.

## Notes

- This code now uses a true RK4 stage evaluation inside `BBGKY.py`.
- For the 1D spin-boson model, the L1 truncation in `run.py` is the active hierarchy choice.
- Filtering can improve speed, but for sensitive long-time dynamics it is worth checking convergence with filtering disabled.
- Convergence should be checked by varying `TIERS`, `dt`, and filtering settings.
