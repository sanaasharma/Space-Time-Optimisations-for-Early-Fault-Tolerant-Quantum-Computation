# Space-Time Optimisations for Early Fault-Tolerant Quantum Computation
 

## Setup

Requires Python 3.12, **Qiskit 2.0.0** (other Qiskit versions are not supported; 2.5.x fails at import) and `pygridsynth` (Clifford+T approximation of rotations).

```bash
python3.12 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Running an experiment

```bash
python script.py --program ising --edge 6 --iterations 5
```

This compiles the 2D transverse-field Ising benchmark on a 6x6 grid (4x4 = 16 data qubits surrounded by a ring of routing cells) five times. It prints the execution time of each run and writes the minimum to `output.txt`. 

```
t_MSF=11 magic_states=40 total_time=480.0 overhead=37.0 lower_bound=440
```

Options (`python script.py --help`):

| option | meaning |
|---|---|
| `--program` | `ising`, `heisenberg` or `fermi-hubbard` (defined in `Hamiltonians.py`) |
| `--edge N` | grid of N x N cells with (N-2)^2 data qubits; the paper uses N = 4 to 12 |
| `--steps K` | number of Trotter steps |
| `--factories F` | number of magic state factories |
| `--add-paths X` | insert extra routing rows and columns 0..X with `add_paths_upto` |
| `--msf-time T ...` | magic state factory processing time(s) used by `generate_time` |
| `--angle A` | rz angle of every benchmark interaction term, in radians (default and paper: pi/4) |
| `--rotation-epsilon E` | accuracy of the Clifford+T approximation of rotations that are not multiples of pi/4 (default 1e-4) |
| `--qasm FILE` | compile your own OpenQASM 2 circuit instead of a benchmark  |
| `--layout FILE` | use your own initial layout instead of the default grid (see below) |
| `--iterations`, `--seed` | number of runs and the first seed |

### Rotations and your own circuits

The router executes H, S, S†, T, T†, X, √X, CNOT and Rz. Every Rz becomes a sequence of Clifford+T operations, and every T or T† consumes one magic state:

 

### Output

`generate_time` (in `Postprocessing.py`) returns `(overheads, totals, dag, n_magic_states)`, with one entry per factory processing time t_MSF:

* `n_magic_states` is the number of magic states consumed, which is the number of `ct` operations. The factory lower bound is `n_magic_states * t_MSF`.
* `totals` is the estimated execution time. The logical qubit waits for the factory, so each stretch between two consecutive magic-state consumptions counts as at least t_MSF, plus the work after the last magic state.
* `overheads` is the time by which those stretches (and the final one) exceed t_MSF.
  

* **Layouts.** A grid is a list of rows. `0` is a free routing cell, a positive integer is a data qubit, and `-1` is an unusable cell. You can change the default layout from `find_initial_layout_2D`, add routing paths with `add_paths_upto`, or block cells.
* **Circuits.** Always build the circuit from the grid it will be compiled on. Circuit qubit `i*cols + j` is grid cell `(i, j)`. `apply_sequence` raises `ValueError` if the circuit and the grid have different sizes.
* **Infeasible layouts.** If the grid has too little free space to route the circuit, the compiler raises `Functions.LayoutInfeasible` with an explanation. When routing stalls, the compiler first moves one waiting qubit to break the deadlock; whether that succeeds can depend on the seed, so use `--iterations` to try several.
 

### Your own initial layouts

Write the layout as a text file, one grid row per line:

```
# 0 = free routing cell, -1 = unusable cell, 1..N = data qubits (each label once)
0  1  2  3  4  0
0  5  6  7  8 -1
0  9 10 11 12  0
0 13 14 15 16  0
0  0  0  0  0  0
```

```bash
python script.py --program heisenberg --layout my_layout.txt --iterations 5
```

The grid can be any rectangle, and the data qubits can sit anywhere in it. The benchmarks need a square number of data qubits (an L x L lattice; data qubit k is lattice site `divmod(k-1, L)`). `--qasm` circuits can use any number up to N. The magic state factory outputs are at the midpoints of the grid edges (`Functions.factory_position`: the left edge for one factory, then right, top and bottom), and those cells must not be -1. Invalid layouts are rejected with a `ValueError` that names the problem. On tight layouts some seeds can stall, so use `--iterations` to try several.

## Tests

```bash
python -m pytest tests
```

The tests check that every gate acts on the logical qubits the circuit asks for, and pin the results of a few fixed-seed compiles. They check that rotations are exact for multiples of π/4 and within ε otherwise, and that your own circuits with unsupported gates compile. They also check custom layouts and their validation, and that reduced-space layouts, added routing paths, odd grid sizes, several factories and T gates terminate.
 
## Citation

If you use this code, please cite:

```bibtex
@inproceedings{inproceedings,
author = {Sharma, Sanaa and Murali, Prakash},
year = {2026},
month = {01},
pages = {281-294},
title = {Space-Time Optimisations for Early Fault-Tolerant Quantum Computation},
doi = {10.1109/CGO68049.2026.11395205}
}
```

 
