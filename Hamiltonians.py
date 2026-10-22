
from qiskit import QuantumCircuit
import numpy as np
import math
from Functions import find_layout, validate_grid


# Every benchmark is built from the grid it will be compiled on. Data qubit k (k = 1..L*L,
# the labels written into the grid by find_initial_layout_2D) is the logical lattice site
# (row, col) = divmod(k-1, L), and its circuit qubit is its grid cell i*cols + j, so the
# benchmarks work on any layout: added routing paths (add_paths_upto), blocked cells (-1)
# or a custom placement of the data qubits.
# `angle` is the rz angle of every interaction term (2 * coupling * time step). The paper uses
# pi/4 (one T gate per term); other angles are approximated in Clifford+T (Rotations.py).

def _lattice(grid):
    """Returns (circuit with one qubit per grid cell, side L, layout) where layout[k-1]
    is the circuit qubit of data qubit k."""
    validate_grid(grid)
    n_data = sum(1 for row in grid for v in row if isinstance(v, int) and v > 0)
    L = math.isqrt(n_data)
    if L * L != n_data:
        raise ValueError(f"the lattice benchmarks need a square number of data qubits, the grid has {n_data}")
    layout = find_layout(grid, n_data)
    return QuantumCircuit(len(grid) * len(grid[0])), L, layout


def _pairs(L, dir, grp):
    """Nearest-neighbour pairs of the L x L lattice, in benchmark order.
    dir=True: horizontal pairs, False: vertical. grp=True: pairs starting at even
    offsets, False: odd."""
    start = 0 if grp else 1
    for i in range(L):
        for j in range(start, L - 1, 2):
            if dir:
                yield (i, j), (i, j + 1)
            else:
                yield (j, i), (j + 1, i)


def fermi_hubbard(grid, angle=np.pi/4):

    qc, _, initial_layout = _lattice(grid)

    def apply_double_x(qc, initial_layout):

        for i in range(0,len(initial_layout)-1,2):

            qc.h([initial_layout[i], initial_layout[i+1]])
            qc.cx(initial_layout[i], initial_layout[i+1])
            qc.rz(angle, initial_layout[i+1])
            qc.cx(initial_layout[i], initial_layout[i+1])
            qc.h([initial_layout[i], initial_layout[i+1]])

        return qc


    def apply_double_z(qc, initial_layout):

        x=int(len(initial_layout)/2)

        for i in range(0,x):
                qc.cx(initial_layout[i], initial_layout[i+x])
                qc.rz(angle, initial_layout[i+x])
                qc.cx(initial_layout[i], initial_layout[i+x])

        return qc

    def apply_double_y(qc,initial_layout):

        for i in range(0,len(initial_layout)-1,2):

            qc.s([initial_layout[i], initial_layout[i+1]])
            qc.h([initial_layout[i], initial_layout[i+1]])
            qc.cx(initial_layout[i], initial_layout[i+1])
            qc.rz(angle, initial_layout[i+1])
            qc.cx(initial_layout[i], initial_layout[i+1])
            qc.h([initial_layout[i], initial_layout[i+1]])
            qc.s([initial_layout[i], initial_layout[i+1]])

        return qc

    qc= apply_double_x(qc,initial_layout)
    qc= apply_double_y(qc,initial_layout)
    qc= apply_double_z(qc,initial_layout)

    return qc

def ising_2D(grid, angle=np.pi/4):

    qc, L, initial_layout = _lattice(grid)

    def q(site):
        return initial_layout[site[0] * L + site[1]]

    def apply_all_rx(qc, initial_layout):

        for q in initial_layout :
            qc.h(q)
            qc.rz(angle, q)
            qc.h(q)

        return qc

    def apply_double_z(qc, dir, grp):

        for s1, s2 in _pairs(L, dir, grp):
                q1, q2 = q(s1), q(s2)

                qc.cx(q1, q2)
                qc.rz(angle, q2)
                qc.cx(q1, q2)

        return qc

    for y in initial_layout:
        qc.h(y)

    for (grp,dir) in [(True,True),(True,False),(False,True),(False,False)]:
        qc= apply_double_z(qc,dir,grp)

    qc= apply_all_rx(qc, initial_layout)


    return qc

def heisenberg(grid, angle=np.pi/4):

    qc, L, initial_layout = _lattice(grid)

    def q(site):
        return initial_layout[site[0] * L + site[1]]

    def apply_double_y(qc, dir, grp):

        for s1, s2 in _pairs(L, dir, grp):
                q1, q2 = q(s1), q(s2)

                qc.sdg([q1, q2])
                qc.h([q1, q2])
                qc.cx(q1, q2)
                qc.rz(angle, q2)
                qc.cx(q1, q2)
                qc.h([q1, q2])
                qc.s([q1, q2])

        return qc

    def apply_double_x(qc, dir, grp):

        for s1, s2 in _pairs(L, dir, grp):
                q1, q2 = q(s1), q(s2)

                qc.h([q1, q2])
                qc.cx(q1, q2)
                qc.rz(angle, q2)
                qc.cx(q1, q2)
                qc.h([q1, q2])

        return qc

    def apply_double_z(qc, dir, grp):

        for s1, s2 in _pairs(L, dir, grp):
                q1, q2 = q(s1), q(s2)

                qc.cx(q1, q2)
                qc.rz(angle, q2)
                qc.cx(q1, q2)

        return qc

    for (grp,dir) in [(True,True),(True,False),(False,True),(False,False)]:
        qc= apply_double_x(qc,dir,grp)
        qc= apply_double_y(qc,dir,grp)
        qc= apply_double_z(qc,dir,grp)

    return qc


BENCHMARKS = {"ising": ising_2D, "heisenberg": heisenberg, "fermi-hubbard": fermi_hubbard}
