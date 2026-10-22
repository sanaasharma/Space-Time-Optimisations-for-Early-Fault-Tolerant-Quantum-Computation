"""Regression tests. Run with `python -m pytest tests` (under a minute).

Every compile must terminate: it either succeeds or raises LayoutInfeasible."""
import collections
import io
import math
import contextlib
import random

import pytest
from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag

from Functions import (LayoutInfeasible, space_search, find_best_path, find_initial_layout_2D,
                       add_global_barrier, apply_sequence)
from levels import level_1, level_2
from Postprocessing import generate_time
from script import build_grid, build_circuit, compile_circuit

pytestmark = pytest.mark.timeout(300)

PINNED = [("ising", 0, 480.0, 40), ("heisenberg", 1, 963.0, 72), ("fermi-hubbard", 2, 321.5, 24)]


def run(program, seed=0, edge=6, **grid_options):
    random.seed(seed)
    grid = build_grid(edge, **grid_options)
    qc = build_circuit(program, grid)
    with contextlib.redirect_stdout(io.StringIO()):
        _, dag = compile_circuit(qc, grid)
    return dag


# --- search helpers -----------------------------------------------------------------------

def test_space_search_raises_when_grid_is_full():
    with pytest.raises(LayoutInfeasible):
        space_search((0, 0), [[1, 2], [3, 4]], "h")


def test_space_search_reaches_far_corner():
    grid = [[9] * 6 for _ in range(6)]
    grid[5][5] = 0  # the far corner, Manhattan distance 10 from (0, 0)
    assert space_search((0, 0), grid, "h") == {10: [(5, 5)]}


def test_find_best_path_raises_when_target_is_walled_in():
    grid = find_initial_layout_2D(6, 16)
    for i, j in [(1, 3), (3, 3), (1, 5), (3, 5)]:  # diagonal neighbours of qubit 8 at (2, 4)
        grid[i][j] = -1
    with pytest.raises(LayoutInfeasible):
        find_best_path(grid, 1, 1, 2, 4, False)


def test_circuit_must_match_grid():
    grid = find_initial_layout_2D(6, 16)
    dag = add_global_barrier(circuit_to_dag(QuantumCircuit(49)))
    with pytest.raises(ValueError, match="grid"):
        apply_sequence(dag, grid, {"level1": level_1, "level2": level_2}, 1)


# --- every gate is executed on the right qubits ------------------------------------------
# The apply_* functions raise RoutingError for physically illegal operations (non-adjacent
# magic state or ancilla, occupied ancilla, wrong CNOT geometry). These tests also check that
# each gate acts on the logical qubits the input circuit asked for.

@pytest.mark.parametrize("program, edge, seed", [
    ("ising", 6, 0), ("heisenberg", 6, 1), ("fermi-hubbard", 8, 2), ("ising", 9, 3),
])
def test_gates_act_on_the_intended_qubits(program, edge, seed, monkeypatch):
    import levels
    random.seed(seed)
    grid = build_grid(edge)
    cols = len(grid[0])
    label = {i * cols + j: v for i, row in enumerate(grid) for j, v in enumerate(row)}
    wrong = []

    def checked_cnot(ancilla, node_qargs, g, dag, node):
        applied = tuple(g[q // cols][q % cols] for q in node_qargs)
        intended = tuple(label[q._index] for q in node.qargs)
        if applied != intended:
            wrong.append(("cx", intended, applied))
        return Gates.apply_CNOT(ancilla, node_qargs, g, dag, node)

    def checked_h(ancilla, qindex, g, dag, node, remove_rotation=None):
        intended = label[(remove_rotation or node).qargs[0]._index]
        if g[qindex // cols][qindex % cols] != intended:
            wrong.append(("h", intended, g[qindex // cols][qindex % cols]))
        return Gates.apply_H(ancilla, qindex, g, dag, node, remove_rotation)

    import Gates, Rotations
    monkeypatch.setattr(levels, "apply_CNOT", checked_cnot)
    monkeypatch.setattr(levels, "apply_H", checked_h)
    monkeypatch.setattr(Rotations, "apply_H", checked_h)
    with contextlib.redirect_stdout(io.StringIO()):
        _, dag = compile_circuit(build_circuit(program, grid), grid)
    assert wrong == []
    assert generate_time(dag)[3] == sum(1 for n in dag.op_nodes() if n.op.name == "ct")


# Fixed-seed results on the default layout (random.seed(seed)).
@pytest.mark.parametrize("program, seed, total, n_magic", PINNED)
def test_default_layout_regression(program, seed, total, n_magic):
    overheads, totals, _, n = generate_time(run(program, seed))
    assert (totals[0], n) == (total, n_magic)


# --- rotations ----------------------------------------------------------------------------

def test_multiples_of_pi_over_4_are_exact():
    from Rotations import decompose_rz
    assert decompose_rz(math.pi / 4, 1e-4) == ("t",)
    assert decompose_rz(-math.pi / 4, 1e-4) == ("tdg",)
    assert decompose_rz(math.pi / 2, 1e-4) == ("s",)
    assert decompose_rz(0.0, 1e-4) == ()


@pytest.mark.parametrize("angle", [0.3, 1.234, -2.0])
@pytest.mark.parametrize("epsilon", [1e-3, 1e-6])
def test_arbitrary_angles_are_approximated_to_epsilon(angle, epsilon):
    import numpy as np
    from qiskit.quantum_info import Operator
    from Rotations import decompose_rz
    gates = decompose_rz(angle, epsilon)
    approx = QuantumCircuit(1)
    for g in gates:
        getattr(approx, g)(0)
    exact = QuantumCircuit(1)
    exact.rz(angle, 0)
    overlap = abs(np.trace(Operator(exact).data.conj().T @ Operator(approx).data)) / 2
    assert math.sqrt(max(0.0, 1 - overlap)) <= epsilon
    assert gates.count("t") + gates.count("tdg") <= 4 * math.log2(1 / epsilon) + 10


def test_arbitrary_angle_benchmark_compiles():
    from Rotations import decompose_rz, ROTATION_EPSILON
    random.seed(0)
    grid = build_grid(5)
    qc = build_circuit("fermi-hubbard", grid, angle=0.3)
    with contextlib.redirect_stdout(io.StringIO()):
        _, dag = compile_circuit(qc, grid)
    t_per_rotation = sum(g in ("t", "tdg") for g in decompose_rz(0.3, ROTATION_EPSILON))
    assert generate_time(dag)[3] == qc.count_ops()["rz"] * t_per_rotation


def test_user_circuit_with_unsupported_gates():
    from script import place_circuit
    qc = QuantumCircuit(5, 5)
    qc.h(0); qc.cz(0, 1); qc.rx(0.7, 2); qc.ry(0.2, 3); qc.swap(3, 4); qc.ccx(0, 1, 2)
    qc.measure(range(5), range(5))
    random.seed(0)
    grid = build_grid(5)
    placed = place_circuit(qc, grid)
    assert set(placed.count_ops()) <= {"h", "cx", "rz", "s", "sdg", "t", "tdg", "x", "sx"}
    with contextlib.redirect_stdout(io.StringIO()):
        _, dag = compile_circuit(placed, grid)
    assert generate_time(dag)[3] > 0


# the whole routing ring of the 6x6 grid except the magic state factory cell (3, 0)
RING = [(i, j) for i in range(6) for j in range(6) if (i in (0, 5) or j in (0, 5)) and (i, j) != (3, 0)]


# --- reduced free space: terminate with success or LayoutInfeasible ------------------------

@pytest.mark.parametrize("program, options", [
    ("ising", {"blocked": [(0, 1)]}),                       # one top-row cell blocked
    ("heisenberg", {"blocked": [(0, 1)]}),
    ("heisenberg", {"blocked": [(1, 5)]}),                  # one right-column cell blocked
    ("fermi-hubbard", {"blocked": [(1, 5), (2, 5), (3, 5)]}),  # reported segmentation fault
    ("ising", {"blocked": RING}),
])
def test_blocked_cells_terminate(program, options):
    try:
        run(program, **options)
    except LayoutInfeasible:
        pass


def test_no_free_space_is_reported_as_infeasible():
    with pytest.raises(LayoutInfeasible):
        run("ising", blocked=RING)


@pytest.mark.parametrize("add_paths", [1, 2, 4])
def test_added_routing_paths_compile(add_paths):
    dag = run("ising", add_paths=add_paths)
    assert generate_time(dag)[3] == 40


# --- other configurations the README advertises ------------------------------------------

@pytest.mark.parametrize("edge", [5, 7])
def test_odd_edge_sizes(edge):
    run("fermi-hubbard", edge=edge)


@pytest.mark.parametrize("factories", [2, 4])
def test_multiple_magic_state_factories(factories):
    random.seed(0)
    grid = build_grid(6)
    with contextlib.redirect_stdout(io.StringIO()):
        _, dag = compile_circuit(build_circuit("fermi-hubbard", grid), grid, factories)
    assert generate_time(dag)[3] == 24


def test_t_gates():
    random.seed(0)
    grid = build_grid(6)
    qc = QuantumCircuit(36)
    for q in (7, 8, 13, 14):
        qc.h(q); qc.t(q); qc.tdg(q); qc.x(q)
    with contextlib.redirect_stdout(io.StringIO()):
        _, dag = compile_circuit(qc, grid)
    assert generate_time(dag)[3] == 8


# --- your own initial layouts --------------------------------------------------------------

def test_layout_file_with_row_removed_and_blocked_cell(tmp_path):
    from script import load_layout
    path = tmp_path / "layout.txt"
    path.write_text("# top row of the default 6x6 grid removed\n"
                    "0 1 2 3 4 0\n0 5 6 7 8 -1\n0 9 10 11 12 0\n0 13 14 15 16 0\n0 0 0 0 0 0\n")
    grid = load_layout(path)
    random.seed(0)
    with contextlib.redirect_stdout(io.StringIO()):
        _, dag = compile_circuit(build_circuit("ising", grid), grid)
    assert generate_time(dag)[3] == 40


def test_freely_placed_data_qubits():
    grid = [[0] * 7 for _ in range(7)]
    for k, (i, j) in enumerate([(1, 1), (1, 3), (1, 5), (3, 1), (3, 3), (3, 5), (5, 1), (5, 3), (5, 5)], 1):
        grid[i][j] = k
    random.seed(0)
    with contextlib.redirect_stdout(io.StringIO()):
        _, dag = compile_circuit(build_circuit("ising", grid), grid)
    assert generate_time(dag)[3] == 21


@pytest.mark.parametrize("change, message", [
    (lambda g: g[2].pop(), "rectangle"),
    (lambda g: g[0].__setitem__(0, 5), "duplicated"),
    (lambda g: g[3].__setitem__(0, -1), "factory"),
    (lambda g: g[0].__setitem__(0, "x"), "cell"),
])
def test_invalid_layouts_are_rejected(change, message):
    from Functions import validate_grid
    grid = build_grid(6)
    change(grid)
    with pytest.raises(ValueError, match=message):
        validate_grid(grid)
