"""Compile a Hamiltonian-simulation benchmark (or your own OpenQASM circuit) onto a 2D
lattice-surgery grid and estimate its execution time.  Run `python script.py --help` for the options."""
import argparse
import math
import random

import numpy as np
from qiskit import QuantumCircuit, transpile, qasm2
from qiskit.converters import circuit_to_dag
from qiskit.circuit import Measure
from qiskit.circuit.library import Barrier

from Functions import (add_global_barrier, apply_sequence, find_initial_layout_2D, add_paths_upto,
                       find_layout, validate_grid, LayoutInfeasible)
from Hamiltonians import BENCHMARKS
from Postprocessing import generate_time
from Rotations import set_rotation_epsilon
from levels import level_1, level_2

# gates the router executes directly; anything else is translated to these first
SUPPORTED_GATES = ["h", "cx", "rz", "s", "sdg", "t", "tdg", "x", "sx"]


def build_grid(edge_size, add_paths=None, blocked=()):
    """Default layout: an edge_size x edge_size grid whose interior holds (edge_size-2)^2 data
    qubits, surrounded by a ring of routing cells (0).
    add_paths: if given, insert routing paths 0..add_paths with add_paths_upto.
    blocked: (row, col) cells to mark unusable (-1); they must be free routing cells."""
    grid = find_initial_layout_2D(edge_size, (edge_size - 2) ** 2)
    if add_paths is not None:
        grid = add_paths_upto(grid, add_paths)
    for i, j in blocked:
        if grid[i][j] != 0:
            raise ValueError(f"cell {(i, j)} holds data qubit {grid[i][j]}; only free cells can be blocked")
        grid[i][j] = -1
    return grid


def load_layout(path):
    """Read a layout file: one grid row per line, cells separated by spaces or commas.
    0 = free routing cell, -1 = unusable cell, 1..N = data qubits. Lines starting with # are ignored."""
    grid = []
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].replace(",", " ").split()
            if line:
                grid.append([int(v) for v in line])
    validate_grid(grid)
    return grid


def build_circuit(program, grid, steps=1, angle=np.pi / 4):
    """The benchmark circuit for `grid`; steps > 1 repeats it (Trotter steps).
    angle: rz angle of every interaction term (the paper uses pi/4)."""
    one_step = BENCHMARKS[program](grid, angle)
    qc = one_step.copy()
    for _ in range(steps - 1):
        qc = qc.compose(one_step)
    return qc


def edge_for(num_qubits):
    """Smallest default grid whose interior holds num_qubits data qubits."""
    return max(4, math.isqrt(num_qubits - 1) + 3)


def place_circuit(qc, grid):
    """Map a user circuit onto the grid: circuit qubit k becomes data qubit k+1 (row-major in the
    grid interior). Measurements and barriers are dropped and unsupported gates are translated
    into SUPPORTED_GATES. Returns a circuit with one qubit per grid cell."""
    layout = find_layout(grid, sum(1 for row in grid for v in row if isinstance(v, int) and v > 0))
    if qc.num_qubits > len(layout):
        raise ValueError(f"the circuit has {qc.num_qubits} qubits but the grid holds only {len(layout)} data qubits")
    qc = qc.remove_final_measurements(inplace=False)
    body = QuantumCircuit(qc.num_qubits)
    for inst in qc.data:
        if inst.operation.name in ("barrier", "measure"):
            continue
        if inst.clbits or inst.operation.name == "reset":
            raise ValueError(f"'{inst.operation.name}' with classical bits or resets is not supported")
        body.append(inst.operation, [qc.find_bit(q).index for q in inst.qubits])
    if set(body.count_ops()) - set(SUPPORTED_GATES):
        body = transpile(body, basis_gates=SUPPORTED_GATES, optimization_level=0)
    placed = QuantumCircuit(len(grid) * len(grid[0]))
    placed.compose(body, qubits=layout[:qc.num_qubits], inplace=True)
    return placed


def compile_circuit(qc, grid, magic_state_factories=1, weights=None):
    """Route qc on grid. Returns (final grid, compiled dag). Raises LayoutInfeasible."""
    dag = add_global_barrier(circuit_to_dag(qc))  # global barrier marks the end of the input circuit
    opt_levels = {'level1': level_1, 'level2': level_2}
    weights = weights or {'level1': 0.1, 'level2': 0.9}
    grid, dag, _ = apply_sequence(dag, grid, opt_levels, magic_state_factories, weights)
    for node in list(dag.op_nodes()):
        if isinstance(node.op, (Measure, Barrier)):
            dag.remove_op_node(node)
    return grid, dag


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--program", choices=sorted(BENCHMARKS), default="ising")
    p.add_argument("--qasm", metavar="FILE", help="compile this OpenQASM 2 circuit instead of a benchmark")
    p.add_argument("--layout", metavar="FILE",
                   help="your own initial layout (see load_layout); replaces --edge and --add-paths")
    p.add_argument("--edge", type=int, default=None,
                   help="qubits on one edge of the grid; (edge-2)^2 data qubits (paper: 4..12). "
                        "Default 6, or the smallest grid that fits the --qasm circuit")
    p.add_argument("--steps", type=int, default=1, help="Trotter steps")
    p.add_argument("--angle", type=float, default=np.pi / 4,
                   help="rz angle of every benchmark interaction term, in radians (paper: pi/4)")
    p.add_argument("--rotation-epsilon", type=float, default=1e-4,
                   help="accuracy of the Clifford+T approximation of rotations that are not multiples of pi/4")
    p.add_argument("--factories", type=int, default=1, choices=[1, 2, 3, 4],
                   help="magic state factories (paper: 1)")
    p.add_argument("--add-paths", type=int, default=None, metavar="X",
                   help="insert routing paths 0..X with add_paths_upto")
    p.add_argument("--msf-time", type=float, nargs="+", default=[11],
                   help="magic state factory processing time(s), in the units of Postprocessing.get_move_time_total")
    p.add_argument("--iterations", type=int, default=1,
                   help="compile this many times (seeds seed..seed+iterations-1) and report the minimum")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--output", default="output.txt")
    args = p.parse_args()

    set_rotation_epsilon(args.rotation_epsilon)
    user_circuit = qasm2.load(args.qasm, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS) if args.qasm else None
    edge = args.edge or (edge_for(user_circuit.num_qubits) if user_circuit else 6)
    name = args.qasm or args.program

    best = None
    for it in range(args.iterations):
        seed = args.seed + it
        random.seed(seed)
        grid = load_layout(args.layout) if args.layout else build_grid(edge, args.add_paths)
        if user_circuit is not None:
            qc = place_circuit(user_circuit, grid)
        else:
            qc = build_circuit(args.program, grid, args.steps, args.angle)
        try:
            _, dag = compile_circuit(qc, grid, args.factories)
        except LayoutInfeasible as e:
            print(f"seed {seed}: layout infeasible: {e}")
            continue
        overheads, totals, _, n_magic = generate_time(dag, args.msf_time)
        print(f"seed {seed}: magic states {n_magic}, total time {totals}, overhead {overheads}")
        if best is None or totals[0] < best[1][0]:
            best = (seed, totals, overheads, n_magic)

    if best is None:
        raise SystemExit("no iteration compiled successfully")
    seed, totals, overheads, n_magic = best
    with open(args.output, "w") as f:
        print(f"program={name} layout={args.layout or f'edge {edge}'} steps={args.steps} angle={args.angle} "
              f"rotation_epsilon={args.rotation_epsilon} factories={args.factories} "
              f"add_paths={args.add_paths} best_seed={seed}", file=f)
        for t_msf, total, overhead in zip(args.msf_time, totals, overheads):
            # paper Eq. 2. Note: generate_time models a single factory, so with more factories
            # the total time does not fall towards this bound.
            print(f"t_MSF={t_msf} magic_states={n_magic} total_time={total} overhead={overhead} "
                  f"lower_bound={n_magic * t_msf / args.factories:g}", file=f)
    print(f"minimum over {args.iterations} iteration(s) written to {args.output}")


if __name__ == "__main__":
    main()
