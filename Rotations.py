from Functions import find_index, find_best_path, find_path, space_search, LayoutInfeasible, adjacent, factory_position
from Gates import check_H, apply_Z, apply_X, apply_H, check_S_new, check_CT, apply_CT, apply_S
from Move_functions import build_ancilla, move_extended, MagicMoveGate, free_a_neighbour, move_to_free_cell
from qiskit.converters import circuit_to_dag
from qiskit import QuantumCircuit
from functools import lru_cache
import math


MAX_ROTATION_ATTEMPTS = 5

# Accuracy of the Clifford+T approximation of rotations that are not multiples of pi/4
# (operator-norm error, up to global phase). Set with set_rotation_epsilon().
ROTATION_EPSILON = 1e-4

# rz(k*pi/4) up to global phase, k = 0..7: exact, no approximation needed
_EXACT_RZ = {0: [], 1: ["t"], 2: ["s"], 3: ["s", "t"], 4: ["z"], 5: ["z", "t"], 6: ["sdg"], 7: ["tdg"]}
_GRIDSYNTH_GATES = {"H": "h", "S": "s", "T": "t", "X": "x"}  # gridsynth also emits W/I (global phase)


def set_rotation_epsilon(epsilon):
    global ROTATION_EPSILON
    ROTATION_EPSILON = float(epsilon)


@lru_cache(maxsize=None)
def decompose_rz(angle, epsilon):
    """Clifford+T gate names (circuit order) for rz(angle). Multiples of pi/4 are exact; any other
    angle is approximated to `epsilon` with gridsynth (Ross & Selinger), which uses about
    3*log2(1/epsilon) T gates."""
    k = angle / (math.pi / 4)
    if abs(k - round(k)) < 1e-9:
        return tuple(_EXACT_RZ[round(k) % 8])
    try:
        from pygridsynth.gridsynth import gridsynth_gates
    except ImportError as e:
        raise ImportError("rotations that are not multiples of pi/4 need pygridsynth "
                          "(pip install -r requirements.txt)") from e
    word = gridsynth_gates(theta=repr(float(angle)), epsilon=repr(float(epsilon)), seed=0, up_to_phase=True)
    # the word is an operator product (rightmost gate acts first), so reverse it for circuit order
    return tuple(_GRIDSYNTH_GATES[g] for g in reversed(word) if g in _GRIDSYNTH_GATES)


def rotation_circuit_for(node_input, sx=False):
    if sx:
        names = ("sdg", "h", "sdg")  # sx up to global phase
    else:
        names = decompose_rz(float(node_input.op.params[0]), ROTATION_EPSILON)
    circuit = QuantumCircuit(1)
    for name in names:
        getattr(circuit, name)(0)
    return circuit_to_dag(circuit)


def _clear_magic_path(grid, dag, start, target, S_presence, path_log):
    """Plan a magic state path from the factory cell `start` to a usable neighbour of `target`
    and clear it. Re-plans around cells that cannot be cleared. Returns (dag, grid, path, end);
    path is None if no clear path was found."""
    forbidden = [target]  # the path must not run through the qubit that will consume the state
    for _attempt in range(4):
        try:
            path,f_i, f_j = find_best_path(grid,start[0],start[1],target[0],target[1],S_presence, forbidden=forbidden)
        except LayoutInfeasible:
            return dag, grid, None, None
        end=f_i,f_j
        found = find_path(grid,start,end,forbidden)
        if found is None:
            return dag, grid, None, None
        path,_ = found
        path_log.append(path)
        # clear the path (not the factory cell itself, which may already hold a waiting magic state);
        # the cleared qubits must not be pushed through the consuming qubit either
        dag,grid= move_extended(path[1:] if grid[start[0]][start[1]] == 'm' else path, grid,dag,path + [target])
        blocked = [c for c in path[1:] if grid[c[0]][c[1]] != 0]
        if not blocked:
            return dag, grid, path, end
        forbidden = forbidden + blocked  # re-plan around the cells that could not be cleared
    return dag, grid, None, None


def _relocate_qubit(grid, dag, target):
    """Last resort: move the qubit at `target` through free cells to the nearest free cell that has
    a free cell above or below it, where a magic state can be delivered. Returns True if moved."""
    rows = len(grid)
    def has_free_vertical_neighbour(c):
        return any(0 <= r < rows and grid[r][c[1]] == 0 for r in (c[0] - 1, c[0] + 1))
    return move_to_free_cell(grid, dag, target, has_free_vertical_neighbour)


def magic_state_routing_first_pass(grid,qubit,dag, number,path_log, S_presence=True):

    start = factory_position(grid, number)
    target = find_index(grid,qubit)

    dag, grid, path, end = _clear_magic_path(grid, dag, start, target, S_presence, path_log)
    if path is None and _relocate_qubit(grid, dag, target):
        target = find_index(grid,qubit)
        dag, grid, path, end = _clear_magic_path(grid, dag, start, target, S_presence, path_log)
    if path is None:
        raise LayoutInfeasible(f"no clear path for a magic state from factory {start} to qubit {qubit} at {target}")
    if grid[start[0]][start[1]] not in (0, 'm'):
        return None,None,None, path_log  # the factory output cell is still occupied

    grid[start[0]][start[1]]='m'  # a new magic state, or the one still waiting from a failed delivery
    custom_magicmove= MagicMoveGate(path, start,end,len(path),operation_time=0)
    dag, _, _, grid = custom_magicmove.update_dag_and_grid( grid,dag,path, start, end)
    path_log.append(grid)
    return grid, dag , end,path_log
 
def convert_index(grid,initial_grid,node):
        node_qargs = [q._index for q in node.qargs]
        i,j= divmod(node_qargs[0],len(grid[0]))
        true_i,true_j=find_index(grid,initial_grid[i][j])
        qindex=len(grid[0])*true_i+true_j
        return qindex , true_i, true_j
  
 

def rotation_consumption(qubit,dag,initial_grid,grid,node_input,number,path_log, SXGATE=None, executed_nodes=None):
    """Execute rz / sx gate `node_input` as a sequence of Clifford+T operations.
    `qubit` is ignored; the qubit is always taken from node_input."""

    i, j = divmod(node_input.qargs[0]._index, len(grid[0]))
    qubit = initial_grid[i][j]

    if executed_nodes is None:
        executed_nodes = set()

    rotation_circuit = rotation_circuit_for(node_input, sx=SXGATE is not None)
    success=False
     

    def rotation_check(rotation_circuit,dag,grid,executed_nodes,qubit,initial_grid,node_input,path_log):
       
        for node in list(rotation_circuit.op_nodes()):
            
            
            if node in executed_nodes:
                continue

            try:

                if node.op.name in ('t', 'tdg'): 

                    
                    qindex,ci,cj=convert_index(grid,initial_grid,node_input)
                    grid,dag,final_index,path_log = magic_state_routing_first_pass(grid,qubit,dag,number,path_log, S_presence=True)
                    if grid is not None:
                         
                        qindex,ci,cj=convert_index(grid,initial_grid,node_input)
                       
                        check, dag,grid=check_CT(qubit,( ci,cj),final_index,grid,dag)
                        
                        if check== True:
                            
                            qindex,ci,cj=convert_index(grid,initial_grid,node_input)
                            dag=apply_CT( ( ci,cj), final_index,grid, dag )
                             
                            executed_nodes.add(node)
                        else:  
                            return None,None,executed_nodes,path_log
                    else: 
                        return None,None,executed_nodes,path_log 
                    
                 
                if node.op.name in ('z', 'x'):  # Pauli frame update, no routing needed
                     
                    qindex,_,_=convert_index(grid,initial_grid,node_input)
                    dag=(apply_Z if node.op.name == 'z' else apply_X)(qindex,dag,node,remove=False)
                    executed_nodes.add(node)
                  

                if node.op.name=='h':
                    
                    qindex,true_i,true_j=convert_index(grid,initial_grid,node_input)
                    
                    if check_H((true_i,true_j),grid):
                            ancilla=check_H((true_i,true_j),grid)
                            dag=apply_H(ancilla, qindex,grid, dag, node,node_input)
                            executed_nodes.add(node)

                    else:
                        dof= space_search((true_i,true_j), grid, "h")
                        dag,grid= build_ancilla(dof,grid,(true_i,true_j),node.op.name,dag)
                        qindex,true_i,true_j=convert_index(grid,initial_grid,node_input)
                        if not check_H((true_i,true_j),grid):
                            dag,grid= free_a_neighbour((true_i,true_j), grid, dag)
                        if check_H((true_i,true_j),grid):
                            ancilla=check_H((true_i,true_j),grid)
                            dag=apply_H(ancilla, qindex,grid, dag, node,node_input)
                            executed_nodes.add(node)
                        else: 
                            return None,None,executed_nodes,path_log

                
                if node.op.name in ('s', 'sdg'):

                    qindex,true_i,true_j=convert_index(grid,initial_grid,node_input)
                    ancilla= check_S_new((true_i,true_j), grid,dag)
                    qindex,true_i,true_j=convert_index(grid,initial_grid,node_input)  # check_S_new can move qubits
                    if ancilla and adjacent((true_i,true_j), ancilla) and grid[ancilla[0]][ancilla[1]] == 0:
                        dag= apply_S(grid, (true_i,true_j),ancilla, dag, node,node_input)
                        executed_nodes.add(node)
                    else: 
                        dof= space_search((true_i,true_j), grid, "s")
                        dag,grid= build_ancilla(dof,grid,(true_i,true_j),node.op.name,dag)
                        qindex,true_i,true_j=convert_index(grid,initial_grid,node_input)  # build_ancilla can move the qubit
                        ancilla= check_S_new((true_i,true_j), grid,dag)
                        qindex,true_i,true_j=convert_index(grid,initial_grid,node_input)
                        if ancilla and adjacent((true_i,true_j), ancilla) and grid[ancilla[0]][ancilla[1]] == 0:
                            dag= apply_S(grid, (true_i,true_j), ancilla, dag, node,node_input)
                            executed_nodes.add(node)

                        else:
                            return None, None ,executed_nodes  ,path_log  

 
                 
            except LayoutInfeasible:
                return  None,None ,executed_nodes,path_log

        return dag,grid,executed_nodes,path_log


    for _attempt in range(MAX_ROTATION_ATTEMPTS):
        new_dag,_,executed_nodes,path_log= rotation_check(rotation_circuit,dag,grid,executed_nodes,qubit,initial_grid,node_input,path_log)
        if new_dag is not None and len(executed_nodes) == len(list(rotation_circuit.op_nodes())):
            dag.remove_op_node(node_input)
            success = True
            break

    if not success:
        if not executed_nodes:
            return None, None, path_log  # nothing applied yet: the caller can retry this rotation later
        # part of the decomposed rotation is already in the DAG, so it cannot be retried later
        raise LayoutInfeasible(
            f"could not route rotation on qubit {qubit} after {MAX_ROTATION_ATTEMPTS} attempts; "
            "the grid has too little free space")

    return dag,grid,path_log
