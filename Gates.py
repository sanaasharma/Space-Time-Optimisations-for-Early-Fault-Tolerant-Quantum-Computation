from qiskit.circuit.library import PhaseGate, SXGate, UGate, CXGate, IGate, TGate, TdgGate, SGate, CCXGate, SXGate, SdgGate
from qiskit.circuit import Instruction, Gate
from qiskit import QuantumCircuit, transpile, QuantumRegister
from Functions import find_index, find_shortest_path, find_path
import random  
 
 
  
 

def check_H(coords, grid):
    
    i, j = coords
    adjacents = [(i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)]
    options = []

    for x, y in adjacents:
        if 0 <= x < len(grid) and 0 <= y < len(grid[0]) and grid[x][y] == 0:
            options.append((x, y))
            return options   

    return False

 

def check_S(coords, grid):
    
    i, j = coords
    adjacents = [(i + 1, j), (i - 1, j)]
    options = []

    for x, y in adjacents:
        if 0 <= x < len(grid) and 0 <= y < len(grid[0]) and grid[x][y] == 0:
            options.append((x, y))
            return options  # Return immediately on finding the first available

    return False

 
def check_CT(qubit, control_coords, target_coords, grid, dag):
     
    from Move_functions import move_along_path_reversed

    qi,qj = find_index(grid, qubit)
    ci,cj = control_coords
    ti,tj = target_coords

    # Move qubit to control if needed
    if (qi, qj) != (ci, cj):
        
        path = find_shortest_path((qi, qj), (ci, cj), grid, true_coords=target_coords)
        temp_dag, temp_grid = move_along_path_reversed(path, grid, dag)
        if temp_dag is not None:
            dag, grid = temp_dag, temp_grid
    
    if grid[ti][tj]=='m':
        grid[ti][tj]=0
        return True, dag,grid
 
    return False, dag, grid

def check_CNOT(control_coords, target_coords, grid, dag):
    from Move_functions import move_ancilla
    
    ci, cj = control_coords
    ti, tj = target_coords
    ancilla_options = [(ti, cj)]
    potential_cells = []
    valid_cells = []

    for x, y in ancilla_options:
        if 0 <= x < len(grid) and 0 <= y < len(grid[0]) and grid[x][y] != -1:
            if abs(ci - ti) == 1 and abs(cj - tj) == 1:
                potential_cells.append((x, y))
                if grid[x][y] == 0:
                    valid_cells.append((x, y))

    if valid_cells:
        return random.choice(valid_cells)
    elif potential_cells:
        grid, dag = move_ancilla(control_coords, target_coords, grid, dag)
        if grid is not None:
            return check_CNOT(control_coords, target_coords, grid, dag)
        else:
            return False
    else:
        return False

 

def check_CNOT_new(control_coords, target_coords, grid, dag):
    from  Move_functions import move_extended
    ci,cj = control_coords
    ti,tj = target_coords
    ancilla_options = [(ti, cj)]
    potential_cells = []
    valid_cells = []

    for x, y in ancilla_options:
        if abs(ci - ti) == 1 and abs(cj - tj) == 1 and 0 <= x < len(grid) and 0 <= y < len(grid[0]) and grid[x][y] != -1:
            potential_cells.append((x, y))
            if grid[x][y] == 0:
                valid_cells.append((x, y))

    if valid_cells:
        return random.choice(valid_cells)
    elif potential_cells:
        true_coords = [control_coords, target_coords]
        for cell in potential_cells:
            dag, grid = move_extended([cell], grid, dag, true_coords)
            if grid is not None:
                return check_CNOT(control_coords, target_coords, grid, dag)
        return False
    else:
        return False


def check_S_new(coords, grid, dag, depth=0, max_depth=1000):
    from  Move_functions import move_extended
    i, j = coords
    ancilla_options = [(i + 1, j), (i - 1, j)]

    if depth > max_depth:
        return False

    potential_cells = []
    valid_cells = []

    for x, y in ancilla_options:
        if 0 <= x < len(grid) and 0 <= y < len(grid[0]) and grid[x][y] != -1:
            potential_cells.append((x, y))
            if grid[x][y] == 0:
                valid_cells.append((x, y))

    if len(valid_cells) > 1:
        return random.choice(valid_cells)
    elif len(valid_cells) == 1:
        return valid_cells[0]
    elif potential_cells:
        forbidden_cells = [coords]  
        for cell in potential_cells:
            if grid[cell[0]][cell[1]] != 0:
                new_dag, new_grid = move_extended([cell], grid, dag, forbidden_cells)
                if new_grid is not None:
                    return check_S_new(coords, new_grid, new_dag, depth + 1, max_depth)
        return False
    else:
        return False


def apply_S(grid, coords, ancilla, dag, node, required_node_removal=None):
    
    q_indices = [coords, ancilla]
    nodes = [i * len(grid[0]) + j for (i, j) in q_indices]
    qargs = [dag.qubits[nodes[0]], dag.qubits[nodes[1]]]
    custom_s = CustomSGate()
    dag.apply_operation_back(custom_s, qargs=qargs)
    if required_node_removal is None:
        dag.remove_op_node(node)

    return dag
   
def apply_H(ancilla,qindex,grid,dag,node, remove_rotation=None):
    
    ancilla_indices = [i * len(grid[0]) + j for (i, j) in ancilla] 
    qargs = [dag.qubits[qindex],dag.qubits[ancilla_indices[0]]]
    custom_h=CustomHGate()
    dag.apply_operation_back(custom_h,qargs=qargs)
    if remove_rotation==None:
        dag.remove_op_node(node)
    
    return dag

def apply_X(qindex,dag,node):
 
    qargs = [dag.qubits[qindex]]
    custom_x=CustomXGate()
    dag.apply_operation_back(custom_x,qargs=qargs)
    dag.remove_op_node(node)

    return dag

def apply_Z(qindex,dag,node):
 
    qargs = [dag.qubits[qindex]]
    custom_z=CustomZGate()
    dag.apply_operation_back(custom_z,qargs=qargs)
    dag.remove_op_node(node)

    return dag
  
 
def apply_CT(node1, node2,grid,dag ):

    qindices=[node1]+[node2]
    nodes=[i * len(grid[0]) + j  for (i,j) in qindices]
    qargs = [dag.qubits[nodes[0]],dag.qubits[nodes[1]]]
    custom_t=CustomTGate()
    dag.apply_operation_back(custom_t,qargs=qargs)

    return dag
 
       
 
def apply_CNOT(ancilla,node_qargs,grid,dag,node):

    i,j=ancilla
    ancilla_indices =  i * len(grid[0]) + j  
    qargs = [dag.qubits[node_qargs[0]],dag.qubits[node_qargs[1]],dag.qubits[ancilla_indices]]
    custom_cx=CustomCXGate()
    dag.apply_operation_back(custom_cx,qargs=qargs)
    dag.remove_op_node(node)
     
    return dag


class CustomHGate(Gate):
    def __init__(self):
        super().__init__("ch", 2, [])

    def _define(self):
       
        q = QuantumRegister(2, "q")
        qc = QuantumCircuit(q, name=self.name)
        qc.h(q[0])
        qc.id(q[1])
        self.definition = qc
    
    

class CustomXGate(Gate):
    def __init__(self):
        super().__init__("customx", 1, [])
    
    def _define(self):
       
        qc = QuantumCircuit(1, name=self.name)
        self.definition = qc
    
class CustomZGate(Gate):
    def __init__(self):
        super().__init__("customz", 1, [])

    def _define(self):
       
        qc = QuantumCircuit(1, name=self.name)
        self.definition = qc
  
class CustomCXGate(Gate):
    def __init__(self):
        super().__init__("ccx",3,[])

    def _define(self):
       
        q = QuantumRegister(3, "q")
        qc = QuantumCircuit(q, name=self.name)
        qc.cx(q[0], q[1])  
        self.definition = qc

class CustomTGate(Gate):
    def __init__(self):
        super().__init__("ct",2,[])

    def _define(self):
       
        q = QuantumRegister(2, "q")
        qc = QuantumCircuit(q, name=self.name)
        qc.append(TGate(), [q[0]])
        self.definition = qc

class CustomSGate(Gate):
   
    def __init__(self):
         
        super().__init__("cs", 2, [])
        
    def _define(self):
        q = QuantumRegister(2, "q")
        qc = QuantumCircuit(q, name=self.name)
        qc.append(SGate(), [q[0]])
        self.definition = qc

