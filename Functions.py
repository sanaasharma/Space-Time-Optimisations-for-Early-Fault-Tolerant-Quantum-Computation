 
from qiskit import QuantumCircuit
from qiskit.circuit import Instruction
from qiskit import QuantumCircuit, transpile
 
from qiskit.transpiler import CouplingMap
#from qiskit.providers.fake_provider.fake_backend import FakeBackend
import numpy as np
# is operations swap, how to figure out where to apply to the CNOT
from qiskit.circuit.library import Barrier
from qiskit.dagcircuit import DAGCircuit

from qiskit.passmanager import BasePassManager, GenericPass, ConditionalController
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.dagcircuit import DAGCircuit, DAGOpNode
from qiskit.circuit.library.standard_gates import SwapGate
import copy
from qiskit.circuit import Gate, QuantumCircuit

from qiskit.providers import BackendV2 as Backend
from qiskit.transpiler import Target
from qiskit.providers import Options
from qiskit.circuit import Parameter, Measure
from qiskit.circuit.library import PhaseGate, SXGate, UGate, CXGate, IGate, TGate, TdgGate, SGate, CCXGate, SXGate, SdgGate, HGate
from qiskit.transpiler.passes.synthesis import SolovayKitaev
from qiskit.quantum_info import Operator
from qiskit.synthesis import generate_basic_approximations
from heapq import heappush, heappop
import random
import copy
import math
from qiskit.converters import circuit_to_dag

from collections import defaultdict
                   
  


def find_index(grid, target):

     
    for i, row in enumerate(grid):
        for j, value in enumerate(row):
            if value == target:
                return (i, j)
    return None

def find_layout(grid,circuit_size):
    initial_layout=[]
    for x in range(1,circuit_size+1):
        i,j=find_index(grid,x)
        initial_layout.append(len(grid[0])*i+j)
    
    return initial_layout

def find_initial_layout_2D(n,circuit_size, bottomsideoff=bool):
    
    if bottomsideoff is not True:
        grid = [[0] * n for _ in range(n)]
        qubit=1
        while circuit_size>=qubit:
            for row in range(1,n-1):
                #if row != int(n/2):
                    for col in range(1,n-1):
                       # if col!=int(n/2):
                            grid[row][col]=qubit
                            qubit+=1
    else: 
        grid = [[0] * (n-1) for _ in range(n-1)]
        qubit=1
        while circuit_size>=qubit:
            for row in range(1,n-1):
                if row != int(n/2):
                    for col in range(1,n-1):
                        if col!=int(n/2):
                            grid[row][col]=qubit
                            qubit+=1
          

    return grid
    
def get_op_layers(dag):
    op_layers = []
    for layer in dag.multigraph_layers():
        ops_only = [node for node in layer if isinstance(node, DAGOpNode)]
        if ops_only:
            op_layers.append(ops_only)
    return op_layers[0]

def apply_sequence(dag, grid, opt_levels, divisor, weights=None):
    initial_grid = copy.deepcopy(grid)
    
    path_log=[]
    assert 'level1' in opt_levels, "'level1' must be present in opt_levels"

    all_levels = list(opt_levels.keys())

    counter=0
    if weights is None:
        level_weights = [1.0 for _ in all_levels]
    else:
        level_weights = [weights[lvl] for lvl in all_levels]

     
    dag, grid, counter,path_log = opt_levels['level1'](dag, initial_grid, grid,counter,divisor,path_log)
    selected_levels = ['level1']

    
    y = list(get_op_layers(dag))
 
    while y[0].op.name != "barrier":
        
        level = random.choices(all_levels, weights=level_weights, k=1)[0]

        new_dag, new_grid,counter,path_log = opt_levels[level](dag, initial_grid, grid,counter,divisor,path_log)
        if new_dag is not None:
            dag=new_dag
            grid=new_grid
            selected_levels.append(level)
            y = list(get_op_layers(dag))
        else: 
            y = list(get_op_layers(dag))
            continue

    return grid, dag, path_log

def choose_best_ancilla(grid, op_type, qindex, target=None):
     
    def get_neighbors(i, j):
        directions = [(0,1), (1,0), (0,-1), (-1,0)]
        return [(i+di, j+dj) for di, dj in directions]

    if op_type == "h":
        i, j = find_index(grid, qindex)

        ancilla_options = [
            (i + 1, j), (i, j + 1) ,
             (i - 1, j,) ,(i, j - 1) 
        ]

        for x,y in ancilla_options:
            if 0<=x<len(grid) and 0<=y<len(grid[0]):
                if grid[x][y]==0:
                    return [(x,y)]
        
        for x,y in ancilla_options:
            if 0<=x<len(grid) and 0<=y<len(grid[0]):
                if grid[x][y]!=-1:
                    for a,b in get_neighbors(x,y):
                        if 0<=a<len(grid) and 0<=b<len(grid[0]):
                            if grid[a][b]==0:
                                return [(x,y),(a,b)]
        
            
        return  None
         
            

        

    elif op_type in ("cx", "t"):
        
        if target is None:
            raise ValueError(f"Operation {op_type} requires a target qubit.")

        ci, cj = find_index(grid,qindex)
        ti, tj = find_index(grid, target)
        
        path,final_i,final_j=find_best_path(grid, ci, cj,ti, tj, False)

        return path,final_i,final_j

def add_global_barrier(dag):
     
    all_qubits = list(dag.qubits)
    barrier = Barrier(len(all_qubits))
    dag.apply_operation_back(barrier, qargs=all_qubits)
    return dag
  
def find_best_path(grid, start_i, start_j, target_i, target_j, S_presence=bool):
    all_neighbors = [
        (target_i - 1, target_j - 1),
        (target_i + 1, target_j - 1),
        (target_i - 1, target_j + 1),
        (target_i + 1, target_j + 1), (target_i  , target_j + 1), (target_i + 1, target_j ), (target_i - 1, target_j ), (target_i  , target_j - 1)
    ]
    if S_presence is False:
        neighbors=all_neighbors[:4]
    else:
        neighbors=all_neighbors[5:7]

    least_cost = float('inf')
    shortest_path = None
    final_i = final_j = None

    for i, j in neighbors:
        if 0 <= i < len(grid) and 0 <= j < len(grid[0]) and grid[i][j]!=-1:
            start = (start_i, start_j)
            end = (i, j)
            #print(start,end)
            path, cost = find_path(grid, start, end)

            if path:
                if cost < least_cost:
                    
                    least_cost = cost
                    shortest_path = path
                    final_i, final_j = i, j

                elif cost == least_cost:
                    if len(path) < len(shortest_path):
                        
                        shortest_path = path
                        final_i, final_j = i, j

                    elif len(path) == len(shortest_path):
                       
                        if grid[i][j] == 0:
                            shortest_path = path
                            
                            final_i, final_j = i, j

    if shortest_path is not None:
        return shortest_path, final_i, final_j
    else:
        return False


def find_path( grid, start, end,forbidden=None):
        if forbidden is not None:
            forbidden= forbidden

        def is_valid(i, j):
            return 0 <= i < len(grid) and 0 <= j < len(grid[0]) and grid[i][j]!=-1

        def manhattan(p1, p2):
            return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])

        def has_zero_neighbor(i, j):
            for di, dj in directions:
                ni, nj = i + di, j + dj
                if is_valid(ni, nj) and grid[ni][nj] == 0:
                    return True
            return False

        directions = [(0,1), (1,0), (0,-1), (-1,0)]
        
        visited = set()
        heap=[]
        initial_count = 0 if grid[start[0]][start[1]] in [0, 'm'] else 1
        initial_dist = manhattan(start, end)
        initial_cost = initial_dist * initial_count

        heappush(heap, (initial_cost, initial_dist, start, [start], initial_count))


        while heap:
            cost, dist_to_end, current, path, count = heappop(heap)
            

            if current == end:
                print(f"Final cost: {cost}, Non-zero count: {count}")
                return path,cost

            if current in visited:
                continue
            visited.add(current)
             
             
            for di, dj in directions:
                ni, nj = current[0] + di, current[1] + dj
                 
                neighbor = (ni, nj)
                apply_penalty = False
                if is_valid(ni, nj) and neighbor not in visited and (forbidden is None or (ni, nj) not in forbidden):

                    if not has_zero_neighbor(ni,nj):
                        apply_penalty=True

                    if grid[ni][nj] in [0, 'm']:
                        is_nonzero = 0  
                        step_count = count + is_nonzero
                    else:
                        is_nonzero = 2
                        step_count = 2*(count + is_nonzero) 

                    
                    dist_to_end = manhattan(neighbor, end)
                    if dist_to_end!=0:
                        new_cost = manhattan(start, (ni,nj))+ step_count * 2 if apply_penalty else step_count + manhattan(start, (ni,nj))

                        #new_cost=dist_to_end*penalty_count
                        heappush(heap, (new_cost, dist_to_end, neighbor, path + [neighbor], step_count))
                    else: 
                        path=path+[neighbor]
                        return(path,cost)
                        
        return None



def space_search(coords, grid, operation):
    i, j = coords
    rows, cols = len(grid), len(grid[0])

 
    Found = False
    max_distance=max(rows,cols)
    gate_ancilla = {
        'h': 1,
        'rz': 1,
        'sdg': 1,
        's': 1,
        'cx':1
    }
    dof = {}
    ancilla_required = gate_ancilla.get(operation, 1)
    
    while not Found:
        
        for dist in range(max_distance):
            
            valid_cells=[]
             
            for di in range(-dist, dist + 1):
                for dj in range(-dist, dist + 1):
                    if abs(di) + abs(dj) == dist:
                        ni, nj =  i + di,  j + dj
                        if 0 <= ni < rows and 0 <= nj < cols and grid[ni][nj]==0:
                            valid_cells.append((ni, nj))
                             
        
            if valid_cells:
                dof[dist] = valid_cells
         

            total_found = sum(len(v) for v in dof.values())
            if total_found >= ancilla_required:
                Found = True
                sorted_dof = dict(sorted(dof.items()))
                return sorted_dof 
 
def find_path_dp(grid, start, end):
    rows, cols = len(grid), len(grid[0])
    directions = [(0,1), (1,0), (0,-1), (-1,0)]

    def manhattan(p1, p2):
        return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])

    def is_valid(i, j):
        return 0 <= i < rows and 0 <= j < cols

    def has_zero_neighbor(i, j):
        for di, dj in directions:
            ni, nj = i + di, j + dj
            if is_valid(ni, nj) and grid[ni][nj] == 0:
                return True
        return False

    dp = [[math.inf] * cols for _ in range(rows)]
    parent = [[None] * cols for _ in range(rows)]

    ei, ej = end
    initial_penalty = 0 if grid[ei][ej] == 0 else 2
    dp[ei][ej] = initial_penalty
    COST_LIMIT = 1000
    heap = []
    heappush(heap, (initial_penalty, ei, ej))

    while heap:
        
        curr_cost, i, j = heappop(heap)
        if curr_cost > COST_LIMIT:
            continue

        apply_penalty = False
        # Skip if we've already found a better path
        if curr_cost > dp[i][j]:
            continue

        for di, dj in directions:
            ni, nj = i + di, j + dj
            if not is_valid(ni, nj):
                continue

            if not has_zero_neighbor(ni,nj):
                apply_penalty=True
            
            if grid[ni][nj] in [0, 'm']:
                is_nonzero = 0  
                step_count = curr_cost + is_nonzero
            else:
                is_nonzero = 2
                step_count = 2*(curr_cost + is_nonzero) 
                    
            step_cost = manhattan(start, (ni,nj))+ step_count * 2 if apply_penalty else step_count + manhattan(start, (ni,nj))
                        
            if step_cost < dp[ni][nj]:
                dp[ni][nj] = step_cost
                parent[ni][nj] = (i, j)
                heappush(heap, (step_cost, ni, nj))

    # If no path found
    si, sj = start
    if dp[si][sj] == math.inf:
        return [], -1

    # Reconstruct path from start to end
    path = []
    i, j = si, sj
    while (i, j) != end:
        path.append((i, j))
        i, j = parent[i][j]
    path.append(end)

    return path, dp[si][sj] 

def find_shortest_path(coords,end,grid,true_coords):
        if true_coords is not None:
            forbidden_coords= true_coords

        i,j=coords
        ei,ej=end
 
        di = ei-i
        dj = ej-j

        new_path = [(i,j)]

         
        if abs(di) !=0:
            step = 1 if di > 0 else -1
            for m in range(i + step, ei + step, step):
                if grid[m][ej]!=-1 and (true_coords is None or (m,ej) not in forbidden_coords):
                    new_path.append((m, ej))

        if abs(dj)!=0:  
            step = 1 if dj > 0 else -1
            for k in range(j + step, ej + step, step):
                if grid[ei][k]!=-1 and (true_coords is None or (ei,k) not in forbidden_coords) :
                    new_path.append((ei, k))

        return new_path
 
def _insert_single_path(grid, t: int):
  
    n_rows, n_cols = len(grid), len(grid[0])
    k = t // 2 + 2 + t//2   

    if t % 2 == 0:
         
        insert_col = k  
         
        new_grid = []
        for r in range(n_rows):
            row = grid[r]
            new_grid.append(row[:insert_col] + [0] + row[insert_col:])
        return new_grid
    else:
        
        insert_row = k  
        
        zero_row = [0] * n_cols
        return grid[:insert_row] + [zero_row] + grid[insert_row:]


def add_paths_upto(grid, x: int):
   
    out = [row[:] for row in grid]
    for t in range(x + 1):
        out = _insert_single_path(out, t)
    return out

 
 


