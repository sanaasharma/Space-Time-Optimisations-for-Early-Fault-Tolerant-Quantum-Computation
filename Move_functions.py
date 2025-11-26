from qiskit.circuit import Instruction, Gate
from Functions import find_path , find_shortest_path


def perform_move(final_i, final_j, x, y, dag, grid, path=None):
    start=(x, y)
    final= (final_i, final_j)
    def is_nearest_neighbor(c1, c2):
            return abs(c1[0] - c2[0]) + abs(c1[1] - c2[1]) == 1
    
    if path is None:
            path, _ = find_path(grid, start, final) \
                if not is_nearest_neighbor(start, final) \
                else ([start, final], 1)
   
    if not all(grid[x][y]==0 for x,y in path[1:]):
        dag,grid= move_extended(path[1:], grid,dag, [(x,y),(final_i,final_j)]) #first key doesn't need to move
        
    if all(grid[x][y]==0 for x,y in path[1:] ):
        start_pos = x * len(grid[0]) + y
        end_pos = final_i * len(grid[0]) + final_j
        move_gate = MoveGate(
                start_pos, end_pos, 2,
                operation_time=0,
                label=f"Move {start_pos} → {end_pos}"
            )
        dag, check, _, grid = move_gate.update_dag_and_grid(
                dag, grid, start_pos, end_pos, path #adding the entire path here because the movegate also removes the start index
            )

     
        return dag, grid, final_i, final_j
    else:
        return dag,grid,None, None
        
        
        

def move_ancilla(control_coords,target_coords,grid,dag):

    
    ci,cj=control_coords
    ti,tj=target_coords
    ancilla_options = [(ti,cj)]

    def get_neighbors(i, j):
        directions = [(0,1), (1,0), (0,-1), (-1,0)]
        return [(i+di, j+dj) for di, dj in directions if 0 <= i+di < len(grid) and 0 <= j+dj < len(grid[0]) and grid[i+di][j+dj]==0  ]
    
    for option in ancilla_options:
        x,y =option
        if grid[x][y]!=-1:
                for ni,nj in get_neighbors(x,y):
                    dag,grid,_,_ = perform_move(ni, nj, x, y, dag, grid)
                    
                    return grid,dag
    
    return None , None

def build_ancilla(dof,grid,coords, operation,dag, target_coords=None):

                
    def manhattan(p1, p2):
                return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])

    def find_best_path_manh(grid, start_i, start_j, target_i, target_j, operation):
            all_neighbors = [
                (target_i - 1, target_j - 1),
                (target_i + 1, target_j - 1),
                (target_i - 1, target_j + 1),
                (target_i + 1, target_j + 1), (target_i  , target_j + 1), (target_i + 1, target_j ), (target_i - 1, target_j ), (target_i  , target_j - 1)
            ]
            
            if operation=="h":
                neighbors= all_neighbors[4:]
            elif operation in {"s","sdg","rz"}:
                neighbors=all_neighbors[5:7]
            else:
                neighbors=all_neighbors[:4]

            max = float('inf')
            shortest_path = None
            final_i = final_j = None

            for i, j in neighbors:
                if 0 <= i < len(grid) and 0 <= j < len(grid[0]):
                    if grid[i][j]!=-1:
                        start = (start_i, start_j)
                        end = (i, j)
                        
                        path = find_shortest_path(  start, end,grid, target_coords)
                        
                        if path:
                            
                            if len(path)< max :
                                shortest_path = path
                                final_i, final_j = i, j

            if shortest_path is not None:
                return shortest_path, final_i, final_j
            else:
                return None, None, None

    def move_along_path_reversed(path, grid, dag):
        
        if path is not None:
            if all(grid[x][y]!=-1 for (x,y) in path):
                num_cols = len(grid[0])

                
                one_move=[path[0]]
                 
                for i in range(len(path)-1):
                        
                        j=i+1
                        src_i, src_j = path[j] 
                        if grid[src_i][src_j]!=0 :
                            dst_i, dst_j = one_move[0]
                            if grid[dst_i][dst_j]==0:

                                qindex = src_i * num_cols + src_j
                                new_position = dst_i * num_cols + dst_j

                                move_gate = MoveGate(
                                    qindex, new_position,2,
                                    operation_time=0,
                                    label=f"Move {qindex} → {new_position}"
                                )
                                dag, operation_time, _, grid = move_gate.update_dag_and_grid(
                                    dag, grid, qindex, new_position, [(src_i, src_j),(dst_i, dst_j)]
                                )
                                one_move=[path[j]]
                                
                        else: 
                            
                            one_move.append((src_i,src_j))
                        
            return dag, grid

        
        return None, None

    def get_neighbors(i, j):
        directions = [(0,1), (1,0), (0,-1), (-1,0), (1,1),(1,-1),(-1,1),(-1,-1)]
    
        return [(i+di, j+dj) for di, dj in directions ]
    
    
    true_i,true_j=coords
    first_key = next(iter(dof))
    final_i,final_j=dof[first_key][0]

    if str(first_key) == "1":
        
        if operation == "cx":
            if manhattan(coords,(final_i,final_j))!=2:
                for cell in get_neighbors(true_i,true_j):
                    x,y=cell
                    if 0<=x<len(grid) and 0<=y<len(grid[0]) and manhattan(cell,(final_i,final_j))==1:
                        dag,grid,final_i,final_j= perform_move(final_i,final_j,x,y,dag,grid)
                        break

        elif operation in {"s", "sdg", "rz"}:
            if manhattan(coords,(final_i,final_j))!=1:
                dag,grid,final_i,final_j= perform_move(final_i,final_j,final_i,true_j,dag,grid)

        else:
            if manhattan(coords,(final_i,final_j))!=1:
                for x,y in get_neighbors(true_i,true_j):
                    if 0<=x<len(grid) and 0<=y<len(grid[0]) and  manhattan((x,y),(final_i,final_j))==1:
                            dag,grid,final_i,final_j= perform_move(final_i,final_j,x,y,dag,grid)
                            break
                
    total_found = sum(len(v) for v in dof.values())
    if total_found > 1:
                all_coords = [key_coord for coords in dof.values() for key_coord in coords]
                sorted_keys = sorted(all_coords, key=lambda key_coord: manhattan(key_coord, (true_i, true_j)))
                top_keys = [sorted_keys[0]]
    else: 
                top_keys=[key_coord for coords in dof.values() for key_coord in coords]
         
    if operation in {"h","s","sdg","rz"}:
        
        start_i,start_j=top_keys[0]
        path, final_i, final_j = find_best_path_manh(grid, start_i, start_j, true_i, true_j, operation)
        final_ancl, final_keys= (final_i,final_j) ,top_keys 

    elif operation=="cx": 
        ti,tj=target_coords
        final_ancl, final_keys= (ti,true_j) ,top_keys 
        
 
        
    
    for i in range(len(final_ancl)):
           
            start = final_keys 
            end=final_ancl 
            a,b= end
            
             
            if grid[a][b]!=0:
                if target_coords is not None:
                    forbidden=[coords,target_coords]
                    if find_path( grid, start[0], end,forbidden ) is not None:
                        path, cost=  find_path( grid, start[0], end,forbidden )
                        if len(path)>1:
                            dag,grid=move_along_path_reversed(path, grid, dag)
                         
                    else:
                        None, None
                else:
                    path, cost=  find_path( grid, start[0], end  )

                    if len(path)>1:
                        dag,grid=move_along_path_reversed(path, grid, dag)
                       


    return dag, grid

def move_ancilla_for_S( coords,grid,dag):

    ci,cj= coords
    directions = [(-1,0), (1,0)]
  
    def get_neighbors(i, j):
        directions = [(0,1), (1,0), (0,-1), (-1,0)]
        return [(i+di, j+dj) for di, dj in directions if 0 <= i+di < len(grid) and 0 <= j+dj < len(grid[0]) and grid[i+di][j+dj]==0  ]
    
    for (di,dj) in directions :
        x,y =ci+di, cj+dj
        if 0 <= x < len(grid) and 0 <= y < len(grid[0]):
            if grid[x][y]!=-1:
                    for ni,nj in get_neighbors(x,y):
                        dag, grid, _, _= perform_move(ni, nj, x, y, dag, grid)
                        return grid,dag
                    
                
            
    return None , None



def cnot_move(path, final_i,final_j,grid,qindex,dag):
         
    def zero_count(path, grid):
        if all(grid[x][y]==0 for x,y in path):
                return True
        else:
            False  # No zero found in the path

 
    if zero_count(path,grid)==False:
        return None, None
    else: 
        new_position= len(grid[0])*final_i+final_j
        custom_move = MoveGate(qindex, new_position,len(path)+1,operation_time=0,  label=f"Move {qindex} → {new_position}")
        dag, _, _, grid = custom_move.update_dag_and_grid(dag, grid, qindex, new_position,path)
        
        return dag,grid
    
  
 
  
 
def move_along_path_reversed(path, grid, dag):
    
    if path is not None :
            if all(grid[x][y]!=-1 for x,y in path):
                for i in range(len(path) - 2, -1, -1):
                    src_i, src_j = path[i]
                    dst_i, dst_j = path[i + 1]
                    if grid[src_i][src_j]!=0:
                       
                        qindex = src_i * len(grid[0]) + src_j
                        new_position = dst_i * len(grid[0]) + dst_j

                        move_gate = MoveGate(
                            qindex, new_position, 2,
                            operation_time=0,
                            label=f"Move {qindex} → {new_position}"
                        )
                        dag, operation_time, _, grid = move_gate.update_dag_and_grid(
                            dag, grid, qindex, new_position,[path[i],path[i + 1]]
                        )
                            
                return dag, grid
            else:
                None, None
    return None, None

 

def move_extended(path, grid,dag, true_coords=None):


    def space_search(coords, grid,path, true_coords):
        if true_coords is not None:
            forbidden_coords= true_coords
        i, j = coords
        rows, cols = len(grid), len(grid[0])
        Found = False
        max_distance=max(rows,cols)
         
       
        while not Found:
            
            for dist in range(max_distance):
                
                 
                # Iterate over all (i, j) such that Manhattan distance from center is exactly 'dist'
                for di in range(-dist, dist + 1):
                    for dj in range(-dist, dist + 1):
                        if abs(di) + abs(dj) == dist:
                            ni, nj =  i + di,  j + dj
                            if 0 <= ni < rows and 0 <= nj < cols and grid[ni][nj]==0 and (ni,nj) not in path and (true_coords is None or (ni,nj) not in forbidden_coords):
                                return (ni,nj)
                                 
    

    if all(grid[x][y]!=-1 for x,y in path):
        for pos in path:
            (x,y)=pos
            if grid[x][y]!=0 :
               
                end=space_search((x,y), grid,path, true_coords=true_coords)
                cell_path= find_shortest_path( (x,y), end ,grid,  true_coords=true_coords)
                
                temp_dag,temp_grid=move_along_path_reversed(cell_path, grid,dag)
                if temp_dag is not None:
                    dag, grid=temp_dag, temp_grid
                    
        
        return dag,grid
 

class MagicMoveGate(Gate):
 
    def __init__(self, path, start_pos, end_pos, num_qubits, operation_time=None, label=None):
        name = "magicmove"
        super().__init__(name, num_qubits, [], label=label or f"{start_pos}->{end_pos}")
        self.path = path 
        self.num_qubits = num_qubits
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.operation_time = operation_time

    def _define(self):
        qc = QuantumCircuit(self.num_qubits)
        self.definition = qc

    def update_dag_and_grid(self, grid, dag, path, start, end):
        num_cols = len(grid[0])
        operation_time = 0
        swap_sequence = []
        start_i, start_j = start
        final_i, final_j = end
        old_qubit = start_i * num_cols + start_j
        new_qubit = final_i * num_cols + final_j

        for (old_i, old_j), (new_i, new_j) in zip(path, path):
            old_q = old_i * num_cols + old_j
            new_q = new_i * num_cols + new_j
            swap_sequence.append((old_q, new_q))
            operation_time += 1

        if len(path) <= 1:
            return dag, operation_time, swap_sequence, grid

        if all(grid[x][y] == 0 for x, y in path[1:]):
            qarg_indices = [i * num_cols + j for i, j in path]
            move_gate = MagicMoveGate(path, start, end, len(path), operation_time=operation_time,label=f"MagicMove {old_qubit } → {new_qubit }")
            dag.apply_operation_back(move_gate, qargs=[dag.qubits[q] for q in qarg_indices])
            final_i, final_j = end
            
            grid[final_i][final_j] = grid[start_i][start_j]
            grid[start_i][start_j] = 0
            return dag, operation_time, swap_sequence, grid
        else:
            return dag, operation_time, swap_sequence, grid

class MoveGate(Gate): 
    def __init__(self, start_pos, end_pos, num_qubits, operation_time=None, label=None):
        name = "move"
        super().__init__(name, num_qubits, [], label=label or f"{start_pos}->{end_pos}")
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.operation_time = operation_time

    def _define(self):
        qc = QuantumCircuit(self.num_qubits)
        self.definition = qc

    def update_dag_and_grid(self, dag, grid, start_pos, final_pos, path=None):
        num_cols = len(grid[0])
        start_i, start_j = divmod(start_pos, num_cols)
        final_i, final_j = divmod(final_pos, num_cols)
 

        operation_time = 0
        swap_sequence = [
            (old_i * num_cols + old_j, new_i * num_cols + new_j)
            for (old_i, old_j), (new_i, new_j) in zip(path, path[1:])
        ]
        operation_time = len(swap_sequence)

        if start_pos == final_pos:
            return dag, operation_time, swap_sequence, grid

        
        if all(grid[x][y]==0 for x,y in path[1:] ):
            qarg_indices = [i * num_cols + j for i, j in path]
            move_gate = MoveGate(start_pos, final_pos,num_qubits=len(qarg_indices), operation_time=operation_time, label=f"Move {start_pos} → {final_pos}")
            dag.apply_operation_back(move_gate, qargs=[dag.qubits[q] for q in qarg_indices])
            grid[final_i][final_j] = grid[start_i][start_j]
            grid[start_i][start_j] = 0
            return dag, operation_time, swap_sequence, grid
        else: 
             return dag, operation_time, swap_sequence, grid



   
   
"""
class MoveGate(Gate): 
    def __init__(self, start_pos, end_pos, num_qubits, operation_time=None, label=None):
        name = "move"
        super().__init__(name, num_qubits, [], label=label or f"{start_pos}->{end_pos}")
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.operation_time = operation_time

    def _define(self):
        qc = QuantumCircuit(self.num_qubits)
        self.definition = qc

    def update_dag_and_grid(self, dag, grid, start_pos, final_pos, path=None):
        num_cols = len(grid[0])
        start_i, start_j = divmod(start_pos, num_cols)
        final_i, final_j = divmod(final_pos, num_cols)
        start = (start_i, start_j)
        final = (final_i, final_j)

        def is_nearest_neighbor(c1, c2):
            return abs(c1[0] - c2[0]) + abs(c1[1] - c2[1]) == 1

        if path is None:
            path, _ = find_path(grid, start, final) \
                if not is_nearest_neighbor(start, final) \
                else ([start, final], 1)

        operation_time = 0
        swap_sequence = [
            (old_i * num_cols + old_j, new_i * num_cols + new_j)
            for (old_i, old_j), (new_i, new_j) in zip(path, path)
        ]
        operation_time = len(swap_sequence)

        if start_pos == final_pos:
            return dag, operation_time, swap_sequence, grid
        
        if not all(grid[x][y]==0 for x,y in path):
            return dag,None,None,grid

        
        if all(grid[x][y]==0 for x,y in path ):
            qarg_indices = [i * num_cols + j for i, j in path]
            move_gate = MoveGate(start_pos, final_pos,num_qubits=len(qarg_indices)+1, operation_time=operation_time, label=f"Move {start_pos} → {final_pos}")
            dag.apply_operation_back(move_gate, qargs=[dag.qubits[q] for q in qarg_indices])
            grid[final_i][final_j] = grid[start_i][start_j]
            grid[start_i][start_j] = 0
            return dag, operation_time, swap_sequence, grid
        else: 
             return dag, operation_time, swap_sequence, grid


"""
         