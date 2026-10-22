
from qiskit.circuit.library import Barrier
from qiskit.dagcircuit import DAGOpNode
from heapq import heappush, heappop
import random
import copy
import math


class LayoutInfeasible(Exception):
    """Raised when the grid does not have enough free (routing) space to execute the circuit."""


class RoutingError(AssertionError):
    """Internal error: the router tried to emit a physically impossible operation."""


def is_data(grid, cell):
    v = grid[cell[0]][cell[1]]
    return isinstance(v, int) and v > 0


def adjacent(a, b):
    return abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1


# gates that come from the input circuit; the compile is finished once none of them remain
INPUT_GATES = {"h", "cx", "rz", "s", "sdg", "t", "tdg", "x", "sx"}


def count_pending_gates(dag):
    return sum(1 for node in dag.op_nodes() if node.op.name in INPUT_GATES)


def factory_position(grid, number):
    """Output cell of magic state factory `number`, at the midpoint of a grid edge
    (left, right, top, bottom). Factory 0 is the one used in the paper."""
    rows, cols = len(grid), len(grid[0])
    positions = [(int(rows/2), 0), (int(rows/2), cols-1), (0, int(cols/2)), (rows-1, int(cols/2))]
    if not 0 <= number < len(positions):
        raise ValueError(f"at most {len(positions)} magic state factories are supported, got factory index {number}")
    return positions[number]


def validate_grid(grid, factories=1):
    """A layout is a rectangular list of rows: 0 = free routing cell, -1 = unusable cell, and the
    data qubits labelled 1..N, each exactly once. The output cells of the magic state factories
    (factory_position) must not be unusable."""
    if not grid or not grid[0] or any(len(row) != len(grid[0]) for row in grid):
        raise ValueError("the layout must be a non-empty rectangle: all rows the same length")
    labels = []
    for i, row in enumerate(grid):
        for j, v in enumerate(row):
            if not isinstance(v, int) or isinstance(v, bool) or v < -1:
                raise ValueError(f"cell {(i, j)} is {v!r}; use 0 (free), -1 (unusable) or a data qubit label 1..N")
            if v > 0:
                labels.append(v)
    if sorted(labels) != list(range(1, len(labels) + 1)):
        dup = sorted({v for v in labels if labels.count(v) > 1})
        raise ValueError(f"data qubit labels must be 1..N, each once (found {len(labels)} labels"
                         + (f", duplicated: {dup}" if dup else ", with gaps") + ")")
    for k in range(factories):
        i, j = factory_position(grid, k)
        if grid[i][j] == -1:
            raise ValueError(f"cell {(i, j)} is the output cell of magic state factory {k} and cannot be unusable (-1)")


def check_grid_matches_dag(grid, dag):
    if any(len(row) != len(grid[0]) for row in grid):
        raise ValueError("grid rows have different lengths")
    if dag.num_qubits() != len(grid) * len(grid[0]):
        raise ValueError(
            f"circuit has {dag.num_qubits()} qubits but the grid has {len(grid)}x{len(grid[0])} cells; "
            "build the circuit from the same grid that is passed to apply_sequence")


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

def find_initial_layout_2D(n,circuit_size, bottomsideoff=False):

    capacity = (n-2)**2 if not bottomsideoff else (n-3)**2
    if circuit_size > capacity:
        raise ValueError(f"{circuit_size} data qubits do not fit in the interior of a {n}x{n} grid (max {capacity})")

    if bottomsideoff is not True:
        grid = [[0] * n for _ in range(n)]
        qubit=1
        while circuit_size>=qubit:
            for row in range(1,n-1):
                    for col in range(1,n-1):
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

def _perturb(dag, grid, initial_grid):
    """Stall breaker: push one qubit of a waiting gate to the nearest reachable free cell
    (shifting other qubits along the way if needed), so that routing can continue from a
    different configuration."""
    from Move_functions import move_extended
    cols = len(grid[0])
    waiting = [n for n in get_op_layers(dag) if n.op.name != "barrier"]
    if not waiting:
        return
    node = random.choice(waiting)
    cells = [find_index(grid, initial_grid[q._index // cols][q._index % cols]) for q in node.qargs]
    cell = random.choice(cells)
    move_extended([cell], grid, dag, [c for c in cells if c != cell])


def apply_sequence(dag, grid, opt_levels, divisor, weights=None, max_stalled_rounds=40, perturb_every=5):
    """Route `dag` on `grid`. After every `perturb_every` rounds without progress one waiting
    qubit is moved (see _perturb). Raises LayoutInfeasible if `max_stalled_rounds` consecutive
    rounds fail to execute any gate of the input circuit."""
    validate_grid(grid, divisor)
    check_grid_matches_dag(grid, dag)
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


    y = list(get_op_layers(dag))
    pending = count_pending_gates(dag)
    stalled_rounds = 0

    while y[0].op.name != "barrier":

        level = random.choices(all_levels, weights=level_weights, k=1)[0]

        new_dag, new_grid,counter,path_log = opt_levels[level](dag, initial_grid, grid,counter,divisor,path_log)
        if new_dag is not None:
            dag=new_dag
            grid=new_grid
        y = list(get_op_layers(dag))

        now_pending = count_pending_gates(dag)
        if now_pending < pending:
            stalled_rounds = 0
        else:
            stalled_rounds += 1
            if stalled_rounds % perturb_every == 0:
                _perturb(dag, grid, initial_grid)
            if stalled_rounds >= max_stalled_rounds:
                raise LayoutInfeasible(
                    f"no gate could be executed in {max_stalled_rounds} consecutive rounds "
                    f"({now_pending} gates left): routing is stuck, the grid most likely has too little free space; "
                    "try more routing paths (add_paths_upto) or fewer blocked cells")
        pending = now_pending

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
  
def find_best_path(grid, start_i, start_j, target_i, target_j, S_presence=False, forbidden=None):
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
            found = find_path(grid, start, end, forbidden)
            if found is None:
                continue
            path, cost = found

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
    raise LayoutInfeasible(
        f"no free path from {(start_i, start_j)} to any usable neighbour of {(target_i, target_j)}")


def find_path( grid, start, end,forbidden=None):
        """Returns (path, cost), or None if `end` cannot be reached."""
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
    """Free cells grouped by Manhattan distance from `coords`, nearest first.
    Raises LayoutInfeasible if the grid has fewer free cells than the operation needs."""
    i, j = coords
    rows, cols = len(grid), len(grid[0])

    max_distance=rows+cols-1  # covers every cell of the grid
    gate_ancilla = {
        'h': 1,
        'rz': 1,
        'sdg': 1,
        's': 1,
        'cx':1
    }
    dof = {}
    ancilla_required = gate_ancilla.get(operation, 1)

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
            sorted_dof = dict(sorted(dof.items()))
            return sorted_dof

    raise LayoutInfeasible(f"no free cell available for '{operation}' at {coords}")

def find_shortest_path(coords,end,grid,true_coords):
    """Contiguous path of adjacent cells from `coords` to `end` that avoids blocked (-1) cells
    and the cells in `true_coords`. Occupied cells may lie on it (move_along_path_reversed
    shifts their qubits along). Prefers an L-shaped path; falls back to breadth-first search.
    Returns None if `end` cannot be reached."""
    forbidden = set(true_coords or ())
    rows, cols = len(grid), len(grid[0])

    def usable(c):
        return c in (coords, end) or (grid[c[0]][c[1]] != -1 and c not in forbidden)

    def l_path(vertical_first):
        (i, j), (ei, ej) = coords, end
        path, (ci, cj) = [(i, j)], (i, j)
        for axis in ((0, 1) if vertical_first else (1, 0)):
            while (ci, cj)[axis] != end[axis]:
                if axis == 0:
                    ci += 1 if ei > ci else -1
                else:
                    cj += 1 if ej > cj else -1
                path.append((ci, cj))
        return path

    for vertical_first in (True, False):
        path = l_path(vertical_first)
        if all(usable(c) for c in path):
            return path

    parent = {coords: None}
    queue = [coords]
    for c in queue:
        if c == end:
            path = []
            while c is not None:
                path.append(c)
                c = parent[c]
            return path[::-1]
        for di, dj in ((0, 1), (1, 0), (0, -1), (-1, 0)):
            n = (c[0] + di, c[1] + dj)
            if 0 <= n[0] < rows and 0 <= n[1] < cols and n not in parent and usable(n):
                parent[n] = c
                queue.append(n)
    return None
 
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

 
 


