from Functions import get_op_layers,find_index,  choose_best_ancilla, space_search, LayoutInfeasible, is_data, adjacent
from Move_functions import  cnot_move, move_ancilla,move_ancilla_for_S, move_extended, build_ancilla, perform_move, free_a_neighbour, move_to_free_cell
from Gates import check_CNOT, apply_CNOT, apply_CT,check_S, apply_S, apply_X, check_H, apply_H, check_S_new
from Rotations import magic_state_routing_first_pass, rotation_consumption, check_CT, convert_index
import random


def cnot_is_legal(grid, c, t, ancilla):
    """Lattice-surgery CNOT: control and target diagonal neighbours, free ancilla at (target row, control column)."""
    return (is_data(grid, c) and is_data(grid, t) and abs(c[0] - t[0]) == 1 and abs(c[1] - t[1]) == 1
            and tuple(ancilla) == (t[0], c[1]) and grid[t[0]][c[1]] == 0)


def process_T(grid, initial_grid, node, dag, counter, divisor, path_log):
    """Shared by level_1 and level_2. Returns (dag, grid, counter, path_log)."""
    _,true_i,true_j =convert_index(grid,initial_grid,node )
    qubit_number=grid[true_i][true_j]
    counter=counter % divisor
    try:
        temp_grid,temp_dag,final_index,path_log = magic_state_routing_first_pass(grid,qubit_number,dag,counter,path_log, S_presence=True)
    except LayoutInfeasible:
        return dag, grid, counter, path_log
    if temp_grid is not None:
        dag,grid=temp_dag,temp_grid
        _,ci,cj =convert_index(grid,initial_grid,node )
        check, dag,grid=check_CT(qubit_number ,( ci,cj),final_index,grid,dag)
        if check is True:
            dag=apply_CT( ( ci,cj), final_index,grid, dag )
            dag.remove_op_node(node)
            counter+=1
            return dag,grid, counter, path_log

        else:
            temp_grid,temp_dag=move_ancilla(( ci,cj),final_index ,grid,dag)
            if temp_grid is not None:
                dag,grid=temp_dag,temp_grid
                check, dag,grid=check_CT(qubit_number,( ci,cj),final_index,grid,dag)
                if check:
                    dag=apply_CT(( ci,cj), final_index,grid, dag )
                    dag.remove_op_node(node)
                    counter+=1
                    return dag,grid, counter, path_log

    return dag,grid, counter, path_log



def level_1(dag,initial_grid,grid,counter ,divisor,path_log):

    
    def process_H(grid,initial_grid,node,dag):

        qindex,true_i,true_j=convert_index(grid,initial_grid,node )

        if check_H((true_i,true_j),grid):
            ancilla=check_H((true_i,true_j),grid)
            dag=apply_H(ancilla, qindex,grid, dag, node)
                        
        else:
                        
            suggested_moves= choose_best_ancilla(grid, "h", grid[true_i][true_j])        
            if suggested_moves:
                if len(suggested_moves)==1:
                    if check_H((true_i,true_j),grid):
                        ancilla=check_H((true_i,true_j),grid)
                        dag=apply_H(ancilla, qindex,grid, dag, node)
                                    
                else:
                    start_i,start_j=suggested_moves[0]
                    new_i,new_j= suggested_moves[1]
                    dag, grid, _, _= perform_move(new_i,new_j,start_i,start_j , dag, grid)
                    qindex,true_i,true_j=convert_index(grid,initial_grid,node )
                    if check_H((true_i,true_j),grid):
                        ancilla=check_H((true_i,true_j),grid)
                        dag=apply_H(ancilla, qindex,grid, dag, node)

        return dag,grid
        
    def process_CNOT(grid, initial_grid, node, dag):

        def level_3_for_CNOT(grid,dag,control,target,node):
            ci,cj= find_index(grid,control)
            ti,tj=find_index(grid,target)
            if grid[ti][cj] not in (0, -1):  # a blocked (-1) ancilla cell cannot be cleared
                dag,grid= move_extended([(ti,cj)],grid,dag,[(ci,cj),(ti,tj)])
                try:
                    path, final_i, final_j = choose_best_ancilla(grid, "cx", control, target)
                except LayoutInfeasible:
                    return dag, grid
                x, y = find_index(grid, control)
                temp_dag, temp_grid = cnot_move(path, final_i, final_j, grid, x*n_cols+y, dag)
                if temp_dag is not None:
                    
                        dag, grid = temp_dag, temp_grid
                        true_ci, true_cj = find_index(grid,control)
                        true_ti, true_tj = find_index(grid, target)
                        dag, grid, applied = try_apply_CNOT(dag, grid, true_ci, true_cj, true_ti, true_tj, node)
                        if applied:
                            return dag, grid
            
            return dag,grid

        
        node_qargs = [q._index for q in node.qargs]
        n_cols = len(grid[0])
        ci, cj = divmod(node_qargs[0], n_cols)
        ti, tj = divmod(node_qargs[1], n_cols)
        true_ci, true_cj = find_index(grid, initial_grid[ci][cj])
        true_ti, true_tj = find_index(grid, initial_grid[ti][tj])
        control,target=initial_grid[ci][cj],initial_grid[ti][tj]

        
        def try_apply_CNOT(dag, grid, ci, cj, ti, tj, node):
            ancilla = check_CNOT((ci, cj), (ti, tj), grid, dag)
            (ci, cj), (ti, tj) = find_index(grid, control), find_index(grid, target)  # the check can move qubits
            if ancilla and cnot_is_legal(grid, (ci, cj), (ti, tj), ancilla):
                new_node_qargs = [n_cols * ci + cj, n_cols * ti + tj]
                dag = apply_CNOT(ancilla, new_node_qargs, grid, dag, node)
                return dag, grid, True
            return dag, grid, False

        
        dag, grid, applied = try_apply_CNOT(dag, grid, true_ci, true_cj, true_ti, true_tj, node)
        if applied:
            return dag, grid

        # Try to move ancilla and retry
        for a,b in [(grid[true_ci][true_cj], grid[true_ti][true_tj]),
                                (grid[true_ti][true_tj], grid[true_ci][true_cj])]:

            try:
                path, final_i, final_j = choose_best_ancilla(grid, "cx", a, b)
            except LayoutInfeasible:
                continue  # try moving the other qubit instead
            x, y = find_index(grid, a)
            temp_dag, temp_grid = cnot_move(path, final_i, final_j, grid, x*n_cols+y, dag)
            if temp_dag is not None:
                dag, grid = temp_dag, temp_grid
                true_ci, true_cj = find_index(grid, control)
                true_ti, true_tj = find_index(grid, target)
                dag, grid, applied = try_apply_CNOT(dag, grid, true_ci, true_cj, true_ti, true_tj, node)
                if applied:
                    return dag, grid
                else:
                    dag,grid= level_3_for_CNOT(grid,dag,control,target,node)
                    return dag,grid


         
        return dag, grid
    

    def process_S(grid, initial_grid,node,dag):

        _,true_i,true_j =convert_index(grid,initial_grid,node )
        if check_S((true_i,true_j), grid):
            ancilla= check_S((true_i,true_j), grid)
            dag= apply_S(grid,(true_i,true_j), random.choice(ancilla), dag,node)
            return dag,grid
        else: 
            temp_grid,temp_dag= move_ancilla_for_S( (true_i,true_j),grid,dag)
            if temp_grid is not None:
                grid,dag=temp_grid,temp_dag
            if check_S((true_i,true_j), grid):
                ancilla= check_S((true_i,true_j), grid)
                dag= apply_S(grid,(true_i,true_j), random.choice(ancilla), dag,node)
                return dag,grid
            
        return dag,grid

    stuck_counter = 0
    max_stuck_iterations = 4# threshold to break or warn
    previous_ops = set() 
    y=list(get_op_layers(dag))
            
            
    while y[0].op.name!="barrier":

        for node in y:

            if node.op.name=="x":
                qindex,_,_=convert_index(grid,initial_grid,node )
                dag=apply_X(qindex,dag,node)
                    
            if node.op.name == "h":
               
                dag,grid= process_H(grid, initial_grid,node,dag)
                        

            if node.op.name== "cx":
                
                dag,grid= process_CNOT(grid, initial_grid,node,dag)
                            
    
            if node.op.name=="rz":
                
                qindex,true_i,true_j =convert_index(grid,initial_grid,node )
                counter=counter % divisor
                temp_dag,temp_grid,path_log =rotation_consumption(None,dag,initial_grid,grid,node,counter,path_log)
                if temp_dag is not None:
                    counter+=1
                    dag,grid=temp_dag,temp_grid
                                

            if node.op.name=="t" or node.op.name=="tdg":
                
                dag,grid, counter, path_log= process_T(grid, initial_grid,node,dag, counter, divisor, path_log)
                            
                
            if node.op.name == "s" or node.op.name=="sdg":
                
                dag,grid= process_S(grid, initial_grid,node,dag )
                        

            if node.op.name=="sx":
                
                qindex,true_i,true_j =convert_index(grid,initial_grid,node )
                counter=counter % divisor
                temp_dag,temp_grid,path_log =rotation_consumption(None,dag,initial_grid,grid,node,counter,path_log, SXGATE=True)
                if temp_dag is not None:
                    counter+=1
                    dag,grid=temp_dag,temp_grid

        
                        
        y=list(get_op_layers(dag))
        if y and set(y) == previous_ops:
            stuck_counter += 1
            if stuck_counter >= max_stuck_iterations:
                break  # hand back to apply_sequence, which raises LayoutInfeasible if no level makes progress
        else:
            stuck_counter = 0  # Reset if progress is made

        previous_ops = set(y)



    return dag, grid, counter, path_log

        
     

def level_2(dag, initial_grid,grid,counter, divisor,path_log):

    def check_CNOT_new(control_coords,target_coords,grid,dag):
        ci,cj=control_coords
        ti,tj=target_coords
        ancilla_options=[ (ti,cj)]
       
        potential_cells = []  # (x, y) that satisfy the coordinate condition
        valid_cells = []      # (x, y) that satisfy both coordinate and grid condition

        for option in ancilla_options:
            x,y=option
            if abs(ci - ti) == 1 and abs(cj - tj) == 1 and 0<=x< len(grid) and 0<=y<len(grid[0]):
                potential_cells.append((x, y))
                if grid[x][y] == 0:
                    valid_cells.append((x, y))
        

        if valid_cells:
            
            return random.choice(valid_cells)
        elif potential_cells:
            try:
                dof= space_search(control_coords, grid, "cx")
            except LayoutInfeasible:
                return False
            dag,grid= build_ancilla(dof,grid,control_coords,"cx",dag,target_coords)
            if  grid is not None:
                ancilla=check_CNOT(control_coords,target_coords,grid,dag)
                return ancilla
             
        else:
            return  False

    def process_H(grid, initial_grid, node, dag):
        
        qindex,true_i,true_j=convert_index(grid, initial_grid, node)
        if check_H((true_i,true_j),grid):
            ancilla=check_H((true_i,true_j),grid)
            dag=apply_H(ancilla, qindex,grid, dag, node)

        else:
            try:
                dof= space_search((true_i,true_j), grid, "h")
            except LayoutInfeasible:
                return dag, grid
            dag,grid= build_ancilla(dof,grid,(true_i,true_j),node.op.name,dag)
            qindex,true_i,true_j=convert_index(grid, initial_grid, node)  # build_ancilla can move the qubit
            if not check_H((true_i,true_j),grid):
                dag,grid= free_a_neighbour((true_i,true_j), grid, dag)
            if check_H((true_i,true_j),grid):
                ancilla=check_H((true_i,true_j),grid)
                dag=apply_H(ancilla, qindex,grid, dag, node)

        return dag,grid   
    
    def process_CNOT(grid, initial_grid, node, dag):
         
        node_qargs = [q._index for q in node.qargs]
        n_cols = len(grid[0])
        ci, cj = divmod(node_qargs[0], n_cols)
        ti, tj = divmod(node_qargs[1], n_cols)
        control = initial_grid[ci][cj]
        target = initial_grid[ti][tj]
        true_ci, true_cj = find_index(grid, control)
        true_ti, true_tj = find_index(grid, target)

        def try_apply(check, dag, grid):
            c, t = find_index(grid, control), find_index(grid, target)
            ancilla = check(c, t, grid, dag)
            c, t = find_index(grid, control), find_index(grid, target)  # the check can move qubits
            if ancilla and cnot_is_legal(grid, c, t, ancilla):
                dag = apply_CNOT(ancilla, [n_cols * c[0] + c[1], n_cols * t[0] + t[1]], grid, dag, node)
                return dag, True
            return dag, False

        # Direct CNOT attempt
        dag, applied = try_apply(check_CNOT, dag, grid)
        if applied:
            return dag,grid
        else:
            # Try moving ancilla and CNOT via routing
            try:
                path, final_i, final_j = choose_best_ancilla(grid, "cx", grid[true_ci][true_cj], grid[true_ti][true_tj])
            except LayoutInfeasible:
                return dag, grid
            dag, grid = move_extended(path[1:], grid, dag, [(true_ci,true_cj),(true_ti, true_tj)])
            true_ci, true_cj = find_index(grid, control)
            temp_dag, temp_grid = cnot_move(path[1:], final_i, final_j, grid, true_ci*n_cols+ true_cj, dag)

            if temp_dag is None:
                # Retry after moving the target if initial move failed
                true_ci, true_cj = find_index(grid, control)
                true_ti, true_tj = find_index(grid, target)
                dag, grid = move_extended(path[1:], grid, dag, [(true_ci,true_cj),(true_ti, true_tj)])
                new_dag, new_grid = cnot_move(path[1:], final_i, final_j, grid, n_cols * true_ci + true_cj, dag)
                if new_grid is not None:
                    dag, grid = new_dag, new_grid

                dag, applied = try_apply(check_CNOT_new, dag, grid)
                if applied:
                    return dag, grid
            else:
                # If ancilla move/CNOT was successful
                dag, grid = temp_dag, temp_grid
                dag, applied = try_apply(check_CNOT_new, dag, grid)
                if applied:
                    return dag, grid

        c, t = find_index(grid, control), find_index(grid, target)
        if abs(c[0] - t[0]) == 1 and abs(c[1] - t[1]) == 1 and grid[t[0]][c[1]] not in (0, -1):
            dag, grid = move_extended([(t[0], c[1])], grid, dag, [c, t])  # push the ancilla cell's qubit away
            dag, applied = try_apply(check_CNOT, dag, grid)
            if applied:
                return dag, grid

        # last resort: move the control, or else the target, through free cells into a diagonal
        # position whose ancilla cell is free
        c, t = find_index(grid, control), find_index(grid, target)
        free = lambda x: 0 <= x[0] < len(grid) and 0 <= x[1] < n_cols and grid[x[0]][x[1]] == 0
        if not move_to_free_cell(grid, dag, c, lambda d: abs(d[0] - t[0]) == 1 and abs(d[1] - t[1]) == 1
                                 and free((t[0], d[1]))):
            move_to_free_cell(grid, dag, t, lambda d: abs(d[0] - c[0]) == 1 and abs(d[1] - c[1]) == 1
                              and free((d[0], c[1])))
        dag, applied = try_apply(check_CNOT, dag, grid)
        return dag, grid


    def process_S(grid, initial_grid, node, dag):
        qindex,true_i,true_j=convert_index(grid, initial_grid, node)
        ancilla= check_S_new((true_i,true_j), grid,dag)
        qindex,true_i,true_j=convert_index(grid, initial_grid, node)  # check_S_new can move qubits
        if ancilla and adjacent((true_i,true_j), ancilla) and grid[ancilla[0]][ancilla[1]] == 0:
            dag= apply_S(grid, (true_i,true_j),ancilla, dag, node)
            return dag,grid
        else:
            try:
                dof= space_search((true_i,true_j), grid, "s")
            except LayoutInfeasible:
                return dag, grid
            dag,grid= build_ancilla(dof,grid,(true_i,true_j),node.op.name,dag)
            qindex,true_i,true_j=convert_index(grid, initial_grid, node)  # build_ancilla can move the qubit
            ancilla= check_S_new((true_i,true_j), grid,dag)
            qindex,true_i,true_j=convert_index(grid, initial_grid, node)
            if ancilla and is_data(grid, (true_i,true_j)) and adjacent((true_i,true_j), ancilla) and grid[ancilla[0]][ancilla[1]] == 0:
                dag= apply_S(grid, (true_i,true_j), ancilla, dag, node)
                return dag,grid
                
            return dag,grid
      
    y=list(get_op_layers(dag))

    for node in y:
            if node.op.name=="x":
                 
                qindex,_,_=convert_index(grid, initial_grid, node)
                dag=apply_X(qindex,dag,node)
                
            if node.op.name == "h":
                
                dag,grid= process_H(grid, initial_grid,node,dag)
                            

            if node.op.name== "cx":
               
                dag,grid= process_CNOT(grid, initial_grid,node,dag)
                                
                                
            if node.op.name=="rz":
                
                qindex,true_i,true_j=convert_index(grid, initial_grid, node)
                try:
                    dof= space_search((true_i,true_j), grid, "rz")
                    dag,grid= build_ancilla(dof,grid,(true_i,true_j),node.op.name,dag)
                except LayoutInfeasible:
                    pass  # no ancilla prepared; rotation_consumption looks for one itself
                counter=counter%divisor
                temp_dag,temp_grid,path_log =rotation_consumption(None,dag,initial_grid,grid,node,counter,path_log)
                if temp_dag is not None:
                    counter+=1
                    dag,grid=temp_dag,temp_grid
    
            if node.op.name=="sx":
                
                qindex,true_i,true_j =convert_index(grid,initial_grid,node )
                counter=counter % divisor
                temp_dag,temp_grid,path_log =rotation_consumption(None,dag,initial_grid,grid,node,counter,path_log, SXGATE=True)
                if temp_dag is not None:
                    counter+=1
                    dag,grid=temp_dag,temp_grid
            
            if node.op.name=="t" or node.op.name=="tdg":
                
                dag,grid, counter, path_log= process_T(grid, initial_grid,node,dag, counter, divisor, path_log)
                                
            if node.op.name == "sdg" or node.op.name=="s":
                 
                dag,grid=process_S(grid, initial_grid, node, dag)
                        
                

    return dag,grid, counter,path_log
 


 