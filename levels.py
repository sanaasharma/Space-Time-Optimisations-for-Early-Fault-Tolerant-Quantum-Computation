from Functions import get_op_layers,find_index,  choose_best_ancilla, space_search
from Move_functions import  cnot_move, move_ancilla,move_ancilla_for_S, move_extended, move_along_path_reversed, build_ancilla, perform_move , MoveGate
from Gates import check_CNOT, apply_CNOT, apply_CT,check_S, apply_S, apply_X, check_H, apply_H, check_S_new
from Rotations import magic_state_routing_first_pass, rotation_consumption, check_CT, convert_index
import random



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
            if grid[ti][cj]!= 0:
                dag,grid= move_extended([(ti,cj)],grid,dag,[(ci,cj)])
                path, final_i, final_j = choose_best_ancilla(grid, "cx", control, target)
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
            if check_CNOT((ci, cj), (ti, tj), grid, dag):
                ancilla = check_CNOT((ci, cj), (ti, tj), grid, dag)
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
             
            path, final_i, final_j = choose_best_ancilla(grid, "cx", a, b)
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
    
    def process_T(grid, initial_grid, node, dag, counter):
        _,true_i,true_j =convert_index(grid,initial_grid,node )
        qubit_number=grid[true_i][true_j]
        counter=counter % divisor
        temp_dag,temp_grid,final_index = magic_state_routing_first_pass(grid,grid[true_i][true_j],dag,counter, S_presence=True)
        if temp_grid is not None:
            dag,grid=temp_dag,temp_grid
            _,ci,cj =convert_index(grid,initial_grid,node )
            check, dag,grid=check_CT(qubit_number ,( ci,cj),final_index,grid,dag)
            if check is True:
                dag=apply_CT( ( ci,cj), final_index,grid, dag )
                dag.remove_op_node(node)
                counter+=1
                return dag,grid, counter
            
            else: 
                temp_grid,temp_dag=move_ancilla(( ci,cj),final_index ,grid,dag)
                if temp_grid is not None:
                    dag,grid=temp_dag,temp_grid
                    check, dag,grid=check_CT(qubit_number,( ci,cj),final_index,grid,dag)
                    if check:
                        dag=apply_CT(( ci,cj), final_index,grid, dag )
                        dag.remove_op_node(node)
                        counter+=1
                        return dag,grid
                    
        return dag,grid
     
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
                temp_dag,temp_grid,path_log =rotation_consumption(grid[true_i][true_j],dag,initial_grid,grid,node,counter,path_log)
                if temp_dag is not None:
                    counter+=1
                    dag,grid=temp_dag,temp_grid
                                

            if node.op.name=="t" or node.op.name=="tdg":
                
                dag,grid, counter= process_T(grid, initial_grid,node,dag, counter)
                            
                
            if node.op.name == "s" or node.op.name=="sdg":
                
                dag,grid= process_S(grid, initial_grid,node,dag )
                        

            if node.op.name=="sx":
                
                qindex,true_i,true_j =convert_index(grid,initial_grid,node )
                counter=counter % divisor
                temp_dag,temp_grid,path_log =rotation_consumption(grid[true_i][true_j],dag,initial_grid,grid,node,counter,path_log, SXGATE=True)
                if temp_dag is not None:
                    counter+=1
                    dag,grid=temp_dag,temp_grid

        
                        
        y=list(get_op_layers(dag))
        if y and set(y) == previous_ops:
            stuck_counter += 1
            print(f"[Warning] Loop appears stuck. Iteration {stuck_counter}.")
            if stuck_counter >= max_stuck_iterations:
                print("Terminating due to potential infinite loop.")
                break
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
            dof= space_search(control_coords, grid, "cx")
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
            dof= space_search((true_i,true_j), grid, "h")
            dag,grid= build_ancilla(dof,grid,(true_i,true_j),node.op.name,dag)
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

        # Direct CNOT attempt
        if check_CNOT((true_ci, true_cj), (true_ti, true_tj), grid, dag):
            ancilla = check_CNOT((true_ci, true_cj), (true_ti, true_tj), grid, dag)
            new_node_qargs = [n_cols * true_ci + true_cj, n_cols * true_ti + true_tj]
            dag = apply_CNOT(ancilla, new_node_qargs, grid, dag, node)
            return dag,grid
        else:
            # Try moving ancilla and CNOT via routing
            path, final_i, final_j = choose_best_ancilla(grid, "cx", grid[true_ci][true_cj], grid[true_ti][true_tj])
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

                # Update locations after potential changes
                true_ci, true_cj = find_index(grid, control)
                true_ti, true_tj = find_index(grid, target)
                if check_CNOT_new((true_ci, true_cj), (true_ti, true_tj), grid, dag):
                    new_node_qargs = [n_cols * true_ci + true_cj, n_cols * true_ti + true_tj]
                    ancilla = check_CNOT_new((true_ci, true_cj), (true_ti, true_tj), grid, dag)
                    dag = apply_CNOT(ancilla, new_node_qargs, grid, dag, node)
                    return dag, grid
            else:
                # If ancilla move/CNOT was successful
                dag, grid = temp_dag, temp_grid
                true_ci, true_cj = find_index(grid, control)
                true_ti, true_tj = find_index(grid, target)
                if check_CNOT_new((true_ci, true_cj), (true_ti, true_tj), grid, dag):
                    new_node_qargs = [n_cols * true_ci + true_cj, n_cols * true_ti + true_tj]
                    ancilla = check_CNOT_new((true_ci, true_cj), (true_ti, true_tj), grid, dag)
                    dag = apply_CNOT(ancilla, new_node_qargs, grid, dag, node)
                    return dag, grid

        return dag, grid

    def process_T(grid, initial_grid, node, dag, counter):
        _,true_i,true_j =convert_index(grid,initial_grid,node )
        qubit_number=grid[true_i][true_j]
        counter=counter % divisor
        temp_dag,temp_grid,final_index = magic_state_routing_first_pass(grid,grid[true_i][true_j],dag,counter, S_presence=True)
        if temp_grid is not None:
            dag,grid=temp_dag,temp_grid
            _,ci,cj =convert_index(grid,initial_grid,node )
            check, dag,grid=check_CT(qubit_number ,( ci,cj),final_index,grid,dag)
            if check is True:
                dag=apply_CT( ( ci,cj), final_index,grid, dag )
                dag.remove_op_node(node)
                counter+=1
                return dag,grid, counter
            
            else: 
                temp_grid,temp_dag=move_ancilla(( ci,cj),final_index ,grid,dag)
                if temp_grid is not None:
                    dag,grid=temp_dag,temp_grid
                    check, dag,grid=check_CT(qubit_number,( ci,cj),final_index,grid,dag)
                    if check:
                        dag=apply_CT(( ci,cj), final_index,grid, dag )
                        dag.remove_op_node(node)
                        counter+=1
                        return dag,grid
                    
        return dag,grid
       
    def process_S(grid, initial_grid, node, dag):
        qindex,true_i,true_j=convert_index(grid, initial_grid, node)
        if check_S_new((true_i,true_j), grid,dag):
            ancilla= check_S_new((true_i,true_j), grid,dag)
            dag= apply_S(grid, (true_i,true_j),ancilla, dag, node)
            return dag,grid
        else: 
            dof= space_search((true_i,true_j), grid, "s")
            dag,grid= build_ancilla(dof,grid,(true_i,true_j),node.op.name,dag)
            if check_S_new((true_i,true_j), grid,dag):
                ancilla= check_S_new((true_i,true_j), grid,dag)
                dag= apply_S(grid, (true_i,true_j), random.choice(ancilla), dag, node) 
                return dag,grid
                
            return dag,grid
      
    y=list(get_op_layers(dag))

    for node in y:
            if node.op.name=="x":
                 
                qindex,_,_=convert_index(grid, initial_grid, node, dag)
                dag=apply_X(qindex,dag,node)
                
            if node.op.name == "h":
                
                dag,grid= process_H(grid, initial_grid,node,dag)
                            

            if node.op.name== "cx":
               
                dag,grid= process_CNOT(grid, initial_grid,node,dag)
                                
                                
            if node.op.name=="rz":
                
                qindex,true_i,true_j=convert_index(grid, initial_grid, node)
                dof= space_search((true_i,true_j), grid, "rz")
                dag,grid= build_ancilla(dof,grid,(true_i,true_j),node.op.name,dag)
                counter=counter%divisor
                temp_dag,temp_grid,path_log =rotation_consumption(grid[true_i][true_j],dag,initial_grid,grid,node,counter,path_log)
                if temp_dag is not None:
                    counter+=1
                    dag,grid=temp_dag,temp_grid
    
            if node.op.name=="sx":
                
                qindex,true_i,true_j =convert_index(grid,initial_grid,node )
                counter=counter % divisor
                temp_dag,temp_grid,path_log =rotation_consumption(grid[true_i][true_j],dag,initial_grid,grid,node,counter,path_log, SXGATE=True)
                if temp_dag is not None:
                    counter+=1
                    dag,grid=temp_dag,temp_grid
            
            if node.op.name=="t" or node.op.name=="tdg":
                
                dag,grid, counter= process_T(grid, initial_grid,node,dag, counter)
                                
            if node.op.name == "sdg" or node.op.name=="s":
                 
                dag,grid=process_S(grid, initial_grid, node, dag)
                        
                

    return dag,grid, counter,path_log
 


 