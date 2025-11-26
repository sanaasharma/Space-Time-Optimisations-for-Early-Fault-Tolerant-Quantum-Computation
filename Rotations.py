from Functions import find_index, find_best_path, find_path, space_search 
from Gates import CustomTGate, CustomHGate, CustomSGate , check_S, check_H, apply_Z, apply_H, check_S_new, check_CT, apply_CT, apply_S
from Move_functions import move_ancilla,build_ancilla, move_extended, perform_move , MagicMoveGate ,MoveGate,move_along_path_reversed
from qiskit.synthesis import generate_basic_approximations
from qiskit.transpiler.passes.synthesis import SolovayKitaev
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.circuit.library import PhaseGate, SXGate, UGate, CXGate, IGate, TGate, TdgGate, SGate, CCXGate, SXGate, SdgGate, HGate
from qiskit import QuantumCircuit, transpile
import random

import numpy as np



 


def magic_state_routing_first_pass(grid,qubit,dag, number,path_log, S_presence=bool):
     
  
    
    if number==0:
        start_i , start_j=int(len(grid)/2) , 0
    "For more magic state factories, add the factory placement as given below:"
    "if number==1:"
    "  start_i , start_j=int(len(grid)/2)-1 , int(len(grid[0])) "
    "add factories around the edges."
 
    start=start_i , start_j
     

    target_i, target_j= find_index(grid,qubit)
     
    path,f_i, f_j = find_best_path(grid,start_i,start_j,target_i,target_j,S_presence)
    end=f_i,f_j
    
    temp_grid,temp_dag= process_nonzero_cells(grid,path,end,dag)
    
    
    if temp_grid is not None:
        
        grid,dag =temp_grid,temp_dag
        grid[start_i][start_j]='m'
        custom_magicmove= MagicMoveGate(path, start,end,len(path),operation_time=0)
        dag, _, _, grid = custom_magicmove.update_dag_and_grid( grid,dag,path, start, end) 
        return grid, dag , (f_i,f_j),path_log
    
    else: 
        path,_ = find_path(grid,start,end)
        path_log.append(path)
        temp_dag,temp_grid= move_extended(path, grid,dag,path)
        if temp_grid is not None:
            
            grid,dag =temp_grid,temp_dag
            grid[start_i][start_j]='m'
            custom_magicmove= MagicMoveGate(path, start,end,len(path),operation_time=0)
            dag, _, _, grid = custom_magicmove.update_dag_and_grid( grid,dag,path, start, end)
            path_log.append(grid)
            return grid, dag , (f_i,f_j),path_log
        
    return None,None,None, path_log
 
def convert_index(grid,initial_grid,node):
        node_qargs = [q._index for q in node.qargs]
        i,j= divmod(node_qargs[0],len(grid[0]))
        true_i,true_j=find_index(grid,initial_grid[i][j])
        qindex=len(grid[0])*true_i+true_j
        return qindex , true_i, true_j
  
 

def rotation_consumption(qubit,dag,initial_grid,grid,node_input,number,path_log, SXGATE=None, executed_nodes=None):

    

    if executed_nodes is None:
        executed_nodes = set()

    if SXGATE==None:
        angle=node_input.op.params[0]
        circuit = QuantumCircuit(1)
        circuit.rz(angle, 0)
    else:
        circuit = QuantumCircuit(1)
        circuit.sx(0)
    
    basis = ["s", "sdg", "t", "tdg", "z", "h"]
    approx = generate_basic_approximations(basis, depth=3) 
    skd = SolovayKitaev(recursion_degree=2, basic_approximations= approx)
    discretized = skd(circuit)
    rotation_circuit=circuit_to_dag(discretized)
    success=False
     

    def rotation_check(rotation_circuit,dag,grid,executed_nodes,qubit,initial_grid,node_input,path_log):
       
        for node in list(rotation_circuit.op_nodes()):
            
            
            if node in executed_nodes:
                continue

            try:

                if node.op.name =='t': 

                    
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
                    
                 
                if node.op.name=='z':
                     
                    qindex,_,_=convert_index(grid,initial_grid,node_input)
                    dag=apply_Z(qindex,dag,node)
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
                        if check_H((true_i,true_j),grid):
                            ancilla=check_H((true_i,true_j),grid)
                            dag=apply_H(ancilla, qindex,grid, dag, node,node_input)
                            executed_nodes.add(node)
                        else: 
                            return None,None,executed_nodes,path_log

                
                if node.op.name in ('s', 'sdg'):

                    qindex,true_i,true_j=convert_index(grid,initial_grid,node_input)
                    if check_S_new((true_i,true_j), grid,dag):
                        ancilla= check_S_new((true_i,true_j), grid,dag)
                        dag= apply_S(grid, (true_i,true_j),ancilla, dag, node,node_input)
                        executed_nodes.add(node)
                    else: 
                        dof= space_search((true_i,true_j), grid, "s")
                        dag,grid= build_ancilla(dof,grid,(true_i,true_j),node.op.name,dag)
                        ancilla= check_S_new((true_i,true_j), grid,dag)
                        if check_S_new((true_i,true_j), grid,dag):
                            
                            ancilla= check_S_new((true_i,true_j), grid,dag)
                            dag= apply_S(grid, (true_i,true_j), random.choice(ancilla), dag, node,node_input)
                            executed_nodes.add(node)

                        else:
                            return None, None ,executed_nodes  ,path_log  

 
                 
            except Exception as e:
                return  None,None ,executed_nodes,path_log
            
        return dag,grid,executed_nodes,path_log
    
    
    while not success:
        new_dag,_,executed_nodes,path_log= rotation_check(rotation_circuit,dag,grid,executed_nodes,qubit,initial_grid,node_input,path_log)
        if new_dag is None:
            continue
        else:  
            if len(executed_nodes) == len(list(rotation_circuit.op_nodes())):
                dag.remove_op_node(node_input)
                success = True


    return dag,grid,path_log 
 
  
def process_nonzero_cells(grid, path, end,dag ):

    def is_valid(i, j, grid):
        return 0 <= i < len(grid) and 0 <= j < len(grid[0]) and grid[i][j] == 0

    def get_neighbors(i, j):
        directions = [(0,1), (1,0), (0,-1), (-1,0)]
        return [(i+di, j+dj) for di, dj in directions]

    def move_cell(grid, cell, path_set,dag):
        i, j = cell
        Found= False
        for ni, nj in get_neighbors(i, j): 
            if (ni, nj) not in path_set  and is_valid(ni, nj, grid):
                Found=True
                return (ni, nj),grid,dag
        
        if Found==False:

            for ni, nj in get_neighbors(i, j):
                    if (ni,nj) not in path_set  and 0 <= ni < len(grid) and 0 <= nj < len(grid[0]):
                        if grid[ni][nj]!=-1: 
                            for nni,nnj in get_neighbors(ni,nj):
                                if (nni, nnj) not in path_set  and is_valid(nni, nnj, grid) and (nni,nnj) != (i,j):
                                    dag, grid, _, _= perform_move(nni, nnj,ni, nj , dag, grid)
                                    return (ni, nj),grid,dag
                    elif (ni,nj) in path_set  and is_valid(ni,nj,grid):
                            for nni,nnj in get_neighbors(ni,nj):
                                if (nni, nnj) not in path_set  and is_valid(nni, nnj, grid) and (nni,nnj) != (i,j):
                                   return (nni, nnj),grid,dag

                        


        return None, None,None   

    path_set = set(path)
    if all(grid[x][y]!=-1 for (x,y) in path):
        return None, None
    
    for (i, j) in path:
        if grid[i][j] != 0:
            end,grid,dag = move_cell(grid, (i, j), path_set,dag)
            if end is not None:
                ti,tj=end[0],end[1]
                dag, grid, _, _= perform_move(ti, tj, i,  j , dag, grid)
                    
        else:
            return None, None

             
    return grid,dag 


 