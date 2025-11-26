from qiskit.converters import circuit_to_dag 
from qiskit.circuit import  Measure
from Functions import find_layout, add_global_barrier, apply_sequence, find_initial_layout_2D
from Hamiltonians import ising_2D, fermi_hubbard, heisenberg
from Postprocessing import generate_time
from levels import level_1, level_2
from qiskit.circuit.library import Barrier
import numpy as np
 


iterations=1 #execution time is different everytime so we normally run a few to determine the minimum overall time. 

 
#val is the same as the number of qubit. The experiment iterates over qubits from 4^{2} to 12^{2} (data qubits from 2^{2} to 10^{2}). 
# The default layout here is of a grid with 3 routing paths. 

edge_size=14 #number qubits on one edge of the grid. 
    
for i in range (0,iterations):

        n=edge_size
        circuit_size=(edge_size-2)**2
        grid=find_initial_layout_2D(n,circuit_size ) 
        #del grid[n-1]
        
        ##change layout by using the add_paths_upto functions
        # 17 applies paths 0,1,2,3...,17 in order
        #grid = add_paths_upto(base,17 )
        
        initial_layout=find_layout(grid,circuit_size)
        qc=ising_2D(n,initial_layout)        
        dag = add_global_barrier(circuit_to_dag(qc)) ## global barrier for distinction between old and new dag
        nedsss=[node for node in list(dag.topological_op_nodes()) if node.op.name=="rz"]
        opt_levels = {
            'level1': level_1,
            'level2': level_2
        }
        weights = {
            'level1':  0.1,
            'level2':0.9,
        }
        magic_state_factories = 1 
        """""grid,dag, path_log = apply_sequence(dag,grid, opt_levels, magic_state_factories, weights)

        #postprocessing the resultant dag
        for node in list(dag.op_nodes()) :
            if isinstance(node.op,  Measure):
                dag.remove_op_node(node) #remove measure operations
        for node in list(dag.op_nodes()):
            if isinstance(node.op, Barrier):
                    dag.remove_op_node(node) #remove unwanted barriers operations
        
        overhead , execution_time, new_dag= generate_time(dag) # processing the output dag and generating the total time 
        ### overheads+=overhead
        ###total_time+=execution_time 

        current_layer=[]
        layers=[]
        operations=[]
        operation_layer=[]
        cnot_operation=[]
        all_cnot=[]
        for node in list(new_dag.topological_op_nodes()) :
                if node.op.name == 'barrier':
                    if current_layer:
                        layers.append(current_layer)
                        current_layer = []
                        operation_layer.append(operations)
                        operations=[]
                        all_cnot.append(cnot_operation)
                        cnot_operation=[]


                else:
                    if node.op.name=="ccx":
                        cnot_operation.append(node)
                    operations.append(node)
                    qubits1 = [q._index for q in node.qargs]
                    for qubits in qubits1:
                        current_layer.append(qubits)

        len_layer= [len(layer) for layer in operation_layer]
        len_cnot=[len(layer) for layer in all_cnot]
        len_layer_nodes=[len(layer) for layer in layers]
        active_qubits=[]
        for i in range(len(len_layer_nodes)):
            active_layers_non_qubits=len_layer_nodes[i]-len_layer[i]-len_cnot[i]
            active_qubits.append(100+active_layers_non_qubits)

        sum_=sum(  (active_qubits[i] * len_layer[i]) for i in range(len(len_layer)))

        total1= sum_* 3* (0.1)**9
        total2= sum_* 3* (0.1)**8
        

average_time= [x / iterations for x in overhead]
average_time_overhead= [x / iterations for x in execution_time]"""""

    
 
with open("output.txt", "w") as f:
        print(len(nedsss), file=f)


 

    
 
        
