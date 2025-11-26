from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.dagcircuit import DAGOpNode
from qiskit.circuit import Instruction, Gate
from qiskit import QuantumCircuit
from copy import deepcopy
import copy
 
def check_commutation(dag_nodes):
    layers = []
    current_layer = []
    
     
    for node in dag_nodes:
        if node.op.name == 'barrier':
            if current_layer:
                layers.append(current_layer)
                current_layer = []
        else:
            current_layer.append(node)
    
     
    if current_layer:
        layers.append(current_layer)
    
    commute_list=[]
    for k,layer in enumerate(layers):
        for i, node1 in enumerate(layer):
            for j, node2 in enumerate(layer):
                if i >= j:
                    continue   
                if not commute(node1, node2):  
                    if node1.op.name=='magicmove' and node2.op.name=='ct':
                        continue
                    else:
                        commute_list.append((node2,k,j))
    if len(commute_list)>0:
        return commute_list
    else:
        return None   

def commute(node1, node2):
    
    qubits1 = [q._index for q in node1.qargs]
    qubits2 = [q._index for q in node2.qargs]
     
    if set(qubits1).intersection(set(qubits2)):
        
        return False

    return True

def move_commuting_nodes(non_commuting_nodes, dag_nodes):

    layers = []
    current_layer = []
    for node in dag_nodes:
        if node.op.name == 'barrier':
            if current_layer:
                layers.append(current_layer)
                current_layer = []
        else:
            current_layer.append(node)
    
    if current_layer:
        layers.append(current_layer)
    offset = 0 
    for k,value in enumerate(non_commuting_nodes):
        node, x, y=value
        if k>0:
            _,temp,_=  non_commuting_nodes[k-1]
            if temp!=x:
                same_layer_offset=0
                offset2=0
        else: 
            same_layer_offset=0
            offset2=0

        next_layer = x + 1+ offset- same_layer_offset   
        if next_layer < len(layers):
            next_layer_nodes = layers[next_layer]
            can_move = True

            for next_layer_node in next_layer_nodes:
                if not commute(node, next_layer_node):  
                    can_move = False
                    break
            if can_move:
                if same_layer_offset>0:
                         
                        layers[x+offset-same_layer_offset].pop(y-offset2)
                else:
                        layers[x+offset].pop(y-offset2)
                next_layer_=x + 1 +offset- same_layer_offset 
                layers[next_layer_].append(node)
                offset2+=1
                
                
            else:
                if same_layer_offset>0:
                        layers[x+offset-same_layer_offset].pop(y-offset2)
                else:
                        layers[x+offset].pop(y-offset2)
                   
                new_layer_index = x + 1 +offset- same_layer_offset  # The new layer will be after the current layer
                new_layer = [node]  # Create a new layer with just this node
                layers.insert(new_layer_index, new_layer)
                offset+=1
                same_layer_offset+=1
                offset2+=1
                         
    return layers
 
def get_time_between_magicmoves(layers ):
     

     
    cumulative_time = 0
    results = []

    for i, layer in enumerate(layers):
        

            x= sum(1 for node in layer if node.op.name=='ct')
            
            if any(node.op.name == 'ct' for node in layer):
                ct_index = next((i for i, node in enumerate(layer) if node.op.name == 'ct'), None)
                cidx= layer[ct_index].qargs[-1]._index 
                results.append((i, cumulative_time))
                cumulative_time = 0  # reset after each magicmove
                if any(node.op.name == 'magicmove' for node in layer):
                    ct_index = next((i for i, node in enumerate(layer) if node.op.name == 'magicmove'), None)
                    midx=layer[ct_index].qargs[-1]._index 
                    if cidx==midx:
                        temp = max(get_move_time_total(node) for node in layer if node.op.name not in {"magicmove" or node.op.name!= "ct" } )
                         
                        if  temp ==0:
                             
                            layer_time=2.5 #2.5
                        else:
                            layer_time= max(temp,2.5) #2.5
                        
                        cumulative_time += layer_time

                    if x ==1:
                        continue
                    else:
                        for tmp in range(x-1):
                            results.append((i, cumulative_time))
                            cumulative_time = 0  # reset after each magicmove
                            cumulative_time += 1
                            continue
                else:  
                    layer_time = max(get_move_time_total(node) for node in layer )
                    cumulative_time += layer_time
                    if x ==1:
                        continue
                    else:
                        for tmp in range(x-1):
                            results.append((i, cumulative_time))
                            cumulative_time = 0  # reset after each magicmove
                            cumulative_time += 1
                            continue
            else:
                if layer:
                    layer_time = max(get_move_time_total(node) for node in layer)
                    cumulative_time += layer_time
                 
    results.append((i, cumulative_time))

    return results

def remove_repetitions(data):
    seen = set()   
    unique_data = []  

    for item in data:
         
        operation_tuple = (item[0], item[1], item[2])
        
        
        if operation_tuple not in seen:
            unique_data.append(item)   
            seen.add(operation_tuple)  

    return unique_data   
        
def tensor_product_instruction(op1, op2):
    
    name = f"{op1.op.name}_{op2.op.name}_tensor"
    qargs = op1.qargs + op2.qargs
    return Instruction(name=name, num_qubits=len(qargs), num_clbits=0, params=[]), qargs
 
def move_magicmoves_through_commuting_ops(flat_nodes):

    def commutes(op1, op2):
        
        
        q1 = set(op1.qargs)
        q2 = set(op2.qargs)
        return q1.isdisjoint(q2)
    
    dag_nodes = copy.deepcopy(flat_nodes)
     
    i = 0
    while i < len(dag_nodes):
        node = dag_nodes[i]
        if node.op.name == 'magicmove':
             
             
          
            if i + 1 < len(dag_nodes) and dag_nodes[i + 1].op.name == 'ct':
                qarg_indices = [q._index for q in dag_nodes[i+1].qargs]
                magic_indices = [node.qargs[0]._index, node.qargs[-1]._index]

                if magic_indices[1] == qarg_indices[1]:
                 

                    ct_node = dag_nodes[i + 1]
                    composite_op, qargs = tensor_product_instruction(node, ct_node)

                    # Create a virtual composite node (only for commutativity testing)
                    virtual_node = deepcopy(node)
                    virtual_node.op = composite_op
                    virtual_node.qargs = qargs
                     
                    j = i + 2  # Start after ct
                    while j < len(dag_nodes):
                        next_node = dag_nodes[j]
                        if next_node.op.name=="magicmove":
                            #if dag_nodes[j+1]=='ct':
                            break
                        elif next_node.op.name == 'barrier':
                            j += 1
                             
                            continue
                        elif commutes(virtual_node, next_node):
                             
                            j += 1
                            
                            continue

                        break
                     

                    if j - 2 > i:
                        # Move both magicmove and ct to index j
                         
                        magic_node = dag_nodes.pop(i)
                        ct_node = dag_nodes.pop(i)  # note: i stays same because list shifted
                        dag_nodes.insert(j - 2, ct_node)
                        dag_nodes.insert(j - 2, magic_node)
                        i+=1
                        continue  # Check again from the same i
                    else:
                        
                       i+=1

                else: 
                    j = i + 1
                    while j < len(dag_nodes):
                        next_node = dag_nodes[j]
                        if next_node.op.name == 'barrier':
                            
                            j += 1
                            continue
                        if commutes(node, next_node):
                            j += 1
                            
                            continue
                        break

                    if j - 1 > i:
                         
                        dag_nodes.insert(j-1, dag_nodes.pop(i))
                        i+=1
                        continue  # Reprocess current index
                    else:
                         
                        i+=1

            else:
                # Regular magicmove case
                j = i + 1
                while j < len(dag_nodes):
                    next_node = dag_nodes[j]
                    if next_node.op.name=="magicmove":
                             
                        break

                    elif next_node.op.name == 'barrier':
                        
                        j += 1
                        continue
                    
                    elif commutes(node, next_node):
                        j += 1
                        
                        continue
                    break
                
                if j==i+2:
                    if dag_nodes[j]=='magicmove':
                        i+=1
                        continue
                    
                if j - 1 > i:
                     
                    dag_nodes.insert(j-1, dag_nodes.pop(i))
                    continue  # Reprocess current index
                else:
                     
                    i+=1

        i += 1
    return dag_nodes
   
def get_move_time_total(node):
    gate_durations = {
        'ch': 3,#3
        'move': 1,
        'ccx':2 ,#2
        'magicmove':1,
        'cs': 1.5  , #1.5
        'ct':2.5 #2.5
    }
      
    if node.op.__class__.__name__ == 'CustomDelayGate':
        print(node.op)
        return node.op.delay_time

    
        
    return gate_durations.get(node.op.name,1 )

def update_grid_from_dag_moves(layer, grid):
   
     
    rows, cols = len(grid), len(grid[0])

    

    for node in layer:
        if node.op.name == "move":
            qarg_indices = [q._index for q in node.qargs]
             

            source, target = qarg_indices[0], qarg_indices[-1]
            

            source_i, source_j = divmod(source,cols)
            target_i, target_j = divmod(target,cols)
           
            if grid[target_i][target_j]==0:
                grid[target_i][target_j]=grid[source_i][source_j]
                grid[source_i][source_j]=0

    return grid

def build_gridprocess(layers, initial_grid):
    grid= copy.deepcopy(initial_grid)
    grid_process=[copy.deepcopy(grid)]
    for layer in layers:
        grid= update_grid_from_dag_moves(layer, initial_grid, grid)
        grid_process.append(copy.deepcopy(grid))
    
    return grid_process

def generate_time(dag):

    op_layers = []
    for layer in dag.multigraph_layers():
            opss_only = [node for node in layer if isinstance(node, DAGOpNode)]
            if opss_only:
                op_layers.append(opss_only)

    new_circuit = QuantumCircuit(*dag.qregs.values())
    for layer in op_layers:
                for node in layer:
                    
                    new_circuit.append(node.op, node.qargs, node.cargs)
                # Insert a barrier on all qubits in this layer
                qubits_in_layer = sorted(
                    {q for node in layer for q in node.qargs},
                    key=lambda q:dag.qubits.index(q)
                )
                new_circuit.barrier(*qubits_in_layer)

    flat_nodes = list(circuit_to_dag(new_circuit).topological_op_nodes())  
    magic_moves_dag=move_magicmoves_through_commuting_ops(flat_nodes)
    check_temp=check_commutation(magic_moves_dag)
    removed_dag=remove_repetitions(check_temp)
    sorted_dag = sorted(removed_dag, key=lambda magic_moves_dag: (magic_moves_dag[1], magic_moves_dag[2]))
    final_dag=move_commuting_nodes(sorted_dag, magic_moves_dag )
    new_circuit = QuantumCircuit(*dag.qregs.values())
    for layer in final_dag:
        for node in layer:
            new_circuit.append(node.op, node.qargs, node.cargs)
                            # Insert a barrier on all qubits in this layer
        qubits_in_layer = sorted(
                                {q for node in layer for q in node.qargs},
                                key=lambda q:dag.qubits.index(q)
                            )
        new_circuit.barrier(*qubits_in_layer )
            
    new_dag= circuit_to_dag(new_circuit)
    result=get_time_between_magicmoves(final_dag )
    processing_time=[11] #for different distillation factories and their processing times
    overheads=[]
    totals=[]
    for pro_time in processing_time:
        overhead=sum(time- pro_time if time> pro_time else 0 for _,time in result)
        _,y=result[-1]
        total= sum(time if time> pro_time else pro_time for _,time in result[:-1]) + y
        overheads.append(overhead) #more than processing time used for distillation
        totals.append(total)
         
    return overheads, totals, new_dag
        #total_time+=time
        
    
 
 
        
