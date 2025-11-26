
from qiskit import QuantumCircuit 
import numpy as np
   
 
def fermi_hubbard(n, initial_layout):

    qc = QuantumCircuit(n**2 )  # Create circuit with 16 qubits

    def apply_double_x(qc, initial_layout):
         
        #start = 1 if grp else 2
        
        for i in range(0,len(initial_layout),2):

            qc.h([initial_layout[i], initial_layout[i+1]])
            qc.cx(initial_layout[i], initial_layout[i+1])
            qc.rz(np.pi/4, initial_layout[i+1])
            #qc.t( initial_layout[i+1])
            qc.cx(initial_layout[i], initial_layout[i+1])
            qc.h([initial_layout[i], initial_layout[i+1]])
   
        return qc


    def apply_double_z(qc, initial_layout):
        
        #start = 1 if grp else 2
        x=int(len(initial_layout)/2)

        for i in range(0,x):
                qc.cx(initial_layout[i], initial_layout[i+x])
                qc.rz(np.pi/4, initial_layout[i+x])
                #qc.t( initial_layout[i+x])
                qc.cx(initial_layout[i], initial_layout[i+x])

        return qc
    
    def apply_double_y(qc,initial_layout):
        
        #start = 1 if grp else 2

        for i in range(0,len(initial_layout),2):

            #qc.sdg([initial_layout[i], initial_layout[i+1]])
            qc.s([initial_layout[i], initial_layout[i+1]])
            qc.h([initial_layout[i], initial_layout[i+1]])
            qc.cx(initial_layout[i], initial_layout[i+1])
            qc.rz(np.pi/4, initial_layout[i+1])
            #qc.t( initial_layout[i+1])
            qc.cx(initial_layout[i], initial_layout[i+1])
            qc.h([initial_layout[i], initial_layout[i+1]])
            qc.s([initial_layout[i], initial_layout[i+1]])
   
        return qc
    
    qc= apply_double_x(qc,initial_layout)
    qc= apply_double_y(qc,initial_layout)
    qc= apply_double_z(qc,initial_layout)

    return qc

def ising_2D(n,initial_layout):

    qc = QuantumCircuit(n**2 )  # Create circuit with 16 qubits
    """blocked_row= [int(n/2)-1,int(n/2),int(n/2)+1]
    blocked_col= [int(n/2)-1,int(n/2)]"""
   
    def get_qubit_index(row, col, n):
         
        return row * n + col
    
    
    def apply_all_rx(qc, initial_layout):
         
        for q in initial_layout :
            qc.h(q)
            qc.rz(np.pi/4, q)
            #qc.t(q)
            qc.h(q)
        
        return qc
    
    def apply_double_z(qc, n, dir=bool, grp=bool):
        
        start = 1 if grp else 2

        for i in range(1,n-1):
            for j in range(start, n - 2, 2):
                row, col = (i, j) if dir else (j, i)
 

                q1 = get_qubit_index(row, col, n    )
                q2 = get_qubit_index(row, col + 1, n    ) if dir else get_qubit_index(row + 1, col, n )

                qc.cx(q1, q2)
                qc.rz(np.pi/4, q2)
                #qc.t(q2)
                qc.cx(q1, q2)

        return qc

    for y in initial_layout:
        qc.h(y) 

    for (grp,dir) in [(True,True),(True,False),(False,True),(False,False)]:
        qc= apply_double_z(qc,n,dir,grp)

    qc= apply_all_rx(qc, initial_layout)
    

    return qc

def heisenberg(n ):

    qc = QuantumCircuit(n**2 )  # Create circuit with 16 qubits
    
   
    def get_qubit_index(row, col, n):
        
        return row * n + col

     
    def apply_double_y(qc, n, dir=bool, grp=bool):
        
        start = 1 if grp else 2

        for i in range(1,n-1):
            for j in range(start, n - 2, 2):
                row, col = (i, j) if dir else (j, i)
 

                q1 = get_qubit_index(row, col, n   )
                q2 = get_qubit_index(row, col + 1, n    ) if dir else get_qubit_index(row + 1, col, n   )
                
                qc.sdg([q1, q2])
                #qc.s([q1, q2])
                qc.h([q1, q2])
                qc.cx(q1, q2)
                qc.rz(np.pi/4, q2)
                #qc.t(q2)
                qc.cx(q1, q2)
                qc.h([q1, q2])
                qc.s([q1, q2])
                
        return qc
    
    def apply_double_x(qc, n, dir=bool, grp=bool):
         
        start = 1 if grp else 2

        for i in range(1,n-1):
            for j in range(start, n - 2, 2):
                row, col = (i, j) if dir else (j, i)
 

                q1 = get_qubit_index(row, col, n  )
                q2 = get_qubit_index(row, col + 1, n   ) if dir else get_qubit_index(row + 1, col, n  )

                
                qc.h([q1, q2])
                qc.cx(q1, q2)
                qc.rz(np.pi/4, q2)
                #qc.t(q2)
                qc.cx(q1, q2)
                qc.h([q1, q2])
                
        return qc
    def apply_double_z(qc, n, dir=bool, grp=bool):
        
        start = 1 if grp else 2

        for i in range(1,n-1):
            for j in range(start, n - 2, 2):
                row, col = (i, j) if dir else (j, i)
 

                q1 = get_qubit_index(row, col, n    )
                q2 = get_qubit_index(row, col + 1, n    ) if dir else get_qubit_index(row + 1, col, n )

                qc.cx(q1, q2)
                qc.rz(np.pi/4, q2)
                #qc.t(q2)
                qc.cx(q1, q2)

        return qc
     
    for (grp,dir) in [(True,True),(True,False),(False,True),(False,False)]:
        qc= apply_double_x(qc,n,dir,grp)
        qc= apply_double_y(qc,n,dir,grp)
        qc= apply_double_z(qc,n,dir,grp)
                
    return qc
 
 
 