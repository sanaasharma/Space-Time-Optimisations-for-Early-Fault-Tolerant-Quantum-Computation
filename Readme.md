Space-Time Optimisations for Early Fault-Tolerant Quantum Computation
Based on work described in https://doi.org/10.48550/arXiv.2511.08848


To run the experiments, use the script.py file. The default hamiltonian is set to ising_model. The experiments are done for data qubits ranging from n=4 to n=100. The initial layouts can be changed accordingly. The output file generates the average execution times. In the paper, benchmarking is done using minimum execution time for a set of iterations. Magic state factories can be changed in the script. To add magic state factories, edit the magic_state_routing_first_pass function in Rotations.py file and add output cells for distillation factories on the edge of the grid. 

- Uses Qiskit 2.0.0
- Ising, Fermi-Hubbard and Heisenberg hamiltonians used in the paper are given in Hamiltonians.py file.
- The initial layout can be changed using the find_initial_layout function in the script.py file. Use "add_paths_upto" to add more routing paths.
-  "generate_time" function can be used to estimate the total execution time assuming different processing times for distillation factories.  

 