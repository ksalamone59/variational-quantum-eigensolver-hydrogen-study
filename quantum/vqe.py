import time 
import numpy as np
from scipy.optimize import minimize
from qiskit.circuit.library import efficient_su2
from qiskit.quantum_info import Operator, SparsePauliOp, Statevector

'''
Cost function for minimization
'''
def cost(theta, ansatz, H_pauli):
    qc = ansatz.assign_parameters(theta)
    sv = Statevector.from_instruction(qc)
    return sv.expectation_value(H_pauli).real

EPSILON=1e-8
R_MAX = 6
hartree_to_ev = 27.2114
truth_ground_state_energy = -0.5 # Ha
'''
num_qubits = ceil(log2(N_STATES)) to represent the state space of the system
'''
N_STATES = 16
'''
    Approximating the Hamiltonian
'''
def create_hamiltonian(dr = 0.1, l = 0) -> np.array:
    H = np.zeros((N_STATES, N_STATES))
    for i in range(N_STATES):
        r = (i+1)*dr + EPSILON
        V = -1/r + l*(l+1)/(2*r**2)
        H[i, i] = 1.0/dr**2 + V
        if i > 0:
            H[i, i-1]= -0.5/dr**2
        if i < N_STATES-1:
            H[i, i+1] = -0.5/dr**2
    return H  

'''
Comparisons to theory
'''
def compare_to_theory(min_energy):
    print(f"Minimum energy from VQE: {min_energy:.6f} Ha")
    print(f"Theory: {truth_ground_state_energy:.6f} Ha. Error: {(min_energy - truth_ground_state_energy):.4f} Ha")
    print(f"Percent error: {100 * abs(min_energy - truth_ground_state_energy) / abs(truth_ground_state_energy):.2f}%")

if __name__ == "__main__":
    np.random.seed(42)
    H = create_hamiltonian(dr=R_MAX/N_STATES)

    H_operator = Operator(H)
    H_pauli = SparsePauliOp.from_operator(H_operator)

    ansatz = efficient_su2(num_qubits=int(np.ceil(np.log2(N_STATES))), reps=3, entanglement='full')

    x0 = np.random.uniform(0, 2*np.pi, ansatz.num_parameters)
    start_time = time.time()
    result = minimize(cost, x0, args=(ansatz, H_pauli), method='L-BFGS-B', bounds=[(0, 2*np.pi)] * ansatz.num_parameters, options={'maxiter': 5000, 'maxfun': 2e4, "ftol" : 1e-6, "gtol" : 1e-5})
    end_time = time.time()
    print(f"VQE optimization for {int(np.ceil(np.log2(N_STATES)))} qubits took {end_time - start_time:.2f} seconds.")
    print(f"Converged: {result.success}. Message: {result.message}")
    print(f"Function evaluations: {result.nfev}, Iterations: {result.nit}")
    compare_to_theory(result.fun)