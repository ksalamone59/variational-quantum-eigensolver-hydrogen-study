'''
Helper file to quantify some facts about out VQE model
Example: what is the "best" attainable ground state energy as a function of max radius,
for a fixed N?
'''

import numpy as np 
import matplotlib.pyplot as plt
import time
from scipy.optimize import minimize
from qiskit.circuit.library import efficient_su2
from qiskit.quantum_info import Operator, SparsePauliOp, Statevector

TRUTH_GS_ENERGY = -0.5 # Ha

'''
Hamiltonian construction and exact diagonalization to find ground state energy for given N and r_max
'''
def find_gs_energy(N_STATES, r_max, l=0):
    dr = r_max/N_STATES 
    rs = np.array([(i+1)*dr for i in range(N_STATES)])
    H = np.zeros((N_STATES, N_STATES))
    for i in range(N_STATES):
        r = rs[i]
        V = -1/r + l*(l+1)/(2*r**2)
        H[i, i] = 1.0/dr**2 + V
        if i > 0:
            H[i, i-1] = -0.5/dr**2
        if i < N_STATES-1:
            H[i, i+1] = -0.5/dr**2
    eigenvalues, _ = np.linalg.eigh(H)
    return eigenvalues[0], H

'''
Fix N and scan across r_max to find "best" attainable energy given that N
'''
def plot_max_energy_vs_r(N_STATES, r_max, l=0):
    energies = [] 
    r = 1.
    dr = 1.
    min_energy = float('inf')
    min_energy_r = float('inf')
    while r < r_max:
        energy, _ = find_gs_energy(N_STATES, r, l)
        energies.append(energy)
        if energy < min_energy:
            min_energy = energy
            min_energy_r = r
        r += dr
    plt.figure(figsize=(8, 5))
    plt.axhline(TRUTH_GS_ENERGY, color='red', linestyle='--', label='Exact Ground State Energy')
    plt.plot(np.arange(1., r_max, dr), energies, marker='o', label='Estimated Ground State Energy')
    plt.xlabel('r_max')
    plt.ylabel('Ground State Energy (Ha)')
    plt.title(f'Energy Estimation vs. r_max for N_STATES={N_STATES}')
    # Put a vertical line at minima 
    plt.axvline(min_energy_r, color='green', linestyle='--', label=f'Min Energy at r_max={min_energy_r}')
    plt.legend()
    plt.tight_layout()
    plt.savefig("../Plots/energy_estimation.png")

'''
Fix r_max and scan across N to find "best" attainable energy given that r
'''
def plot_max_energy_vs_N(N_STATES, r_max, l=0):
    energies = [] 
    errors = []
    N_values = [2**i for i in range (1, N_STATES)]
    for N in N_values:
        energy, _ = find_gs_energy(N, r_max, l)
        energies.append(energy)
        errors.append(abs(energy - TRUTH_GS_ENERGY))
    plt.figure(figsize=(8, 5))
    plt.axhline(TRUTH_GS_ENERGY, color='red', linestyle='--', label='Exact Ground State Energy')
    plt.plot(N_values, energies, marker='o', label='Estimated Ground State Energy')
    plt.xlabel('Number of States (N)')
    plt.ylabel('Ground State Energy (Ha)')
    plt.title(f'Energy Estimation vs. Number of States for r_max={r_max}')
    plt.legend()
    plt.tight_layout()
    plt.savefig("../Plots/energy_estimation_vs_N.png")

    plt.figure(figsize=(8, 5))
    plt.plot(N_values, errors, marker='o', label='Absolute Error')
    plt.xlabel('Number of States (N)')
    plt.ylabel('Absolute Error (Ha)')
    plt.title(f'Error vs. Number of States for r_max={r_max}')
    plt.legend()
    plt.tight_layout()
    plt.savefig("../Plots/error_vs_N.png")

'''
Discretize phase space of N and r_max; for each N, find "best" attainable error across r_max
'''
def optimize_error_for_fixed_N(N_STATES, r_max, l=0):
    n_values = [2**i for i in range (1, N_STATES)]
    r_range = np.arange(1., r_max, 1.)
    best = []
    for n in n_values:
        best_err = float('inf')
        for r in r_range:
            energy, _ = find_gs_energy(n, r, l)
            err = abs(energy - TRUTH_GS_ENERGY)
            if err < best_err:
                best_err = err
        best.append(best_err)
    plt.figure(figsize=(8, 5))
    plt.plot(n_values, best, marker='o', label='Best Error for Each N')
    plt.xlabel('Number of States (N)')
    plt.ylabel('Best Absolute Error (Ha)')
    plt.title('Best Error vs. Number of States')
    plt.legend()
    plt.tight_layout() 
    plt.savefig("../Plots/best_error_vs_N.png")

    plt.figure(figsize=(8, 5))
    plt.plot(np.log2(n_values), best, marker='o', label='Best Error vs Number of Qubits')
    plt.xlabel('log2(N)')
    plt.ylabel('Best Absolute Error (Ha)')
    plt.title('Best Error vs. Number of Qubits')
    plt.legend()
    plt.tight_layout()
    plt.savefig("../Plots/best_error_vs_N_qubits.png")

'''
Duplicated code from vqe.py for separation of concerns
Cost function for minimization
'''
def cost(theta, ansatz, H):
    qc = ansatz.assign_parameters(theta)
    sv = Statevector.from_instruction(qc).data
    return np.real(np.vdot(sv, H @ sv))
def run_vqe(N_states, H):
    scale = np.max(np.abs(H))
    scale = scale if scale > 0 else 1.0
    Hs = H / scale

    n_qubits = int(np.ceil(np.log2(N_states)))
    n_reps = 1
    
    ansatz = efficient_su2(num_qubits=n_qubits, reps=n_reps, entanglement="full")
    rng = np.random.default_rng()
    count = 0
    best_res = None
    best_energy = float('inf')
    times = []
    while count < 10:
        x0 = rng.uniform(0, 2*np.pi, ansatz.num_parameters)
        t_start = time.time()
        result = minimize(cost, x0, args=(ansatz, Hs), method='L-BFGS-B', bounds=[(0, 2*np.pi)] * ansatz.num_parameters, options={'maxiter': 500, 'maxfun': 1e4, "ftol" : 1e-4, "gtol" : 1e-4})
        t_end = time.time()
        if t_end > t_start:
            times.append(t_end - t_start)
        vqe_energy = result.fun * scale
        if vqe_energy < best_energy: # Idea is minimum energy is our GS energy
            best_energy = vqe_energy
            best_res = result
        count += 1
    avg_time = sum(times) / len(times)
    stddev_time = np.std(times)
    print(f"[VQE] iters={best_res.nit}, params={len(best_res.x)}, energy={best_energy:.6f}, time={avg_time:.2f} +/- {stddev_time:.2f} seconds")
    return best_energy, avg_time
'''
Scan across phase space of N and r_max to find "best" attainable energy and, thus, "best" parameters
Also, creates heatmap of ratio discretization / variational error across phase space
'''
def plot_error_heatmap(N_STATES, r_max, l=0):
    print("Generating error heatmap and ratios; this may take a while...")
    N_values = np.array([2**i for i in range(1, N_STATES)])
    r_values = np.arange(1., r_max, 1.)
    error_matrix = np.zeros((len(N_values), len(r_values)))
    vqe_err_matrix = np.zeros((len(N_values), len(r_values)))
    times = np.zeros((len(N_values), len(r_values)))
    min_error = float('inf')
    min_N, min_r = None, None

    max_n_for_vqe = len(N_values) # To keep runtime manageable; increased to max for plots in README
    for i, N in enumerate(N_values):
        # For speed, we try to reuse VQE parameters across r for a given N
        for j, r in enumerate(r_values):
            print(f"Evaluating i={i}, N={N}, r_max={r}...")
            energy, H = find_gs_energy(N, r, l)
            disc_err = abs(energy - TRUTH_GS_ENERGY)
            error_matrix[i, j] = disc_err

            if i < max_n_for_vqe: # Skip VQE for large N to save time
                H_operator = Operator(H)
                vqe_energy, vqe_time = run_vqe(N, H_operator)
                print(f"Disc energy = {energy:.6f}")
                var_err = abs(vqe_energy - energy)
                vqe_err_matrix[i, j] = var_err
                times[i, j] = vqe_time

            if disc_err < min_error:
                min_error = disc_err
                min_N, min_r = N, r
    log_error_matrix = np.log10(error_matrix + 1e-12)
    vqe_err_matrix = np.log10(vqe_err_matrix + 1e-12)
    _, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(
        log_error_matrix,
        origin='lower',
        aspect='auto',
        cmap='viridis',
        extent=[
            r_values[0], r_values[-1],
            np.log2(N_values[0]), np.log2(N_values[-1])
        ]
    )
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("log10 Absolute Error (Ha)")
    ax.set_xlabel("r_max")
    ax.set_ylabel("Number of Qubits (log2 N)")
    ax.scatter(
        min_r,
        np.log2(min_N),
        color='red',
        edgecolors='white',
        s=80,
        label=f"min err = {min_error:.2e}\n(qubits={int(np.log2(min_N))}, r={min_r})"
    )
    ax.legend()
    ax.set_title("Discretization Error Landscape (log scale)")
    plt.tight_layout()
    plt.savefig("../Plots/error_heatmap.png", dpi=300)
    with open("../Plotting/heatmap/heatmap.dat", "w") as f:
        for i, N in enumerate(N_values):
            for j, r in enumerate(r_values):
                f.write(f"{np.log2(N)} {r} {error_matrix[i, j]}\n")
    with open("../Plotting/vqe_err/vqe_err.dat", "w") as f:
        for i, N in enumerate(N_values):
            for j, r in enumerate(r_values):
                f.write(f"{np.log2(N)} {r} {vqe_err_matrix[i, j]}\n")
    with open("../Plotting/time_scan/time_scan.dat", "w") as f:
        for i, N in enumerate(N_values):
            for j, r in enumerate(r_values):
                f.write(f"{np.log2(N)} {r} {times[i, j]}\n")
    print("\n=== BEST PARAMS FROM HEATMAP ===")
    print(f"log2(N) = {int(np.log2(min_N))}")
    print(f"N       = {min_N}")
    print(f"r_max   = {min_r}")
    print(f"error   = {min_error:.6e}")

if __name__ == "__main__":
    N_STATES = 16 # To match our simulation
    plot_max_energy_vs_r(N_STATES=N_STATES, r_max=100)
    plot_max_energy_vs_N(N_STATES=10, r_max=20)
    optimize_error_for_fixed_N(N_STATES=10, r_max=20)
    plot_error_heatmap(N_STATES=6, r_max=20)