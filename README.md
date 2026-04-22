# Overview
A project that benchmarks and compares two approaches to handling the energy levels of the hydrogen atom:

1. A pure C++ Eigen based linear decomposition of the Hamiltonian
2. Using qiskit in Python to perform a VQE on the Hamiltonian 

The analytic hamiltonian for the radial component of the wavefunction $u(r)=rR(r)$ may be represented as 

$$
 \hat{H}=\left(-\dfrac{\hbar^{2}}{2m_{e}}\dfrac{d^{2}}{dr^{2}} - \dfrac{e^{2}}{4\pi\varepsilon_{0}r} + \dfrac{\hbar^{2}\ell\left(\ell+1\right)}{2m_{e}r^{2}}\right)
$$

in SI units and 

$$
\hat{H}=\left(-\dfrac{1}{2}\dfrac{d^{2}}{dr^{2}} - \dfrac{1}{r} + \dfrac{\ell\left(\ell+1\right)}{2r^{2}}\right)
$$

in atomic units. A brief note that the C++ eigensolver uses SI units, and the python VQE uses atomic units (for consistency in Qiskit implementations). The ground state energy is 

$$
E_{gs}=-13.6 \text{ eV (SI Units)} = -0.5 \text{ Ha (Atomic Units)}.
$$

## Methodology
While exact diagonalization (C++ approach) scales as O($N^{3}$) where N is the size of the discretized Hamiltonian basis, VQE embeds this problem into exponentially fewer qubits. The number of qubits chosen takes the form $\lceil\log_{2}\left(N_{states}\right)\rceil $. As written the code only accepts exact powers of 2 as the number of states. This can easily be expanded upon by padding the basis to get to the next power of 2. Optimization is performed classically on the expectation value of the Hamiltonian.

## Key Results
- Classical Solver ($N=1000$)
    - Error: $~0.02\%$
    - Runtime: $~0.3$ seconds
    - Deterministic result 
- VQE (4 qubits, $N=16$)
    - Mean error: $4.13\pm0.95\%$
    - Best error: $3.40\%$
    - Runtime: $13.83\pm7.29$ seconds
    - Best runtime: $5.11$ seconds

Variability in the VQE approach arises from sensitivity to initial parameter choice, as well as the non-convex optimization landscape. In a more expressive ansatz, the barren plateau effect may become more significant due to the expressibility and depth of the ansatz.

## Interpretation 
The classical approach uses a significantly larger basis ($N$=1000) than the quantum implementation ($N$=16). This is intentional: the goal is not a like-for-like comparison, but to evaluate how well a low-qubit VQE can approximate the ground state energy under constrained Hilbert space representations and a limited ansatz.

The key distinction between these approaches is as $N$ grows, the classical approach becomes increasingly computationally expensive, and VQE compresses the Hilbert space. VQE encodes the system into a logarithmic number of qubits with respect to the discretized basis size, but shifts computational complexity into repeated circuit evaluations and classical optimization. For example: for $N=2048$, classically you have to use eigendecomposition on a $2048\times2048$ matrix. However, using VQE, one requires 11 qubits to represent the discretized basis under binary encoding.

## Future Work
Future extensions could incorporate shot-based estimation and noise models to study the robustness of VQE under realistic NISQ hardware constraints.

## How to Run the Code
After installing all dependencies, you can simply run `make` from the main directory of the repository. If you wish to run only the classical approach, run `make classic`. For just the quantum portion: `make quantum`. For just the studying of the phase space of the VQE: `make evaluation`. 

# Repository Layout
```plaintext 
├── quantum/
├── classical/
├── data/
├── Plots/
├── Plotting/
├── Makefile
├── README.md
```

- `quantum/`: Directory that performs a VQE in Python and qiskit to solve for the ground state energy of the Hydrogen atom. Contains two files:
    - `vqe.py`: file that handles the VQE solving, as well as creating the Hamiltonian for the VQE
    - `error_estimation.py`: file that benchmarks and evaluates top achievable VQE model performance for varying number of qubits, N and max radius
- `classical/`: Directory that houses C++ code to perform classical eigendecomposition of Hamiltonian to solve for energy levels of Hydrogen atom. This file also outputs the data for the reconstructed wavefunctions of the Hydrogen atom and stores them in `data/wavefunctions.dat`
    - `src/` Houses the source code to run the eigensolver 
    - `include/` Houses the header file to define the eigenSolver class, various constants and other methods
- `data/`: Directory that contains data for output wavefunctions from classical approach
- `Makefile`: Top-level Makefile that, upon running "make", will run all the code and produce the output comparisons
- `Plots/`: Directory that contains output of `error_estimation.py`, including various quantitative evaluations of the VQE project on Hydrogen. As well, output wave functions from the classical approach (done in the `Plotting/` directory)
- `Plotting/`: usage of [gnuplot-latex-utils](https://github.com/ksalamone59/gnuplot_latex_utils) as a submodule to create publication-quality plots on the fly through gnuplot and LaTeX. Please see original documentation for more information.
    - `wave_functions/`: Directory that plots the extracted wave functions from the classical approach vs the analytic curves.
    - `pdfs/`: Where the output pdf is stored/

# Requirements 
- Classical:
    - C++ >= C++17
    - CMake
    - Eigen3
- Quantum:
    - numpy
    - matplotlib
    - scipy (scipy.optimize)
    - qiskit
- Dependencies outlined in [gnuplot-latex-utils](https://github.com/ksalamone59/gnuplot_latex_utils) including 
    - Gnuplot 
    - LaTeX

# Results 
1. Heatmap of minimum achievable VQE error as a function of qubit count and maximum radial cutoff r_max. The optimal configuration (9 qubits, r_max = 9.0) is highlighted.
![](heatmap.png)
2. Comparison of numerically computed eigenfunctions (C++ diagonalization) with analytic hydrogen wavefunctions for the ground and first excited states.
![](waveFunctions.png)