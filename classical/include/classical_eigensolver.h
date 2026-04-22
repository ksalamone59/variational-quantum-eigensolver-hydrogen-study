#ifndef GENERAL_EIGENSOLVER_H
#define GENERAL_EIGENSOLVER_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <eigen3/Eigen/Dense>

// Physical constants in SI units
constexpr double hbar = 1.054571817e-34;   // J·s
constexpr double m_e  = 9.10938356e-31;    // kg
constexpr double e_charge    = 1.602176634e-19;   // C
constexpr double eps0 = 8.8541878128e-12;   // F/m
constexpr double a0   = 5.29177210903e-11; // Bohr radius in meters
constexpr double hartree_to_ev = 27.2114; // Conversion factor from Hartree to eV

class eigenSolver
{
    private:
        Eigen::MatrixXd eigenvectors; 
        Eigen::VectorXd eigenvalues;
    public:
        eigenSolver() = default;
        ~eigenSolver() = default;
        eigenSolver(const Eigen::MatrixXd& hamiltonian) 
        {
            if(hamiltonian.rows() != hamiltonian.cols()) throw std::invalid_argument("Hamiltonian must be a square matrix.");
            if(!hamiltonian.isApprox(hamiltonian.adjoint())) throw std::invalid_argument("Hamiltonian must be Hermitian.");
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(hamiltonian);
            if(es.info() != Eigen::Success) throw std::runtime_error("Eigenvalue decomposition failed.");
            eigenvectors = es.eigenvectors();
            eigenvalues = es.eigenvalues();
            std::transform(eigenvalues.begin(), eigenvalues.end(), eigenvalues.begin(), [](double val) { return val / e_charge; }); // Convert to eV
        }
        const Eigen::VectorXd &get_eigenvalues() const
        {
            return eigenvalues;
        }
        const Eigen::MatrixXd &get_eigenvectors() const
        {
            return eigenvectors;
        }
};

void output_wavefunctions_to_gnu(const eigenSolver &es, const int N, const int l, const double dr);
Eigen::MatrixXd construct_hydrogen_hamiltonian(const int N, const double r_max, const double dr, const int l);

#endif 