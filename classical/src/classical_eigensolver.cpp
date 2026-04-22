#include "classical_eigensolver.h"

// H = T + V = -0.5 hbar^2/2m_e * d^2/dr^2 -e^2/(4pie0r) + hbar^2l(l+1)/(2m_e r^2)
Eigen::MatrixXd construct_hydrogen_hamiltonian(const int N, const double r_max, const double dr, const int l)
{
    if(l < 0) throw std::invalid_argument("Angular momentum quantum number l must be non-negative.");
    if(N <= 0) throw std::invalid_argument("Number of grid points N must be positive.");
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(N, N);
    const double kinematic_term =  (hbar * hbar) / (2.0 * m_e);
    const double coulomb_term = (e_charge * e_charge) / (4.0 * M_PI * eps0);
    for(int i=1;i<N-1;i++)
    {
        double r = i*dr;
        double V = -coulomb_term/r + (hbar*hbar*l*(l+1))/(2.0*m_e*r*r);
        H(i,i-1) = -kinematic_term / (dr * dr);
        H(i,i)   =  2.0 * kinematic_term / (dr * dr) + V;
        H(i,i+1) = -kinematic_term / (dr * dr); 
    }
    H(0,0) = 1.0;   H(0,1) = 0.0;   H(1,0) = 0.0;      // r=0
    H(N-1,N-1) = 1.0; H(N-1,N-2) = 0.0; H(N-2,N-1) = 0.0; // r=r_max
    return H;
}

void output_wavefunctions_to_gnu(const eigenSolver &es, const int N, const int l, const double dr)
{
    const auto &eigenvectors = es.get_eigenvectors();
    std::ofstream file("../../data/wavefunctions.dat", std::ios_base::trunc);
    if (!file.is_open()) throw std::runtime_error("Failed to open output file.");
    file << "# r(m) u_n(r) (radial wavefunctions, normalized)\n";
    const int n_states = std::min(5, static_cast<int>(eigenvectors.cols()));
    std::vector<double> norms(n_states);
    for (int n = 0; n < n_states; n++)
    {
        norms[n] = std::sqrt(eigenvectors.col(n).squaredNorm() * dr);
        if (norms[n] == 0.0) throw std::runtime_error("Zero normalization encountered.");
    }
    for (int i = 0; i < N; i++)
    {
        double r = i * dr;
        file << r;
        for (int n = 0; n < n_states; n++)
        {
            double u = eigenvectors(i, n) / norms[n];
            file << " " << u;
        }
        file << "\n";
    }
    file.close();
}

int main() 
{
    const int N = 1000; // Number of grid points
    const double r_max = 50.0 * a0;
    const double dr = r_max / (N-1);
    const int l = 0;  // angular momentum, >= 0

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    eigenSolver es(construct_hydrogen_hamiltonian(N, r_max, dr, l));
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Eigenvalue decomposition took " << elapsed.count() << " seconds.\n";
    const Eigen::VectorXd &eigenvalues = es.get_eigenvalues();

    std::cout << "First energy level: \n";
    // Make this loop longer to show more energy levels/where it breaks down based on our N and r_max
    for(int i=0;i<1;i++)
    {
        int n = i + l + 1; // For fixed l, allowed n starts at l (since l = 0,1,...,n-1)
        double exact = -13.6 / (n * n); // Exact energy levels for hydrogen in eV
        std::cout << "Energy level " << i << ". E = " << eigenvalues(i) << " eV. Exact: " << exact << " eV. Error = " << std::fabs(eigenvalues(i) - exact) << " eV\n";
        std::cout << "In Ha: " << eigenvalues(i) / hartree_to_ev << " Ha. Exact: " << exact / hartree_to_ev << " Ha. Error = " << std::fabs(eigenvalues(i) - exact) / hartree_to_ev << " Ha";
        std::cout << " Percent error: " << 100 * std::fabs(eigenvalues(i) - exact) / std::fabs(exact) << "%\n";
    }
    output_wavefunctions_to_gnu(es, N, l, dr);

    return 0;
}