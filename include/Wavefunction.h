#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "../include/Eigen/Sparse"

class Wavefunction
{
    public:
        Wavefunction(
                const unsigned int& n,
                const unsigned int& l, 
                const unsigned int& m, 
                const unsigned int& Z, 
                const Eigen::VectorXd& grid
                );

        Eigen::VectorXd grid;
        Eigen::VectorXd values;

        unsigned int n; // principal quantum number
        unsigned int l; // angular momentum quantum number
        unsigned int m; // angular momentum z-projection
        unsigned int Z; // charge of nucleus for the problem of interest

        bool converged; // flag indicating whether this wavefunction's solution has converged
};

#endif /* WAVEFUNCTION_H */
