#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include <vector>

#include "../include/Eigen/Sparse"
#include "../include/Wavefunction.h"

Eigen::VectorXd calculateKineticEnergyTerm(const Wavefunction& wf);
Eigen::VectorXd calculateSpinOrbitTerm(const Wavefunction& wf);
Eigen::VectorXd calculateExternalPotentialTerm(const Wavefunction& wf);
Eigen::VectorXd calculateHartreeTerm(const Wavefunction& wf, const std::vector<Wavefunction>& wavefunctions);
Eigen::VectorXd calculateFockTerm(const Wavefunction& wf, const std::vector<Wavefunction>& wavefunctions);
Eigen::VectorXd solveTridiagonalMatrix(const Wavefunction& wf, const std::vector<Wavefunction>& wavefunctions, const Eigen::VectorXd& Y, const double eigenvalue);

#endif /* MATH_FUNCTIONS_H */
