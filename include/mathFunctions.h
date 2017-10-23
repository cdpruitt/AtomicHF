#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include <vector>

#include "../include/Eigen/Sparse"

Eigen::VectorXd calculateKineticEnergyTerm(const std::vector<double>& grid, const Eigen::VectorXd& wf, const unsigned int& l);
Eigen::VectorXd calculateSpinOrbitTerm(const Eigen::VectorXd& wf);
Eigen::VectorXd calculateExternalPotentialTerm(const Eigen::VectorXd& wf);
Eigen::VectorXd calculateHartreeTerm(const Eigen::VectorXd& wf);

#endif /* MATH_FUNCTIONS_H */
