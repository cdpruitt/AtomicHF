#ifndef BUILD_HYDROGENIC_WF
#define BUILD_HYDROGENIC_WF

#include <vector>

#include "../include/Eigen/Sparse"

Eigen::VectorXd buildHydrogenicWF(
        const unsigned int& n,
        const unsigned int& l,
        const int& m,
        const unsigned int& Z,
        const std::vector<double>& grid);

#endif /* BUILD_HYDROGENIC_WF */
