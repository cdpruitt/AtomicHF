#ifndef BUILD_HYDROGENIC_WF
#define BUILD_HYDROGENIC_WF

#include "../include/Eigen/Sparse"

Eigen::VectorXd buildHydrogenicWF(
        const unsigned int& n,
        const unsigned int& l,
        const int& m,
        const unsigned int& Z,
        const double& minimum,
        const double& maximum,
        const unsigned int& numberOfPoints);

#endif /* BUILD_HYDROGENIC_WF */
