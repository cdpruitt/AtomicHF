#include "../include/Wavefunction.h"
#include "../include/buildHydrogenicWF.h"

using namespace std;

Wavefunction::Wavefunction(
        const unsigned int& _n,
        const unsigned int& _l,
        const unsigned int& _m,
        const unsigned int& _Z,
        const Eigen::VectorXd& _grid
        )
    : n(_n), l(_l), m(_m), Z(_Z), grid(_grid)
{
    values = buildHydrogenicWF(n, l, m, Z, grid);
    converged = false;
}
