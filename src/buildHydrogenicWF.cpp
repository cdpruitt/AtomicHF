#include <iostream>
#include <vector>

#include "../include/Eigen/Sparse"

#include "../include/PhysicalConstants.h"
#include "../include/MathConstants.h"

#include <math.h>

using namespace std;

Eigen::VectorXd buildHydrogenicWF(
        const unsigned int& n,
        const unsigned int& l,
        const int& m,
        const unsigned int& Z,
        const Eigen::VectorXd& grid)
{
    Eigen::VectorXd wavefunction = Eigen::VectorXd(grid.size(),1);

    if(n==1 && l==0)
    {
        for(unsigned int i=0; i<grid.size(); i++)
        {
            double r = exp(grid[i]);
            wavefunction(i) = pow(r, 1.5)*exp(-r*(Z/A_0))
                *2*pow(Z/A_0,1.5);
        }
    }

    if(n==2)
    {
        if(l==0)
        {
            for(unsigned int i=0; i<grid.size(); i++)
            {
                double r = exp(grid[i]);
                wavefunction(i) = pow(r, 1.5)*exp(-r*(Z/2*A_0))
                    *(2-r*(Z/A_0))
                    *pow(Z/A_0,1.5)
                    /(2*sqrt(2));
            }
        }

        if(l==1)
        {
            for(unsigned int i=0; i<grid.size(); i++)
            {
                double r = exp(grid[i]);
                wavefunction(i) = pow(r, 1.5)*exp(-r*(Z/2*A_0))
                    *r*(Z/A_0)
                    *pow(Z/A_0,1.5)
                    /(2*sqrt(6));
            }
        }
    }

    return wavefunction;
}
