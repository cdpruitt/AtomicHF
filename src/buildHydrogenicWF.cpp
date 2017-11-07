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
            wavefunction(i) = exp(-(Z/A_0)*grid(i))
                *(pow(Z/A_0,1.5)/(sqrt(PI)));
        }
    }

    if(n==2)
    {
        if(l==0)
        {
            for(unsigned int i=0; i<grid.size(); i++)
            {
                wavefunction(i) = exp(-(Z/(2*A_0))*grid(i))
                    *(2-(Z/A_0)*grid(i))
                    *pow(Z/A_0,1.5)
                    /(4*sqrt(2*PI));
            }
        }

        if(l==1)
        {
            if(m==1)
            {
            }

            if(m==0)
            {

            }

            if(m==-1)
            {

            }
        }
    }

    return wavefunction;
}
