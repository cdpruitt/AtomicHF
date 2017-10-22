#include <iostream>

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
        const double& minimum,
        const double& maximum,
        const unsigned int& numberOfPoints)
{
    double stepSize = log(maximum/minimum)/(numberOfPoints-1);
    double radius = 0;

    Eigen::VectorXd wavefunction = Eigen::VectorXd(numberOfPoints,1);

    if(n==1 && l==0)
    {
        for(unsigned int i=0; i<numberOfPoints; i++)
        {
            radius = log(minimum) + i*stepSize;
            wavefunction(i) = exp(-(Z/A_0)*radius)
                *(pow(Z/A_0,1.5)/(sqrt(PI)));
        }
    }

    if(n==2)
    {
        if(l==0)
        {
            for(unsigned int i=0; i<numberOfPoints; i++)
            {
                radius = log(minimum) + i*stepSize;
                wavefunction(i) = exp(-(Z/(2*A_0))*radius)
                    *(2-(Z/A_0)*radius)
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
