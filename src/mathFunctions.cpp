#include <iostream>

#include "../include/Eigen/Sparse"
#include "../include/Wavefunction.h"

using namespace std;

/*Eigen::MatrixXd constructHamiltonian(Wavefunction wf, Eigen::VectorXd grid)
{
    Eigen::MatrixXd Hamiltonian = Eigen::MatrixXd(grid.size(), grid.size());

    for(unsigned int i=0; i<grid.size(); i++)
    {
        Hamiltonian(i,i) += wf.Z/exp(grid(i)); // external potential term
        Hamiltonian(i,i) +=  ; // Hartree term

        Hamiltonian(i,i-1) += ; // kinetic term
        Hamiltonian(i,i) += ; // kinetic term
        Hamiltonian(i,i+1) += ; // kinetic term
    }
}*/

Eigen::VectorXd calculateKineticEnergyTerm(const Wavefunction& wf)
{
    Eigen::VectorXd KETerm(wf.grid.size(), 1);
    Eigen::VectorXd grid_r(wf.grid.size(), 1);

    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        grid_r(i) = exp(wf.grid(i));
    }

    double stepSize = wf.grid(1)-wf.grid(0);

    double r_prev, r_0, r_next;
    double wf_prev, wf_0, wf_next;

    for(unsigned int i=1; i<grid_r.size()-1; i++)
    {
        r_prev = grid_r(i-1);
        wf_prev = wf.values(i-1);

        r_0 = grid_r(i);
        wf_0 = wf.values(i);

        r_next = grid_r(i+1);
        wf_next = wf.values(i+1);

        //cout << "i = " << i << ", wf_prev = " << wf_prev
        //    << ", wf_0 = " << wf_0 << ", wf_next = " << wf_next << endl;

        KETerm(i) = -2*wf_0/pow(r_0,2) +
            wf_next/(r_0*r_next) +
            wf_prev/(r_0*r_prev);

        KETerm(i) /= pow(stepSize,2);

        KETerm(i) -= wf_0*pow(wf.l+0.5,2)/pow(r_0,2);
        KETerm(i) *= -0.5;
    }

    KETerm(0) = KETerm(1);
    KETerm(grid_r.size()-1) = KETerm(grid_r.size()-2);

    return KETerm;
}

Eigen::VectorXd calculateExternalPotentialTerm(const Wavefunction& wf)
{
    Eigen::VectorXd EPTerm(wf.grid.size(), 1);

    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        double r = exp(wf.grid(i));
        EPTerm(i) = wf.values(i)*(-(wf.Z/r));
    }

    return EPTerm;
}

Eigen::VectorXd calculateHartreeTerm(const Wavefunction& wf, const vector<Wavefunction>& wavefunctions)
{
    Eigen::VectorXd HTerm(wf.values.size(), 1);
    Eigen::VectorXd r(wf.grid.size(), 1);

    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        r(i) = exp(wf.grid(i));
    }

    double stepSize = wf.grid(1)-wf.grid(0);

    vector<double> V;
    vector<double> W;

    V.push_back(0);
    W.push_back(0);

    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        double G = 0;

        for(auto& wf_prime : wavefunctions)
        {
            G += pow(wf_prime.values(i),2);
        }

        V.push_back(V.back() + G);
        W.push_back(W.back() + G/r(i));
    }

    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        HTerm(i) = stepSize*(V[i+1]/r(i) + W[wf.grid.size()]-W[i+1]);
        HTerm(i) = HTerm(i)*wf.values(i);
    }

    return HTerm;
}

Eigen::VectorXd calculateFockTerm(const Wavefunction& wf, const vector<Wavefunction>& wavefunctions)
{
    Eigen::VectorXd FTerm(wf.values.size(), 1);
    Eigen::VectorXd r(wf.grid.size(), 1);

    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        r(i) = exp(wf.grid(i));
    }

    double stepSize = wf.grid(1)-wf.grid(0);

    vector<double> V;
    vector<double> W;

    V.push_back(0);
    W.push_back(0);

    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        double G = 0;

        for(auto& wf_prime : wavefunctions)
        {
            G += wf_prime.values(i)*wf.values(i);
        }

        V.push_back(V.back() + G);
        //V.push_back(V.back() + G*r(i));
        W.push_back(W.back() + G/r(i));
    }

    //for(unsigned int L=0; L<
    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        FTerm(i) = stepSize*(V[i+1]/r(i) + (W[wf.grid.size()]-W[i+1]));
        //FTerm(i) = stepSize*(V[i+1]/r(i) + r(i)*(W[wf.grid.size()]-W[i+1]));
        FTerm(i) = FTerm(i)*wf.values(i);
    }

    return FTerm;
}
