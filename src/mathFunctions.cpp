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

double calculateKineticEnergyTerm(const Wavefunction& wf, unsigned int i)
{
    double result = 0;

    // identify relevant grid values
    double r_prev;
    double r_0;
    double r_next;

    double wf_prev;
    double wf_0;
    double wf_next;

    if(i==0)
    {
        r_prev = exp(wf.grid(i+1));
        wf_prev = exp(wf.values(i+1));
    }

    else
    {
        r_prev = exp(wf.grid(i-1));
        wf_prev = exp(wf.values(i-1));
    }

    if(i==wf.grid.size()-1)
    {
        r_next = exp(wf.grid(i-1));
        wf_next = exp(wf.values(i-1));
    }

    else
    {
        r_next = exp(wf.grid(i+1));
        wf_next = exp(wf.values(i+1));
    }

    r_0 = exp(wf.grid(i));
    wf_0 = exp(wf.values(i));

    double stepSize = wf.grid(1)-wf.grid(0);

    result = -2*wf_0/pow(r_0,2) +
        wf_next/(r_0*r_next) +
        wf_prev/(r_0*r_prev);

    result /= pow(stepSize,2);

    result += pow(wf.l+0.5,2)/pow(r_0,2);
    result *= -0.5;

    return result;
}

Eigen::VectorXd calculateKineticEnergyTerm(const Wavefunction& wf)
{
    Eigen::VectorXd KETerm(wf.grid.size(), 1);
    Eigen::VectorXd r(wf.grid.size(), 1);

    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        r(i) = exp(wf.grid(i));
    }

    double stepSize = log((r(r.size()-1))/r(0))/(r.size()-1);

    KETerm(0) = -2*wf.values(0)/pow(r(0),2) +
        2*wf.values(0)/(r(1)*r(0));

    KETerm(0) /= pow(stepSize,2);

    KETerm(0) += pow(wf.l+0.5,2)/pow(r(0),2);

    KETerm(0) *= -0.5;

    for(unsigned int i=1; i<r.size()-1; i++)
    {
        KETerm(i) = -2*wf.values(i)/pow((r(i)),2) +
                wf.values(i+1)/(r(i)*r(i+1)) +
                wf.values(i-1)/(r(i)*r(i-1));
        KETerm(i) /= pow(stepSize,2);

        KETerm(i) += pow(wf.l+0.5,2)/pow(r(i),2);
        KETerm(i) *= -0.5;

    }

    unsigned int n = r.size()-1;

    KETerm(n) = -2*wf.values(n)/pow(r(n),2) +
        2*wf.values(n)/(r(n)*r(n-1));
    KETerm(n) /= pow(stepSize,2);

    KETerm(n) += pow(wf.l+0.5,2)/pow(r(n),2);
    KETerm(n) *= -0.5;

    return KETerm;
}

Eigen::VectorXd calculateExternalPotentialTerm(const Wavefunction& wf)
{
    Eigen::VectorXd EPTerm(wf.grid.size(), 1);

    Eigen::VectorXd r(wf.grid.size(), 1);

    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        r(i) = exp(wf.grid(i));
    }

    for(unsigned int i=0; i<r.size(); i++)
    {
        EPTerm(i) = wf.Z/r[i];
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

    double stepSize = log((r(r.size()-1))/r(0))/(r.size()-1);

    vector<double> V;
    vector<double> W;

    V.push_back(0);
    W.push_back(0);

    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        double G = 0;

        for(auto& wf_prime : wavefunctions)
        {
            G += 2*(2*wf_prime.l+1)*pow(wf.values(i),2);
        }

        V.push_back(V.back() + G);
        W.push_back(W.back() + G/r(i));
    }

    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        HTerm(i) = stepSize*(V[i+1]/r(i) + W[wf.grid.size()]-W[i+1]);
    }

    return HTerm;
}

/*double calculateFockTerm(const Wavefunction& wf, const vector<Wavefunctions> wavefunctions)
{
}*/
