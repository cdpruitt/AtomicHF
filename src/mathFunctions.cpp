#include "../include/Eigen/Sparse"
#include "../include/Wavefunction.h"

using namespace std;

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
        V.push_back(V.back() + pow(wf.values(i),2));
        W.push_back(W.back() + pow(wf.values(i),2)/r(i));
    }

    for(unsigned int i=0; i<wf.grid.size(); i++)
    {
        HTerm(i) = stepSize*(V[i+1]/r(i) + W[wf.grid.size()]-W[i+1]);
    }

    for(unsigned int n_prime=0; n_prime<=wf.n; n_prime++)
    {
        for(unsigned int l_prime=0; l_prime<=wf.l; l_prime++)
        {
            //HTerm(i) *= 2*(2*l_prime+1)
        }
    }

    return HTerm;
}

/*double calculateFockTerm(const Wavefunction& wf, const vector<Wavefunctions> wavefunctions)
{
}*/
