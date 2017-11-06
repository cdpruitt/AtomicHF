#include "../include/Eigen/Sparse"

using namespace std;

double calculateKineticEnergyTerm(const Wavefunction& wf)
{
    Eigen::VectorXd KETerm(wf.size(), 1);

    vector<double> r;

    for(auto& p : grid)
    {
        r.push_back(exp(p));
    }

    double stepSize = log((r[r.size()-1])/r[0])/(r.size()-1);

    KETerm(0) = -2*wf[0]/pow(r[0],2) +
        2*wf[0]/(r[1]*r[0]);
    KETerm(0) /= pow(stepSize,2);
    KETerm(0) += pow(l+0.5,2)/pow(r[0],2);
    KETerm(0) *= -0.5;

    for(unsigned int i=1; i<r.size()-1; i++)
    {
        KETerm(i) = -2*wf[i]/pow((r[i]),2) +
                wf[i+1]/(r[i]*r[i+1]) +
                wf[i-1]/(r[i]*r[i-1]);
        KETerm(i) /= pow(stepSize,2);

        KETerm(i) += pow(l+0.5,2)/pow(r[i],2);
        KETerm(i) *= -0.5;

    }

    unsigned int n = r.size()-1;

    KETerm(n) = -2*wf[n]/pow(r[n],2) +
        2*wf[n]/(r[n]*r[n-1]);
    KETerm(n) /= pow(stepSize,2);

    KETerm(n) += pow(l+0.5,2)/pow(r[n],2);
    KETerm(n) *= -0.5;

    return result;
}

double calculateExternalPotentialTerm(const Wavefunction& wf)
{
    Eigen::VectorXd EPTerm(wf.size(), 1);

    vector<double> r;

    for(auto& p : grid)
    {
        r.push_back(exp(p));
    }

    for(unsigned int i=0; i<r.size(); i++)
    {
        EPTerm(i) = Z/r[i];
    }

    return result;
}

double calculateHartreeTerm(const Wavefunction& wf, const vector<Wavefunctions> wavefunctions)
{
    Eigen::VectorXd HTerm(wf.size(), 1);

    vector<double> r;

    for(auto& p : grid)
    {
        r.push_back(exp(p));
    }

    double stepSize = log((r[r.size()-1])/r[0])/(r.size()-1);

    vector<double> V;
    vector<double> W;

    V.push_back(0);
    W.push_back(0);

    for(unsigned int i=0; i<grid.size(); i++)
    {
        V.push_back(V.back() + pow(wf(i),2));
        W.push_back(W.back() + pow(wf(i),2)/r[i]);
    }

    for(unsigned int i=0; i<grid.size(); i++)
    {
        HTerm(i) = stepSize*(V[i+1]/r[i] + W[grid.size()]-W[i+1]);
    }

    for(unsigned int n_prime=0; n_prime<=n; n_prime++)
    {
        for(unsigned int l_prime=0; l_prime<=l; l_prime++)
        {
            2*(2*l_prime+1)*
        }
    }

    return result;
}

double calculateFockTerm(const Wavefunction& wf, const vector<Wavefunctions> wavefunctions)
{
    return result;
}
