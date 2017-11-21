#include <iostream>
#include <string>
#include <vector>

#include "../include/Eigen/Sparse"
#include "../include/Eigen/Dense"
#include "../include/Eigen/Eigenvalues"
#include "../include/buildHydrogenicWF.h"
#include "../include/mathFunctions.h"
#include "../include/Wavefunction.h"

#include "TFile.h"
#include "TGraph.h"

using namespace std;

const double CONVERGENCE_LIMIT = 0.9999;

void plot(Eigen::VectorXd grid, Eigen::VectorXd values, string name)
{
    // extract info into std::vectors for plotting
    vector<double> dummyGrid;
    vector<double> dummyValues;

    for(unsigned int i=0; i<values.size(); i++)
    {
        dummyGrid.push_back(exp(grid(i)));
        dummyValues.push_back(values(i));
    }

    TGraph* graph = new TGraph(dummyGrid.size(), &dummyGrid[0], &dummyValues[0]);
    graph->SetNameTitle(name.c_str(), name.c_str());
    graph->Write();
}

double integrate(Eigen::VectorXd vector1, Eigen::VectorXd vector2, double stepSize)
{
    if(vector1.size()!=vector2.size())
    {
        cerr << "Error: can't integrate vectors of different sizes."
            << "Vector 1 had size " << vector1.size() << "; vector 2 had size "
            << vector2.size() << endl;
        return 0;
    }

    double result = 0;

    for(unsigned int i=0; i<vector1.size(); i++)
    {
        result += vector1(i)*vector2(i)*stepSize;
    }

    return result;
}

Eigen::VectorXd normalize(Eigen::VectorXd values, double stepSize)
{
    Eigen::VectorXd result(values.size(),1);

    double norm = 0;

    for(unsigned int i=0; i<values.size(); i++)
    {
        norm += values(i)*values(i)*stepSize;
    }

    if(norm<0)
    {
        cerr << "Error: norm of wavefunction cannot be negative." << endl;
        exit(1);
    }

    norm = sqrt(norm);

    for(unsigned int i=0; i<result.size(); i++)
    {
        result(i) = values(i)/norm;
    }

    return result;
}

bool testConvergence(Eigen::VectorXd wf1, Eigen::VectorXd wf2, double stepSize)
{
    if(wf1.size()!=wf2.size())
    {
        cerr << "Error: can't test convergence of two vectors of different sizes."
            << "Vector 1 had size " << wf1.size() << "; vector 2 had size "
            << wf2.size() << endl;
        exit(1);
    }

    double result = integrate(wf1, wf2, stepSize);

    return (result>CONVERGENCE_LIMIT && result<(2-CONVERGENCE_LIMIT));
}

int main(int argc, char** argv)
{
    if(argc<2)
    {
        cerr << "Error: incorrect number of input parameters." << endl;
        cerr << "Expected: \"calculateHF [Z of system]\"" << endl;
        return 1;
    }

    TFile* outputFile = new TFile("output.root", "RECREATE");

    // define constants
    unsigned int Z = stoi(argv[1]); // Z of nucleus

    const double gridMinimum = log(pow(10,-7)/Z);
    const double gridMaximum = log(25);
    const unsigned int numberOfPoints = 1000;

    const double stepSize = (gridMaximum-gridMinimum)/(numberOfPoints-1);

    // initialize grid
    Eigen::VectorXd grid = Eigen::VectorXd(numberOfPoints,1);

    for(unsigned int i=0; i<numberOfPoints; i++)
    {
        grid(i) = gridMinimum+i*stepSize;
    }

    // initialize wavefunctions
    vector<Wavefunction> wavefunctions;
    vector<Wavefunction> nlWavefunctions;

    if(Z==2)
    {
        wavefunctions.push_back(Wavefunction(1, 0, 0, Z, grid));
        wavefunctions.push_back(Wavefunction(1, 0, 0, Z, grid));

        nlWavefunctions.push_back(Wavefunction(1, 0, 0, Z, grid));
    }

    if(Z==10)
    {
        wavefunctions.push_back(Wavefunction(1, 0, 0, Z, grid));
        wavefunctions.push_back(Wavefunction(1, 0, 0, Z, grid));

        wavefunctions.push_back(Wavefunction(2, 0, 0, Z, grid));
        wavefunctions.push_back(Wavefunction(2, 0, 0, Z, grid));

        wavefunctions.push_back(Wavefunction(2, 1, -1, Z, grid));
        wavefunctions.push_back(Wavefunction(2, 1, -1, Z, grid));

        wavefunctions.push_back(Wavefunction(2, 1, 0, Z, grid));
        wavefunctions.push_back(Wavefunction(2, 1, 0, Z, grid));

        wavefunctions.push_back(Wavefunction(2, 1, 1, Z, grid));
        wavefunctions.push_back(Wavefunction(2, 1, 1, Z, grid));

        nlWavefunctions.push_back(Wavefunction(1, 0, 0, Z, grid));
        nlWavefunctions.push_back(Wavefunction(2, 0, 0, Z, grid));
        nlWavefunctions.push_back(Wavefunction(2, 1, -1, Z, grid));
        nlWavefunctions.push_back(Wavefunction(2, 1, 0, Z, grid));
        nlWavefunctions.push_back(Wavefunction(2, 1, 1, Z, grid));
    }

    unsigned int counter = 0;
    bool allConverged = false;

    // iterate wavefunctions until they converge
    while(counter < 10)
    {
        // plot initial wavefunctions
        for(unsigned int i=0; i<wavefunctions.size(); i += 2)
        {
            string name = "Electron" + to_string(i);
            plot(wavefunctions[i].grid, wavefunctions[i].values, name);
        }

        vector<Eigen::VectorXd> newWavefunctions;

        for(auto& wf : wavefunctions)
        {
            Eigen::VectorXd KETerm = calculateKineticEnergyTerm(wf);
            Eigen::VectorXd EPTerm = calculateExternalPotentialTerm(wf);
            Eigen::VectorXd HartreeTerm = calculateHartreeTerm(wf, wavefunctions);
            Eigen::VectorXd FockTerm = calculateFockTerm(wf, nlWavefunctions);

            //plot(wf.grid, KETerm, "KETerm");
            //plot(wf.grid, EPTerm, "EPTerm");
            //plot(wf.grid, HartreeTerm, "HartreeTerm");
            //plot(wf.grid, FockTerm, "FockTerm");

            double KEValue = integrate(KETerm, wf.values, stepSize);
            double EPValue = integrate(EPTerm, wf.values, stepSize);
            double HValue = integrate(HartreeTerm, wf.values, stepSize);
            double FValue = integrate(FockTerm, wf.values, stepSize);

            wf.eigenvalue = KEValue + EPValue + HValue - FValue;

            //cout << "KEValue = " << KEValue << ", EPValue = " << EPValue
            //     << ", HValue = " << HValue << ", FValue = " << FValue << endl;

            if(wf.eigenvalue>0)
            {
                wf.eigenvalue = 0;
            }

            Eigen::VectorXd unnormalizedNewWF = solveTridiagonalMatrix(wf, wavefunctions, FockTerm, wf.eigenvalue);

            newWavefunctions.push_back(normalize(unnormalizedNewWF, stepSize));

            //wf.converged = testConvergence(wf.values, newWavefunctions.back(), stepSize);
        }

        allConverged = true;

        for(auto&wf : wavefunctions)
        {
            allConverged = (allConverged && wf.converged);

        }

        for(unsigned int i=0; i<wavefunctions.size(); i++)
        {
            wavefunctions[i].values = newWavefunctions[i];
        }


        cout << endl << "Iteration " << counter << " results:" << endl;

        for(unsigned int i=0; i<wavefunctions.size(); i++)
        {
            if(i%2==0)
            {
                cout << "n = " << wavefunctions[i].n
                    << ", l = " << wavefunctions[i].l
                    << " eigenvalue = " << wavefunctions[i].eigenvalue << endl;
            }
        }

        counter++;
    }

    outputFile->Close();

    return 0;
}
