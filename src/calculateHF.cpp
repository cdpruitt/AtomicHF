#include <iostream>
#include <string>
#include <vector>

#include "../include/Eigen/Sparse"
#include "../include/buildHydrogenicWF.h"
#include "../include/mathFunctions.h"
#include "../include/Wavefunction.h"

#include "TFile.h"
#include "TGraph.h"

using namespace std;

const double CONVERGENCE_LIMIT = 0.99;

void plot(Eigen::VectorXd grid, Eigen::VectorXd values, string name)
{
    // extract info into std::vectors for plotting
    vector<double> dummyGrid;
    vector<double> dummyValues;

    for(unsigned int i=0; i<values.size(); i++)
    {
        dummyGrid.push_back(grid(i));
        dummyValues.push_back(values(i));
    }

    TGraph* graph = new TGraph(dummyGrid.size(), &dummyGrid[0], &dummyValues[0]);
    graph->SetNameTitle(name.c_str(), name.c_str());
    graph->Write();
}

double integrate(Eigen::VectorXd values, double stepSize)
{
    double result = 0;

    for(unsigned int i=0; i<values.size(); i++)
    {
        result += values(i)*stepSize;
    }

    return result;
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

    unsigned int Z = stoi(argv[1]); // Z of nucleus

    const double gridMinimum = pow(10,-7)/Z;
    const double gridMaximum = 25;
    const unsigned int numberOfPoints = 1000;

    const double stepSize = log(gridMaximum/gridMinimum)/(numberOfPoints-1);

    // initialize grid
    Eigen::VectorXd grid = Eigen::VectorXd(numberOfPoints,1);

    for(int i=0; i<numberOfPoints; i++)
    {
        grid(i) = gridMinimum+i*stepSize;
    }

    vector<Wavefunction> wavefunctions; // starting point for calculation of all electron wavefunctions

    // hardcoding initial conditions for helium
    wavefunctions.push_back(Wavefunction(1, 0, 0, Z, grid));
    wavefunctions.push_back(Wavefunction(1, 0, 0, Z, grid));

    for(auto& wf : wavefunctions)
    {
        Eigen::VectorXd KETerm = calculateKineticEnergyTerm(wf);
        Eigen::VectorXd EPTerm = calculateExternalPotentialTerm(wf);
        Eigen::VectorXd HartreeTerm = calculateHartreeTerm(wf, wavefunctions);
        //double FockTerm = calculateFockTerm(wf, wavefunctions);

        plot(wf.grid, KETerm, "KETerm");
        plot(wf.grid, EPTerm, "EPTerm");
        plot(wf.grid, HartreeTerm, "HartreeTerm");

        double eigenvalue = integrate(KETerm, stepSize) - integrate(EPTerm, stepSize);// - integrate(HartreeTerm, grid);
        cout << eigenvalue << endl;
    }

    // plot initial wavefunctions
    for(unsigned int i=0; i<wavefunctions.size(); i++)
    {
        string name = "Electron" + to_string(i);
        plot(wavefunctions[i].grid, wavefunctions[i].values, name);
    }

    outputFile->Close();

    return 0;
}
