#include <iostream>
#include <string>
#include <vector>

#include "../include/Eigen/Sparse"
#include "../include/buildHydrogenicWF.h"
#include "../include/mathFunctions.h"

#include "TFile.h"
#include "TGraph.h"

using namespace std;

const double CONVERGENCE_LIMIT = 0.99;

int main(int argc, char** argv)
{
    if(argc<2)
    {
        cerr << "Error: incorrect number of input parameters." << endl;
        cerr << "Expected: \"calculateHF [Z of system]\"" << endl;
        return 1;
    }

    unsigned int Z = stoi(argv[1]); // Z of nucleus

    const double gridMinimum = pow(10,-7)/Z;
    const double gridMaximum = 25;
    const unsigned int numberOfPoints = 1000;

    const double stepSize = log(gridMaximum/gridMinimum)/(numberOfPoints-1);

    vector<double> grid;

    for(int i=0; i<numberOfPoints; i++)
    {
        grid.push_back(gridMinimum+i*stepSize);
    }

    vector<Wavefunction> wavefunctions; // starting point for calculation of all electron wavefunctions

    // hardcoding initial conditions for helium
    wavefunctions.push_back(Wavefunction(1, 0, 0, Z, grid));
    wavefunctions.push_back(Wavefunction(1, 0, 0, Z, grid));
    wavefunctions.push_back(Wavefunction(2, 0, 0, Z, grid));
    wavefunctions.push_back(Wavefunction(2, 0, 0, Z, grid));

    // perform HF calculation until wavefunctions show convergence within some prescribed limit

    double eigenvalue = 0;

    for(auto& wf : wavefunctions)
    {
        double KETerm = calculateKineticEnergyTerm(wf);
        double EPTerm = calculateExternalPotentialTerm(wf);
        double HartreeTerm = calculateHartreeTerm(wf, wavefunctions);
        //double FockTerm = calculateFockTerm(wf, wavefunctions);
    }

    TFile* outputFile = new TFile("output.root", "RECREATE");

    // plot initial wavefunctions
    for(unsigned int i=0; i<wavefunctions.size(); i++)
    {
        // extract info into vector
        vector<double> value;

        for(unsigned int j=0; j<wavefunctions[i].rows(); j++)
        {
            value.push_back(wavefunctions[i](j));
        }

        TGraph* graph = new TGraph(grid.size(), &grid[0], &value[0]);
        string name = "Electron " + to_string(i);
        graph->SetNameTitle(name.c_str(), name.c_str());
        graph->Write();
    }

    // plot KETerm
    for(unsigned int i=0; i<KETerms.size(); i++)
    {
        // extract info into vector
        vector<double> value;

        for(unsigned int j=0; j<KETerms[i].rows(); j++)
        {
            value.push_back(KETerms[i](j));
        }

        TGraph* graph = new TGraph(grid.size(), &grid[0], &value[0]);
        string name = "KETerm " + to_string(i);
        graph->SetNameTitle(name.c_str(), name.c_str());
        graph->Write();
    }

    outputFile->Close();

    return 0;
}
