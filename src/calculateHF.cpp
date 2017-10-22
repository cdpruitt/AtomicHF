#include <iostream>
#include <string>
#include <vector>

#include "../include/Eigen/Sparse"
#include "../include/buildHydrogenicWF.h"

#include "TFile.h"
#include "TGraph.h"

using namespace std;

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

    vector<Eigen::VectorXd> initialWFs; // starting point for calculation of all electron wavefunctions

    // hardcoding initial conditions for helium
    initialWFs.push_back(buildHydrogenicWF(1, 0, 0, Z, gridMinimum, gridMaximum, numberOfPoints));
    initialWFs.push_back(buildHydrogenicWF(1, 0, 0, Z, gridMinimum, gridMaximum, numberOfPoints));
    initialWFs.push_back(buildHydrogenicWF(2, 0, 0, Z, gridMinimum, gridMaximum, numberOfPoints));
    initialWFs.push_back(buildHydrogenicWF(2, 0, 0, Z, gridMinimum, gridMaximum, numberOfPoints));

    TFile* outputFile = new TFile("output.root", "RECREATE");

    // plot initial wavefunctions
    for(unsigned int i=0; i<initialWFs.size(); i++)
    {
        // extract info into vector
        vector<double> position;
        vector<double> value;

        for(unsigned int j=0; j<initialWFs[i].rows(); j++)
        {
            position.push_back(log(gridMinimum) + j*stepSize);
            value.push_back(initialWFs[i](j));
        }

        TGraph* graph = new TGraph(numberOfPoints, &position[0], &value[0]);
        string name = "Electron " + to_string(i);
        graph->SetNameTitle(name.c_str(), name.c_str());
        graph->Write();
    }

    outputFile->Close();

    return 0;
}
