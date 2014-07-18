#ifndef FDTD_INPUTS
#define FDTD_INPUTS

#include "Obj.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
class programInputs
{
public:
	std::string filename;
	//Computational cell parameters and run parameters
	int procs;
	double x_size;
	double y_size;
    double z_size;
	double res;
	double courant;
	double t_pml;
    //vector<Source<double>> srcArr();
    vector<Obj> objArr;
	//Source Parameters Start with just Gaussian complicate later
    
	//Geometry Will add parameters as I write
	int n_struct;

	//Ouputs Will add parameters as I write
	
	//Imports the parameters from json file, converts to atomic units
	programInputs(std::string infile);

    vector<double> getDielectricParams(string mat);

	//Removes comment lines from the json file
	void stripComments();

	//Convert parameters to FDTD units
	void toFDTD();
};
#endif

