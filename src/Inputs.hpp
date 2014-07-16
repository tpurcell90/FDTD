#ifndef FDTD_INPUTS
#define FDTD_INPUTS

#include <string>
#include <iostream>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
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
	double res;
	double courant;
	double t_lim;

	//Source Parameters Start with just Gaussian complicate later
	double f_cen;
	double f_width;
	double t_start;
	double cut_off;
	string pol;

	//Geometry Will add parameters as I write
	int n_struct;

	//Ouputs Will add parameters as I write
	
	//Imports the parameters from json file, converts to atomic units
	programInputs(std::string infile);

	//Removes comment lines from the json file
	void stripComments();

	//Convert parameters to FDTD units
	void toFDTD();
};
#endif

