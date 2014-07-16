#include "Inputs.hpp"

using namespace std;

programInputs::programInputs(std::string fn) : filename(fn)
{
	stripComments();
	boost::property_tree::ptree IP;
	boost::property_tree::json_parser::read_json(filename,IP);

	procs   = IP.get<int>("CompCell.Procs",1);
        x_size	= IP.get<double>("CompCell.Xsize", 10.00);
        y_size  = IP.get<double>("CompCell.Ysize", 10.00);
        res  	= IP.get<double>("CompCell.Res", 30.00);
        courant	= IP.get<double>("CompCell.Courant", 0.50);
        t_lim	= IP.get<double>("CompCell.T_lim", 10.00);

        f_cen	= IP.get<double>("SourceParam.F_cen", 0.200);
        f_width	= IP.get<double>("SourceParam.F_width", 0.100);
        t_start	= IP.get<double>("SourceParam.T_start", 0.00);
        cut_off	= IP.get<double>("SourceParam.Cut_off", 10.00);
	pol	= IP.get<string>("SourceParam.Pol", "EZ");

	n_struct= IP.get<int>("Geometry.N_struct",1);
	toFDTD();
  //Copies the json data to a file to check for debugging
  ofstream ss;
  ss.open("output_data/inputs_check.json");
  boost::property_tree::write_json(ss,IP);
  ss.close();
}

void programInputs::stripComments()
{
	//Open input and output file
	string newfn = "output_data/" + filename;
	fstream inputfile;
	inputfile.open(filename);
	ofstream inputcopy;
	inputcopy.open(newfn);

	//search for '//', delete everything following, print remainder to new file
	string line;
	int found, found2;
	while (getline(inputfile,line))
	{
		found  = line.find('/');
		found2 = line.find('/', found+1);
		if (found != line.npos && found2 == found+1)
			inputcopy << line.erase(found, line.length()) << endl;
		else
			inputcopy << line << endl;
	}
	inputcopy.close();

	//update filename;
	filename = newfn;
}

void programInputs::toFDTD()
{ /*
	epsilon     /= 27.21;
	eaffin      /= 27.21;
	workfxn     /= 27.21;

	omega       = sqrt(72.0 * epsilon/pow(requil,2)*(1.0/m_C60));
	lbarrieritn = 2*M_PI*M_PI/pow(lbarrier,2);
	rbarrieritn = 2*M_PI*M_PI/pow(rbarrier,2);
   */
   // Will add after Geometry
}

