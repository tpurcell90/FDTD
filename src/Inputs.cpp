#include "Inputs.hpp"
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <iterator>
#include <vector>
using namespace std;

programInputs::programInputs(std::string fn) : filename(fn)
{
	stripComments();
    boost::property_tree::ptree IP;
    boost::property_tree::json_parser::read_json(filename,IP);
    // Basic FDTD params
    procs   = IP.get<int>("CompCell.Procs",1);
    x_size  = IP.get<double>("CompCell.x_size",10.1);
    y_size  = IP.get<double>("CompCell.y_size",10.1);
    z_size  = IP.get<double>("CompCell.z_size",10.1);
    res     = IP.get<int>("CompCell.res", 10);
    t_pml   = IP.get<double>("CompCell.PML_Thick", 1.1);
    courant = IP.get<double>("CompCell.courant", 0.5);
    // Crating the srcArr
    //for(int ii = 0; ii < IP.get_child("SourceList").size(); ii ++)
    for (auto& iter : IP.get_child("SourceList"))
    {
        string prof = iter.second.get<string>("profile");
        /*boost::property_tree::ptree& size = iter.second.get_child("size");
        if (size.size() != 2)
          throw runtime_error("\"size\" needs two elements");
        vector<double> tmpsize;
        for (auto& i : size)
            tmpsize.push_back(i.second);
        */   //This will be fixed for vectors but for now hard coding in the vectors
        vector<double> size(2,0.0);
        vector<double> loc(2,0.0);
        vector<double> dir(2,0.0);
        size.push_back(iter.second.get<double>("size_x"));
        size.push_back(iter.second.get<double>("size_y"));
        loc.push_back(iter.second.get<double>("loc_x"));
        loc.push_back(iter.second.get<double>("loc_y"));
        dir.push_back(iter.second.get<double>("dir_x"));
        dir.push_back(iter.second.get<double>("dir_y"));
        //I need to fix this I know it's not good
        pol = iter.second.get<string>("pol");
        vector<double> fxn(4,100.8);
        fxn.push_back(iter.second.get<double>("fcen"));
        fxn.push_back(iter.second.get<double>("fwidth"));
        fxn.push_back(iter.second.get<double>("cutoff"));
        fxn.push_back(iter.second.get<double>("t_start"));
        // Now I need to write the funciton class for the source to take in
        // srcArr.push_back(Source(.....));
    }
    for (auto& iter : IP.get_child("ObjectList"))
    {
        string s(iter.second.get<string>("shape"));
        vector<double> loc(2,0.0);
        vector<double> size(2,0.0);
        size.push_back(iter.second.get<double>("size_x"));
        size.push_back(iter.second.get<double>("size_y"));
        loc.push_back(iter.second.get<double>("loc_x"));
        loc.push_back(iter.second.get<double>("loc_y"));
        string mater = iter.second.get<string>("material");
        objArr.push_back(Obj(sphere,getDielectricParams(mater),size,loc));
    }
    string pol	= IP.get<string>("SourceParam.Pol", "EZ");
    //Copies the json data to a file to check for debugging
    ofstream ss;
    ss.open("output_data/inputs_check.json");
    boost::property_tree::write_json(ss,IP);
    ss.close();
}

vector<double> programInputs::getDielectricParams(string mat)
{
    vector<double> dielec(1,1.0);
    if(mat.compare("Ag"))
    {
        dielec.push_back(1.00); // so on and so on I'll add all this later
    }
    return dielec;
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

