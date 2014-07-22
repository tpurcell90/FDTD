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
        string p = iter.second.get<string>("profile");
        ProfType prof = string2prof(p);
        /*boost::property_tree::ptree& size = iter.second.get_child("size");
        if (size.size() != 2)
          throw runtime_error("\"size\" needs two elements");
        vector<double> tmpsize;
        for (auto& i : size)
            tmpsize.push_back(i.second);
        */   //This will be fixed for vectors but for now hard coding in the vectors
        pol = iter.second.get<string>("pol");
        Polarization pols = string2pol(pol);
        vector<double> fxn ={};
        fxn.push_back(iter.second.get<double>("fcen"));
        fxn.push_back(iter.second.get<double>("fwidth"));
        fxn.push_back(iter.second.get<double>("cutoff"));
        fxn.push_back(iter.second.get<double>("t_start"));
        // Now I need to write the funciton class for the source to take in
        double sz_x  = iter.second.get<double>("size_x");
        double sz_y  = iter.second.get<double>("size_y");
        double loc_x = iter.second.get<double>("loc_x");
        double loc_y = iter.second.get<double>("loc_y");
        // Make proper rounding function
        int x_min = find_pt(loc_x-sz_x/2.0+x_size/2.0);
        int x_max = find_pt(loc_x+sz_x/2.0+x_size/2.0);
        int y_min = find_pt(loc_y-sz_y/2.0+y_size/2.0);
        int y_max = find_pt(loc_y+sz_y/2.0+y_size/2.0);
        for(int x = x_min; x <= x_max; x ++)
            for(int y = y_min; y <= y_max; y++)
            {
                vector<int> loc(2,x);
                loc[1] = y;
                srcArr.push_back(Source<double>(Pulse<double>(fxn,prof), pols, loc));
            }
    }
    vector<double> loc(2,0.0);
    vector<double> size = {x_size,y_size};
    vector<double> mat(1,1.0);
    objArr.push_back(Obj(block,mat,size,loc));
    for (auto& iter : IP.get_child("ObjectList"))
    {
        string sStr(iter.second.get<string>("shape"));
        Shape s = string2shape(sStr);
        vector<double> loc = {};
        vector<double> size = {};
        size.push_back(iter.second.get<double>("size_x"));
        size.push_back(iter.second.get<double>("size_y"));
        loc.push_back(iter.second.get<double>("loc_x"));
        loc.push_back(iter.second.get<double>("loc_y"));
        string mater = iter.second.get<string>("material");
        objArr.push_back(Obj(s,getDielectricParams(mater),size,loc));
    }
    for (auto& iter : IP.get_child("DetectorList"))
    {
        string tt(iter.second.get<string>("type"));
        OupuptsData t = string2out(tt);
        double sz_x = iter.second.get<double>("size_x");
        double sz_y = iter.second.get<double>("size_y");
        double loc_x = iter.second.get<double>("loc_x");
        double loc_y = iter.second.get<double>("loc_y");
        // Make proper rounding function
        int x_min = find_pt(loc_x-sz_x/2.0+x_size/2.0);
        int x_max = find_pt(loc_x+sz_x/2.0+x_size/2.0);
        int y_min = find_pt(loc_y-sz_y/2.0+y_size/2.0);
        int y_max = find_pt(loc_y+sz_y/2.0+y_size/2.0);
        for(int x = x_min; x <= x_max; x ++)
            for(int y = y_min; y <= y_max; y++)
            {
                vector<int> loc = {x,y};
                dctArr.push_back(Detector<double>(loc,t));
            }
    }
    string pol	= IP.get<string>("SourceParam.Pol", "EZ");
    //Copies the json data to a file to check for debugging
    //ofstream ss;
    //ss.open("output_data/inputs_check.json");
    //boost::property_tree::write_json(ss,IP);
    //ss.close();
}

vector<double> programInputs::getDielectricParams(string mat)
{
    vector<double> dielec = {};
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
Polarization programInputs::string2pol(string p)
{
    if(p.compare("Ez") == 0)
        return EZ;
    else if(p.compare("Ey") == 0)
        return EY;
    else if(p.compare("Ex") == 0)
        return EX;
    else if(p.compare("Hx") == 0)
        return HX;
    else if(p.compare("Hz") == 0)
        return HZ;
    else if(p.compare("Hx") == 0)
        return HY;
    else
        return EZ; //Throw an error but first need to look up how to do that
}
Shape programInputs::string2shape(string s)
{
    if(s.compare("sphere") == 0)
        return sphere;
    else if (s.compare("block") == 0)
        return block;
    else if (s.compare("cone") == 0)
        return cone;
    else if (s.compare("ellipse") == 0)
        return ellipse;
    else if (s.compare("cylinder") == 0)
        return cylinder;
    else
        return block; //Throw an error I know
}
OupuptsData programInputs::string2out(string t)
{
    if(t.compare("field") == 0)
        return field;
    else if (t.compare("flux") == 0)
        return flux;
    else
        return field; //Yes, yes error, I know
}
ProfType programInputs::string2prof(string p)
{
    if(p.compare("gaussian") == 0)
        return gaussian;
    else if(p.compare("continuous") == 0)
        return continuous;
    else
        return gaussian; // I should be an error again
}

int programInputs::find_pt(double pt)
{
    return int(pt*res); // make it better
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

