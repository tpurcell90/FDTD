#include "Inputs.hpp"
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <iterator>
#include <vector>
using namespace std;

programInputs::programInputs(std::string fn) : filename_(fn)
{
    boost::property_tree::ptree IP;
    boost::property_tree::json_parser::read_json(filename_,IP);
    // Basic FDTD params
    procs_       = IP.get<int>("CompCell.Procs",1);
    x_size_      = IP.get<double>("CompCell.x_size",10.1);
    y_size_      = IP.get<double>("CompCell.y_size",10.1);
    z_size_      = IP.get<double>("CompCell.z_size",10.1);
    res_         = IP.get<int>("CompCell.res", 10);
    t_pml_       = IP.get<double>("CompCell.PML_Thick", 0);
    courant_     = IP.get<double>("CompCell.courant", 0.5);
    output_base_ = IP.get<string>("CompCell.output", "dtc_out");
    
    for (auto& iter : IP.get_child("SourceList"))
    {
        string p = iter.second.get<string>("profile");
        plsShape prof = string2prof(p);
        /*boost::property_tree::ptree& size = iter.second.get_child("size");
        if (size.size() != 2)
          throw runtime_error("\"size\" needs two elements");
        vector<double> tmpsize;
        for (auto& i : size)
            tmpsize.push_back(i.second);
        */   //This will be fixed for vectors but for now hard coding in the vectors
        pol_ = iter.second.get<string>("pol");
        Polarization pols = string2pol(pol_);
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
        int x_min = find_pt(loc_x-sz_x/2.0+x_size_/2.0);
        int x_max = find_pt(loc_x+sz_x/2.0+x_size_/2.0);
        int y_min = find_pt(loc_y-sz_y/2.0+y_size_/2.0);
        int y_max = find_pt(loc_y+sz_y/2.0+y_size_/2.0);
        for(int x = x_min; x <= x_max; x ++)
            for(int y = y_min; y <= y_max; y++)
            {
                vector<int> loc = {x,y};
                srcArr_.push_back(Source<double>(Pulse<double>(fxn,prof), pols, loc));
            }
    }
    vector<double> loc(2,0.0);
    vector<double> size = {x_size_,y_size_};
    vector<double> mat(1,1.0);
    objArr_.push_back(Obj(block,mat,size,loc));
    for (auto& iter : IP.get_child("PML"))
    {
        string dd = iter.second.get<string>("direction");
        Direction d = string2dir(dd);
        int thickness = iter.second.get<int>("thickness");
        if (d == X)
        {
            pmlArr_.push_back(UPML<double>(thickness,d, 2.0, 1.0e-7, int(res_*y_size_),1.0/res_,1.0/res_,string2pol(pol_)));
        }
        else if(d == Y)
        {
            pmlArr_.push_back(UPML<double>(thickness,d, 2.0, 1.0e-7, int(res_*x_size_),1.0/res_,1.0/res_,string2pol(pol_)));
        }
        else if (d == Z)
            throw logic_error("While yes we could have a thrid dimension to run, I have yet to be implimented to do such a thing. So please accept this error as my sincerest appology.");
        else
            throw logic_error("I would appricate it if you stick to the standard X,Y,Z directions. While it's fun to invent new ones, it is very hard to do calculations if I don't even understand what dimension I am in. Please try again!");

    }
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
        objArr_.push_back(Obj(s,getDielectricParams(mater),size,loc));
    }
    int ii = 0;
    int jj = 0;
    for (auto& iter : IP.get_child("DetectorList"))
    {
        string tt(iter.second.get<string>("type"));
        dtcOutType t = string2out(tt);
        string pp(iter.second.get<string>("pol"));
        Polarization p = string2pol(pp);
        string out_name;
        if(t == field)
        {
            out_name = "output_data/" + output_base_ + "_field_" + to_string(ii) + ".dat";
            ii ++;
        }
        else
        {
            out_name = "output_data/" + output_base_ + "_flux_" + to_string(jj) + ".dat";
            jj++;
        }
        double sz_x = iter.second.get<double>("size_x");
        double sz_y = iter.second.get<double>("size_y");
        double loc_x = iter.second.get<double>("loc_x");
        double loc_y = iter.second.get<double>("loc_y");
        // Make proper rounding function
        int x_min = find_pt(loc_x-sz_x/2.0+x_size_/2.0);
        int x_max = find_pt(loc_x+sz_x/2.0+x_size_/2.0);
        int y_min = find_pt(loc_y-sz_y/2.0+y_size_/2.0);
        int y_max = find_pt(loc_y+sz_y/2.0+y_size_/2.0);
        for(int x = x_min; x <= x_max; x ++)
            for(int y = y_min; y <= y_max; y++)
            {
                vector<int> loc = {x,y};
                dctArr_.push_back(Detector<double>(loc,t,out_name,p));
            }
    }
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
	string newfn = "output_data/" + filename_;
	fstream inputfile;
	inputfile.open(filename_);
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
	filename_ = newfn;
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
    else if(p.compare("Hy") == 0)
        return HY;
    else
        throw logic_error("Polarization undefined");

}
Direction programInputs::string2dir(string dir)
{
    if((dir.compare("x") == 0) || (dir.compare("X") == 0))
        return X;
    else if((dir.compare("y") == 0) || (dir.compare("Y") == 0))
        return Y;
    else if ((dir.compare("z") == 0) || (dir.compare("Z") == 0))
        return Z;
    else
        throw logic_error("I would appricate it if you stick to the standard X,Y,Z directions. While it's fun to invent new ones, it is very hard to do calculations if I don't even understand what dimension I am in. Please try again!");

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
        throw logic_error("Shape undefined");

}
dtcOutType programInputs::string2out(string t)
{
    if(t.compare("field") == 0)
        return field;
    else if (t.compare("flux") == 0)
        return flux;
    else
        throw logic_error("Output type undefined");
}
plsShape programInputs::string2prof(string p)
{
    if(p.compare("gaussian") == 0)
        return gaussian;
    else if(p.compare("continuous") == 0)
        return continuous;
    else
        throw logic_error("Pulse sahpe undefined");
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

