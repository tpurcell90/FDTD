#include "Inputs.hpp"
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iterator>
#include <vector>

using namespace std;

/**
 * @brief Constructs a Input file manager
 * @details Takes in an input file name and converts real space parameters into FDTD space parameters
 *
 * @param fn Input File name
 */
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
    xPml_        = 0;
    yPml_        = 0;
    courant_     = IP.get<double>("CompCell.courant", 0.5);
    output_base_ = IP.get<string>("CompCell.output", "dtc_out");
    periodic_    = IP.get<bool>("CompCell.PBC", false);
    tMax_        = IP.get<double>("CompCell.tLim",100.0);
    pmlCalc_     = IP.get<bool>("CompCell.precalcPML",true);
    for (auto& iter : IP.get_child("SourceList"))
    {
        string p = iter.second.get<string>("profile");
        plsShape prof = string2prof(p);
        pol_ = iter.second.get<string>("pol");
        Polarization pols = string2pol(pol_);
        vector<double> fxn ={};
        if(prof == gaussian)
        {
            fxn.push_back(iter.second.get<double>("fcen") * courant_/res_);
            fxn.push_back(iter.second.get<double>("fwidth")*res_/courant_);
            fxn.push_back(iter.second.get<double>("cutoff"));
        }
        else if(prof == continuous)
        {
            fxn.push_back(iter.second.get<double>("fcen") * courant_/res_);
            fxn.push_back(iter.second.get<double>("cutoff"));
            fxn.push_back(courant_/res_);
        }
        else if(prof == ricker)
        {
            fxn.push_back(iter.second.get<double>("fcen") * courant_/res_);
            fxn.push_back(iter.second.get<double>("fwidth"));
            fxn.push_back(iter.second.get<double>("cutoff"));
        }
        // Now I need to write the funciton class for the source to take in
        double sz_x  = iter.second.get<double>("size_x");
        double sz_y  = iter.second.get<double>("size_y");
        double loc_x = iter.second.get<double>("loc_x");
        double loc_y = iter.second.get<double>("loc_y");
        // Make proper rounding function
        int x_min = find_pt(loc_x-sz_x/2.0+x_size_/2.0)+1;
        int x_max = find_pt(loc_x+sz_x/2.0+x_size_/2.0);
        int y_min = find_pt(loc_y-sz_y/2.0+y_size_/2.0);
        int y_max = find_pt(loc_y+sz_y/2.0+y_size_/2.0);
        for(int x = x_min; x <= x_max; x ++)
            for(int y = y_min; y <= y_max; y++)
            {
                vector<int> loc = {x,y};
                srcArr_.push_back(Source<complex<double>>(Pulse<complex<double>>(fxn,prof), pols, loc));
            }
    }
    for (auto& iter : IP.get_child("PML"))
    {
        string dd = iter.second.get<string>("direction");
        Direction d = string2dir(dd);
        int thickness = iter.second.get<double>("thickness");
        if (d == X)
        {
            xPml_ = find_pt(thickness);
            //pmlArr_.push_back(UPML<double>(thickness,d, 4.0, exp(-16), find_pt(y_size_),1.0/res_,1.0/res_,string2pol(pol_)));
            pmlArr_.push_back(UPML<complex<double>>(xPml_,d, 4.0, exp(-16), find_pt(y_size_) + 1,1.0/res_,1.0/res_,string2pol(pol_),pmlCalc_));
        }
        else if(d == Y)
        {
            yPml_ = find_pt(thickness);
            //pmlArr_.push_back(UPML<double>(thickness,d, 4.0, exp(-16), find_pt(x_size_),1.0/res_,1.0/res_,string2pol(pol_)));
            pmlArr_.push_back(UPML<complex<double>>(yPml_,d, 4.0, exp(-16), find_pt(x_size_) + 1,1.0/res_,1.0/res_,string2pol(pol_),pmlCalc_));
        }
        else if (d == Z)
            throw logic_error("While yes we could have a thrid dimension to run, I have yet to be implimented to do such a thing. So please accept this error as my sincerest appology.");
        else
            throw logic_error("I would appricate it if you stick to the standard X,Y,Z directions. While it's fun to invent new ones, it is very hard to do calculations if I don't even understand what dimension I am in. Please try again!");

    }
    vector<double> loc(2,0.0);
    vector<double> size = {x_size_,y_size_};
    vector<double> mat = {1.0};
    objArr_.push_back(Obj(block,mat,size,loc));
    for (auto& iter : IP.get_child("ObjectList"))
    {
        string sStr(iter.second.get<string>("shape"));
        Shape s = string2shape(sStr);
        vector<double> loc = {};
        vector<double> size = {};
        if(s == block)
        {
            size.push_back(iter.second.get<double>("size_x",0.0));
            size.push_back(iter.second.get<double>("size_y",0.0));
        }
        else if(s == sphere)
            size.push_back(iter.second.get<double>("radius",0.0));

        loc.push_back(iter.second.get<double>("loc_x",0.0));
        loc.push_back(iter.second.get<double>("loc_y",0.0));
        std::vector<double> mater = {iter.second.get<double>("eps",1.00)};

        boost::property_tree::ptree& pols = iter.second.get_child("pols");
        for(auto& iter2 : pols)
        {
            mater.push_back(iter2.second.get<double>("sigma"));
            mater.push_back(iter2.second.get<double>("gamma"));
            mater.push_back(iter2.second.get<double>("kappa"));
        }
        objArr_.push_back(Obj(s,mater,size,loc));
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
                dctArr_.push_back(Detector<complex<double>>(loc,t,out_name,p));
            }
    }
    //Copies the json data to a file to check for debugging
    //ofstream ss;
    //ss.open("output_data/inputs_check.json");
    //boost::property_tree::write_json(ss,IP);
    //ss.close();
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
/**
 * @brief Converts the string from an input file into a polarization
 *
 * @param p String type of polarization
 * @return The polarization
 */
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
/**
 * @brief Converts the string from an input file into a direction
 *
 * @param p String type of direction
 * @return The direction
 */
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
/**
 * @brief Converts the string from an input file into a shape
 *
 * @param p String type of shape
 * @return The shape
 */
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
/**
 * @brief Converts the string from an input file into a Detector Type
 *
 * @param p String type of Detector Type
 * @return The Detector Type
 */
dtcOutType programInputs::string2out(string t)
{
    if(t.compare("field") == 0)
        return field;
    else if (t.compare("flux") == 0)
        return flux;
    else
        throw logic_error("Output type undefined");
}
/**
 * @brief Converts the string from an input file into a Pulse Shape
 *
 * @param p String type of Pulse Shape
 * @return The Pulse Shape
 */
plsShape programInputs::string2prof(string p)
{
    if(p.compare("gaussian") == 0)
        return gaussian;
    else if(p.compare("continuous") == 0)
        return continuous;
    else if(p.compare("ricker") == 0)
        return ricker;
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

