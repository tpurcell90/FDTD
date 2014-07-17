#include "FDTDField.hpp"

// #include <assert.h>
// #include <iomanip>
// #include <iostream>
#include <memory>
// #include <random>
// #include <stdexcept>
// #include <string>
// #include <vector>
// #include <complex>

using namespace std;

FDTDField::FDTDField(programInputs *IP)
{
    //Define cell parameters

    if (IP)
    {
        res = IP->res;
        dx = 1.0/res;
        dy = 1.0/res;
        dt = 1.0/(IP->courant*res);
        nx = IP->x_size/dx;
        ny = IP->y_size/dy;
        nt = IP->t_lim/dt;

        // The source array should go here, I need to read up on Vectors before I do that
        physGrid = make_shared<Grid2D<double>>(nx,ny,dx,dy);

        // Create the Grids. Do I need a null constructor for the set that I disregard?
        if(IP->pol.compare("Hz") == 0)
        {
            Ex = make_shared<Grid2D<double>>(nx-1,ny,dx,dy);
            Ey = make_shared<Grid2D<double>>(nx,ny-1,dx,dy);
            Hz = make_shared<Grid2D<double>>(nx-1,ny-1,dx,dy);
           // Hx = null;
           // Hy = null;
           // Ez = null;
        }
        else
        {
            Hx = make_shared<Grid2D<double>>(nx-1,ny,dx,dy);
            Hy = make_shared<Grid2D<double>>(nx,ny-1,dx,dy);
            Ez = make_shared<Grid2D<double>>(nx-1,ny-1,dx,dy);
            //Ex = null;
            //Ey = null;
            //Hz = null;
        }
    }
}

void FDTDField::initializeGrid(programInputs *IP)
{
    // This will plot the geometry structures onto the physical grid still thinking about this one I was thinking about how best to do this
}

void FDTDField::ouputField()
{
    // Work on a plan for this one
    // This again will be bade on how I set up the geometry of the cell, but basically will just be storing the data in an output file with points specified 
    // Need to think of the best format to do this in
}

void FDTDField::step()
{

}
Obj makeSphere(vector<double> mater, double rad)
{
    vector<double> geo(1,rad);
    return Obj(sphere, mater, geo);
}

