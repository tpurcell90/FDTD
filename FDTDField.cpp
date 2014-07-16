#include "FDTDField.hpp"

using namespace std;

FDTDField::(programInputs &IP)
{
    //Define cell parameters
    res = IP.res;
    dx = 1.0/res;
    dy = 1.0/res;
    dt = 1.0/(IP.courant*res);
    nx = IP.size_x/dx;
    ny = IP.size_y/dy;
    nt = IP.t_lim/dt;
    // The source array should go here, I am not sure how I am going to do the source input yet... From this I can get the polarizations
    srcArray[1] = new Source(IP.prof,IP.pol);
    physGrid = new Grid2D(nx,ny,dx,dy)
    // Create the Grids. Do I need a null constructor for the set that I disregard?
    // If people do calculations with mulptiple different polarizations this will need to get more complicated
    // I also need to change the grid stuff (They should not be the same but will look up later)
    if(srcArray[0].polarization.compare("Hz") == 0)
    {
        Ex = new Grid2D(nx,ny,dx,dy)
        Ey = new Grid2D(nx,ny,dx,dy)
        Hz = new Grid2D(nx,ny,dx,dy)
        Hx = null;
        Hy = null;
        Ez = null;
    }
    else
    {
        Hx = new Grid2D(nx,ny,dx,dy);
        Hy = new Grid2D(nx,ny,dx,dy);
        Ez = new Grid2D(nx,ny,dx,dy);
        Ex = null;
        Ey = null;
        Hz = null;
    }
}

void FDTDField::initializeGrid(programInputs &IP)
{
    // This will plot the geometry structures onto the physical grid still thinking about this one I was thinking about how best to do this
}

void FDTDField::ouputField()
{
    // Work on a plan for this one
    // This again will be bade on how I set up the geometry of the cell, but basically will just be storing the data in an output file with points specified 
    // Need to think of the best format to do this in
}



