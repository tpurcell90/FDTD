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
        // The source array should go here, I need to read up on Vectors before I do that
        vector<Source<double>> srcArr();
        // The Object array will be inserted here map int to objects based on order in input file
        map<int,Obj> objArr;
        physGrid = make_shared<Grid2D<int>>(nx,ny,dx,dy);
        
        // Create the Grids. Do I need a null constructor for the set that I disregard?
        if(IP->pol.compare("Hz") == 0)
        {
            Ex = make_shared<Grid2D<double>>(nx-1,ny,dx,dy);
            Ey = make_shared<Grid2D<double>>(nx,ny-1,dx,dy);
            Hz = make_shared<Grid2D<double>>(nx-1,ny-1,dx,dy);
            phys_Ex = make_shared<Grid2D<int>>(nx-1,ny,dx,dy);
            phys_Ey = make_shared<Grid2D<int>>(nx,ny-1,dx,dy);
            phys_Hz = make_shared<Grid2D<int>>(nx-1,ny-1,dx,dy);

           // Hx = null;
           // Hy = null;
           // Ez = null;
        }
        else
        {
            Hx = make_shared<Grid2D<double>>(nx-1,ny,dx,dy);
            Hy = make_shared<Grid2D<double>>(nx,ny-1,dx,dy);
            Ez = make_shared<Grid2D<double>>(nx-1,ny-1,dx,dy);
            phys_Hx = make_shared<Grid2D<int>>(nx-1,ny,dx,dy);
            phys_Hy = make_shared<Grid2D<int>>(nx,ny-1,dx,dy);
            phys_Ez = make_shared<Grid2D<int>>(nx-1,ny-1,dx,dy);
            //Ex = null;
            //Ey = null;
            //Hz = null;
        }
    }
}

double dist(vector<double> pt1, vector<double> pt2)
{
    double sum = 0;
    for(int cc = 0; cc < pt1.size(); cc ++)
        sum += pow((pt1[cc]-pt2[cc]),2);
    return sqrt(sum);
}

void FDTDField::initializeGrid(programInputs *IP)
{
    for(int ii = 0; ii < objArr.size(); ii ++) 
    {
        if(IP ->pol.compare("Hz") == 0)
        {
            if(objArr.at(ii).s() == sphere)
            {
                double rad = objArr.at(ii).geo()[0];
                for(int jj = 0; jj < nx-1;jj ++)
                {
                    for(int kk = 0; kk < ny-1; kk ++)
                    {
                        vector<double> pt({jj*dx,kk*dy});
                        if(dist(pt,objArr.at(ii).loc()) <= rad)
                        {
                            phys_Ex->point(jj,kk) = ii;
                            phys_Ey->point(jj,kk) = ii;
                            phys_Hz->point(jj,kk) = ii;
                        }
                    }
                }
            }
            else if(objArr.at(ii).s() == block)
            {
                for(int jj = 0; jj < nx-1;jj ++)
                {
                    for(int kk = 0; kk < ny-1; kk ++)
                    {
                        if((jj*dx < (objArr.at(ii).loc()[0] + objArr.at(ii).geo()[0])) && (jj*dx >= (objArr.at(ii).loc()[0] - objArr.at(ii).geo()[0])) && (kk*dy <= (objArr.at(ii).loc()[1] + objArr.at(ii).geo()[1])) && (kk*dy >= (objArr.at(ii).loc()[1] - objArr.at(ii).geo()[1])) )
                        {
                            phys_Ex->point(jj,kk) = ii;
                            phys_Ey->point(jj,kk) = ii;
                            phys_Hz->point(jj,kk) = ii;
                        }
                    }
                }
            }
        }
        else
        {
            if(objArr.at(ii).s() == sphere)
            {
                double rad = objArr.at(ii).geo()[0];
                for(int jj = 0; jj < nx-1;jj ++)
                {
                    for(int kk = 0; kk < ny-1; kk ++)
                    {
                        vector<double> pt({jj*dx,kk*dy});
                        if(dist(pt,objArr.at(ii).loc()) <= rad)
                        {
                            phys_Hx->point(jj,kk) = ii;
                            phys_Hy->point(jj,kk) = ii;
                            phys_Ez->point(jj,kk) = ii;
                        }
                    }
                }
            }
            else if(objArr.at(ii).s() == block)
            {
                for(int jj = 0; jj < nx-1;jj ++)
                {
                    for(int kk = 0; kk < ny-1; kk ++)
                    {
                        if((jj*dx < (objArr.at(ii).loc()[0] + objArr.at(ii).geo()[0])) && (jj*dx >= (objArr.at(ii).loc()[0] - objArr.at(ii).geo()[0])) && (kk*dy <= (objArr.at(ii).loc()[1] + objArr.at(ii).geo()[1])) && (kk*dy >= (objArr.at(ii).loc()[1] - objArr.at(ii).geo()[1])) )
                        {
                            phys_Hx->point(jj,kk) = ii;
                            phys_Hy->point(jj,kk) = ii;
                            phys_Ez->point(jj,kk) = ii;
                        }
                    }
                }
            }
        }
    }
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
Obj makeSphere(vector<double> mater, double rad, vector<double> loc)
{
    vector<double> geo(1,rad);
    return Obj(sphere, mater, geo,loc);
}

Obj makeBlock(vector<double> mater, vector<double> dims, vector<double> loc)
{
    return Obj(block, mater, dims,loc);
}
