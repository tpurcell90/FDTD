#include "FDTDField.hpp"

// #include <assert.h>
// #include <iomanip>
#include <iostream>
#include <memory>
#include <fstream>
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
        t_cur = 0;
        res = IP->res;
        dx = 1.0/res;
        dy = 1.0/res;
        dt = 1.0/(IP->courant*res);
        nx = int(IP->x_size*dx); //Better way here
        ny = int(IP->y_size*dy); //Better way here
        physGrid = make_shared<Grid2D<int>>(nx,ny,dx,dy);
        vector<Source<double>> srcArr(IP -> srcArr);
        vector<Obj> objArr(IP -> objArr); 
        vector<Detector<double>> dtcArr (IP -> dctArr);
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
void FDTDField::inc_t()
{
    t_cur += dt;
}


void FDTDField::initializeGrid(programInputs *IP)
{
    for(int ii = 0; ii < srcArr.size(); ii++)
    {
        //Take the size and location of the source array convert that into the pulse in the step fxn
    }
    for(int ii = 0; ii < objArr.size(); ii ++) 
    {
        if(IP ->pol.compare("Hz") == 0)
        {
            if(objArr.at(ii).s() == sphere)
            {
                for(int jj = 0; jj < nx-1;jj ++)
                {
                    for(int kk = 0; kk < ny-1; kk ++)
                    {
                        vector<double> pt({(jj+0.5)*dx,kk*dy});
                        if(objArr.at(ii).isObj(pt)) 
                            phys_Ex->point(jj,kk) = ii+1;
                        pt[1] += 0.5*dy;
                        if(objArr.at(ii).isObj(pt))
                            phys_Hz->point(jj,kk) = ii+1;
			            pt[0] -= 0.5*dx;
                        if(objArr.at(ii).isObj(pt))
                            phys_Ey->point(jj,kk) = ii+1;
                    }
                }
            }
            else if(objArr.at(ii).s() == block)
            {
                for(int jj = 0; jj < nx-1;jj ++)
                {
                    for(int kk = 0; kk < ny-1; kk ++)
                    {
                        vector<double> pt({(jj+0.5)*dx,kk*dy});
                        if(objArr.at(ii).isObj(pt))
                            phys_Ex->point(jj,kk) = ii+1;
                        pt[1] += 0.5*dy;
                        if(objArr.at(ii).isObj(pt))
                            phys_Hz->point(jj,kk) = ii+1;
                        pt[0] -= 0.5*dx;
                        if(objArr.at(ii).isObj(pt))
                            phys_Ey->point(jj,kk) = ii+1;
                    }
                }
            }
        }
        else
        {
            if(objArr.at(ii).s() == sphere)
            {
                for(int jj = 0; jj < nx-1;jj ++)
                {
                    for(int kk = 0; kk < ny-1; kk ++)
                    {
                        vector<double> pt({(jj+0.5)*dx,kk*dy});
                        if(objArr.at(ii).isObj(pt))
                            phys_Ex->point(jj,kk) = ii+1;
                        pt[1] += 0.5*dy;
                        if(objArr.at(ii).isObj(pt))
                            phys_Hz->point(jj,kk) = ii+1;
                        pt[0] -= 0.5*dx;
                        if(objArr.at(ii).isObj(pt))
                            phys_Ey->point(jj,kk) = ii+1;
                    }
                }
            }
            else if(objArr.at(ii).s() == block)
            {
                for(int jj = 0; jj < nx-1;jj ++)
                {
                    for(int kk = 0; kk < ny-1; kk ++)
                    {
                        vector<double> pt({(jj+0.5)*dx,kk*dy});
                        if(objArr.at(ii).isObj(pt))
                            phys_Ex->point(jj,kk) = ii+1;
                        pt[1] += 0.5*dy;
                        if(objArr.at(ii).isObj(pt))
                            phys_Hz->point(jj,kk) = ii+1;
                        pt[0] -= 0.5*dx;
                        if(objArr.at(ii).isObj(pt))
                            phys_Ey->point(jj,kk) = ii+1;
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
    ofstream outFile;
    outFile.open("Out.dat");
    for(int ii = 0; ii < dtcArr.size(); ii ++)
    {
        outFile << t_cur << "\t" << dtcArr[ii].loc()[0] << "\t" << dtcArr[ii].loc()[1] << "\t" << dtcArr[ii].output(Ez); 
    }
    
    
}

void FDTDField::step()
{
    // disregard PML's to start
    inc_t();
    // Time step and other fdtd constants. When Other things introduced these will change.
    double dt_mu0 = dt/4.0e-7;
    double den_ex = 1.0/dx;
    double den_hx = 1.0/dx;
    double den_ey = 1.0/dy;
    double den_hy = 1.0/dy;
    // PML and source stuff will change this
    for(int ii = 1; ii < nx - 1; ii ++)
    {
        for(int jj = 1; jj < ny -1; jj ++)
        {
            // Update source
            for(int kk = 0; kk < srcArr.size(); kk ++)
            {
                if(srcArr[kk].pol() == EZ)
                    Ez -> point(ii,jj) += srcArr[kk].prof().pulse(t_cur);
                else if(srcArr[kk].pol() == EX)
                    Ex -> point(ii,jj) += srcArr[kk].prof().pulse(t_cur);
                else if(srcArr[kk].pol() == EY)
                    Ey -> point(ii,jj) += srcArr[kk].prof().pulse(t_cur);
                else if(srcArr[kk].pol() == HX)
                    Hx -> point(ii,jj) += srcArr[kk].prof().pulse(t_cur);
                else if(srcArr[kk].pol() == HY)
                    Hy -> point(ii,jj) += srcArr[kk].prof().pulse(t_cur);
                else if(srcArr[kk].pol() == HZ)
                    Hz -> point(ii,jj) += srcArr[kk].prof().pulse(t_cur);
            }
            Hz->point(ii,jj) += dt_mu0*((Ey->point(ii,jj)-Ey->point(ii+1,jj))*den_hx + (Ex->point(ii,jj+1)-Ex->point(ii,jj))*den_hy);
            // Look up how D-L models affects the Electric fields, Something about Jx/Jy in Maxim's code also
            // look up how to deal with the frequency dependence of eps in the time domain
            double freq = 1.0/t_cur; // I know this is not right
            double dt_eps = dt/(1.0/(4.0e-7*pow(299792458.0,2)) * real(objArr[phys_Ey->point(ii,jj)].dielectric(freq)));
            Ey->point(ii,jj) += dt_eps*(Hz->point(ii-1,jj)-Hz->point(1,jj))*den_ex;
            dt_eps = dt/(1.0/(4.0e-7*pow(299792458.0,2)) * real(objArr[phys_Ex->point(ii,jj)].dielectric(freq)));
            Ex->point(ii,jj) += dt_eps*(Hz->point(ii,jj)-Hz->point(ii,jj-1))*den_ey;
        }
    }
    ouputField();
}

/*std::vector<double> FDTDField::pml(int npml, int m, int ma)
{
    
    double mu0 = 4.0e-7;
    double eps0 = 1.0/(4.0e-7*pow(299792458.0,2));
    double sigmaCPML = 0.75*(0.8*(m+1)/(dx*sqrt(mu0/eps0*eps_delectric)));
    double alphaCPML = 0.24;
    double kappaCPML = 20.0;
    double precision psi_Hzy_1[Nx-1,npml-1],psi_Exy_1[Nx-1,npml]                              
    double precision psi_Hzy_2(Nx-1,npml-1),psi_Exy_2(Nx-1,npml)
    double precision be_y(npml),ce_y(npml),alphae_y(npml),sige_y(npml),kappae_y(npml)
    double precision bh_y(npml-1),ch_y(npml-1),alphah_y(npml-1),sigh_y(npml-1),kappah_y(npml-1)
    double precision den_ex(Nx),den_hx(Nx),den_ey(N_loc),den_hy(N_loc)
    double precision psi_Hzy_1_inc(Nx-1,npml-1),psi_Exy_1_inc(Nx-1,npml)                              
    double precision psi_Hzy_2_inc(Nx-1,npml-1),psi_Exy_2_inc(Nx-1,npml)


    for(int jj = 0; jj < npml; jj++)
    {
        sige_y[jj] = sigmaCPML* pow((npml-jj) / (npml-1.0),m);
        alphae_y[jj] = alphaCPML * pow((jj-1) / (npml-1.0),m);    
        kappae_y[jj] = 1.0 + (kappaCPML - 1.0) * pow((npml- jj) / (npml - 1.0),m);
        be_y[jj] = exp(-(sige_y[jj] / kappae_y[jj] + alphae_y[jj]) * dt/eps0);
        if((sige_y[jj] == 0.0) && (alphae_y[jj] == 0.0) && (j == npml))
            ce_y[jj] = 0.0;
        else
            ce_y[jj] = sige_y[jj] * (be_y[jj]-1.0) / (sige_y[jj] + kappae_y[jj] * alphae_y[jj]) / kappae_y[jj];
    }
   
    vector<double> test(3,0.0); 
    return test;
    */
//}

Obj makeSphere(vector<double> mater, double rad, vector<double> loc)
{
    vector<double> geo(1,rad);
    return Obj(sphere, mater, geo,loc);
}

Obj makeBlock(vector<double> mater, vector<double> dims, vector<double> loc)
{
    return Obj(block, mater, dims,loc);
}
