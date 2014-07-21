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

FDTDField::FDTDField(programInputs IP)
{
    //Define cell parameters
    t_cur = 0;
    res = IP.res;
    dx = 1.0/res;
    dy = 1.0/res;
    dt = IP.courant/(res);
    nx = int(IP.x_size*res); //Better way here
    ny = int(IP.y_size*res); //Better way here
    srcArr = IP.srcArr;
    objArr = IP.objArr; 
    dtcArr = IP.dctArr;
    // Create the Grids. Do I need a null constructor for the set that I disregard?
    if(IP.pol.compare("Hz") == 0 || IP.pol.compare("Ey") == 0 || IP.pol.compare("Ex") == 0)
    {
        Ex = make_shared<Grid2D<double>>(nx-1,ny,dx,dy);
        Ey = make_shared<Grid2D<double>>(nx,ny-1,dx,dy);
        Hz = make_shared<Grid2D<double>>(nx-1,ny-1,dx,dy);
        phys_Ex = make_shared<Grid2D<int>>(nx-1,ny,dx,dy);
        phys_Ey = make_shared<Grid2D<int>>(nx,ny-1,dx,dy);
        phys_Hz = make_shared<Grid2D<int>>(nx-1,ny-1,dx,dy);
        //These are never used in the TE mode
        Hx = NULL;
        Hy = NULL;
        Ez = NULL;
        phys_Hx = NULL;
        phys_Hy = NULL;
        phys_Ez = NULL;
    }
    else
    {
        Hx = make_shared<Grid2D<double>>(nx-1,ny,dx,dy);
        Hy = make_shared<Grid2D<double>>(nx,ny-1,dx,dy);
        Ez = make_shared<Grid2D<double>>(nx-1,ny-1,dx,dy);
        phys_Hx = make_shared<Grid2D<int>>(nx-1,ny,dx,dy);
        phys_Hy = make_shared<Grid2D<int>>(nx,ny-1,dx,dy);
        phys_Ez = make_shared<Grid2D<int>>(nx-1,ny-1,dx,dy);
        // These are never used in the TM mode
        Ex = NULL;
        Ey = NULL;
        Hz = NULL;
        phys_Ex = NULL;
        phys_Ey = NULL;
        phys_Hz = NULL;
    }

}
// Access Functions
size_t FDTDField::n_x() {return nx;}
size_t FDTDField::n_y() {return ny;}
double FDTDField::d_x() {return dx;}
double FDTDField::d_y() {return dy;}
double FDTDField::d_t() {return dt;}
double FDTDField::getTime() {return t_cur;}
int FDTDField::getRes() {return res;}
std::vector<Detector<double>> FDTDField::getDtcArr() {return dtcArr;}
shared_ptr<Grid2D<int>> FDTDField::getPhysEz() {return phys_Ez;}

void FDTDField::inc_t()
{
    t_cur += dt;
}


void FDTDField::initializeGrid(programInputs IP)
{
    for(int ii = 0; ii < objArr.size(); ii ++) 
    {
        if(Hz)
        {
            if(objArr[ii].s() == sphere)
            {
                for(int jj = 0; jj < nx-1;jj ++)
                {
                    for(int kk = 0; kk < ny-1; kk ++)
                    {
                        vector<double> pt({(jj+0.5)*dx,kk*dy});
                        if(objArr[ii].isObj(pt)) 
                            phys_Ex->point(jj,kk) = ii;
                        pt[1] += 0.5*dy;
                        if(objArr[ii].isObj(pt))
                            phys_Hz->point(jj,kk) = ii;
			            pt[0] -= 0.5*dx;
                        if(objArr[ii].isObj(pt))
                            phys_Ey->point(jj,kk) = ii;
                    }
                }
            }
            else if(objArr[ii].s() == block)
            {
                for(int jj = 0; jj < nx-1;jj ++)
                {
                    for(int kk = 0; kk < ny-1; kk ++)
                    {
                        vector<double> pt({(jj+0.5)*dx,kk*dy});
                        if(objArr[ii].isObj(pt))
                            phys_Ex->point(jj,kk) = ii;
                        pt[1] += 0.5*dy;
                        if(objArr[ii].isObj(pt))
                            phys_Hz->point(jj,kk) = ii;
                        pt[0] -= 0.5*dx;
                        if(objArr[ii].isObj(pt))
                            phys_Ey->point(jj,kk) = ii;
                    }
                }
            }
        }
        else
        {
            if(objArr[ii].s() == sphere)
            {
                for(int jj = 0; jj < nx-1;jj ++)
                {
                    for(int kk = 0; kk < ny-1; kk ++)
                    {
                        vector<double> pt({(jj+0.5)*dx,kk*dy});
                        if(objArr[ii].isObj(pt))
                            phys_Hx->point(jj,kk) = ii;
                        pt[1] += 0.5*dy;
                        if(objArr[ii].isObj(pt))
                            phys_Ez->point(jj,kk) = ii;
                        pt[0] -= 0.5*dx;
                        if(objArr[ii].isObj(pt))
                            phys_Hy->point(jj,kk) = ii;
                    }
                }
            }
            else if(objArr[ii].s() == block)
            {
                for(int jj = 0; jj < nx-1;jj ++)
                {
                    for(int kk = 0; kk < ny-1; kk ++)
                    {
                        vector<double> pt({(jj+0.5)*dx,kk*dy});
                        if(objArr[ii].isObj(pt))
                        {
                            phys_Hx->point(jj,kk) = ii;
                        }
                        pt[1] += 0.5*dy;
                        if(objArr[ii].isObj(pt))
                            phys_Ez->point(jj,kk) = ii;
                        pt[0] -= 0.5*dx;
                        if(objArr[ii].isObj(pt))
                            phys_Hy->point(jj,kk) = ii;
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
    outFile.open("Out.dat", ios_base::app);
    for(int ii = 0; ii < dtcArr.size(); ii ++)
    {
        outFile << t_cur << "\t" << dtcArr[ii].loc()[0] << "\t" << dtcArr[ii].loc()[1] << "\t" << dtcArr[ii].output(Ez) << "\t"<< srcArr[0].prof().pulse(t_cur) << "\n";
        cout << t_cur << "\t" << dtcArr[ii].loc()[0] << "\t" << dtcArr[ii].loc()[1] << "\t" << dtcArr[ii].output(Ez) << "\t" <<  srcArr[0].prof().pulse(t_cur) << "\n";
    }
    outFile.close();
    
    
}

void FDTDField::step()
{
    // disregard PML's to start
    inc_t();
    // Source
    for(int kk = 0; kk < srcArr.size(); kk ++)
    {
        int ii = srcArr[kk].loc()[0];
        int jj = srcArr[kk].loc()[1];
        if(srcArr[kk].pol() == EZ)
        {
            Ez -> point(ii,jj) += srcArr[kk].prof().pulse(t_cur);
        }
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
    // Time step and other fdtd constants. When Other things introduced these will change.
    double dt_mu0 = dt;
    double den_ex = dx;
    double den_hx = dx;
    double den_ey = dy;
    double den_hy = dy;
    double dt_eps = dt;    
    for(int ii = 2; ii < nx - 2; ii ++)
    {
        for(int jj = 2; jj < ny -2; jj ++)
        {
            if(Ez)
            {
                
                Ez->point(ii,jj) += dt_mu0*((Hy->point(ii,jj)-Hy->point(ii+1,jj))*den_hx + (Hx->point(ii,jj+1)-Hx->point(ii,jj))*den_hy);
                // Look up how D-L models affects the Electric fields, Something about Jx/Jy in Maxim's code also
                // look up how to deal with the frequency dependence of eps in the time domain
                //double dt_eps = 1.0;
                Hy->point(ii,jj) += dt_eps*(Ez->point(ii-1,jj)-Ez->point(1,jj))*den_ex;
                Hx->point(ii,jj) += dt_eps*(Ez->point(ii,jj)-Ez->point(ii,jj-1))*den_ey;
            }
            else
            {
                //Ez->point(ii,jj) += dt_mu0*((Hy->point(ii,jj)-Hy->point(ii+1,jj))*den_hx + (Hx->point(ii,jj+1)-Hx->point(ii,jj))*den_hy);
                // Look up how D-L models affects the Electric fields, Something about Jx/Jy in Maxim's code also
                // look up how to deal with the frequency dependence of eps in the time domain
                //double freq = 1.0/t_cur; // I know this is not right
                //double dt_eps = dt/(1.0/(4.0e-7*pow(299792458.0,2)) * real(objArr[phys_Hy->point(ii,jj)].dielectric(freq)));
                //Hy->point(ii,jj) += dt_eps*(Ez->point(ii-1,jj)-Ez->point(1,jj))*den_ex;
                //dt_eps = dt/(1.0/(4.0e-7*pow(299792458.0,2)) * real(objArr[phys_Hx->point(ii,jj)].dielectric(freq)));
                //Hx->point(ii,jj) += dt_eps*(Ez->point(ii,jj)-Ez->point(ii,jj-1))*den_ey;
            }
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
