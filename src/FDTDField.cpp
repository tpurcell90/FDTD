#include "FDTDField.hpp"

// #include <assert.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <fstream>
// #include <random>
// #include <stdexcept>
// #include <string>
// #include <vector>
// #include <complex>



using namespace std;

FDTDField::FDTDField(programInputs &IP)
{
    //Define cell parameters
    t_cur = 0;
    res = IP.res;
    dx = 1.0/res;
    dy = 1.0/res;
    dt = 0.99/(sqrt(pow(dx,-2.0) + pow(dy,-2.0)));
    nx = int((IP.x_size)*res) + 1; //Better way here; + 1 to include the 0 point
    ny = int((IP.y_size)*res) + 1; //Better way here; + 1 to include the 0 point
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
        Hx = nullptr;
        Hy = nullptr;
        Ez = nullptr;
        phys_Hx = nullptr;
        phys_Hy = nullptr;
        phys_Ez = nullptr;
    }
    else
    {
        Hx = make_shared<Grid2D<double>>(nx-1,ny,dx,dy);
        Hy = make_shared<Grid2D<double>>(nx,ny-1,dx,dy);
        Ez = make_shared<Grid2D<double>>(nx,ny,dx,dy);
        phys_Hx = make_shared<Grid2D<int>>(nx-1,ny,dx,dy);
        phys_Hy = make_shared<Grid2D<int>>(nx,ny-1,dx,dy);
        phys_Ez = make_shared<Grid2D<int>>(nx-1,ny-1,dx,dy);
        // These are never used in the TM mode
        Ex = nullptr;
        Ey = nullptr;
        Hz = nullptr;
        phys_Ex = nullptr;
        phys_Ey = nullptr;
        phys_Hz = nullptr;
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

void FDTDField::inc_t() {t_cur += dt;}

void FDTDField::initializeGrid(programInputs &IP)
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

void FDTDField::ouputField() //iostream as input parameter?
{
    // Work on a plan for this one
    // This again will be bade on how I set up the geometry of the cell, but basically will just be storing the data in an output file with points specified
    // Need to think of the best format to do this in
    ofstream outFile;
    outFile.open("Out.dat", ios_base::app);
    for(int ii = 0; ii < dtcArr.size(); ii ++)
    {
        outFile << t_cur << "\t" << dtcArr[ii].loc()[0] << "\t" << dtcArr[ii].loc()[1] << "\t" << setw(10) << dtcArr[ii].output(Ez) << "\t"<< srcArr[0].prof().pulse(t_cur) << "\n";
        cout << setw(9) << t_cur << "\t" << dtcArr[ii].loc()[0] << "\t" << dtcArr[ii].loc()[1] << "\t" << setw(10) << dtcArr[ii].output(Ez)<< "\t" <<  srcArr[0].prof().pulse(t_cur) << "\n";
    }
    outFile.close();


}

void FDTDField::step()
{
    // disregard PML's to start
    // Source
    for(int kk = 0; kk < srcArr.size(); kk ++)
    {
        int ii = srcArr[kk].loc()[0];
        int jj = srcArr[kk].loc()[1];
        switch ( srcArr[kk].pol() )
        {
            case EZ: //if(srcArr[kk].pol() == EZ)
                Ez -> point(ii,jj) = srcArr[kk].prof().pulse(t_cur);
                Hx -> point(ii,jj) = 0;
                Hy -> point(ii,jj) = 0;
                break;
            case HX: //else if(srcArr[kk].pol() == HX)
                Hx -> point(ii,jj) = srcArr[kk].prof().pulse(t_cur);
                Ez -> point(ii,jj) = 0;
                Hy -> point(ii,jj) = 0;
                break;
            case HY: //else if(srcArr[kk].pol() == HY)
                Hy -> point(ii,jj) = srcArr[kk].prof().pulse(t_cur);
                Ez -> point(ii,jj) = 0;
                Hx -> point(ii,jj) = 0;
                break;
            case HZ: //else if(srcArr[kk].pol() == HZ)
                Hz -> point(ii,jj) = srcArr[kk].prof().pulse(t_cur);
                Ex -> point(ii,jj) = 0;
                Ey -> point(ii,jj) = 0;
                break;
            case EX: //else if(srcArr[kk].pol() == EX)
                Ex -> point(ii,jj) = srcArr[kk].prof().pulse(t_cur);
                Hz -> point(ii,jj) = 0;
                Ey -> point(ii,jj) = 0;
                break;
            case EY: //else if(srcArr[kk].pol() == EY)
                Ey -> point(ii,jj) = srcArr[kk].prof().pulse(t_cur);
                Hz -> point(ii,jj) = 0;
                Ex -> point(ii,jj) = 0;
                break;
            default:
                throw logic_error("reached a default case in a switch state that should never happen!");
                break;
        }
    }
    for(int ii = 1; ii < nx - 1; ii ++)
    {
        for(int jj = 1; jj < ny - 1; jj ++)
        {
            if(Ez)
            {
                //Always true in a magnetic lossless enviroment, with mu0 = 1, sigma_m = 0;
                double c_hxh = 1.0;
                double c_hxe = 1.0 * dt/dx;
                double c_hyh = 1.0;
                double c_hye = 1.0 * dt/dy;
                Hx->point(ii,jj) = c_hxh * Hx->point(ii,jj) - c_hxe * (Ez->point(ii,jj+1)-Ez->point(ii,jj));
                Hy->point(ii,jj) = c_hyh * Hy->point(ii,jj) + c_hye * (Ez->point(ii+1,jj)-Ez->point(ii,jj));
            }
            else
            {
                //Always true in a magnetic lossless enviroment, with mu0 = 1, sigma_m = 0;
                double c_hzh = 1.0;
                double c_hze = 1.0 * dt/dx;
                Hz->point(ii,jj) = c_hzh * Hz->point(ii,jj) + c_hze * ((Ex->point(ii,jj+1) - Ex->point(ii,jj)) - (Ey->point(ii+1,jj)-Ey->point(ii,jj)));

            }
        }
    }
    // Code for prefect reflectors, which we don't ever really want
    /*for(int ii = 0; ii < nx; ii ++)
    {
        double c_hxh = 1.0;
        double c_hxe = 1.0 * dt/dx;
        double c_hyh = 1.0;
        double c_hye = 1.0 * dt/dy;
        Hx->point(ii,0) = c_hxh * Hx->point(ii,0) - c_hxe * (Ez->point(ii,0+1)-Ez->point(ii,0));
        Hy->point(0,ii) = c_hyh * Hy->point(0,ii) + c_hye * (Ez->point(0+1,ii)-Ez->point(0,ii));
    }*/
    for(int ii = 1; ii < nx - 1; ii ++)
    {
        for(int jj = 1; jj < ny - 1; jj ++)
        {
            if(Ez)
            {
                //Only True in Vac Once Mat/PML introduced vectorize it
                double c_eze = 1.0;
                double c_ezh = 1.0 * dt/dx;
                Ez->point(ii,jj) = c_eze * Ez->point(ii,jj) + c_ezh * ((Hy->point(ii,jj)-Hy->point(ii-1,jj)) - (Hx->point(ii,jj)-Hx->point(ii,jj-1)));
            }
            else
            {
                //Only True in Vac Once Mat/PML introduced vectorize it
                //
                double c_exe = 1.0;
                double c_exh = 1.0 * dt/dx;
                double c_eye = 1.0;
                double c_eyh = 1.0 * dt/dy;
                Ey->point(ii,jj) = c_eye * Ey->point(ii,jj) - c_eyh * (Hz->point(ii,jj) - Hz->point(ii-1,jj));
                Ex->point(ii,jj) = c_exe * Ex->point(ii,jj) + c_exh * (Hz->point(ii,jj) - Hz->point(ii,jj-1));
            }
        }
    }
    ouputField();
    inc_t();
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
