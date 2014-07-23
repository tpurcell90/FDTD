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
    tcur_ = 0;
    res_ = IP.res_;
    dx_ = 1.0/res_;
    dy_ = 1.0/res_;
    dt_ = 0.99/(sqrt(pow(dx_,-2.0) + pow(dy_,-2.0)));
    nx_ = int((IP.x_size_)*res_) + 1; //Better way here; + 1 to include the 0 point
    ny_ = int((IP.y_size_)*res_) + 1; //Better way here; + 1 to include the 0 point
    srcArr_ = IP.srcArr_;
    objArr_ = IP.objArr_;
    dtcArr_ = IP.dctArr_;
    for(int ii =0; ii < dtcArr_.size(); ii++)
    {
        ofstream outFile;
        outFile.open(dtcArr_[ii].outfile());
        outFile << "Output for the " << to_string(ii) << "th detector" << endl;
        outFile.close();
    } 
    // Create the Grids. Do I need a null constructor for the set that I disregard?
    if(IP.pol_.compare("Hz") == 0 || IP.pol_.compare("Ey") == 0 || IP.pol_.compare("Ex") == 0)
    {
        Ex_ = make_shared<Grid2D<double>>(nx_-1,ny_,dx_,dy_);
        Ey_ = make_shared<Grid2D<double>>(nx_,ny_-1,dx_,dy_);
        Hz_ = make_shared<Grid2D<double>>(nx_-1,ny_-1,dx_,dy_);
        phys_Ex_ = make_shared<Grid2D<int>>(nx_-1,ny_,dx_,dy_);
        phys_Ey_ = make_shared<Grid2D<int>>(nx_,ny_-1,dx_,dy_);
        phys_Hz_ = make_shared<Grid2D<int>>(nx_-1,ny_-1,dx_,dy_);
        //These are never used in the TE mode
        Hx_ = nullptr;
        Hy_ = nullptr;
        Ez_ = nullptr;
        phys_Hx_ = nullptr;
        phys_Hy_ = nullptr;
        phys_Ez_ = nullptr;
    }
    else
    {
        Hx_ = make_shared<Grid2D<double>>(nx_-1,ny_,dx_,dy_);
        Hy_ = make_shared<Grid2D<double>>(nx_,ny_-1,dx_,dy_);
        Ez_ = make_shared<Grid2D<double>>(nx_,ny_,dx_,dy_);
        phys_Hx_ = make_shared<Grid2D<int>>(nx_-1,ny_,dx_,dy_);
        phys_Hy_ = make_shared<Grid2D<int>>(nx_,ny_-1,dx_,dy_);
        phys_Ez_ = make_shared<Grid2D<int>>(nx_-1,ny_-1,dx_,dy_);
        // These are never used in the TM mode
        Ex_ = nullptr;
        Ey_ = nullptr;
        Hz_ = nullptr;
        phys_Ex_ = nullptr;
        phys_Ey_ = nullptr;
        phys_Hz_ = nullptr;
    }

}

void FDTDField::initializeGrid(programInputs &IP)
{
    for(int ii = 0; ii < objArr_.size(); ii ++)
    {
        if(Hz_)
        {
            if(objArr_[ii].s() == sphere)
            {
                for(int jj = 0; jj < nx_-1;jj ++)
                {
                    for(int kk = 0; kk < ny_-1; kk ++)
                    {
                        vector<double> pt({(jj+0.5)*dx_,kk*dy_});
                        if(objArr_[ii].isObj(pt))
                            phys_Ex_->point(jj,kk) = ii;
                        pt[1] += 0.5*dy_;
                        if(objArr_[ii].isObj(pt))
                            phys_Hz_->point(jj,kk) = ii;
			            pt[0] -= 0.5*dx_;
                        if(objArr_[ii].isObj(pt))
                            phys_Ey_->point(jj,kk) = ii;
                    }
                }
            }
            else if(objArr_[ii].s() == block)
            {
                for(int jj = 0; jj < nx_-1;jj ++)
                {
                    for(int kk = 0; kk < ny_-1; kk ++)
                    {
                        vector<double> pt({(jj+0.5)*dx_,kk*dy_});
                        if(objArr_[ii].isObj(pt))
                            phys_Ex_->point(jj,kk) = ii;
                        pt[1] += 0.5*dy_;
                        if(objArr_[ii].isObj(pt))
                            phys_Hz_->point(jj,kk) = ii;
                        pt[0] -= 0.5*dx_;
                        if(objArr_[ii].isObj(pt))
                            phys_Ey_->point(jj,kk) = ii;
                    }
                }
            }
        }
        else
        {
            if(objArr_[ii].s() == sphere)
            {
                for(int jj = 0; jj < nx_-1;jj ++)
                {
                    for(int kk = 0; kk < ny_-1; kk ++)
                    {
                        vector<double> pt({(jj+0.5)*dx_,kk*dy_});
                        if(objArr_[ii].isObj(pt))
                            phys_Hx_->point(jj,kk) = ii;
                        pt[1] += 0.5*dy_;
                        if(objArr_[ii].isObj(pt))
                            phys_Ez_->point(jj,kk) = ii;
                        pt[0] -= 0.5*dx_;
                        if(objArr_[ii].isObj(pt))
                            phys_Hy_->point(jj,kk) = ii;
                    }
                }
            }
            else if(objArr_[ii].s() == block)
            {
                for(int jj = 0; jj < nx_-1;jj ++)
                {
                    for(int kk = 0; kk < ny_-1; kk ++)
                    {
                        vector<double> pt({(jj+0.5)*dx_,kk*dy_});
                        if(objArr_[ii].isObj(pt))
                            phys_Hx_->point(jj,kk) = ii;
                        pt[1] += 0.5*dy_;
                        if(objArr_[ii].isObj(pt))
                            phys_Ez_->point(jj,kk) = ii;
                        pt[0] -= 0.5*dx_;
                        if(objArr_[ii].isObj(pt))
                            phys_Hy_->point(jj,kk) = ii;
                    }
                }
            }
        }
    }
}

void FDTDField::ouputField(Detector<double> d) //iostream as input parameter?
{
    ofstream outFile;
    outFile.open(d.outfile(), ios_base::app);
    double eps = 1.00; // Change this but for a test it will work
    switch ( d.pol() )
    {
        case EZ: 
            outFile << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Ez_,eps)<< "\t" << setw(10) << srcArr_[0].prof().pulse(tcur_) << endl;
            cout << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Ez_,eps)<< "\t" << setw(10) << srcArr_[0].prof().pulse(tcur_) << endl;
            break;
        case HX:
            outFile << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Hx_,eps)<< endl;
            cout << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Hx_,eps)<< endl;
            break;
        case HY:
            outFile << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Hy_,eps)<< endl;
            cout << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Hy_,eps)<< endl;
            break;
        case HZ:
            outFile << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Hz_,eps)<< endl;
            cout << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Hz_,eps)<< endl;
            break;
        case EX:
            outFile << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Ex_,eps)<< endl;
            cout << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Ex_,eps)<< endl;
            break;
        case EY:
            outFile << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Ey_,eps)<< endl;
            cout << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Ey_,eps)<< endl;
            break;
        default:
            throw logic_error("reached a default case in a switch state that should never happen!");
            break;
    }
    outFile.close();
}

void FDTDField::step()
{
    // disregard PML's to start
    // Source
    tcur_ += dt_;
    for(int kk = 0; kk < srcArr_.size(); kk ++)
    {
        int ii = srcArr_[kk].loc()[0];
        int jj = srcArr_[kk].loc()[1];
        switch ( srcArr_[kk].pol() )
        {
            case EZ: //if(srcArr[kk].pol() == EZ)
                Ez_ -> point(ii,jj) = srcArr_[kk].prof().pulse(tcur_);
                break;
            case HX: //else if(srcArr[kk].pol() == HX)
                Hx_ -> point(ii,jj) = srcArr_[kk].prof().pulse(tcur_);
                break;
            case HY: //else if(srcArr[kk].pol() == HY)
                Hy_ -> point(ii,jj) = srcArr_[kk].prof().pulse(tcur_);
                break;
            case HZ: //else if(srcArr[kk].pol() == HZ)
                Hz_ -> point(ii,jj) = srcArr_[kk].prof().pulse(tcur_);
                break;
            case EX: //else if(srcArr[kk].pol() == EX)
                Ex_ -> point(ii,jj) = srcArr_[kk].prof().pulse(tcur_);
                break;
            case EY: //else if(srcArr[kk].pol() == EY)
                Ey_ -> point(ii,jj) = srcArr_[kk].prof().pulse(tcur_);
                break;
            default:
                throw logic_error("reached a default case in a switch state that should never happen!");
                break;
        }
    }
    for(int ii = 1; ii < nx_ - 1; ii ++)
    {
        for(int jj = 1; jj < ny_ - 1; jj ++)
        {
            if(Ez_)
            {
                //Always true in a magnetic lossless enviroment, with mu0 = 1, sigma_m = 0;
                double c_hxh = 1.0;
                double c_hxe = 1.0 * dt_/dx_;
                double c_hyh = 1.0;
                double c_hye = 1.0 * dt_/dy_;
                Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) - c_hxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
            }
            else
            {
                //Always true in a magnetic lossless enviroment, with mu0 = 1, sigma_m = 0;
                double c_hzh = 1.0;
                double c_hze = 1.0 * dt_/dx_;
                Hz_->point(ii,jj) = c_hzh * Hz_->point(ii,jj) + c_hze * ((Ex_->point(ii,jj+1) - Ex_->point(ii,jj)) - (Ey_->point(ii+1,jj)-Ey_->point(ii,jj)));

            }
        }
    }
    // Code for prefect reflectors, which we don't ever really want
    for(int ii = 0; ii < nx_; ii ++)
    {
        double c_hxh = 1.0;
        double c_hxe = 1.0 * dt_/dx_;
        double c_hyh = 1.0;
        double c_hye = 1.0 * dt_/dy_;
        Hx_->point(ii,0) = c_hxh * Hx_->point(ii,0) - c_hxe * (Ez_->point(ii,0+1)-Ez_->point(ii,0));
        Hy_->point(0,ii) = c_hyh * Hy_->point(0,ii) + c_hye * (Ez_->point(0+1,ii)-Ez_->point(0,ii));
    }
    for(int ii = 1; ii < nx_-1; ii ++)
    {
        double c_hxh = 1.0;
        double c_hxe = 1.0 * dt_/dx_;
        double c_hyh = 1.0;
        double c_hye = 1.0 * dt_/dy_;
        Hx_->point(nx_-1,ii) = c_hxh * Hx_->point(ii,0) - c_hxe * (Ez_->point(ii,0+1)-Ez_->point(ii,0));
        Hy_->point(ii,nx_-1) = c_hyh * Hy_->point(0,ii) + c_hye * (Ez_->point(0+1,ii)-Ez_->point(0,ii));
    }
    for(int ii = 1; ii < nx_ - 1; ii ++)
    {
        for(int jj = 1; jj < ny_ - 1; jj ++)
        {
            if(Ez_)
            {
                //Only True in Vac Once Mat/PML introduced vectorize it
                double c_eze = 1.0;
                double c_ezh = 1.0 * dt_/dx_;
                Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
            }
            else
            {
                //Only True in Vac Once Mat/PML introduced vectorize it
                //
                double c_exe = 1.0;
                double c_exh = 1.0 * dt_/dx_;
                double c_eye = 1.0;
                double c_eyh = 1.0 * dt_/dy_;
                Ey_->point(ii,jj) = c_eye * Ey_->point(ii,jj) - c_eyh * (Hz_->point(ii,jj) - Hz_->point(ii-1,jj));
                Ex_->point(ii,jj) = c_exe * Ex_->point(ii,jj) + c_exh * (Hz_->point(ii,jj) - Hz_->point(ii,jj-1));
            }
        }
    }
    for(int ii = 0; ii < dtcArr_.size(); ii ++)
        ouputField(dtcArr_[ii]);
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
