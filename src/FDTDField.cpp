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
    tcur_     = 0;
    res_      = IP.res_;
    dx_       = 1.0/res_;
    dy_       = 1.0/res_;
    dt_       = 0.99/(sqrt(pow(dx_,-2.0) + pow(dy_,-2.0)));
    nx_       = int((IP.x_size_)*res_) + 1; //Better way here; + 1 to include the 0 point
    ny_       = int((IP.y_size_)*res_) + 1; //Better way here; + 1 to include the 0 point
    srcArr_   = IP.srcArr_;
    objArr_   = IP.objArr_;
    dtcArr_   = IP.dctArr_;
    pmlArr_   = IP.pmlArr_;
    xPML_     = IP.xPml_;
    yPML_     = IP.yPml_;
    periodic_ = IP.periodic_;
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
        if(periodic_)
        {

        }
        else
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

    // Code for prefect reflectors, which we don't ever really want
    /*for(int ii = 0; ii < nx_; ii ++)
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
    }*/
    // make this more accurate

    for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
    {
        for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
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
    
    // PML part goes here
    for(int kk = 0; kk < pmlArr_.size(); kk++)
    {
        for(int ii = 2; ii < pmlArr_[kk].thickness(); ii ++)
        {
            switch (pmlArr_[kk].d())
            {
                case X:
                    for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
                    {
                        if(Ez_)
                        {
                            double eps = 1.0;
                            double kapy = 1.0; double kapz = 1.0; double sigy = 0.0; double sigz= 0.0;
                            double kapx = pmlArr_[kk].kappa(ii);
                            double sigxx = pmlArr_[kk].sigma(ii,HX);
                            double sigxy = pmlArr_[kk].sigma(ii,HY);
                            //Kappas change throughout
                            double bxstore = pmlArr_[kk].Bx_->point(ii,jj);
                            double bystore = pmlArr_[kk].By_->point(ii,jj);

                            double c_bxb = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            double c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigy*dt_)) ;
                            double c_hxh = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            double c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            double c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            double c_byb = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            double c_bye = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_)) ;
                            double c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            double c_hyb0 = (2*eps*kapy - sigy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            double c_hyb1 = (2*eps*kapy + sigy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                            pmlArr_[kk].Bx_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                            pmlArr_[kk].By_->point(ii,jj) = c_byb * pmlArr_[kk].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));

                            Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                            Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[kk].By_->point(ii,jj) - c_hyb0 * bystore;

                            Hx_->point((nx_-1) - ii,jj) = c_hxh * Hx_->point((nx_-1) - ii,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                            Hy_->point((nx_-1) - ii,jj) = c_hyh * Hy_->point((nx_-1) - ii,jj) + c_hyb1 * pmlArr_[kk].By_->point(ii,jj) - c_hyb0 * bystore;
                        }
                        else
                        {
                            //Fix
                            //Same constants as Ez so not changing them
                            /*double sigma = pmlArr_[kk].sigma(ii);
                            double bstore = pmlArr_[kk].Bz_->point(ii,jj);
                            double c_dzd = (2*pmlArr_[kk].kappa(ii) - dt_ * sigma) / (2*pmlArr_[kk].kappa(ii) + dt_ * sigma);
                            double c_dzh = 2 * dt_ / (2*pmlArr_[kk].kappa(ii) + dt_*sigma) / dx_;
                            double c_eze = c_dzd; // kappa_x v. kappa_y. I don't think there is a difference here if step size is the same
                            double c_ezd1 = 1.0; // kappa difference again with eps = 1
                            double c_ezd0 = c_dzd; // kappa difference with epsislon =1
                            pmlArr_[kk].Bz_->point(ii,jj) = c_dzd * pmlArr_[kk].Bz_->point(ii,jj) + c_dzh * ((Ey_->point(ii,jj)-Ey_->point(ii-1,jj)) - (Ex_->point(ii,jj)-Ex_->point(ii,jj-1)));
                            Hz_->point(ii,jj) = c_eze * Hz_->point(ii,jj) + c_ezd1 * pmlArr_[kk].Bz_->point(ii,jj) - c_ezd0 * bstore;
                            Hz_->point((nx_ - 1) - ii,jj) = c_eze * Hz_->point((nx_ - 1) - ii,jj) + c_ezd1 * pmlArr_[kk].Bz_->point(ii,jj) - c_ezd0 * bstore;*/
                        }
                    }
                    break;
                case Y:
                    //cout <<2 * dt_ / (2*pmlArr_[kk].kappa(ii) + dt_* pmlArr_[kk].sigma(ii)) << "\t\t" << ii << endl;
                    for(int jj = xPML_; jj < nx_ - xPML_; jj ++)
                    {
                        if(Ez_)
                        {
                            double eps = 1.0;
                            double kapx = 1.0; double kapz = 1.0; double sigx = 0.0; double sigz= 0.0;
                            double kapy = pmlArr_[kk].kappa(ii);
                            double sigyy = pmlArr_[kk].sigma(ii, HY);
                            double sigyx = pmlArr_[kk].sigma(ii, HX);
                            //Kappas change throughout
                            double bxstore = pmlArr_[kk].Bx_->point(jj,ii);
                            double bystore = pmlArr_[kk].By_->point(jj,ii);

                            double c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            double c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;
                            double c_hxh = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            double c_hxb0 = (2*eps*kapx - sigx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            double c_hxb1 = (2*eps*kapx + sigx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            double c_byb = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            double c_bye = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_)) ;
                            double c_hyh = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            double c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigx*dt_) / eps;
                            double c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigx*dt_) / eps;

                            pmlArr_[kk].Bx_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_->point(jj,ii) - c_bxe * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                            pmlArr_[kk].By_->point(jj,ii) = c_byb * pmlArr_[kk].By_->point(jj,ii) + c_bye * (Ez_->point(jj+1,ii)-Ez_->point(jj,ii));

                            Hx_->point(jj,ii) = c_hxh * Hx_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,ii) - c_hxb0 * bxstore;
                            Hy_->point(jj,ii) = c_hyh * Hy_->point(jj,ii) + c_hyb1 * pmlArr_[kk].By_->point(jj,ii) - c_hyb0 * bystore;

                            Hx_->point(jj,ny_ - 1 - ii) = c_hxh * Hx_->point(jj,ny_ - 1 - ii) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,ii) - c_hxb0 * bxstore;
                            Hy_->point(jj,ny_ - 1 - ii) = c_hyh * Hy_->point(jj,ny_ - 1 - ii) + c_hyb1 * pmlArr_[kk].By_->point(jj,ii) - c_hyb0 * bystore;
                        }
                        else
                        {
                            /*double sigma = pmlArr_[kk].sigma(ii);
                            double bstore = pmlArr_[kk].Bz_->point(jj,ii);
                            double c_dzd = (2*pmlArr_[kk].kappa(ii) - dt_ * sigma) / (2*pmlArr_[kk].kappa(ii) + dt_ * sigma);
                            double c_dzh = 2 * dt_ / (2*pmlArr_[kk].kappa(ii) + dt_*sigma) / dx_;
                            double c_eze = c_dzd; // kappa_x v. kappa_y. I don't think there is a difference here if step size is the same
                            double c_ezd1 = 1.0; // kappa difference again with eps = 1
                            double c_ezd0 = c_dzd; // kappa difference with epsislon =1

                            pmlArr_[kk].Bz_->point(jj,ii) = c_dzd * pmlArr_[kk].Bz_->point(jj,ii) + c_dzh * ((Ey_->point(jj,ii)-Ey_->point(jj-1,ii)) - (Ex_->point(jj,ii)-Ex_->point(jj,ii-1)));
                            Hz_->point(jj,(ny_ - 1) - ii) = c_eze * Hz_->point(jj,(ny_ - 1) - ii) + c_ezd1 * pmlArr_[kk].Bz_->point(jj,ii) - c_ezd0 * bstore;
                            */
                        }  
                    }
                    break;
                case Z:
                    throw logic_error("While yes we could have a thrid dimension to run, I have yet to be implimented to do such a thing. So please accept this error as my sincerest appology.");
                    break;
                default:
                    throw logic_error("I would appricate it if you stick to the standard X,Y,Z directions. While it's fun to invent new ones, it is very hard to do calculations if I don't even understand what dimension I am in. Please try again!");
                    break;
            }
        }
    }
    //cout<< << endl;
    // Better conditions will be added
    for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
    {
        for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
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
    //cout << pmlArr_[0].thickness()<< endl;
    for(int kk = 0; kk < pmlArr_.size(); kk++)
    {
        for(int ii = 0; ii < pmlArr_[kk].thickness(); ii ++)
        {
            /*double eps = 1.0;
            double kapx = 1.0; double kapz = 1.0; double sigx = 0.0; double sigz= 0.0;
            double kapy = pmlArr_[kk].kappa(ii);
            double sigy = pmlArr_[kk].sigma(ii);
            //Kappas change throughout
            double c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
            double c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_)) ;
            double c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
            double c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
            double c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
            cout << "kapy \t" << kapy <<endl;
            cout << "sigy \t" << sigy <<endl;
            cout << "c_dzd \t" << c_dzd <<endl;
            cout << "c_dzh \t" << c_dzh <<endl;
            cout << "c_eze \t" << c_eze <<endl;
            cout << "c_ezd0 \t" << c_ezd0 <<endl;
            cout << "c_ezd1 \t" << c_ezd1 <<endl;
            cout << "--------------------------------------"<< endl;*/
            switch (pmlArr_[kk].d())
            {
                case X:
                    for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
                    {
                        if(Ez_)
                        {
                            double eps = 1.0;
                            double kapx = 1.0; double kapy = 1.0; double kapz = 1.0; 
                            double sigx = pmlArr_[kk].sigma(ii,EZ); double sigy = 0.0; double sigz = 0.0;
                            //double kapx = pmlArr_[kk].kappa(ii);             
                            //Kappas change throughout
                            double dzstore = pmlArr_[kk].Dz_->point(ii,jj);

                            double c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            double c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_)) ;
                            double c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            double c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            double c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            
                            pmlArr_[kk].Dz_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                            
                            Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                            
                            Ez_->point(nx_-1-ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                        }
                        else
                        {
                            //Kappas change throughout
                            /*double sigma = pmlArr_[kk].sigma(ii);
                            double dxstore = pmlArr_[kk].Dx_->point(ii,jj);
                            double dystore = pmlArr_[kk].Dy_->point(ii,jj);
                            double c_bxb = (2*pmlArr_[kk].kappa(ii) - dt_ * sigma ) / (2*pmlArr_[kk].kappa(ii) + dt_* sigma);
                            double c_bxe = 2 * dt_ / (2*pmlArr_[kk].kappa(ii) + dt_* sigma) /dx_;
                            double c_hxh = c_bxb; // kappa_x v. kappa_y. I don't think there is a difference here if step size is the same
                            double c_hxb1 = 1.0; // kappa difference again with eps = 1
                            double c_hxb0 = c_bxb; // kappa difference with epsislon =1
                            pmlArr_[kk].Dx_->point(ii,jj) = c_bxb * pmlArr_[kk].Dx_->point(ii,jj) - c_bxe * (Hz_->point(ii,jj+1)-Hz_->point(ii,jj));
                            pmlArr_[kk].Dy_->point(ii,jj) = c_bxb * pmlArr_[kk].Dy_->point(ii,jj) - c_bxe * (Hz_->point(ii+1,jj)-Hz_->point(ii,jj));
                            Ex_->point(ii,jj) = c_hxh * Ex_->point(ii,jj) + c_hxb1 * pmlArr_[kk].Dx_->point(ii,jj) - c_hxb0 * dxstore;
                            Ey_->point(ii,jj) = c_hxh * Ey_->point(ii,jj) + c_hxb1 * pmlArr_[kk].Dy_->point(ii,jj) - c_hxb0 * dystore;
                            Ex_->point((nx_ - 1) - ii,jj) = c_hxh * Ex_->point((nx_ - 1) - ii,jj) + c_hxb1 * pmlArr_[kk].Dx_->point(ii,jj) - c_hxb0 * dxstore;
                            Ey_->point((nx_ - 1) - ii,jj) = c_hxh * Ey_->point((nx_ - 1) - ii,jj) + c_hxb1 * pmlArr_[kk].Dy_->point(ii,jj) - c_hxb0 * dystore;*/
                        }
                    }
                    break;
                case Y:
                    for(int jj = xPML_; jj < nx_ - xPML_; jj ++)
                    {
                        if(Ez_)
                        {
                            double eps = 1.0;
                            double kapx = 1.0; double kapz = 1.0; double sigx = 0.0; double sigz= 0.0;
                            double kapy = pmlArr_[kk].kappa(ii);
                            double sigy = pmlArr_[kk].sigma(ii,EZ);
                            //Kappas change throughout
                            double dzstore = pmlArr_[kk].Dz_->point(jj,ii);

                            double c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            double c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_)) ;
                            double c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            double c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            double c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            
                            pmlArr_[kk].Dz_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_->point(jj,ii) + c_dzh * ((Hy_->point(jj,ii)-Hy_->point(jj-1,ii)) - (Hx_->point(jj,ii)-Hx_->point(jj,ii-1)));
                            
                            Ez_->point(jj,ii) = c_eze * Ez_->point(jj,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,ii) - c_ezd0 * dzstore;
                            
                            Ez_->point(jj,ny_ - 1 - ii) = c_eze * Ez_->point(jj,ny_ - 1 - ii) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,ii) - c_ezd0 * dzstore;
                        }
                        else
                        {
                            //Kappas change throughout
                            /*double sigma = pmlArr_[kk].sigma(ii);
                            double dxstore = pmlArr_[kk].Dx_->point(jj,ii);
                            double dystore = pmlArr_[kk].Dy_->point(jj,ii);
                            double c_bxb = (2*pmlArr_[kk].kappa(ii) - dt_ * sigma ) / (2*pmlArr_[kk].kappa(ii) + dt_* sigma);
                            double c_bxe = 2 * dt_ / (2*pmlArr_[kk].kappa(ii) + dt_* sigma) / dx_;
                            double c_hxh = c_bxb; // kappa_x v. kappa_y. I don't think there is a difference here if step size is the same
                            double c_hxb1 = 1.0; // kappa difference again with eps = 1
                            double c_hxb0 = c_bxb; // kappa difference with epsislon =1
                            pmlArr_[kk].Dx_->point(jj,ii) = c_bxb * pmlArr_[kk].Dx_->point(jj,ii) + c_bxe * (Hz_->point(jj,ii+1)-Hz_->point(jj,ii));
                            pmlArr_[kk].Dy_->point(jj,ii) = c_bxb * pmlArr_[kk].Dy_->point(jj,ii) - c_bxe * (Hz_->point(jj+1,ii)-Hz_->point(jj,ii));
                            Ex_->point(jj,ii) = c_hxh * Ex_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Dx_->point(jj,ii) - c_hxb0 * dxstore;
                            Ey_->point(jj,ii) = c_hxh * Ey_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Dy_->point(jj,ii) - c_hxb0 * dystore;
                            Ex_->point(jj,(ny_ - 1) - ii) = c_hxh * Ex_->point(jj,(ny_ - 1) - ii) + c_hxb1 * pmlArr_[kk].Dx_->point(jj,ii) - c_hxb0 * dxstore;
                            Ey_->point(jj,(ny_ - 1) - ii) = c_hxh * Ey_->point(jj,(ny_ - 1) - ii) + c_hxb1 * pmlArr_[kk].Dy_->point(jj,ii) - c_hxb0 * dystore;*/
                        }  
                    }
                    break;
                default:
                    throw logic_error("While yes we could have a thrid dimension to run, I have yet to be implimented to do such a thing. So please accept this error as my sincerest appology.");

            }
        }
    }

    for(int ii = 0; ii < dtcArr_.size(); ii ++)
        ouputField(dtcArr_[ii]);
}