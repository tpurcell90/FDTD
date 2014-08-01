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
    nx_       = floor(static_cast<double>(res_) * static_cast<double>(IP.y_size_)+ 0.5) + 1; //Better way here; + 1 to include the 0 point
    ny_       = floor(static_cast<double>(res_) * static_cast<double>(IP.y_size_)+ 0.5) + 1; //Better way here; + 1 to include the 0 point
    srcArr_   = IP.srcArr_;
    objArr_   = IP.objArr_;
    dtcArr_   = IP.dctArr_;
    pmlArr_   = IP.pmlArr_;
    xPML_     = IP.xPml_;
    yPML_     = IP.yPml_;
    cout << xPML_ << "\t\t" << yPML_ <<endl;
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
        /*Hx_ = make_shared<Grid2D<double>>(nx_,ny_-1,dx_,dy_);
        Hy_ = make_shared<Grid2D<double>>(nx_-1,ny_,dx_,dy_);
        Ez_ = make_shared<Grid2D<double>>(nx_,ny_,dx_,dy_);
        phys_Hx_ = make_shared<Grid2D<int>>(nx_,ny_-1,dx_,dy_);
        phys_Hy_ = make_shared<Grid2D<int>>(nx_-1,ny_,dx_,dy_);
        phys_Ez_ = make_shared<Grid2D<int>>(nx_,ny_,dx_,dy_);*/
        Hx_ = make_shared<Grid2D<double>>(nx_,ny_-1,dx_,dy_);
        Hy_= make_shared<Grid2D<double>>(nx_-1,ny_,dx_,dy_);
        Ez_ = make_shared<Grid2D<double>>(nx_,ny_,dx_,dy_);

        phys_Hx_ = make_shared<Grid2D<int>>(nx_,ny_-1,dx_,dy_);
        phys_Hy_ = make_shared<Grid2D<int>>(nx_-1,ny_,dx_,dy_);
        phys_Ez_ = make_shared<Grid2D<int>>(nx_,ny_,dx_,dy_);
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

void FDTDField::step()
{
    // disregard PML's to start
    // Source
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
                
                if(ii!=nx_-1 && jj != ny_-1)
                {
                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) - c_hxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                }
                else if (ii!=nx_-1)
                {
                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                }
                else if(jj!=ny_-1)
                {
                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) - c_hxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                }
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
    //cout << 295 << endl;
    for(int kk = 0; kk < pmlArr_.size(); kk++)
    {
        //cout << pmlArr_[kk].thickness()<<endl;
        for(int ii = 0; ii < pmlArr_[kk].thickness(); ii ++)
        {
            //cout<< ii <<  "\t" << pmlArr_[kk].sigma(static_cast<double>(ii)) << "\t\t" <<  pmlArr_[kk].kappa(ii) <<endl;
            switch (pmlArr_[kk].d())
            {
                case X:

                    for(int jj =  yPML_ ; jj < ny_ - yPML_; jj ++)
                    {
                        if(Ez_)
                        {
                            double eps = 1.0;
                            double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                            double sigz = 0.0;
                            double sigxx = pmlArr_[kk].sigma(static_cast<double>(ii));
                            double sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5);
                            double sigyx = 0.0;
                            double sigyy = 0.0;
                            //Kappas change throughout
                            //cout << "store1" <<endl;
                            //cout << jj << "\t" << ny_+1<< "\t"<< pmlArr_[kk].Bx_->y()  << "\t" << pmlArr_[kk].By_->y() <<endl;
                            double bxstore = pmlArr_[kk].Bx_->point(ii,jj);
                            double bystore = pmlArr_[kk].By_->point(ii,jj);
                            double c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            double c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;
                            double c_hxh = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            double c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            double c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            double c_byb = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            double c_bye = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_)) ;
                            double c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            double c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            double c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            
                            if(jj != ny_-1 && ii != 0)
                            {
                                pmlArr_[kk].Bx_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                pmlArr_[kk].By_->point(ii,jj) = c_byb * pmlArr_[kk].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                                
                                Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                                Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[kk].By_->point(ii,jj) - c_hyb0 * bystore;
                                
                                sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) - 0.5);
                                c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                
                                bxstore = pmlArr_[kk].Bx_end_->point(ii,jj);
                                bystore = pmlArr_[kk].By_end_->point(ii,jj);

                                pmlArr_[kk].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point((nx_-1) - ii, jj+1)-Ez_->point((nx_-1) - ii,jj));
                                pmlArr_[kk].By_end_->point(ii,jj) = c_byb * pmlArr_[kk].By_end_->point(ii,jj) + c_bye * (Ez_->point((nx_-1) - ii+1, jj)-Ez_->point((nx_-1) - ii,jj));

                                Hx_->point((nx_-1) - ii,jj) = c_hxh * Hx_->point((nx_-1) - ii,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                                Hy_->point((nx_-1) - ii,jj) = c_hyh * Hy_->point((nx_-1) - ii,jj) + c_hyb1 * pmlArr_[kk].By_end_->point(ii,jj) - c_hyb0 * bystore;
                                
                            }
                            else if(ii !=0)
                            {
                                pmlArr_[kk].By_->point(ii,jj) = c_byb * pmlArr_[kk].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                                Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[kk].By_->point(ii,jj) - c_hyb0 * bystore;
                                
                                sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) - 0.5);
                                c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                bystore = pmlArr_[kk].By_end_->point(ii,jj);
                                pmlArr_[kk].By_end_->point(ii,jj) = c_byb * pmlArr_[kk].By_end_->point(ii,jj) + c_bye * (Ez_->point((nx_-1) - ii+1, jj)-Ez_->point((nx_-1) - ii,jj));
                                Hy_->point((nx_-1) - ii,jj) = c_hyh * Hy_->point((nx_-1) - ii,jj) + c_hyb1 * pmlArr_[kk].By_end_->point(ii,jj) - c_hyb0 * bystore;
                            }
                            else if (jj != ny_-1)
                            {
                                pmlArr_[kk].Bx_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                pmlArr_[kk].By_->point(ii,jj) = c_byb * pmlArr_[kk].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                                Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                                Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[kk].By_->point(ii,jj) - c_hyb0 * bystore;
                                bxstore = pmlArr_[kk].Bx_end_->point(ii,jj);
                                pmlArr_[kk].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point((nx_-1) - ii, jj+1)-Ez_->point((nx_-1) - ii,jj));
                                Hx_->point((nx_-1) - ii,jj) = c_hxh * Hx_->point((nx_-1) - ii,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                            }
                                
                        }
                        else
                        {
                            
                        }
                    }
                    break;
                case Y:
                    for(int jj = xPML_; jj < nx_ - xPML_; jj ++)
                    {
                        if(Ez_)
                        {
                            double eps = 1.0;
                            double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                            double sigz = 0.0;
                            double sigxx = 0.0;
                            double sigxy = 0.0;
                            double sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) + 0.5));
                            double sigyy = pmlArr_[kk].sigma(static_cast<double>((ii)));
                            //Kappas change throughout
                            //cout << "store1" <<endl;
                            //cout << jj << "\t" << ny_+1<< "\t"<< pmlArr_[kk].Bx_->y()  << "\t" << pmlArr_[kk].By_->y() <<endl;
                            double bxstore = pmlArr_[kk].Bx_->point(jj,ii);
                            double bystore = pmlArr_[kk].By_->point(jj,ii);
                            double c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            double c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;
                            double c_hxh = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            double c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            double c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            double c_byb = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            double c_bye = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_)) ;
                            double c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            double c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            double c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            
                            if(ii!=0 && jj != nx_-1)
                            {
                                pmlArr_[kk].Bx_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_->point(jj,ii) - c_bxe * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                                pmlArr_[kk].By_->point(jj,ii) = c_byb * pmlArr_[kk].By_->point(jj,ii) + c_bye * (Ez_->point(jj+1,ii)-Ez_->point(jj,ii));
                                //cout <<"H1"<<endl;
                                Hx_->point(jj,ii) = c_hxh * Hx_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,ii) - c_hxb0 * bxstore;
                                Hy_->point(jj,ii) = c_hyh * Hy_->point(jj,ii) + c_hyb1 * pmlArr_[kk].By_->point(jj,ii) - c_hyb0 * bystore;

                                sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) - 0.5));
                                c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;

                                bxstore = pmlArr_[kk].Bx_end_->point(jj,ii);
                                bystore = pmlArr_[kk].By_end_->point(jj,ii);
                                pmlArr_[kk].Bx_end_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_end_->point(jj,ii) - c_bxe * (Ez_->point(jj, (ny_-1)-ii+1)-Ez_->point(jj,(ny_-1)-ii));
                                pmlArr_[kk].By_end_->point(jj,ii) = c_byb * pmlArr_[kk].By_end_->point(jj,ii) + c_bye * (Ez_->point(jj+1, (ny_-1)-ii)-Ez_->point(jj,(ny_-1)-ii));
                                //cout <<"H2"<<endl;
                                Hx_->point(jj,(ny_-1)-ii) = c_hxh * Hx_->point(jj,(ny_-1)-ii) + c_hxb1 * pmlArr_[kk].Bx_end_->point(jj,ii) - c_hxb0 * bxstore;
                                Hy_->point(jj,(ny_-1)-ii) = c_hyh * Hy_->point(jj,(ny_-1)-ii) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,ii) - c_hyb0 * bystore;

                            }
                            else if (jj != nx_-1)
                            {
                                pmlArr_[kk].Bx_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_->point(jj,ii) - c_bxe * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                                pmlArr_[kk].By_->point(jj,ii) = c_byb * pmlArr_[kk].By_->point(jj,ii) + c_bye * (Ez_->point(jj+1,ii)-Ez_->point(jj,ii));
                                //cout <<"H1"<<endl;
                                Hx_->point(jj,ii) = c_hxh * Hx_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,ii) - c_hxb0 * bxstore;
                                Hy_->point(jj,ii) = c_hyh * Hy_->point(jj,ii) + c_hyb1 * pmlArr_[kk].By_->point(jj,ii) - c_hyb0 * bystore;

                                bystore = pmlArr_[kk].By_end_->point(jj,ii);
                                pmlArr_[kk].By_end_->point(jj,ii) = c_byb * pmlArr_[kk].By_end_->point(jj,ii) + c_bye * (Ez_->point(jj+1, (ny_-1)-ii)-Ez_->point(jj,(ny_-1)-ii));
                                Hy_->point(jj,(ny_-1)-ii) = c_hyh * Hy_->point(jj,(ny_-1)-ii) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,ii) - c_hyb0 * bystore;

                            }
                            else if(ii != 0)
                            {
                                pmlArr_[kk].Bx_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_->point(jj,ii) - c_bxe * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                                //cout <<"H1"<<endl;
                                Hx_->point(jj,ii) = c_hxh * Hx_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,ii) - c_hxb0 * bxstore;

                                sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) - 0.5));
                                c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;
                                
                                bxstore = pmlArr_[kk].Bx_end_->point(jj,ii);
                                pmlArr_[kk].Bx_end_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_end_->point(jj,ii) - c_bxe * (Ez_->point(jj, (ny_-1)-ii+1)-Ez_->point(jj,(ny_-1)-ii));
                                //cout <<"H2"<<endl;
                                Hx_->point(jj,(ny_-1)-ii) = c_hxh * Hx_->point(jj,(ny_-1)-ii) + c_hxb1 * pmlArr_[kk].Bx_end_->point(jj,ii) - c_hxb0 * bxstore;
                            }
                            
                            
                        }
                        else
                        {
                           
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
    //Corners: Top Left -> X0, Top Right ->Xend, botLeft ->Y0, Bot Right yEnd
    //fix to match ii being x always
    if(pmlArr_.size() == 2)
    {
        for(int ii = 0; ii < pmlArr_[0].thickness(); ii++)
        {
            for(int jj = 0; jj < pmlArr_[1].thickness(); jj++)
            {
                double kapz = 1.0; double sigz = 0.0; double eps = 1.0;
                switch(pmlArr_[0].d())
                {
                    case X:
                        if(Ez_)
                        {
                            double kapx = 1.0; double kapy = 1.0; 
                            double sigxx = pmlArr_[0].sigma(static_cast<double>(ii));
                            double sigxy = pmlArr_[0].sigma(static_cast<double>(ii) + 0.5);
                            double sigyx = pmlArr_[1].sigma(static_cast<double>(jj) + 0.5);
                            double sigyy = pmlArr_[1].sigma(static_cast<double>(jj));
                            //double sigyx = 0.0;
                            //double sigyy = 0.0;                            
                            // set all the constants
                            double c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            double c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;
                            double c_hxh = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            double c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            double c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            double c_byb = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            double c_bye = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_)) ;
                            double c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            double c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            double c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                            // Bot Left
                            double bxstore = pmlArr_[0].Bx_->point(ii,jj);
                            double bystore = pmlArr_[0].By_->point(ii,jj);
                            pmlArr_[0].Bx_->point(ii,jj) = c_bxb * pmlArr_[0].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                            pmlArr_[0].By_->point(ii,jj) = c_byb * pmlArr_[0].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                            
                            Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[0].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                            Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[0].By_->point(ii,jj) - c_hyb0 * bystore;
                            // Bot Right
                            sigxy = pmlArr_[0].sigma(static_cast<double>(ii) - 0.5);
                            c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            if(ii != 0)
                            {
                                bxstore = pmlArr_[0].Bx_end_->point(ii,jj);
                                bystore = pmlArr_[0].By_end_->point(ii,jj);
                                
                                pmlArr_[0].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[0].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point(nx_-1-ii,jj+1)-Ez_->point(nx_-1-ii,jj));
                                pmlArr_[0].By_end_->point(ii,jj) = c_byb * pmlArr_[0].By_end_->point(ii,jj) + c_bye * (Ez_->point(nx_-1-ii+1,jj)-Ez_->point(nx_-1-ii,jj));
                                
                                Hx_->point(nx_-1-ii,jj) = c_hxh * Hx_->point(nx_-1-ii,jj) + c_hxb1 * pmlArr_[0].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                                Hy_->point(nx_-1-ii,jj) = c_hyh * Hy_->point(nx_-1-ii,jj) + c_hyb1 * pmlArr_[0].By_end_->point(ii,jj) - c_hyb0 * bystore;
                            }
                            else
                            {
                                bxstore = pmlArr_[0].Bx_end_->point(ii,jj);
                                pmlArr_[0].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[0].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point(nx_-1-ii,jj+1)-Ez_->point(nx_-1-ii,jj));
                                Hx_->point(nx_-1-ii,jj) = c_hxh * Hx_->point(nx_-1-ii,jj) + c_hxb1 * pmlArr_[0].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                            }
                            //Top Right
                            sigyx = pmlArr_[1].sigma(static_cast<double>(jj) - 0.5);
                            c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            if(ii != 0 && jj !=0)
                            {
                                bxstore = pmlArr_[0].Bx_end_->point(ii,ny_ - 1 - jj);
                                bystore = pmlArr_[0].By_end_->point(ii,ny_ - 1 - jj);
                                pmlArr_[0].Bx_end_->point(ii,ny_ - 1 - jj) = c_bxb * pmlArr_[0].Bx_end_->point(ii,ny_ - 1 - jj) - c_bxe * (Ez_->point(nx_-1-ii,ny_ - 1 - jj+1)-Ez_->point(nx_-1-ii,ny_ - 1 - jj));
                                pmlArr_[0].By_end_->point(ii,ny_ - 1 - jj) = c_byb * pmlArr_[0].By_end_->point(ii,ny_ - 1 - jj) + c_bye * (Ez_->point(nx_-1-ii+1,ny_ - 1 - jj)-Ez_->point(nx_-1-ii,ny_ - 1 - jj));
                                Hx_->point(nx_-1-ii,ny_ - 1 - jj) = c_hxh * Hx_->point(nx_-1-ii,ny_ - 1 - jj) + c_hxb1 * pmlArr_[0].Bx_end_->point(ii,ny_ - 1 - jj) - c_hxb0 * bxstore;
                                Hy_->point(nx_-1-ii,ny_ - 1 - jj) = c_hyh * Hy_->point(nx_-1-ii,ny_ - 1 - jj) + c_hyb1 * pmlArr_[0].By_end_->point(ii,ny_ - 1 - jj) - c_hyb0 * bystore;
                            }
                            else if(jj!=0)
                            {
                                bxstore = pmlArr_[0].Bx_end_->point(ii,jj);
                                pmlArr_[0].Bx_end_->point(ii,ny_ - 1 - jj) = c_bxb * pmlArr_[0].Bx_end_->point(ii,ny_ - 1 - jj) - c_bxe * (Ez_->point(nx_-1-ii,ny_ - 1 - jj+1)-Ez_->point(nx_-1-ii,ny_ - 1 - jj));
                                Hx_->point(nx_-1-ii,ny_ - 1 - jj) = c_hxh * Hx_->point(nx_-1-ii,ny_ - 1 - jj) + c_hxb1 * pmlArr_[0].Bx_end_->point(ii,ny_ - 1 - jj) - c_hxb0 * bxstore;
                            }
                            else if (ii!=0)
                            {
                                bystore = pmlArr_[0].By_end_->point(ii,ny_ - 1 - jj);
                                pmlArr_[0].By_end_->point(ii,ny_ - 1 - jj) = c_byb * pmlArr_[0].By_end_->point(ii,ny_ - 1 - jj) + c_bye * (Ez_->point(nx_-1-ii+1,ny_ - 1 - jj)-Ez_->point(nx_-1-ii,ny_ - 1 - jj));
                                Hy_->point(nx_-1-ii,ny_ - 1 - jj) = c_hyh * Hy_->point(nx_-1-ii,ny_ - 1 - jj) + c_hyb1 * pmlArr_[0].By_end_->point(ii,ny_ - 1 - jj) - c_hyb0 * bystore;
                            }
                            //Top Left
                            sigxy = pmlArr_[0].sigma(static_cast<double>(ii) + 0.5);
                            c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            if(jj !=0)
                            {
                                bxstore = pmlArr_[0].Bx_->point(ii,ny_ - 1 - jj);
                                bystore = pmlArr_[0].By_->point(ii,ny_ - 1 - jj);
                                pmlArr_[0].Bx_->point(ii,ny_ - 1 - jj) = c_bxb * pmlArr_[0].Bx_->point(ii,ny_ - 1 - jj) - c_bxe * (Ez_->point(ii,ny_ - 1 - jj+1)-Ez_->point(ii,ny_ - 1 - jj));
                                pmlArr_[0].By_->point(ii,ny_ - 1 - jj) = c_byb * pmlArr_[0].By_->point(ii,ny_ - 1 - jj) + c_bye * (Ez_->point(ii+1,ny_ - 1 - jj)-Ez_->point(ii,ny_ - 1 - jj));
                                Hx_->point(ii,ny_ - 1 - jj) = c_hxh * Hx_->point(ii,ny_ - 1 - jj) + c_hxb1 * pmlArr_[0].Bx_->point(ii,ny_ - 1 - jj) - c_hxb0 * bxstore;
                                Hy_->point(ii,ny_ - 1 - jj) = c_hyh * Hy_->point(ii,ny_ - 1 - jj) + c_hyb1 * pmlArr_[0].By_->point(ii,ny_ - 1 - jj) - c_hyb0 * bystore;
                            }
                            else
                            {
                                bystore = pmlArr_[0].By_->point(ii,ny_ - 1 - jj);
                                pmlArr_[0].By_->point(ii,ny_ - 1 - jj) = c_byb * pmlArr_[0].By_->point(ii,ny_ - 1 - jj) + c_bye * (Ez_->point(ii+1,ny_ - 1 - jj)-Ez_->point(ii,ny_ - 1 - jj));
                                Hy_->point(ii,ny_ - 1 - jj) = c_hyh * Hy_->point(ii,ny_ - 1 - jj) + c_hyb1 * pmlArr_[0].By_->point(ii,ny_ - 1 - jj) - c_hyb0 * bystore;

                            }
                        }
                        break;
                    case Y:
                        if(Ez_)
                        {
                            double kapx = 1.0; double kapy = 1.0; 
                            double sigxx = pmlArr_[1].sigma(static_cast<double>(ii));
                            double sigxy = pmlArr_[1].sigma(static_cast<double>(ii) + 0.5);
                            double sigyx = pmlArr_[0].sigma(static_cast<double>(jj) + 0.5);
                            double sigyy = pmlArr_[0].sigma(static_cast<double>(jj));
                            // set all the constants
                            double c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            double c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;
                            double c_hxh = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            double c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            double c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            double c_byb = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            double c_bye = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_)) ;
                            double c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            double c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            double c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                            // Do the Bot left Corner
                            double bxstore = pmlArr_[1].Bx_->point(ii,jj);
                            double bystore = pmlArr_[1].By_->point(ii,jj);
                            pmlArr_[1].Bx_->point(ii,jj) = c_bxb * pmlArr_[1].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                            pmlArr_[1].By_->point(ii,jj) = c_byb * pmlArr_[1].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                            Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[1].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                            Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[1].By_->point(ii,jj) - c_hyb0 * bystore;

                            // Do the Bot Right Corner
                            sigxy = pmlArr_[1].sigma(static_cast<double>(ii) - 0.5);
                            c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            if(ii != 0)
                            {
                                bxstore = pmlArr_[1].Bx_->point((nx_-1) - ii,jj);
                                bystore = pmlArr_[1].By_->point((nx_-1) - ii,jj);
                                pmlArr_[1].Bx_->point((nx_-1) - ii,jj) = c_bxb * pmlArr_[1].Bx_->point((nx_-1) - ii,jj) - c_bxe * (Ez_->point((nx_-1)-ii,jj+1)-Ez_->point((nx_-1)-ii,jj));
                                pmlArr_[1].By_->point((nx_-1) - ii,jj) = c_byb * pmlArr_[1].By_->point((nx_-1) - ii,jj) + c_bye * (Ez_->point((nx_-1)-ii+1,jj)-Ez_->point((nx_-1)-ii,jj));
                                Hx_->point((nx_-1)-ii,jj) = c_hxh * Hx_->point((nx_-1)-ii,jj) + c_hxb1 * pmlArr_[1].Bx_->point((nx_-1) - ii,jj) - c_hxb0 * bxstore;
                                Hy_->point((nx_-1)-ii,jj) = c_hyh * Hy_->point((nx_-1)-ii,jj) + c_hyb1 * pmlArr_[1].By_->point((nx_-1) - ii,jj) - c_hyb0 * bystore;
                            }
                            else
                            {
                                bxstore = pmlArr_[1].Bx_->point((nx_-1) - ii,jj);
                                pmlArr_[1].Bx_->point((nx_-1) - ii,jj) = c_bxb * pmlArr_[1].Bx_->point((nx_-1) - ii,jj) - c_bxe * (Ez_->point((nx_-1)-ii,jj+1)-Ez_->point((nx_-1)-ii,jj));
                                Hx_->point((nx_-1)-ii,jj) = c_hxh * Hx_->point((nx_-1)-ii,jj) + c_hxb1 * pmlArr_[1].Bx_->point((nx_-1) - ii,jj) - c_hxb0 * bxstore;
                            }
                            // Do the Top Right Corner
                            sigyx = pmlArr_[0].sigma(static_cast<double>(jj) - 0.5);
                            c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;
                            if(ii!=0 && jj !=0)
                            {
                                bxstore = pmlArr_[1].Bx_end_->point((nx_-1) - ii,jj);
                                bystore = pmlArr_[1].By_end_->point((nx_-1) - ii,jj);
                                pmlArr_[1].Bx_end_->point((nx_-1) - ii,jj) = c_bxb * pmlArr_[1].Bx_end_->point((nx_-1) - ii,jj) - c_bxe * (Ez_->point((nx_-1)-ii,(ny_-1)-jj+1)-Ez_->point((nx_-1)-ii,(ny_-1)-jj));
                                pmlArr_[1].By_end_->point((nx_-1) - ii,jj) = c_byb * pmlArr_[1].By_end_->point((nx_-1) - ii,jj) + c_bye * (Ez_->point((nx_-1)-ii+1,(ny_-1)-jj)-Ez_->point((nx_-1)-ii,(ny_-1)-jj));
                                Hx_->point((nx_-1)-ii,(ny_-1)-jj) = c_hxh * Hx_->point((nx_-1)-ii,(ny_-1)-jj) + c_hxb1 * pmlArr_[1].Bx_end_->point((nx_-1) - ii,jj) - c_hxb0 * bxstore;
                                Hy_->point((nx_-1)-ii,(ny_-1)-jj) = c_hyh * Hy_->point((nx_-1)-ii,(ny_-1)-jj) + c_hyb1 * pmlArr_[1].By_end_->point((nx_-1) - ii,jj) - c_hyb0 * bystore;
                            }
                            else if(jj!=0)
                            {
                                bxstore = pmlArr_[1].Bx_end_->point((nx_-1) - ii,jj);
                                pmlArr_[1].Bx_end_->point((nx_-1) - ii,jj) = c_bxb * pmlArr_[1].Bx_end_->point((nx_-1) - ii,jj) - c_bxe * (Ez_->point((nx_-1)-ii,(ny_-1)-jj+1)-Ez_->point((nx_-1)-ii,(ny_-1)-jj));
                                Hx_->point((nx_-1)-ii,(ny_-1)-jj) = c_hxh * Hx_->point((nx_-1)-ii,(ny_-1)-jj) + c_hxb1 * pmlArr_[1].Bx_end_->point((nx_-1) - ii,jj) - c_hxb0 * bxstore;
                            }
                            else if(ii!=0)
                            {
                                bystore = pmlArr_[1].By_end_->point((nx_-1) - ii,jj);
                                pmlArr_[1].By_end_->point((nx_-1) - ii,jj) = c_byb * pmlArr_[1].By_end_->point((nx_-1) - ii,jj) + c_bye * (Ez_->point((nx_-1)-ii+1,(ny_-1)-jj)-Ez_->point((nx_-1)-ii,(ny_-1)-jj));
                                Hy_->point((nx_-1)-ii,(ny_-1)-jj) = c_hyh * Hy_->point((nx_-1)-ii,(ny_-1)-jj) + c_hyb1 * pmlArr_[1].By_end_->point((nx_-1) - ii,jj) - c_hyb0 * bystore;
                            }
                            // Do the Top Left Corner
                            sigxy = pmlArr_[1].sigma(static_cast<double>(ii) + 0.5);
                            c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            if(jj!=0)
                            {
                                bxstore = pmlArr_[1].Bx_end_->point(ii,jj);
                                bystore = pmlArr_[1].By_end_->point(ii,jj);
                                pmlArr_[1].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[1].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point(ii,(ny_-1)-jj+1)-Ez_->point(ii,(ny_-1)-jj));
                                pmlArr_[1].By_end_->point(ii,jj) = c_byb * pmlArr_[1].By_end_->point(ii,jj) + c_bye * (Ez_->point(ii+1,(ny_-1)-jj)-Ez_->point(ii,(ny_-1)-jj));
                                Hx_->point(ii,(ny_-1)-jj) = c_hxh * Hx_->point(ii,(ny_-1)-jj) + c_hxb1 * pmlArr_[1].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                                Hy_->point(ii,(ny_-1)-jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[1].By_end_->point(ii,jj) - c_hyb0 * bystore;
                            }
                            else
                            {
                                bystore = pmlArr_[1].By_end_->point(ii,jj);
                                pmlArr_[1].By_end_->point(ii,jj) = c_byb * pmlArr_[1].By_end_->point(ii,jj) + c_bye * (Ez_->point(ii+1,(ny_-1)-jj)-Ez_->point(ii,(ny_-1)-jj));
                                Hy_->point(ii,(ny_-1)-jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[1].By_end_->point(ii,jj) - c_hyb0 * bystore;
                            }

                        }
                        break;
                    case Z:
                        throw logic_error("Z is not yet defined");
                        break;
                    default:
                        throw logic_error("hit switch default");
                        break;
                }
            }
        }
    }
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
                if(ii != nx_-1 && ii != 0)
                {
                    if(jj!=0 && jj!=ny_-1)
                    {
                        Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                    }
                    else if (jj != ny_-1)
                    {
                        Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-0));                        
                    }
                    else if (jj != 0)
                    {
                        Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (0-Hx_->point(ii,jj-1)));
                    }
                }
                else if(ii != nx_-1)
                {
                    if(jj!=0 && jj!=ny_-1)
                    {
                        Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-0) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                    }
                    else if (jj != ny_-1)
                    {
                        Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-0) - (Hx_->point(ii,jj)-0));                        
                    }
                    else if (jj != 0)
                    {
                        Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-0) - (0-Hx_->point(ii,jj-1)));
                    }
                }
                else if (ii !=0)
                {
                    if(jj!=0 && jj!=ny_-1)
                    {
                        Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((0-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                    }
                    else if (jj != ny_-1)
                    {
                        Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((0-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-0));                        
                    }
                    else if (jj != 0)
                    {
                        Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((0-Hy_->point(ii-1,jj)) - (0-Hx_->point(ii,jj-1)));
                    }
                }

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
    for(int kk = 0; kk < pmlArr_.size(); kk++)
    {
        for(int ii = 0; ii < pmlArr_[kk].thickness(); ii ++)
        {
            switch (pmlArr_[kk].d())
            {
                case X:
                    for(int jj = yPML_ ; jj < ny_ -yPML_; jj ++)
                    {
                        if(Ez_)
                        {
                            double eps = 1.0;
                            double kapx = 1.0; double kapy = 1.0; double kapz = 1.0; 
                            double sigx = pmlArr_[kk].sigma(static_cast<double>(ii));
                            double sigy = 0.0; double sigz = 0.0;
                            //double kapx = pmlArr_[kk].kappa(ii);
                            //Kappas change throughout
                            double dzstore = pmlArr_[kk].Dz_->point(ii,jj);

                            double c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            double c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_)) ;
                            double c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            double c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            double c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;

                            if(ii !=0)
                            {
                                if(jj !=0 && jj != ny_-1)
                                {
                                    pmlArr_[kk].Dz_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    dzstore = pmlArr_[kk].Dz_end_->point(ii,jj);
                                    pmlArr_[kk].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point((nx_-1) - ii,jj)-Hy_->point((nx_-1) - ii-1,jj)) - (Hx_->point((nx_-1) - ii,jj)-Hx_->point((nx_-1) - ii,jj-1)));
                                    Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                else if(jj!= 0)
                                {
                                    pmlArr_[kk].Dz_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (0-Hx_->point(ii,jj-1)));
                                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    dzstore = pmlArr_[kk].Dz_end_->point(ii,jj);
                                    pmlArr_[kk].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point((nx_-1) - ii,jj)-Hy_->point((nx_-1) - ii-1,jj)) - (0-Hx_->point((nx_-1) - ii,jj-1)));
                                    Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                else if (jj != ny_-1)
                                {
                                    pmlArr_[kk].Dz_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-0));
                                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    dzstore = pmlArr_[kk].Dz_end_->point(ii,jj);
                                    pmlArr_[kk].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point((nx_-1) - ii,jj)-Hy_->point((nx_-1) - ii-1,jj)) - (Hx_->point((nx_-1) - ii,jj)-0));
                                    Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                
                            }
                            else
                            {
                                if(jj !=0 && jj != ny_-1)
                                {
                                    pmlArr_[kk].Dz_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj) - 0) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1))); // 0 is boundary condtion
                                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    dzstore = pmlArr_[kk].Dz_end_->point(ii,jj);
                                    pmlArr_[kk].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,jj) + c_dzh * ((0 -Hy_->point((nx_-1) - ii-1,jj)) - (Hx_->point((nx_-1) - ii,jj)-Hx_->point((nx_-1) - ii,jj-1))); // 0 is boundary condition
                                    Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                else if(jj != 0)
                                {
                                    pmlArr_[kk].Dz_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj) - 0) - (0-Hx_->point(ii,jj-1))); // 0 is boundary condtion
                                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    dzstore = pmlArr_[kk].Dz_end_->point(ii,jj);
                                    pmlArr_[kk].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,jj) + c_dzh * ((0 -Hy_->point((nx_-1) - ii-1,jj)) - (0-Hx_->point((nx_-1) - ii,jj-1))); // 0 is boundary condition
                                    Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                else if(jj != ny_-1)
                                {
                                    pmlArr_[kk].Dz_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj) - 0) - (Hx_->point(ii,jj)-0)); // 0 is boundary condtion
                                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    dzstore = pmlArr_[kk].Dz_end_->point(ii,jj);
                                    pmlArr_[kk].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,jj) + c_dzh * ((0 -Hy_->point((nx_-1) - ii-1,jj)) - (Hx_->point((nx_-1) - ii,jj)-0)); // 0 is boundary condition
                                    Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                
                            }
                        }
                        else
                        {
                            
                        }
                    }
                    break;
                case Y:
                    for(int jj = xPML_; jj < nx_ - xPML_; jj ++)
                    {
                        if(Ez_)
                        {
                            double eps = 1.0;
                            double kapx = 1.0; double kapy = 1.0; double kapz = 1.0; 
                            double sigx = 0.0;
                            double sigy = pmlArr_[kk].sigma(static_cast<double>(ii)); double sigz = 0.0;
                            //double kapx = pmlArr_[kk].kappa(ii);
                            //Kappas change throughout
                            double dzstore = pmlArr_[kk].Dz_->point(jj,ii);

                            double c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            double c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_)) ;
                            double c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            double c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            double c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            if(ii != 0)
                            {
                                if(jj !=0 && jj != nx_-1)
                                {
                                    pmlArr_[kk].Dz_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_->point(jj,ii) + c_dzh * ((Hy_->point(jj,ii)-Hy_->point(jj-1,ii)) - (Hx_->point(jj,ii)-Hx_->point(jj,ii-1)));
                                    Ez_->point(jj,ii) = c_eze * Ez_->point(jj,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,ii) - c_ezd0 * dzstore;
                                    
                                    dzstore = pmlArr_[kk].Dz_end_->point(jj,ii);
                                    pmlArr_[kk].Dz_end_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_end_->point(jj,ii) + c_dzh * ((Hy_->point(jj,(ny_-1)-ii) - Hy_->point(jj - 1,(ny_-1)-ii)) - (Hx_->point(jj,(ny_-1)-ii) - Hx_->point(jj,(ny_-1)-ii-1)));
                                    Ez_->point(jj,(ny_-1)-ii) = c_eze * Ez_->point(jj,(ny_-1)-ii) + c_ezd1 * pmlArr_[kk].Dz_end_->point(jj,ii) - c_ezd0 * dzstore;
                                }
                                else if(jj != 0)
                                {
                                    pmlArr_[kk].Dz_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_->point(jj,ii) + c_dzh * ((0-Hy_->point(jj-1,ii)) - (Hx_->point(jj,ii)-Hx_->point(jj,ii-1)));
                                    Ez_->point(jj,ii) = c_eze * Ez_->point(jj,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,ii) - c_ezd0 * dzstore;
                                    
                                    dzstore = pmlArr_[kk].Dz_end_->point(jj,ii);
                                    pmlArr_[kk].Dz_end_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_end_->point(jj,ii) + c_dzh * ((0 - Hy_->point(jj - 1,(ny_-1)-ii)) - (Hx_->point(jj,(ny_-1)-ii) - Hx_->point(jj,(ny_-1)-ii-1)));
                                    Ez_->point(jj,(ny_-1)-ii) = c_eze * Ez_->point(jj,(ny_-1)-ii) + c_ezd1 * pmlArr_[kk].Dz_end_->point(jj,ii) - c_ezd0 * dzstore;
                                }
                                else if(jj != nx_-1)
                                {
                                    pmlArr_[kk].Dz_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_->point(jj,ii) + c_dzh * ((Hy_->point(jj,ii)-0) - (Hx_->point(jj,ii)-Hx_->point(jj,ii-1)));
                                    Ez_->point(jj,ii) = c_eze * Ez_->point(jj,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,ii) - c_ezd0 * dzstore;
                                    
                                    dzstore = pmlArr_[kk].Dz_end_->point(jj,ii);
                                    pmlArr_[kk].Dz_end_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_end_->point(jj,ii) + c_dzh * ((Hy_->point(jj,(ny_-1)-ii) - 0) - (Hx_->point(jj,(ny_-1)-ii) - Hx_->point(jj,(ny_-1)-ii-1)));
                                    Ez_->point(jj,(ny_-1)-ii) = c_eze * Ez_->point(jj,(ny_-1)-ii) + c_ezd1 * pmlArr_[kk].Dz_end_->point(jj,ii) - c_ezd0 * dzstore;
                                }
                            }
                            else
                            {
                                if(jj!= 0 && jj != nx_-1)
                                {
                                    pmlArr_[kk].Dz_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_->point(jj,ii) + c_dzh * ((Hy_->point(jj,ii)-Hy_->point(jj-1,ii)) - (Hx_->point(jj,ii)- 0)); // 0 is boundary condition
                                    Ez_->point(jj,ii) = c_eze * Ez_->point(jj,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,ii) - c_ezd0 * dzstore;
                                    
                                    dzstore = pmlArr_[kk].Dz_end_->point(jj,ii);
                                    pmlArr_[kk].Dz_end_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_end_->point(jj,ii) + c_dzh * ((Hy_->point(jj,(ny_-1)-ii) - Hy_->point(jj - 1,(ny_-1)-ii)) - (0 - Hx_->point(jj,(ny_-1)-ii-1))); // 0 is boundary condtion
                                    Ez_->point(jj,(ny_-1)-ii) = c_eze * Ez_->point(jj,(ny_-1)-ii) + c_ezd1 * pmlArr_[kk].Dz_end_->point(jj,ii) - c_ezd0 * dzstore;
                                }
                                else if(jj != 0)
                                {
                                    pmlArr_[kk].Dz_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_->point(jj,ii) + c_dzh * ((0-Hy_->point(jj-1,ii)) - (Hx_->point(jj,ii)- 0)); // 0 is boundary condition
                                    Ez_->point(jj,ii) = c_eze * Ez_->point(jj,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,ii) - c_ezd0 * dzstore;
                                    
                                    dzstore = pmlArr_[kk].Dz_end_->point(jj,ii);
                                    pmlArr_[kk].Dz_end_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_end_->point(jj,ii) + c_dzh * ((0 - Hy_->point(jj - 1,(ny_-1)-ii)) - (0 - Hx_->point(jj,(ny_-1)-ii-1))); // 0 is boundary condtion
                                    Ez_->point(jj,(ny_-1)-ii) = c_eze * Ez_->point(jj,(ny_-1)-ii) + c_ezd1 * pmlArr_[kk].Dz_end_->point(jj,ii) - c_ezd0 * dzstore;
                                }
                                else if(jj != nx_-1)
                                {
                                    pmlArr_[kk].Dz_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_->point(jj,ii) + c_dzh * ((Hy_->point(jj,ii)-0) - (Hx_->point(jj,ii)- 0)); // 0 is boundary condition
                                    Ez_->point(jj,ii) = c_eze * Ez_->point(jj,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,ii) - c_ezd0 * dzstore;
                                    
                                    dzstore = pmlArr_[kk].Dz_end_->point(jj,ii);
                                    pmlArr_[kk].Dz_end_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_end_->point(jj,ii) + c_dzh * ((Hy_->point(jj,(ny_-1)-ii) - 0) - (0 - Hx_->point(jj,(ny_-1)-ii-1))); // 0 is boundary condtion
                                    Ez_->point(jj,(ny_-1)-ii) = c_eze * Ez_->point(jj,(ny_-1)-ii) + c_ezd1 * pmlArr_[kk].Dz_end_->point(jj,ii) - c_ezd0 * dzstore;
                                }
                            }
                        }
                        else
                        {
                            
                        }  
                    }
                    break;
                default:
                    throw logic_error("While yes we could have a thrid dimension to run, I have yet to be implimented to do such a thing. So please accept this error as my sincerest appology.");

            }
        }
    }
    // Corners
    if(pmlArr_.size() == 2)
    {
        for(int ii = 0; ii < pmlArr_[0].thickness(); ii++)
        {
            for(int jj = 0; jj < pmlArr_[1].thickness(); jj++)
            {
                double eps = 1.0;
                switch(pmlArr_[0].d())
                {
                    case X:
                        if(Ez_)
                        {
                            double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                            double sigx = pmlArr_[0].sigma(static_cast<double>(ii));
                            double sigy = pmlArr_[1].sigma(static_cast<double>(jj));                            
                            //double sigy=0.0;
                            double sigz = 0.0;

                            double c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            double c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_)) ;
                            double c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            double c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            double c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            double dzstore = pmlArr_[0].Dz_->point(ii,jj);
                            //Bot Left
                            if(ii !=0)
                            {
                                if(jj !=0)
                                {
                                    pmlArr_[0].Dz_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,jj);
                                    //pmlArr_[0].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point((nx_-1) - ii,jj)-Hy_->point((nx_-1) - ii-1,jj)) - (Hx_->point((nx_-1) - ii,jj)-Hx_->point((nx_-1) - ii,jj-1)));
                                    //Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                else
                                {
                                    pmlArr_[0].Dz_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-0));
                                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,jj);
                                    //pmlArr_[0].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point((nx_-1) - ii,jj)-Hy_->point((nx_-1) - ii-1,jj)) - (Hx_->point((nx_-1) - ii,jj)-0));
                                    //Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                
                            }
                            else
                            {
                                if(jj !=0 )
                                {
                                    pmlArr_[0].Dz_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj) - 0) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1))); // 0 is boundary condtion
                                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,jj);
                                    //pmlArr_[0].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,jj) + c_dzh * ((0 -Hy_->point((nx_-1) - ii-1,jj)) - (Hx_->point((nx_-1) - ii,jj)-Hx_->point((nx_-1) - ii,jj-1))); // 0 is boundary condition
                                    //Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                else
                                {
                                    pmlArr_[0].Dz_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj) - 0) - (Hx_->point(ii,jj)-0)); // 0 is boundary condtion
                                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,jj);
                                    //pmlArr_[0].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,jj) + c_dzh * ((0 -Hy_->point((nx_-1) - ii-1,jj)) - (Hx_->point((nx_-1) - ii,jj)-0)); // 0 is boundary condition
                                    //Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                
                            }
                            // Bot Right
                            dzstore = pmlArr_[0].Dz_end_->point(ii,jj);
                            if(ii !=0)
                            {
                                if(jj !=0)
                                {
                                    //pmlArr_[0].Dz_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                    //Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,jj);
                                    pmlArr_[0].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point((nx_-1) - ii,jj)-Hy_->point((nx_-1) - ii-1,jj)) - (Hx_->point((nx_-1) - ii,jj)-Hx_->point((nx_-1) - ii,jj-1)));
                                    Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                else
                                {
                                    //pmlArr_[0].Dz_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-0));
                                    //Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,jj);
                                    pmlArr_[0].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point((nx_-1) - ii,jj)-Hy_->point((nx_-1) - ii-1,jj)) - (Hx_->point((nx_-1) - ii,jj)-0));
                                    Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                
                            }
                            else
                            {
                                if(jj !=0 )
                                {
                                    //pmlArr_[0].Dz_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj) - 0) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1))); // 0 is boundary condtion
                                    //Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,jj);
                                    pmlArr_[0].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,jj) + c_dzh * ((0 -Hy_->point((nx_-1) - ii-1,jj)) - (Hx_->point((nx_-1) - ii,jj)-Hx_->point((nx_-1) - ii,jj-1))); // 0 is boundary condition
                                    Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                else
                                {
                                    //pmlArr_[0].Dz_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj) - 0) - (Hx_->point(ii,jj)-0)); // 0 is boundary condtion
                                    //Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,jj);
                                    pmlArr_[0].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,jj) + c_dzh * ((0 -Hy_->point((nx_-1) - ii-1,jj)) - (Hx_->point((nx_-1) - ii,jj)-0)); // 0 is boundary condition
                                    Ez_->point((nx_-1) - ii,jj) = c_eze * Ez_->point((nx_-1) - ii,jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                
                            }
                            //Top Right
                            dzstore = pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj);
                            if(ii !=0)
                            {
                                if(jj!=0)
                                {
                                    //pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) + c_dzh * ((Hy_->point(ii,ny_ -1 - jj)-Hy_->point(ii-1,ny_ -1 - jj)) - (Hx_->point(ii,ny_ -1 - jj)-Hx_->point(ii,ny_ -1 - jj-1)));
                                    //Ez_->point(ii,ny_ -1 - jj) = c_eze * Ez_->point(ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj);
                                    pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) + c_dzh * ((Hy_->point((nx_-1) - ii,ny_ -1 - jj)-Hy_->point((nx_-1) - ii-1,ny_ -1 - jj)) - (Hx_->point((nx_-1) - ii,ny_ -1 - jj)-Hx_->point((nx_-1) - ii,ny_ -1 - jj-1)));
                                    Ez_->point((nx_-1) - ii,ny_ -1 - jj) = c_eze * Ez_->point((nx_-1) - ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                }
                                else
                                {
                                    //pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) + c_dzh * ((Hy_->point(ii,ny_ -1 - jj)-Hy_->point(ii-1,ny_ -1 - jj)) - (Hx_->point(ii,ny_ -1 - jj)-0));
                                    //Ez_->point(ii,ny_ -1 - jj) = c_eze * Ez_->point(ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj);
                                    pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) + c_dzh * ((Hy_->point((nx_-1) - ii,ny_ -1 - jj)-Hy_->point((nx_-1) - ii-1,ny_ -1 - jj)) - (0-Hx_->point((nx_-1) - ii,ny_ -1 - jj-1)));
                                    Ez_->point((nx_-1) - ii,ny_ -1 - jj) = c_eze * Ez_->point((nx_-1) - ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                }
                                
                            }
                            else
                            {
                                if(jj !=0 )
                                {
                                    //pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) + c_dzh * ((Hy_->point(ii,ny_ -1 - jj) - 0) - (Hx_->point(ii,ny_ -1 - jj)-Hx_->point(ii,ny_ -1 - jj-1))); // 0 is boundary condtion
                                    //Ez_->point(ii,ny_ -1 - jj) = c_eze * Ez_->point(ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj);
                                    pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) + c_dzh * ((0 -Hy_->point((nx_-1) - ii-1,ny_ -1 - jj)) - (Hx_->point((nx_-1) - ii,ny_ -1 - jj)-Hx_->point((nx_-1) - ii,ny_ -1 - jj-1))); // 0 is boundary condition
                                    Ez_->point((nx_-1) - ii,ny_ -1 - jj) = c_eze * Ez_->point((nx_-1) - ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                }
                                else
                                {
                                    //pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) + c_dzh * ((Hy_->point(ii,ny_ -1 - jj) - 0) - (Hx_->point(ii,ny_ -1 - jj)-0)); // 0 is boundary condtion
                                    //Ez_->point(ii,ny_ -1 - jj) = c_eze * Ez_->point(ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj);
                                    pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) + c_dzh * ((0 -Hy_->point((nx_-1) - ii-1,ny_ -1 - jj)) - (0 - Hx_->point((nx_-1) - ii,ny_ -1 - jj-1))); // 0 is boundary condition
                                    Ez_->point((nx_-1) - ii,ny_ -1 - jj) = c_eze * Ez_->point((nx_-1) - ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                }
                                
                            }
                            // Top Left
                            dzstore = pmlArr_[0].Dz_->point(ii,ny_ -1 - jj);
                            if(ii !=0)
                            {
                                if(jj!=0)
                                {
                                    pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) + c_dzh * ((Hy_->point(ii,ny_ -1 - jj)-Hy_->point(ii-1,ny_ -1 - jj)) - (Hx_->point(ii,ny_ -1 - jj)-Hx_->point(ii,ny_ -1 - jj-1)));
                                    Ez_->point(ii,ny_ -1 - jj) = c_eze * Ez_->point(ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj);
                                    //pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) + c_dzh * ((Hy_->point((nx_-1) - ii,ny_ -1 - jj)-Hy_->point((nx_-1) - ii-1,ny_ -1 - jj)) - (Hx_->point((nx_-1) - ii,ny_ -1 - jj)-Hx_->point((nx_-1) - ii,ny_ -1 - jj-1)));
                                    //Ez_->point((nx_-1) - ii,ny_ -1 - jj) = c_eze * Ez_->point((nx_-1) - ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                }
                                else
                                {
                                    pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) + c_dzh * ((Hy_->point(ii,ny_ -1 - jj)-Hy_->point(ii-1,ny_ -1 - jj)) - (0-Hx_->point(ii,ny_ -1 - jj-1)));
                                    Ez_->point(ii,ny_ -1 - jj) = c_eze * Ez_->point(ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj);
                                    //pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) + c_dzh * ((Hy_->point((nx_-1) - ii,ny_ -1 - jj)-Hy_->point((nx_-1) - ii-1,ny_ -1 - jj)) - (0-Hx_->point((nx_-1) - ii,ny_ -1 - jj-1)));
                                    //Ez_->point((nx_-1) - ii,ny_ -1 - jj) = c_eze * Ez_->point((nx_-1) - ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                }
                                
                            }
                            else
                            {
                                if(jj !=0 )
                                {
                                    pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) + c_dzh * ((Hy_->point(ii,ny_ -1 - jj) - 0) - (Hx_->point(ii,ny_ -1 - jj)-Hx_->point(ii,ny_ -1 - jj-1))); // 0 is boundary condtion
                                    Ez_->point(ii,ny_ -1 - jj) = c_eze * Ez_->point(ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj);
                                    //pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) + c_dzh * ((0 -Hy_->point((nx_-1) - ii-1,ny_ -1 - jj)) - (Hx_->point((nx_-1) - ii,ny_ -1 - jj)-Hx_->point((nx_-1) - ii,ny_ -1 - jj-1))); // 0 is boundary condition
                                    //Ez_->point((nx_-1) - ii,ny_ -1 - jj) = c_eze * Ez_->point((nx_-1) - ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                }
                                else
                                {
                                    pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) + c_dzh * ((Hy_->point(ii,ny_ -1 - jj) - 0) - (0-Hx_->point(ii,ny_ -1 - jj-1))); // 0 is boundary condtion
                                    Ez_->point(ii,ny_ -1 - jj) = c_eze * Ez_->point(ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                    
                                    //dzstore = pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj);
                                    //pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) + c_dzh * ((0 -Hy_->point((nx_-1) - ii-1,ny_ -1 - jj)) - (0 - Hx_->point((nx_-1) - ii,ny_ -1 - jj-1))); // 0 is boundary condition
                                    //Ez_->point((nx_-1) - ii,ny_ -1 - jj) = c_eze * Ez_->point((nx_-1) - ii,ny_ -1 - jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,ny_ -1 - jj) - c_ezd0 * dzstore;
                                }
                                
                            }
                        }
                        break;
                    case Y:
                        if(Ez_)
                        {
                            double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                            double sigx = pmlArr_[1].sigma(static_cast<double>(ii));
                            double sigy = pmlArr_[0].sigma(static_cast<double>(jj));                            
                            double sigz = 0.0;

                            double c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            double c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_)) ;
                            double c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            double c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            double c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            double dzstore = pmlArr_[1].Dz_->point(ii,jj); 

                            // Do the Bot left Corner
                            if(ii!=0 && jj!=0)
                            {   
                                pmlArr_[1].Dz_->point(ii,jj) = c_dzd * pmlArr_[1].Dz_->point(ii,jj) - c_dzh * ((Hy_->point(ii,jj) - Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj) - Hx_->point(ii,jj-1)));
                                Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[1].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                            }
                            else if(ii!=0)
                            {
                                pmlArr_[1].Dz_->point(ii,jj) = c_dzd * pmlArr_[1].Dz_->point(ii,jj) - c_dzh * ((Hy_->point(ii,jj) - Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj) - 0));
                                Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[1].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                            }
                            else if(jj !=0)
                            {
                                pmlArr_[1].Dz_->point(ii,jj) = c_dzd * pmlArr_[1].Dz_->point(ii,jj) - c_dzh * ((Hy_->point(ii,jj) - 0) - (Hx_->point(ii,jj) - Hx_->point(ii,jj-1)));
                                Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[1].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                            }
                            else
                            {
                                pmlArr_[1].Dz_->point(ii,jj) = c_dzd * pmlArr_[1].Dz_->point(ii,jj) - c_dzh * ((Hy_->point(ii,jj) - 0) - (Hx_->point(ii,jj) - 0));
                                
                                Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[1].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                            }
                            // Do the Bot Right Corner
                            if(ii != 0 && jj !=0)
                            {
                                dzstore = pmlArr_[1].Dz_->point((nx_-1)-ii,jj);
                                pmlArr_[1].Dz_->point((nx_-1)-ii,jj) = c_dzd * pmlArr_[1].Dz_->point((nx_-1)-ii,jj) - c_dzh * ((Hy_->point((nx_-1)-ii,jj) - Hy_->point((nx_-1)-ii-1,jj)) - (Hx_->point((nx_-1)-ii,jj) - Hx_->point((nx_-1)-ii,jj-1)));;
                                Ez_->point((nx_-1)-ii,jj) = c_eze * Ez_->point((nx_-1)-ii,jj) + c_ezd1 * pmlArr_[1].Dz_->point((nx_-1)-ii,jj) - c_ezd0 * dzstore;
                            }
                            else if(jj!=0)
                            {
                                dzstore = pmlArr_[1].Dz_->point((nx_-1)-ii,jj);
                                pmlArr_[1].Dz_->point((nx_-1)-ii,jj) = c_dzd * pmlArr_[1].Dz_->point((nx_-1)-ii,jj) - c_dzh * ((0 - Hy_->point((nx_-1)-ii-1,jj)) - (Hx_->point((nx_-1)-ii,jj) - Hx_->point((nx_-1)-ii,jj-1)));//boundray conditions
                                Ez_->point((nx_-1)-ii,jj) = c_eze * Ez_->point((nx_-1)-ii,jj) + c_ezd1 * pmlArr_[1].Dz_->point((nx_-1)-ii,jj) - c_ezd0 * dzstore;
                            }
                            else if(ii!=0)
                            {
                                dzstore = pmlArr_[1].Dz_->point((nx_-1)-ii,jj);
                                pmlArr_[1].Dz_->point((nx_-1)-ii,jj) = c_dzd * pmlArr_[1].Dz_->point((nx_-1)-ii,jj) - c_dzh * ((Hy_->point((nx_-1)-ii,jj) - Hy_->point((nx_-1)-ii-1,jj)) - (Hx_->point((nx_-1)-ii,jj) - 0));;
                                Ez_->point((nx_-1)-ii,jj) = c_eze * Ez_->point((nx_-1)-ii,jj) + c_ezd1 * pmlArr_[1].Dz_->point((nx_-1)-ii,jj) - c_ezd0 * dzstore;
                            }
                            else
                            {
                                dzstore = pmlArr_[1].Dz_->point((nx_-1)-ii,jj);
                                pmlArr_[1].Dz_->point((nx_-1)-ii,jj) = c_dzd * pmlArr_[1].Dz_->point((nx_-1)-ii,jj) - c_dzh * ((0 - Hy_->point((nx_-1)-ii-1,jj)) - (Hx_->point((nx_-1)-ii,jj) - 0));;
                                Ez_->point((nx_-1)-ii,jj) = c_eze * Ez_->point((nx_-1)-ii,jj) + c_ezd1 * pmlArr_[1].Dz_->point((nx_-1)-ii,jj) - c_ezd0 * dzstore;
                            }
                            // Do the Top Left Corner
                            if(jj != 0 && ii !=0)
                            {
                                dzstore = pmlArr_[1].Dz_end_->point(ii,jj);
                                pmlArr_[1].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[1].Dz_end_->point(ii,jj) - c_dzh * ((Hy_->point(ii,(ny_-1)-jj) - Hy_->point(ii-1,(ny_-1)-jj)) - (Hx_->point(ii,(ny_-1)-jj) - Hx_->point(ii,(ny_-1)-jj-1)));
                                Ez_->point(ii,(ny_-1)-jj) = c_eze * Ez_->point(ii,(ny_-1)-jj) + c_ezd1 * pmlArr_[1].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                            }
                            else if(jj!=0)
                            {
                                dzstore = pmlArr_[1].Dz_end_->point(ii,jj);
                                pmlArr_[1].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[1].Dz_end_->point(ii,jj) - c_dzh * ((Hy_->point(ii,(ny_-1)-jj) - 0) - (Hx_->point(ii,(ny_-1)-jj) - Hx_->point(ii,(ny_-1)-jj-1)));
                                Ez_->point(ii,(ny_-1)-jj) = c_eze * Ez_->point(ii,(ny_-1)-jj) + c_ezd1 * pmlArr_[1].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                            }
                            else if(ii!=0)
                            {
                                dzstore = pmlArr_[1].Dz_end_->point(ii,jj);
                                pmlArr_[1].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[1].Dz_end_->point(ii,jj) - c_dzh * ((Hy_->point(ii,(ny_-1)-jj) - Hy_->point(ii-1,(ny_-1)-jj)) - (0 - Hx_->point(ii,(ny_-1)-jj-1))); // boudrays 
                                Ez_->point(ii,(ny_-1)-jj) = c_eze * Ez_->point(ii,(ny_-1)-jj) + c_ezd1 * pmlArr_[1].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;   
                            }
                            else
                            {
                                dzstore = pmlArr_[1].Dz_end_->point(ii,jj);
                                pmlArr_[1].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[1].Dz_end_->point(ii,jj) - c_dzh * ((Hy_->point(ii,(ny_-1)-jj) - 0) - (0 - Hx_->point(ii,(ny_-1)-jj-1))); // boudrays 
                                Ez_->point(ii,(ny_-1)-jj) = c_eze * Ez_->point(ii,(ny_-1)-jj) + c_ezd1 * pmlArr_[1].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                            }
                            // Do the Top Right Corner
                            if(ii!=0 && jj!=0)
                            {
                                dzstore = pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj);
                                pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj) = c_dzd * pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj) - c_dzh * ((Hy_->point((nx_-1)-ii,(ny_-1)-jj) - Hy_->point((nx_-1)-ii-1,(ny_-1)-jj)) - (Hx_->point((nx_-1)-ii,(ny_-1)-jj) - Hx_->point((nx_-1)-ii,(ny_-1)-jj-1)));;
                                Ez_->point((nx_-1)-ii,(ny_-1)-jj) = c_eze * Ez_->point((nx_-1)-ii,(ny_-1)-jj) + c_ezd1 * pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj) - c_ezd0 * dzstore;
                            }
                            else if(jj!=0)
                            {
                                
                                dzstore = pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj);
                                pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj) = c_dzd * pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj) - c_dzh * ((0 - Hy_->point((nx_-1)-ii-1,(ny_-1)-jj)) - (Hx_->point((nx_-1)-ii,(ny_-1)-jj) - Hx_->point((nx_-1)-ii,(ny_-1)-jj-1)));;
                                Ez_->point((nx_-1)-ii,(ny_-1)-jj) = c_eze * Ez_->point((nx_-1)-ii,(ny_-1)-jj) + c_ezd1 * pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj) - c_ezd0 * dzstore;
                            }
                            else if(ii!=0)
                            {
                                dzstore = pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj);
                                pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj) = c_dzd * pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj) - c_dzh * ((Hy_->point((nx_-1)-ii,(ny_-1)-jj) - Hy_->point((nx_-1)-ii-1,(ny_-1)-jj)) - (0 - Hx_->point((nx_-1)-ii,(ny_-1)-jj-1)));;
                                Ez_->point((nx_-1)-ii,(ny_-1)-jj) = c_eze * Ez_->point((nx_-1)-ii,(ny_-1)-jj) + c_ezd1 * pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj) - c_ezd0 * dzstore;
                            }
                            else
                            {
                                dzstore = pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj);
                                pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj) = c_dzd * pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj) - c_dzh * ((0 - Hy_->point((nx_-1)-ii-1,(ny_-1)-jj)) - (0 - Hx_->point((nx_-1)-ii,(ny_-1)-jj-1)));;
                                Ez_->point((nx_-1)-ii,(ny_-1)-jj) = c_eze * Ez_->point((nx_-1)-ii,(ny_-1)-jj) + c_ezd1 * pmlArr_[1].Dz_end_->point((nx_-1)-ii,jj) - c_ezd0 * dzstore;
                            }
                        }
                        break;
                    case Z:
                        throw logic_error("Z is not yet defined");
                        break;
                    default:
                        throw logic_error("hit switch default");
                        break;
                }
            }
        }
    }
    for(int ii = 0; ii < dtcArr_.size(); ii ++)
        ouputField(dtcArr_[ii]);
    tcur_ += dt_;
}