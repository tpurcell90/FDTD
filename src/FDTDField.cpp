#include "FDTDField.hpp"

// #include <assert.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <fstream>
#include <algorithm>
#include "utilities.hpp"
// #include <random>
// #include <stdexcept>
// #include <string>
// #include <vector>
// #include <complex>



using namespace std;
typedef  double (UPML<complex<double>>::*PMLMemFn)(double x, double y);

/**
 * @brief Constructs a FDTD field from a programInputs object
 * @details Constructs the FDTD field manger using the information from the inputs parameter. Sets dx, and dy to 1/res, dt_ is set to S/res, The unused fields are set to null, nx and ny are calculated by rouding the product of the physical cell size with the
 *
 * @param IP Input object created from an input file
 */
FDTDField::FDTDField(programInputs &IP)
{
    //Define cell parameters
    tcur_       = 0;
    t_step_     = 0;
    res_        = IP.res_;
    dx_         = 1.0/res_;
    dy_         = 1.0/res_;
    dt_         = IP.courant_ * dx_;
    nx_         = floor(static_cast<double>(res_) * static_cast<double>(IP.x_size_)+ 0.5) + 1; //Better way here; + 1 to include the 0 point
    ny_         = floor(static_cast<double>(res_) * static_cast<double>(IP.y_size_)+ 0.5) + 1; //Better way here; + 1 to include the 0 point
    srcArr_     = IP.srcArr_;
    objArr_     = IP.objArr_;
    dtcArr_     = IP.dctArr_;
    pmlArr_     = IP.pmlArr_;
    xPML_       = IP.xPml_;
    yPML_       = IP.yPml_;
    precalcPML_ = IP.pmlCalc_;
    periodic_   = IP.periodic_;
    zaxEzList_  = {};
    y0EdgeInd_   = 0;
    ynEdgeInd_   = 0;
    x0EdgeInd_   = 0;
    xnEdgeInd_   = 0;

    if(IP.pol_.compare("Hz") == 0 || IP.pol_.compare("Ey") == 0 || IP.pol_.compare("Ex") == 0)
    {
        if(periodic_)
        {

        }
        else
        {
            Ex_ = make_shared<Grid2D<complex<double>>>(nx_-1,ny_,dx_,dy_);
            Ey_ = make_shared<Grid2D<complex<double>>>(nx_,ny_-1,dx_,dy_);
            Hz_ = make_shared<Grid2D<complex<double>>>(nx_-1,ny_-1,dx_,dy_);
            phys_Ex_ = make_shared<Grid2D<int>>(nx_-1,ny_,dx_,dy_);
            phys_Ey_ = make_shared<Grid2D<int>>(nx_,ny_-1,dx_,dy_);
            //phys_Hz_ = make_shared<Grid2D<int>>(nx_-1,ny_-1,dx_,dy_);
            //These are never used in the TE mode
            Hx_ = nullptr;
            Hy_ = nullptr;
            Ez_ = nullptr;
            //phys_Hx_ = nullptr;
            //phys_Hy_ = nullptr;
            phys_Ez_ = nullptr;
        }
    }
    else
    {
        // Hx_ = make_shared<Grid2D<complex<double>>>(nx_,ny_-1,dx_,dy_);
        // Hy_ = make_shared<Grid2D<complex<double>>>(nx_-1,ny_,dx_,dy_);
        Hx_ = make_shared<Grid2D<complex<double>>>(nx_,ny_,dx_,dy_);
        Hy_ = make_shared<Grid2D<complex<double>>>(nx_,ny_,dx_,dy_);
        Ez_ = make_shared<Grid2D<complex<double>>>(nx_,ny_,dx_,dy_);

        //phys_Hx_ = make_shared<Grid2D<int>>(nx_,ny_-1,dx_,dy_);
        //phys_Hy_ = make_shared<Grid2D<int>>(nx_-1,ny_,dx_,dy_);
        phys_Ez_ = make_shared<Grid2D<int>>(nx_,ny_,dx_,dy_);
        // These are never used in the TM mode
        Ex_ = nullptr;
        Ey_ = nullptr;
        Hz_ = nullptr;
        phys_Ex_ = nullptr;
        phys_Ey_ = nullptr;
        //phys_Hz_ = nullptr;
    }
}

/**
 * @brief Initializes the physical grid for materials look up
 * @details Initializes the physical grid for materials look up, sets the object to the number in the object array will overwrite to the last one if multiple objects exist at the same point
 *
 */
void FDTDField::initializeGrid()
{
    for(int kk = 0; kk < objArr_.size(); kk++)
    {
        vector<double> pt(2,0.0);
        if(Hz_)
        {
            if(objArr_[kk].s() == sphere)
            {
                for(int ii = 0; ii < nx_-1;ii ++)
                {
                    for(int jj = 0; jj < ny_-1; jj ++)
                    {

                        pt = {(ii+0.5)*dx_,jj*dy_};
                        if(objArr_[kk].isObj(pt))
                            phys_Ex_->point(ii,jj) = kk;
                        pt[1] += 0.5*dy_;
                        pt[0] -= 0.5*dx_;
                        if(objArr_[kk].isObj(pt))
                            phys_Ey_->point(ii,jj) = kk;
                    }
                }
            }
            else if(objArr_[kk].s() == block)
            {
                for(int ii = 0; ii < nx_-1;ii ++)
                {
                    for(int jj = 0; jj < ny_-1; jj ++)
                    {
                        pt = {(ii+0.5)*dx_,jj*dy_};
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Ex_->point(ii,jj) = kk;
                        pt[1] += 0.5*dy_;
                        pt[0] -= 0.5*dx_;
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Ey_->point(ii,jj) = kk;
                    }
                }
            }
        }
        else
        {
            if(objArr_[kk].s() == sphere)
            {
                for(int ii = 0; ii < nx_-1;ii ++)
                {
                    for(int jj = 0; jj < ny_-1; jj ++)
                    {
                        pt= {(ii-(nx_-1)/2.0)*dx_,(jj-static_cast<double>(ny_-1)/2.0)*dy_};
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Ez_->point(ii,jj) = kk;
                    }
                }
                for(int ii = 0; ii < nx_-1;ii ++)
                {
                    pt={(ii-(nx_-1)/2.0)*dx_,(ny_-1-(ny_-1)/2.0)*dy_};
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ez_->point(ii,ny_-1) = kk;
                }
                for(int jj = 0; jj < ny_-1; jj ++)
                {
                    pt = {(nx_-1-(nx_-1)/2.0)*dx_,(jj-(ny_-1)/2.0)*dy_};
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ez_->point(nx_-1,jj) = kk;
                }
                pt={(nx_-1-(nx_-1)/2.0)*dx_,(ny_-1-(ny_-1)/2.0)*dy_};
                if(objArr_[kk].isObj(pt)==true)
                    phys_Ez_->point(nx_-1,ny_-1) = kk;
            }
            else if(objArr_[kk].s() == block)
            {
                if(yPML_ != 0 && xPML_ != 0)
                {
                    for(int ii = xPML_; ii < nx_-xPML_;ii ++)
                    {
                        for(int jj = yPML_; jj < ny_-yPML_; jj ++)
                        {
                            pt[0] = (ii-(nx_-1)/2.0)*dx_;
                            pt[1] = (jj-static_cast<double>(ny_-1)/2.0)*dy_;
                            if(objArr_[kk].isObj(pt)==true)
                                phys_Ez_->point(ii,jj) = kk;
                        }
                    }
                    pt[0]=(nx_-1-(nx_-1)/2.0)*dx_;
                    pt[1]=(ny_-1-(ny_-1)/2.0)*dy_;
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ez_->point(nx_-1,ny_-1) = kk;
                }
                else if(yPML_ != 0)
                {
                    for(int ii = 0; ii < nx_-1;ii ++)
                    {
                        for(int jj = yPML_; jj < ny_-yPML_; jj ++)
                        {
                            pt[0] = (ii-(nx_-1)/2.0)*dx_;
                            pt[1] = (jj-static_cast<double>(ny_-1)/2.0)*dy_;
                            if(objArr_[kk].isObj(pt)==true)
                                phys_Ez_->point(ii,jj) = kk;
                        }
                    }
                    for(int jj = 0; jj < ny_-1; jj ++)
                    {
                        pt[0]=(nx_-1-(nx_-1)/2.0)*dx_;
                        pt[1]=(jj-(ny_-1)/2.0)*dy_;
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Ez_->point(nx_-1,jj) = kk;
                    }
                    pt[0]=(nx_-1-(nx_-1)/2.0)*dx_;
                    pt[1]=(ny_-1-(ny_-1)/2.0)*dy_;
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ez_->point(nx_-1,ny_-1) = kk;
                }
                else if(xPML_ != 0)
                {
                    for(int ii = xPML_; ii < nx_-xPML_;ii ++)
                    {
                        for(int jj = 0; jj < ny_-1; jj ++)
                        {
                            pt[0] = (ii-(nx_-1)/2.0)*dx_;
                            pt[1] = (jj-static_cast<double>(ny_-1)/2.0)*dy_;
                            if(objArr_[kk].isObj(pt)==true)
                                phys_Ez_->point(ii,jj) = kk;
                        }
                    }
                    for(int ii = 0; ii < nx_-1;ii ++)
                    {
                        pt[0]=(ii-(nx_-1)/2.0)*dx_;
                        pt[1]=(ny_-1-(ny_-1)/2.0)*dy_;
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Ez_->point(ii,ny_-1) = kk;
                    }
                    pt[0]=(nx_-1-(nx_-1)/2.0)*dx_;
                    pt[1]=(ny_-1-(ny_-1)/2.0)*dy_;
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ez_->point(nx_-1,ny_-1) = kk;
                }
                else
                {
                    for(int ii = 0; ii < nx_-1;ii ++)
                    {
                        for(int jj = 0; jj < ny_-1; jj ++)
                        {
                            pt[0] = (ii-(nx_-1)/2.0)*dx_;
                            pt[1] = (jj-static_cast<double>(ny_-1)/2.0)*dy_;
                            if(objArr_[kk].isObj(pt)==true)
                                phys_Ez_->point(ii,jj) = kk;
                        }
                    }
                    for(int ii = 0; ii < nx_-1;ii ++)
                    {
                        pt[0]=(ii-(nx_-1)/2.0)*dx_;
                        pt[1]=(ny_-1-(ny_-1)/2.0)*dy_;
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Ez_->point(ii,ny_-1) = kk;
                    }
                    for(int jj = 0; jj < ny_-1; jj ++)
                    {
                        pt[0]=(nx_-1-(nx_-1)/2.0)*dx_;
                        pt[1]=(jj-(ny_-1)/2.0)*dy_;
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Ez_->point(nx_-1,jj) = kk;
                    }
                    pt[0]=(nx_-1-(nx_-1)/2.0)*dx_;
                    pt[1]=(ny_-1-(ny_-1)/2.0)*dy_;
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ez_->point(nx_-1,ny_-1) = kk;
                }
            }
        }
    }
    if(xPML_ != 0 && yPML_ != 0)
    {
        for(int jj = yPML_; jj < ny_ - yPML_; jj++)
        {
            int ii = xPML_;
            while(ii < nx_-xPML_)
            {
                int iistore = ii;
                while(ii < nx_-xPML_-1 && phys_Ez_ -> point(ii,jj) == phys_Ez_ -> point(ii+1,jj) )
                    ii ++;
                array<int,4> tempArr = { iistore,jj,ii-iistore+1,phys_Ez_->point(iistore,jj)};
                zaxEzList_.push_back(tempArr);
                ii++;
            }
        }
    }
    else if(xPML_ != 0)
    {
        for(int jj = 1; jj < ny_ - 1; jj++)
        {
            int ii = xPML_;
            while(ii < nx_-xPML_)
            {
                int iistore = ii;
                while(ii < nx_-xPML_-1 && phys_Ez_ -> point(ii,jj) == phys_Ez_ -> point(ii+1,jj) )
                    ii ++;
                array<int,4> tempArr = { iistore,jj,ii-iistore+1,phys_Ez_->point(iistore,jj)};
                zaxEzList_.push_back(tempArr);
                ii++;
            }
        }
        y0EdgeInd_= zaxEzList_.size();
        int ii = xPML_;
        while(ii < nx_-xPML_)
        {
            int iistore = ii;
            while(ii < nx_-xPML_-1 && phys_Ez_ -> point(ii,0) == phys_Ez_ -> point(ii+1,0) )
                ii ++;
            array<int,4> tempArr = { iistore,0,ii-iistore+1,phys_Ez_->point(iistore,0)};
            zaxEzList_.push_back(tempArr);
            ii++;
        }
        ynEdgeInd_= zaxEzList_.size();
        ii = xPML_;
        while(ii < nx_-xPML_)
        {
            int iistore = ii;
            while(ii < nx_-xPML_-1 && phys_Ez_ -> point(ii,ny_-1) == phys_Ez_ -> point(ii+1,ny_-1) )
                ii ++;
            int n = ny_-1;
            array<int,4> tempArr = { iistore,n,ii-iistore+1,phys_Ez_->point(iistore,ny_-1)};
            zaxEzList_.push_back(tempArr);
            ii++;
        }
    }
    else if(yPML_ != 0)
    {
        for(int jj = yPML_; jj < ny_ - yPML_; jj++)
        {
            int ii = 1;
            while(ii < nx_-1)
            {
                int iistore = ii;
                while(ii < nx_-1-1 && phys_Ez_ -> point(ii,jj) == phys_Ez_ -> point(ii+1,jj) )
                    ii ++;
                array<int,4> tempArr = { iistore,jj,ii-iistore+1,phys_Ez_->point(iistore,jj)};
                zaxEzList_.push_back(tempArr);
                ii++;
            }
        }
        x0EdgeInd_= zaxEzList_.size();
        int ii = yPML_;
        while(ii < ny_-yPML_)
        {
            int iistore = ii;
            while(ii < nx_-xPML_-1 && phys_Ez_ -> point(0,ii) == phys_Ez_ -> point(0,ii+1) )
                ii ++;
            array<int,4> tempArr = {0,iistore,ii-iistore+1,phys_Ez_->point(iistore,0)};
            zaxEzList_.push_back(tempArr);
            ii++;
        }
        xnEdgeInd_= zaxEzList_.size();
        ii = yPML_;
        while(ii < nx_-xPML_)
        {
            int iistore = ii;
            while(ii < nx_-xPML_-1 && phys_Ez_ -> point(0,ii) == phys_Ez_ -> point(0,ii+1) )
                ii ++;
            int n = nx_-1;
            array<int,4> tempArr = {n,iistore,ii-iistore+1,phys_Ez_->point(iistore,0)};
            zaxEzList_.push_back(tempArr);
            ii++;
        }
    }
    else
    {
        for(int jj = 1; jj < ny_ - 1; jj++)
        {
            int ii = 1;
            while(ii < nx_-1)
            {
                int iistore = ii;
                while(ii < nx_-1-1 && phys_Ez_ -> point(ii,jj) == phys_Ez_ -> point(ii+1,jj) )
                    ii ++;
                array<int,4> tempArr = { iistore,jj,ii-iistore+1,phys_Ez_->point(iistore,jj)};
                zaxEzList_.push_back(tempArr);
                ii++;
            }
        }
        y0EdgeInd_= zaxEzList_.size();
        int ii = 1;
        while(ii < nx_-1)
        {
            int iistore = ii;
            while(ii < nx_-1-1 && phys_Ez_ -> point(ii,0) == phys_Ez_ -> point(ii+1,0) )
                ii ++;
            array<int,4> tempArr = { iistore,0,ii-iistore+1,phys_Ez_->point(iistore,0)};
            zaxEzList_.push_back(tempArr);
            ii++;
        }
        ynEdgeInd_= zaxEzList_.size();
        ii = 1;
        while(ii < nx_-1)
        {
            int iistore = ii;
            while(ii < nx_-1-1 && phys_Ez_ -> point(ii,ny_-1) == phys_Ez_ -> point(ii+1,ny_-1) )
                ii ++;
            int n = ny_-1;
            array<int,4> tempArr = { iistore,n,ii-iistore+1,phys_Ez_->point(iistore,ny_-1)};
            zaxEzList_.push_back(tempArr);
            ii++;
        }
        x0EdgeInd_= zaxEzList_.size();
        ii = 1;
        while(ii < ny_-1)
        {
            int iistore = ii;
            while(ii < ny_-1-1 && phys_Ez_ -> point(0,ii) == phys_Ez_ -> point(0,ii+1) )
                ii ++;
            array<int,4> tempArr = { 0,iistore,ii-iistore+1,phys_Ez_->point(0,iistore)};
            zaxEzList_.push_back(tempArr);
            ii++;
        }
        xnEdgeInd_= zaxEzList_.size();
        ii = 1;
        while(ii < ny_-1)
        {
            int iistore = ii;
            while(ii < ny_-1-1 && phys_Ez_ -> point(nx_-1,ii) == phys_Ez_ -> point(nx_-1,ii+1) )
                ii ++;
            int n = nx_-1;
            array<int,4> tempArr = { n,iistore,ii-iistore+1,phys_Ez_->point(nx_-1,iistore)};
            zaxEzList_.push_back(tempArr);
            ii++;
        }
    }
    if(pmlArr_.size() > 1)
    {
        for(int kk = 0; kk < pmlArr_.size(); kk++)
        {
            UPML<complex<double>>* opp = &pmlArr_[abs(kk-1)];
            pmlArr_[kk].initializeUPML(objArr_, nx_,ny_,dx_,dy_,dt_, yPML_, xPML_, opp);
        }
    }
    else
        pmlArr_[0].initializeUPML(objArr_, nx_,ny_,dx_,dy_,dt_, yPML_, xPML_, nullptr);
}
/**
 * @brief Outputs the relevant field information to an output file specified by the detector
 * @details Outputs the relevant field information to an output file specified by the detector
 *
 * @param d The Detector used for the output
 */
void FDTDField::ouputField(Detector<complex<double>> d) //iostream as input parameter?
{
    ofstream outFile;
    outFile.open(d.outfile(), ios_base::app);
    double eps = 1.00;
    switch ( d.pol() )
    {
        case EZ:
            outFile << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Ez_,eps).real()<< "\t" << setw(10) << srcArr_[0].prof().pulse(t_step_).real() << endl;
            cout << setw(9) << tcur_ << "\t" << d.loc()[0] << "\t" << d.loc()[1] << "\t" << setw(10) << d.output(Ez_,eps).real()<< "\t" << setw(10) << srcArr_[0].prof().pulse(t_step_).real() << endl;
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

/**
 * @brief Updates the H fields to the next time step throughout all space
 * @details Updates the H fields uding the FDTD update equations
 * For Normal Space TM modes
 * \f$H_x^{q+\frac{3}{2}}\left[i,j+\frac{1}{2}\right] =  \frac{1 - \frac{\sigma^*\left[i,j+\frac{1}{2}\right] d_t}{2*mu\left[i,j+\frac{1}{2}\right]d_t}}{2\mu\left[i,j+\frac{1}{2}\right]}} H_x^{q+\frac{3}{2}}\left[i,j+\frac{1}{2}\right] + \frac{\frac{dt}{mu\left[i,j+\frac{1}{2}\right]}}{\sigma^*\left[i,j+\frac{1}{2}\right] d_t}{2*mu\left[i,j+\frac{1}{2}\right]d_t}{2\mu\left[i,j+\frac{1}{2}\right]}} \left\{ \frac{1}{dy} \left( E_z^{q+1}\left[i,j+1\right] - E_z^{q+1}\left[i,j\right]\right)  \right\}\f$
 * \f$H_y^{q+\frac{3}{2}}\left[i+\frac{1}{2},j\right] =  \frac{1 - \frac{\sigma^*\left[i+\frac{1}{2},j\right] d_t}{2*mu\left[i+\frac{1}{2},j\right]d_t}}{2\mu\left[i+\frac{1}{2},j\right]}} H_y^{q+\frac{3}{2}}\left[i+\frac{1}{2},j\right] + \frac{\frac{dt}{mu\left[i+\frac{1}{2},j\right]}}{\sigma^*\left[i+\frac{1}{2},j\right] d_t}{2*mu\left[i+\frac{1}{2},j\right]d_t}{2\mu\left[i+\frac{1}{2},j\right]}} \left\{ \frac{1}{dy} \left( E_z^{q+1}\left[i+1,j\right] - E_z^{q+1}\left[i,j\right]\right)  \right\}\f$
 *
 * For UPML space TM mode
 * \f$B_x^{q+\frac{3}{2}}\left[i,j+\frac{1}{2}\right] = \frac{2 \epsilon \kappa_y - \sigma_y dt}{2 \epsilon  \kappa_y + \sigma_y dt} B_x^{q+\frac{3}{2}}\left[i,j+\frac{1}{2}\right] - \frac{2  \epsilon  dt}{2  \epsilon  \kappa_y + \sigma_y  dt} \left\{ \frac{1}{dy} \left( E_z^{q+1}\left[i,j+1\right] - E_z^{q+1}\left[i,j\right]\right)  \right\} \f$
 *
 *\f$H_x^{q+\frac{3}{2}}\left[i,j+\frac{1}{2}\right] =  \frac{2 \epsilon \kappa_z - \sigma_z dt}{2 \epsilon  \kappa_z + \sigma_z dt} H_x^{q+\frac{1}{2}}\left[i,j+\frac{1}{2}\right] + \frac{1}{\left(2\epsilon \kappa_z +\sigma_z dt\right)\epsilon} \left\{ \left(2\epsilon\kappa_x + \sigma_x dt\right) B_x^{q+\frac{3}{2}}\left[i,j+\frac{1}{2}\right] -\left(2\epsilon\kappa_x - \sigma_x dt\right) B_x^{q+\frac{1}{2}}\left[i,j+\frac{1}{2}\right]\right\} \f$
 *
 *\f$B_y^{q+\frac{3}{2}}\left[i+\frac{1}{2},j\right] = \frac{2 \epsilon \kappa_z - \sigma_z dt}{2 \epsilon  \kappa_z + \sigma_z dt} B_y^{q+\frac{3}{2}}\left[i+\frac{1}{2},j\right] +  \frac{2  \epsilon  dt}{2  \epsilon  \kappa_z + \sigma_z  dt} \left\{ \frac{1}{dz} \frac{1}{dx} \left(E_z^{q+1}\left[i+1,j\right] - E_z^{q+1}\left[i,j\right)\right]     \right\} \f$
 *
 *\f$H_y^{q+\frac{3}{2}}\left[i+\frac{1}{2},j,\right]  = \frac{2 \epsilon \kappa_x - \sigma_x dt}{2 \epsilon  \kappa_x + \sigma_x dt} H_y^{q+\frac{1}{2}}\left[i+\frac{1}{2},j\right] + \frac{1}{\left(2\epsilon \kappa_x +\sigma_x dt\right)\epsilon} \left\{ \left(2\epsilon\kappa_y + \sigma_y dt\right) B_y^{q+\frac{3}{2}}\left[i+\frac{1}{2},j\right]  -\left(2\epsilon\kappa_y - \sigma_y dt\right) B_y^{q+\frac{1}{2}}\left[i+\frac{1}{2},j\right] \right\}  \f$
 *
 */
void FDTDField::updateH()
{

    if(Ez_)
    {
        complex<double> c_hxh(1.0,0.0);
        double c_hxe = 1.0 * dt_/dx_;
        complex<double> c_hyh(1.0,0.0);
        double c_hye = 1.0 * dt_/dy_;
        if(xPML_ != 0 && yPML_ !=0)
        {
            /*for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
            {
                for(int ii = xPML_; ii < nx_-xPML_; ii++)
                {
                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) - c_hxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                }
            }*/
            //vector<complex<double>> hxstore(nx_-(2*xPML_),0.0);
            for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
            {
                zscal_(nx_-2*xPML_, c_hxh, &Hx_ ->point(xPML_,jj),1);
                zaxpy_(nx_-2*xPML_, -1.0*c_hxe, &Ez_->point(xPML_,jj+1), 1, &Hx_ ->point(xPML_,jj),1);
                zaxpy_(nx_-2*xPML_, c_hxe, &Ez_->point(xPML_,jj), 1, &Hx_ ->point(xPML_,jj),1);

                zscal_(nx_-2*xPML_, c_hyh, &Hy_ ->point(xPML_,jj),1);
                zaxpy_(nx_-2*xPML_, c_hye, &Ez_->point(xPML_+1,jj), 1, &Hy_ ->point(xPML_,jj),1);
                zaxpy_(nx_-2*xPML_, -1.0*c_hye, &Ez_->point(xPML_,jj), 1, &Hy_ ->point(xPML_,jj),1);

            }
        }

        else if(xPML_ != 0)
        {
            for(int jj = 0; jj < ny_-1; jj++)
            {
                zscal_(nx_-2*xPML_, c_hxh, &Hx_ ->point(xPML_,jj),1);
                zaxpy_(nx_-2*xPML_, -1.0*c_hxe, &Ez_->point(xPML_,jj+1), 1, &Hx_ ->point(xPML_,jj),1);
                zaxpy_(nx_-2*xPML_, c_hxe, &Ez_->point(xPML_,jj), 1, &Hx_ ->point(xPML_,jj),1);

                zscal_(nx_-2*xPML_, c_hyh, &Hy_ ->point(xPML_,jj),1);
                zaxpy_(nx_-2*xPML_, c_hye, &Ez_->point(xPML_+1,jj), 1, &Hy_ ->point(xPML_,jj),1);
                zaxpy_(nx_-2*xPML_, -1.0*c_hye, &Ez_->point(xPML_,jj), 1, &Hy_ ->point(xPML_,jj),1);
            }
            zscal_(nx_-2*xPML_, c_hyh, &Hy_ ->point(xPML_,ny_-1),1);
            zaxpy_(nx_-2*xPML_, c_hye, &Ez_->point(xPML_+1,ny_-1), 1, &Hy_ ->point(xPML_,ny_-1),1);
            zaxpy_(nx_-2*xPML_, -1.0*c_hye, &Ez_->point(xPML_,ny_-1), 1, &Hy_ ->point(xPML_,ny_-1),1);

            /*for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                for(int jj = 0; jj < ny_ - 1; jj ++)
                {
                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) - c_hxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                }
            }
            for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                Hy_->point(ii,ny_-1) = c_hyh * Hy_->point(ii,ny_-1) + c_hye * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
            }*/
        }
        else if(yPML_ != 0)
        {
            //vector<complex<double>> hxstore(nx_,0.0);
            for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
            {
                zscal_(nx_, c_hxh, &Hx_ ->point(0,jj),1);
                zaxpy_(nx_, -1.0*c_hxe, &Ez_->point(0,jj+1), 1, &Hx_ ->point(0,jj),1);
                zaxpy_(nx_, c_hxe, &Ez_->point(0,jj), 1, &Hx_ ->point(0,jj),1);

                zscal_(nx_-1, c_hyh, &Hy_ ->point(0,jj),1);
                zaxpy_(nx_-1, c_hye, &Ez_->point(1,jj), 1, &Hy_ ->point(0,jj),1);
                zaxpy_(nx_-1, -1.0*c_hye, &Ez_->point(0,jj), 1, &Hy_ ->point(0,jj),1);
            }
            /*for(int ii = 0; ii < nx_ - 1; ii ++)
            {
                for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
                {
                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) - c_hxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                }
            }
            for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
            {
                Hx_->point(nx_-1,jj) = c_hxh * Hx_->point(nx_-1,jj) - c_hxe * (Ez_->point(nx_-1,jj+1)-Ez_->point(nx_-1,jj));
            }*/
        }
        else
        {
            for(int jj = 0; jj < ny_-1; jj++)
            {
                zscal_(nx_, c_hxh, &Hx_ ->point(0,jj),1);
                zaxpy_(nx_, -1.0*c_hxe, &Ez_->point(0,jj+1), 1, &Hx_ ->point(0,jj),1);
                zaxpy_(nx_, c_hxe, &Ez_->point(0,jj), 1, &Hx_ ->point(0,jj),1);

                zscal_(nx_-1, c_hyh, &Hy_ ->point(0,jj),1);
                zaxpy_(nx_-1, c_hye, &Ez_->point(1,jj), 1, &Hy_ ->point(0,jj),1);
                zaxpy_(nx_-1, -1.0*c_hye, &Ez_->point(0,jj), 1, &Hy_ ->point(0,jj),1);
            }
            zscal_(nx_-1, c_hyh, &Hy_ ->point(0,ny_-1),1);
            zaxpy_(nx_-1, c_hye, &Ez_->point(1,ny_-1), 1, &Hy_ ->point(0,ny_-1),1);
            zaxpy_(nx_-1, -1.0*c_hye, &Ez_->point(0,ny_-1), 1, &Hy_ ->point(0,ny_-1),1);
            /*for(int ii = 0; ii < nx_ - 1; ii ++)
            {
                for(int jj = 0; jj < ny_ - 1; jj ++)
                {

                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) - c_hxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                }
            }
            for(int jj = 0; jj < nx_ - 1; jj ++)
            {
                Hx_->point(nx_-1,jj) = c_hxh * Hx_->point(nx_-1,jj) - c_hxe * (Ez_->point(nx_-1,jj+1)-Ez_->point(nx_-1,jj));
            }
            for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                Hy_->point(ii,ny_-1) = c_hyh * Hy_->point(ii,ny_-1) + c_hye * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
            }*/
        }
    }
    /*if(precalcPML_ ==false)
    {
        if(Ez_)
        {
            for(int kk = 0; kk < pmlArr_.size(); kk++)
            {
                double eps=0.0;
                double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                double sigz = 0.0;
                switch(pmlArr_[kk].d())
                {
                    case X:
                    {
                        double sigxx=0.0;
                        double sigxy=0.0;
                        double sigyx = 0.0;
                        double sigyy = 0.0;

                        double c_bxb=0.0; double c_bxe=0.0; double c_hxh=0.0; double c_hxb0=0.0; double c_hxb1=0.0;
                        double c_byb=0.0; double c_bye=0.0; double c_hyh=0.0;double c_hyb0=0.0; double c_hyb1=0.0;

                        complex<double> bxstore(0.0,0.0); complex<double> bystore(0.0,0.0);
                        if (yPML_== 0)
                        {
                            for (int jj = 0; jj < ny_-1; jj++)
                            {
                                eps    = objArr_[pmlArr_[kk].phys_Hx_->point(0,jj)].dielectric(1.0);
                                sigxx  = pmlArr_[kk].sigma(0.0,eps);
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                eps    = objArr_[pmlArr_[kk].phys_Hy_->point(0,jj)].dielectric(1.0);
                                sigxy = pmlArr_[kk].sigma(0.5,eps);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                //Update both on the left side and Hx on the right
                                bxstore = pmlArr_[kk].Bx_->point(0,jj);
                                bystore = pmlArr_[kk].By_->point(0,jj);

                                pmlArr_[kk].Bx_->point(0,jj) = c_bxb * pmlArr_[kk].Bx_->point(0,jj) - c_bxe * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                                // pmlArr_[kk].By_->point(0,jj) = c_byb * pmlArr_[kk].By_->point(0,jj) + c_bye * (Ez_->point(0+1,jj)-Ez_->point(0,jj));

                                Hx_->point(0,jj) = c_hxh * Hx_->point(0,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(0,jj) - c_hxb0 * bxstore;
                                // Hy_->point(0,jj) = c_hyh * Hy_->point(0,jj) + c_hyb1 * pmlArr_[kk].By_->point(0,jj) - c_hyb0 * bystore;

                                eps    = objArr_[pmlArr_[kk].phys_Hx_end_->point(0,jj)].dielectric(1.0);
                                sigxx  = pmlArr_[kk].sigma(0.0,eps);
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                bxstore = pmlArr_[kk].Bx_end_->point(0,jj);
                                pmlArr_[kk].Bx_end_->point(0,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(0,jj) - c_bxe * (Ez_->point(nx_-1, jj+1)-Ez_->point(nx_-1,jj));
                                Hx_->point(nx_-1,jj) = c_hxh * Hx_->point(nx_-1,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(0,jj) - c_hxb0 * bxstore;
                            }
                            //Top left corner
                            eps    = objArr_[pmlArr_[kk].phys_Hy_->point(0,ny_-1)].dielectric(1.0);
                            sigxy = pmlArr_[kk].sigma(0.5,eps);
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                            bystore = pmlArr_[kk].By_->point(0,ny_-1);
                            // pmlArr_[kk].By_->point(0,ny_-1) = c_byb * pmlArr_[kk].By_->point(0,ny_-1) + c_bye * (Ez_->point(0+1,ny_-1)-Ez_->point(0,ny_-1));
                            // Hy_->point(0,ny_-1) = c_hyh * Hy_->point(0,ny_-1) + c_hyb1 * pmlArr_[kk].By_->point(0,ny_-1) - c_hyb0 * bystore;
                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                eps    = objArr_[pmlArr_[kk].phys_Hy_->point(ii,ny_-1)].dielectric(1.0);
                                sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5,eps);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                //top row
                                bystore = pmlArr_[kk].By_->point(ii,ny_-1);
                                // pmlArr_[kk].By_->point(ii,ny_-1) = c_byb * pmlArr_[kk].By_->point(ii,ny_-1) + c_bye * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
                                // Hy_->point(ii,ny_-1) = c_hyh * Hy_->point(ii,ny_-1) + c_hyb1 * pmlArr_[kk].By_->point(ii,ny_-1) - c_hyb0 * bystore;

                                for(int jj = 0; jj < ny_-1; jj ++)
                                {
                                    eps    = objArr_[pmlArr_[kk].phys_Hx_->point(ii,jj)].dielectric(1.0);
                                    sigxx = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[kk].phys_Hy_->point(ii,jj)].dielectric(1.0);
                                    sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5,eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    //Update everything

                                    bxstore = pmlArr_[kk].Bx_->point(ii,jj);
                                    bystore = pmlArr_[kk].By_->point(ii,jj);

                                    pmlArr_[kk].Bx_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                    // pmlArr_[kk].By_->point(ii,jj) = c_byb * pmlArr_[kk].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));

                                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                                    // Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[kk].By_->point(ii,jj) - c_hyb0 * bystore;

                                    eps    = objArr_[pmlArr_[kk].phys_Hx_end_->point(ii,jj)].dielectric(1.0);
                                    sigxx = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[kk].phys_Hy_end_->point(ii,jj)].dielectric(1.0);
                                    sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) - 0.5,eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    bxstore = pmlArr_[kk].Bx_end_->point(ii,jj);
                                    bystore = pmlArr_[kk].By_end_->point(ii,jj);

                                    pmlArr_[kk].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point(nx_-1-ii, jj+1)-Ez_->point(nx_-1-ii,jj));
                                    pmlArr_[kk].By_end_->point(ii,jj) = c_byb * pmlArr_[kk].By_end_->point(ii,jj) + c_bye * (Ez_->point(nx_-1-ii+1, jj)-Ez_->point(nx_-1-ii,jj));

                                    Hx_->point(nx_-1-ii,jj) = c_hxh * Hx_->point((nx_-1) - ii,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                                    Hy_->point(nx_-1-ii,jj) = c_hyh * Hy_->point((nx_-1) - ii,jj) + c_hyb1 * pmlArr_[kk].By_end_->point(ii,jj) - c_hyb0 * bystore;
                                }
                                // for the top row only upadate Hy
                                eps    = objArr_[pmlArr_[kk].phys_Hy_end_->point(ii,ny_-1)].dielectric(1.0);
                                sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) - 0.5,eps);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                bystore = pmlArr_[kk].By_end_->point(ii,ny_-1);
                                pmlArr_[kk].By_end_->point(ii,ny_-1) = c_byb * pmlArr_[kk].By_end_->point(ii,ny_-1) + c_bye * (Ez_->point(nx_-1-ii+1, ny_-1)-Ez_->point(nx_-1-ii,ny_-1));
                                Hy_->point(nx_-1-ii,ny_-1) = c_hyh * Hy_->point(nx_-1-ii,ny_-1) + c_hyb1 * pmlArr_[kk].By_end_->point(ii,ny_-1) - c_hyb0 * bystore;
                            }
                        }
                        else
                        {
                            for (int jj = yPML_; jj < ny_- yPML_; jj++)
                            {
                                eps    = objArr_[pmlArr_[kk].phys_Hx_->point(0,jj)].dielectric(1.0);
                                sigxx  = pmlArr_[kk].sigma(0.0,eps);
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                eps    = objArr_[pmlArr_[kk].phys_Hy_->point(0,jj)].dielectric(1.0);
                                sigxy  = pmlArr_[kk].sigma(0.5,eps);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                //Update both on the left side and Hx on the right
                                bxstore = pmlArr_[kk].Bx_->point(0,jj);
                                bystore = pmlArr_[kk].By_->point(0,jj);

                                pmlArr_[kk].Bx_->point(0,jj) = c_bxb * pmlArr_[kk].Bx_->point(0,jj) - c_bxe * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                                pmlArr_[kk].By_->point(0,jj) = c_byb * pmlArr_[kk].By_->point(0,jj) + c_bye * (Ez_->point(0+1,jj)-Ez_->point(0,jj));

                                Hx_->point(0,jj) = c_hxh * Hx_->point(0,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(0,jj) - c_hxb0 * bxstore;
                                Hy_->point(0,jj) = c_hyh * Hy_->point(0,jj) + c_hyb1 * pmlArr_[kk].By_->point(0,jj) - c_hyb0 * bystore;

                                eps    = objArr_[pmlArr_[kk].phys_Hx_end_->point(0,jj)].dielectric(1.0);
                                sigxx = pmlArr_[kk].sigma(0.0,eps);
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                bxstore = pmlArr_[kk].Bx_end_->point(0,jj);
                                pmlArr_[kk].Bx_end_->point(0,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(0,jj) - c_bxe * (Ez_->point(nx_-1, jj+1)-Ez_->point(nx_-1,jj));
                                Hx_->point(nx_-1,jj) = c_hxh * Hx_->point(nx_-1,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(0,jj) - c_hxb0 * bxstore;
                            }
                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                for(int jj = yPML_; jj < ny_- yPML_; jj ++)
                                {
                                    //Update everything
                                    eps    = objArr_[pmlArr_[kk].phys_Hx_->point(ii,jj)].dielectric(1.0);
                                    sigxx = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[kk].phys_Hy_->point(ii,jj)].dielectric(1.0);
                                    sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5,eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    bxstore = pmlArr_[kk].Bx_->point(ii,jj);
                                    bystore = pmlArr_[kk].By_->point(ii,jj);

                                    pmlArr_[kk].Bx_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                    pmlArr_[kk].By_->point(ii,jj) = c_byb * pmlArr_[kk].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));

                                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[kk].By_->point(ii,jj) - c_hyb0 * bystore;

                                    eps    = objArr_[pmlArr_[kk].phys_Hx_end_->point(ii,jj)].dielectric(1.0);
                                    sigxx  = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[kk].phys_Hy_end_->point(ii,jj)].dielectric(1.0);
                                    sigxy  = pmlArr_[kk].sigma(static_cast<double>(ii) - 0.5,eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    bxstore = pmlArr_[kk].Bx_end_->point(ii,jj);
                                    bystore = pmlArr_[kk].By_end_->point(ii,jj);

                                    pmlArr_[kk].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point(nx_-1-ii, jj+1)-Ez_->point(nx_-1-ii,jj));
                                    pmlArr_[kk].By_end_->point(ii,jj) = c_byb * pmlArr_[kk].By_end_->point(ii,jj) + c_bye * (Ez_->point(nx_-1-ii+1, jj)-Ez_->point(nx_-1-ii,jj));

                                    Hx_->point(nx_-1-ii,jj) = c_hxh * Hx_->point(nx_-1-ii,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                                    Hy_->point(nx_-1-ii,jj) = c_hyh * Hy_->point(nx_-1-ii,jj) + c_hyb1 * pmlArr_[kk].By_end_->point(ii,jj) - c_hyb0 * bystore;
                                }
                            }
                            if(kk==0)
                            {
                                kapz = 1.0; sigz = 0.0; eps = 1.0;
                                kapx = 1.0; kapy = 1.0;

                                eps    = objArr_[pmlArr_[0].phys_Hx_->point(0,0)].dielectric(1.0);
                                sigxx = pmlArr_[0].sigma(0.0,eps);
                                sigyx = pmlArr_[1].sigma(0.5,eps);
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                eps    = objArr_[pmlArr_[0].phys_Hy_->point(0,0)].dielectric(1.0);
                                sigxy = pmlArr_[0].sigma(0.5,eps);
                                sigyy = pmlArr_[1].sigma(0.0,eps);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                // Bot Left
                                bxstore = pmlArr_[0].Bx_->point(0,0);
                                bystore = pmlArr_[0].By_->point(0,0);
                                pmlArr_[0].Bx_->point(0,0) = c_bxb * pmlArr_[0].Bx_->point(0,0) - c_bxe * (Ez_->point(0,0+1)-Ez_->point(0,0));
                                pmlArr_[0].By_->point(0,0) = c_byb * pmlArr_[0].By_->point(0,0) + c_bye * (Ez_->point(0+1,0)-Ez_->point(0,0));
                                Hx_->point(0,0) = c_hxh * Hx_->point(0,0) + c_hxb1 * pmlArr_[0].Bx_->point(0,0) - c_hxb0 * bxstore;
                                Hy_->point(0,0) = c_hyh * Hy_->point(0,0) + c_hyb1 * pmlArr_[0].By_->point(0,0) - c_hyb0 * bystore;
                                //Bot Right
                                eps    = objArr_[pmlArr_[0].phys_Hx_end_->point(0,0)].dielectric(1.0);
                                sigxx = pmlArr_[0].sigma(0.0,eps);
                                sigyx = pmlArr_[1].sigma(0.5,eps);
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                bxstore = pmlArr_[0].Bx_end_->point(0,0);
                                pmlArr_[0].Bx_end_->point(0,0) = c_bxb * pmlArr_[0].Bx_end_->point(0,0) - c_bxe * (Ez_->point(nx_-1,0+1)-Ez_->point(nx_-1,0));
                                Hx_->point(nx_-1,0) = c_hxh * Hx_->point(nx_-1,0) + c_hxb1 * pmlArr_[0].Bx_end_->point(0,0) - c_hxb0 * bxstore;
                                //Top Left

                                eps    = objArr_[pmlArr_[0].phys_Hy_->point(0,ny_-1)].dielectric(1.0);
                                sigxy = pmlArr_[0].sigma(0.5,eps);
                                sigyy = pmlArr_[1].sigma(0.0,eps);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                bystore = pmlArr_[0].By_->point(0,ny_-1);
                                pmlArr_[0].By_->point(0,ny_-1) = c_byb * pmlArr_[0].By_->point(0,ny_-1) + c_bye * (Ez_->point(0+1,ny_-1)-Ez_->point(0,ny_-1));
                                Hy_->point(0,ny_-1) = c_hyh * Hy_->point(0,ny_-1) + c_hyb1 * pmlArr_[0].By_->point(0,ny_-1) - c_hyb0 * bystore;
                                for(int ii = 1; ii < pmlArr_[0].thickness(); ii++)
                                {
                                    kapz = 1.0; sigz = 0.0; eps = 1.0;
                                    kapx = 1.0; kapy = 1.0;

                                    eps    = objArr_[pmlArr_[0].phys_Hx_->point(ii,0)].dielectric(1.0);
                                    sigxx  = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                    sigyx  = pmlArr_[0].sigma(0.5,eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[0].phys_Hy_->point(ii,0)].dielectric(1.0);
                                    sigxy  = pmlArr_[0].sigma(static_cast<double>(ii) + 0.5,eps);
                                    sigyy  = pmlArr_[0].sigma(0.0,eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    // Bot Left
                                    bxstore = pmlArr_[0].Bx_->point(ii,0);
                                    bystore = pmlArr_[0].By_->point(ii,0);
                                    pmlArr_[0].Bx_->point(ii,0) = c_bxb * pmlArr_[0].Bx_->point(ii,0) - c_bxe * (Ez_->point(ii,0+1)-Ez_->point(ii,0));
                                    pmlArr_[0].By_->point(ii,0) = c_byb * pmlArr_[0].By_->point(ii,0) + c_bye * (Ez_->point(ii+1,0)-Ez_->point(ii,0));
                                    Hx_->point(ii,0) = c_hxh * Hx_->point(ii,0) + c_hxb1 * pmlArr_[0].Bx_->point(ii,0) - c_hxb0 * bxstore;
                                    Hy_->point(ii,0) = c_hyh * Hy_->point(ii,0) + c_hyb1 * pmlArr_[0].By_->point(ii,0) - c_hyb0 * bystore;
                                    //Bot Right
                                    eps    = objArr_[pmlArr_[0].phys_Hx_end_->point(ii,0)].dielectric(1.0);
                                    sigxx  = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                    sigyx  = pmlArr_[0].sigma(0.5,eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[0].phys_Hy_end_->point(ii,0)].dielectric(1.0);
                                    sigxy = pmlArr_[0].sigma(static_cast<double>(ii) - 0.5,eps);
                                    sigyy = pmlArr_[1].sigma(0.0,eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    bxstore = pmlArr_[0].Bx_end_->point(ii,0);
                                    bystore = pmlArr_[0].By_end_->point(ii,0);
                                    pmlArr_[0].Bx_end_->point(ii,0) = c_bxb * pmlArr_[0].Bx_end_->point(ii,0) - c_bxe * (Ez_->point(nx_-1-ii,0+1)-Ez_->point(nx_-1-ii,0));
                                    pmlArr_[0].By_end_->point(ii,0) = c_byb * pmlArr_[0].By_end_->point(ii,0) + c_bye * (Ez_->point(nx_-1-ii+1,0)-Ez_->point(nx_-1-ii,0));
                                    Hx_->point(nx_-1-ii,0) = c_hxh * Hx_->point(nx_-1-ii,0) + c_hxb1 * pmlArr_[0].Bx_end_->point(ii,0) - c_hxb0 * bxstore;
                                    Hy_->point(nx_-1-ii,0) = c_hyh * Hy_->point(nx_-1-ii,0) + c_hyb1 * pmlArr_[0].By_end_->point(ii,0) - c_hyb0 * bystore;
                                    //Top Right
                                    eps    = objArr_[pmlArr_[0].phys_Hy_end_->point(ii,ny_-1)].dielectric(1.0);
                                    sigxy = pmlArr_[0].sigma(static_cast<double>(ii) - 0.5,eps);
                                    sigyy = pmlArr_[1].sigma(0.0,eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    bystore = pmlArr_[0].By_end_->point(ii,ny_-1);
                                    pmlArr_[0].By_end_->point(ii,ny_-1) = c_byb * pmlArr_[0].By_end_->point(ii,ny_-1) + c_bye * (Ez_->point(nx_-1-ii+1,ny_-1)-Ez_->point(nx_-1-ii,ny_-1));
                                    Hy_->point(nx_-1-ii,ny_-1) = c_hyh * Hy_->point(nx_-1-ii,ny_-1) + c_hyb1 * pmlArr_[0].By_end_->point(ii,ny_-1) - c_hyb0 * bystore;
                                    //Top Left
                                    eps    = objArr_[pmlArr_[0].phys_Hy_->point(ii,ny_-1)].dielectric(1.0);
                                    sigxy = pmlArr_[0].sigma(static_cast<double>(ii) + 0.5,eps);
                                    sigyy = pmlArr_[1].sigma(0.0,eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    bystore = pmlArr_[0].By_->point(ii,ny_-1);
                                    pmlArr_[0].By_->point(ii,ny_-1) = c_byb * pmlArr_[0].By_->point(ii,ny_-1) + c_bye * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
                                    Hy_->point(ii,ny_-1) = c_hyh * Hy_->point(ii,ny_-1) + c_hyb1 * pmlArr_[0].By_->point(ii,ny_ - 1) - c_hyb0 * bystore;
                                }
                                for(int jj = 1; jj < pmlArr_[0].thickness(); jj++)
                                {
                                    kapz = 1.0; sigz = 0.0; eps = 1.0;
                                    kapx = 1.0; kapy = 1.0;

                                    eps    = objArr_[pmlArr_[0].phys_Hx_->point(0,jj)].dielectric(1.0);
                                    sigxx = pmlArr_[0].sigma(0.0,eps);
                                    sigyx = pmlArr_[1].sigma(static_cast<double>(jj) + 0.5,eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[0].phys_Hy_->point(0,jj)].dielectric(1.0);
                                    sigxy = pmlArr_[0].sigma(0.5,eps);
                                    sigyy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    // Bot Left
                                    bxstore = pmlArr_[0].Bx_->point(0,jj);
                                    bystore = pmlArr_[0].By_->point(0,jj);
                                    pmlArr_[0].Bx_->point(0,jj) = c_bxb * pmlArr_[0].Bx_->point(0,jj) - c_bxe * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                                    pmlArr_[0].By_->point(0,jj) = c_byb * pmlArr_[0].By_->point(0,jj) + c_bye * (Ez_->point(0+1,jj)-Ez_->point(0,jj));
                                    Hx_->point(0,jj) = c_hxh * Hx_->point(0,jj) + c_hxb1 * pmlArr_[0].Bx_->point(0,jj) - c_hxb0 * bxstore;
                                    Hy_->point(0,jj) = c_hyh * Hy_->point(0,jj) + c_hyb1 * pmlArr_[0].By_->point(0,jj) - c_hyb0 * bystore;
                                    //Bot Right
                                    eps    = objArr_[pmlArr_[0].phys_Hx_end_->point(0,jj)].dielectric(1.0);
                                    sigxx = pmlArr_[0].sigma(0.0,eps);
                                    sigyx = pmlArr_[1].sigma(static_cast<double>(jj) + 0.5,eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    bxstore = pmlArr_[0].Bx_end_->point(0,jj);
                                    pmlArr_[0].Bx_end_->point(0,jj) = c_bxb * pmlArr_[0].Bx_end_->point(0,jj) - c_bxe * (Ez_->point(nx_-1,jj+1)-Ez_->point(nx_-1,jj));
                                    Hx_->point(nx_-1,jj) = c_hxh * Hx_->point(nx_-1,jj) + c_hxb1 * pmlArr_[0].Bx_end_->point(0,jj) - c_hxb0 * bxstore;
                                    //Top Rigt
                                    eps    = objArr_[pmlArr_[0].phys_Hx_end_->point(0,ny_-1-jj)].dielectric(1.0);
                                    sigxx = pmlArr_[0].sigma(0.0,eps);
                                    sigyx = pmlArr_[1].sigma(static_cast<double>(jj) - 0.5,eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    bxstore = pmlArr_[0].Bx_end_->point(0,ny_-1-jj);
                                    pmlArr_[0].Bx_end_->point(0,ny_-1-jj) = c_bxb * pmlArr_[0].Bx_end_->point(0,ny_-1-jj) - c_bxe * (Ez_->point(nx_-1,ny_-1-jj+1)-Ez_->point(nx_-1,ny_-1-jj));
                                    Hx_->point(nx_-1,ny_-1-jj) = c_hxh * Hx_->point(nx_-1,ny_-1-jj) + c_hxb1 * pmlArr_[0].Bx_end_->point(0,ny_-1-jj) - c_hxb0 * bxstore;
                                    //Top Left
                                    eps    = objArr_[pmlArr_[0].phys_Hx_->point(0,ny_-1-jj)].dielectric(1.0);
                                    sigxx = pmlArr_[0].sigma(0.0,eps);
                                    sigyx = pmlArr_[1].sigma(static_cast<double>(jj) - 0.5,eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[0].phys_Hy_->point(0,ny_-1-jj)].dielectric(1.0);
                                    sigxy = pmlArr_[0].sigma(0.5,eps);
                                    sigyy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    bxstore = pmlArr_[0].Bx_->point(0,ny_-1-jj);
                                    bystore = pmlArr_[0].By_->point(0,ny_-1-jj);
                                    pmlArr_[0].Bx_->point(0,ny_-1-jj) = c_bxb * pmlArr_[0].Bx_->point(0,ny_-1-jj) - c_bxe * (Ez_->point(0,ny_-1-jj+1)-Ez_->point(0,ny_-1-jj));
                                    pmlArr_[0].By_->point(0,ny_-1-jj) = c_byb * pmlArr_[0].By_->point(0,ny_-1-jj) + c_bye * (Ez_->point(0+1,ny_-1-jj)-Ez_->point(0,ny_-1-jj));
                                    Hx_->point(0,ny_-1-jj) = c_hxh * Hx_->point(0,ny_-1-jj) + c_hxb1 * pmlArr_[0].Bx_->point(0,ny_-1-jj) - c_hxb0 * bxstore;
                                    Hy_->point(0,ny_-1-jj) = c_hyh * Hy_->point(0,ny_-1-jj) + c_hyb1 * pmlArr_[0].By_->point(0,ny_-1-jj) - c_hyb0 * bystore;
                                }
                                for(int ii = 1; ii < pmlArr_[0].thickness(); ii++)
                                {
                                    for(int jj = 1; jj < pmlArr_[0].thickness(); jj++)
                                    {
                                        kapz = 1.0; sigz = 0.0; eps = 1.0;
                                        kapx = 1.0; kapy = 1.0;

                                        eps    = objArr_[pmlArr_[0].phys_Hx_->point(ii,jj)].dielectric(1.0);
                                        sigxx = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                        sigyx = pmlArr_[1].sigma(static_cast<double>(jj) + 0.5,eps);
                                        c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                        c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                        c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                        eps    = objArr_[pmlArr_[0].phys_Hy_->point(ii,jj)].dielectric(1.0);
                                        sigxy = pmlArr_[0].sigma(static_cast<double>(ii) + 0.5,eps);
                                        sigyy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                        c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                        c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                        // Bot Left
                                        bxstore = pmlArr_[0].Bx_->point(ii,jj);
                                        bystore = pmlArr_[0].By_->point(ii,jj);
                                        pmlArr_[0].Bx_->point(ii,jj) = c_bxb * pmlArr_[0].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                        pmlArr_[0].By_->point(ii,jj) = c_byb * pmlArr_[0].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                                        Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[0].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                                        Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[0].By_->point(ii,jj) - c_hyb0 * bystore;
                                        //Bot Right
                                        eps    = objArr_[pmlArr_[0].phys_Hx_end_->point(ii,jj)].dielectric(1.0);
                                        sigxx = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                        sigyx = pmlArr_[1].sigma(static_cast<double>(jj) + 0.5,eps);
                                        c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                        c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                        c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                        eps    = objArr_[pmlArr_[0].phys_Hy_end_->point(ii,jj)].dielectric(1.0);
                                        sigxy = pmlArr_[0].sigma(static_cast<double>(ii) - 0.5,eps);
                                        sigyy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                        c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                        c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                        bxstore = pmlArr_[0].Bx_end_->point(ii,jj);
                                        bystore = pmlArr_[0].By_end_->point(ii,jj);
                                        pmlArr_[0].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[0].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point(nx_-1-ii,jj+1)-Ez_->point(nx_-1-ii,jj));
                                        pmlArr_[0].By_end_->point(ii,jj) = c_byb * pmlArr_[0].By_end_->point(ii,jj) + c_bye * (Ez_->point(nx_-1-ii+1,jj)-Ez_->point(nx_-1-ii,jj));
                                        Hx_->point(nx_-1-ii,jj) = c_hxh * Hx_->point(nx_-1-ii,jj) + c_hxb1 * pmlArr_[0].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                                        Hy_->point(nx_-1-ii,jj) = c_hyh * Hy_->point(nx_-1-ii,jj) + c_hyb1 * pmlArr_[0].By_end_->point(ii,jj) - c_hyb0 * bystore;
                                        //Top Right
                                        eps    = objArr_[pmlArr_[0].phys_Hx_end_->point(ii,ny_-1-jj)].dielectric(1.0);
                                        sigyx = pmlArr_[1].sigma(static_cast<double>(jj) - 0.5,eps);
                                        sigxx = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                        c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                        c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                        c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                        eps    = objArr_[pmlArr_[0].phys_Hy_end_->point(ii,ny_-1-jj)].dielectric(1.0);
                                        sigxy = pmlArr_[0].sigma(static_cast<double>(ii) - 0.5,eps);
                                        sigyy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                        c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                        c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        bxstore = pmlArr_[0].Bx_end_->point(ii,ny_-1-jj);
                                        bystore = pmlArr_[0].By_end_->point(ii,ny_-1-jj);
                                        pmlArr_[0].Bx_end_->point(ii,ny_-1-jj) = c_bxb * pmlArr_[0].Bx_end_->point(ii,ny_-1-jj) - c_bxe * (Ez_->point(nx_-1-ii,ny_-1-jj+1)-Ez_->point(nx_-1-ii,ny_-1-jj));
                                        pmlArr_[0].By_end_->point(ii,ny_-1-jj) = c_byb * pmlArr_[0].By_end_->point(ii,ny_-1-jj) + c_bye * (Ez_->point(nx_-1-ii+1,ny_-1-jj)-Ez_->point(nx_-1-ii,ny_-1-jj));
                                        Hx_->point(nx_-1-ii,ny_-1-jj) = c_hxh * Hx_->point(nx_-1-ii,ny_-1-jj) + c_hxb1 * pmlArr_[0].Bx_end_->point(ii,ny_-1-jj) - c_hxb0 * bxstore;
                                        Hy_->point(nx_-1-ii,ny_-1-jj) = c_hyh * Hy_->point(nx_-1-ii,ny_-1-jj) + c_hyb1 * pmlArr_[0].By_end_->point(ii,ny_-1-jj) - c_hyb0 * bystore;
                                        //Top Left
                                        eps    = objArr_[pmlArr_[0].phys_Hx_->point(ii,ny_-1-jj)].dielectric(1.0);
                                        sigyx = pmlArr_[1].sigma(static_cast<double>(jj) - 0.5,eps);
                                        sigxx = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                        c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                        c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                        c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                        eps    = objArr_[pmlArr_[0].phys_Hy_->point(ii,ny_-1-jj)].dielectric(1.0);
                                        sigxy = pmlArr_[0].sigma(static_cast<double>(ii) + 0.5,eps);
                                        sigyy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                        c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                        c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        bxstore = pmlArr_[0].Bx_->point(ii,ny_-1-jj);
                                        bystore = pmlArr_[0].By_->point(ii,ny_-1-jj);
                                        pmlArr_[0].Bx_->point(ii,ny_-1-jj) = c_bxb * pmlArr_[0].Bx_->point(ii,ny_-1-jj) - c_bxe * (Ez_->point(ii,ny_-1-jj+1)-Ez_->point(ii,ny_-1-jj));
                                        pmlArr_[0].By_->point(ii,ny_-1-jj) = c_byb * pmlArr_[0].By_->point(ii,ny_-1-jj) + c_bye * (Ez_->point(ii+1,ny_-1-jj)-Ez_->point(ii,ny_-1-jj));
                                        Hx_->point(ii,ny_-1-jj) = c_hxh * Hx_->point(ii,ny_-1-jj) + c_hxb1 * pmlArr_[0].Bx_->point(ii,ny_-1-jj) - c_hxb0 * bxstore;
                                        Hy_->point(ii,ny_-1-jj) = c_hyh * Hy_->point(ii,ny_-1-jj) + c_hyb1 * pmlArr_[0].By_->point(ii,ny_-1-jj) - c_hyb0 * bystore;
                                    }
                                }
                            }
                        }
                        break;
                    }
                    case Y:
                    {
                        double eps = 1.0;
                        double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                        double sigz = 0.0;
                        double sigxx = 0.0;
                        double sigxy = 0.0;
                        double sigyx = 0.0;
                        double sigyy = 0.0;

                        double c_bxb = 0.0; double c_bxe = 0.0; double c_hxh = 0.0; double c_hxb0 = 0.0; double c_hxb1 = 0.0;
                        double c_byb = 0.0; double c_bye = 0.0; double c_hyh = 0.0; double c_hyb0 = 0.0; double c_hyb1 = 0.0;

                        complex<double>bxstore(0.0,0.0); complex<double>bystore(0.0,0.0);
                        if (xPML_== 0)
                        {
                            for (int jj = 0; jj < nx_-1; jj++)
                            {
                                eps    = objArr_[pmlArr_[kk].phys_Hx_->point(jj,0)].dielectric(1.0);
                                sigyx  = pmlArr_[kk].sigma(0.5,eps);
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                eps    = objArr_[pmlArr_[kk].phys_Hy_->point(jj,0)].dielectric(1.0);
                                sigyy  = pmlArr_[kk].sigma(0.0,eps);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                //Update both on the bottom side and Hy on the top
                                bxstore = pmlArr_[kk].Bx_->point(jj,0);
                                bystore = pmlArr_[kk].By_->point(jj,0);
                                pmlArr_[kk].Bx_->point(jj,0) = c_bxb * pmlArr_[kk].Bx_->point(jj,0) - c_bxe * (Ez_->point(jj,0+1)-Ez_->point(jj,0));
                                pmlArr_[kk].By_->point(jj,0) = c_byb * pmlArr_[kk].By_->point(jj,0) + c_bye * (Ez_->point(jj+1,0)-Ez_->point(jj,0));
                                //cout <<"H1"<<endl;
                                Hx_->point(jj,0) = c_hxh * Hx_->point(jj,0) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,0) - c_hxb0 * bxstore;
                                Hy_->point(jj,0) = c_hyh * Hy_->point(jj,0) + c_hyb1 * pmlArr_[kk].By_->point(jj,0) - c_hyb0 * bystore;

                                eps    = objArr_[pmlArr_[kk].phys_Hy_end_->point(jj,0)].dielectric(1.0);
                                sigyy  = pmlArr_[kk].sigma(0.0,eps);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                bystore = pmlArr_[kk].By_end_->point(jj,0);
                                pmlArr_[kk].By_end_->point(jj,0) = c_byb * pmlArr_[kk].By_end_->point(jj,0) + c_bye * (Ez_->point(jj+1, ny_-1)-Ez_->point(jj,ny_-1));
                                Hy_->point(jj,ny_-1) = c_hyh * Hy_->point(jj,ny_-1) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,0) - c_hyb0 * bystore;
                            }
                            //Top left corner
                            eps    = objArr_[pmlArr_[kk].phys_Hx_->point(nx_-1,0)].dielectric(1.0);
                            sigyx  = pmlArr_[kk].sigma(0.5,eps);
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                            bxstore = pmlArr_[kk].Bx_->point(nx_-1,0);
                            pmlArr_[kk].Bx_->point(nx_-1,0) = c_bxb * pmlArr_[kk].Bx_->point(nx_-1,0) - c_bxe * (Ez_->point(nx_-1,0+1)-Ez_->point(nx_-1,0));
                            Hx_->point(nx_-1,0) = c_hxh * Hx_->point(nx_-1,0) + c_hxb1 * pmlArr_[kk].Bx_->point(nx_-1,0) - c_hxb0 * bxstore;
                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                eps    = objArr_[pmlArr_[kk].phys_Hx_->point(nx_-1,ii)].dielectric(1.0);
                                sigyx = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5,eps);
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                //right col
                                bxstore = pmlArr_[kk].Bx_->point(nx_-1,ii);
                                pmlArr_[kk].Bx_->point(nx_-1,ii) = c_bxb * pmlArr_[kk].Bx_->point(nx_-1,ii) - c_bxe * (Ez_->point(nx_-1,ii+1)-Ez_->point(nx_-1,ii));
                                Hx_->point(nx_-1,ii) = c_hxh * Hx_->point(nx_-1,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(nx_-1,ii) - c_hxb0 * bxstore;
                                for(int jj = 0; jj < nx_-1; jj ++)
                                {
                                    //Update everything
                                    eps    = objArr_[pmlArr_[kk].phys_Hx_->point(jj,ii)].dielectric(1.0);
                                    sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) + 0.5),eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[kk].phys_Hy_->point(jj,ii)].dielectric(1.0);
                                    sigyy = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    bxstore = pmlArr_[kk].Bx_->point(jj,ii);
                                    bystore = pmlArr_[kk].By_->point(jj,ii);
                                    pmlArr_[kk].Bx_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_->point(jj,ii) - c_bxe * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                                    pmlArr_[kk].By_->point(jj,ii) = c_byb * pmlArr_[kk].By_->point(jj,ii) + c_bye * (Ez_->point(jj+1,ii)-Ez_->point(jj,ii));
                                    //cout <<"H1"<<endl;
                                    Hx_->point(jj,ii) = c_hxh * Hx_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,ii) - c_hxb0 * bxstore;
                                    Hy_->point(jj,ii) = c_hyh * Hy_->point(jj,ii) + c_hyb1 * pmlArr_[kk].By_->point(jj,ii) - c_hyb0 * bystore;

                                    eps    = objArr_[pmlArr_[kk].phys_Hx_end_->point(jj,ii)].dielectric(1.0);
                                    sigyx = pmlArr_[kk].sigma(static_cast<double>((ii)-0.5),eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[kk].phys_Hy_end_->point(jj,ii)].dielectric(1.0);
                                    sigyy = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    bxstore = pmlArr_[kk].Bx_end_->point(jj,ii);
                                    bystore = pmlArr_[kk].By_end_->point(jj,ii);
                                    pmlArr_[kk].Bx_end_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_end_->point(jj,ii) - c_bxe * (Ez_->point(jj,ny_-1-ii+1)-Ez_->point(jj,ny_-1-ii));
                                    pmlArr_[kk].By_end_->point(jj,ii) = c_byb * pmlArr_[kk].By_end_->point(jj,ii) + c_bye * (Ez_->point(jj+1,ny_-1-ii)-Ez_->point(jj,ny_-1-ii));
                                    //cout <<"H2"<<endl;
                                    Hx_->point(jj,ny_-1-ii) = c_hxh * Hx_->point(jj,ny_-1-ii) + c_hxb1 * pmlArr_[kk].Bx_end_->point(jj,ii) - c_hxb0 * bxstore;
                                    Hy_->point(jj,ny_-1-ii) = c_hyh * Hy_->point(jj,ny_-1-ii) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,ii) - c_hyb0 * bystore;
                                }
                                // for the right only upadate Hx
                                eps    = objArr_[pmlArr_[kk].phys_Hx_end_->point(nx_-1,ii)].dielectric(1.0);
                                sigyx = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5,eps);
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                bxstore = pmlArr_[kk].Bx_end_->point(nx_-1,ii);
                                pmlArr_[kk].Bx_end_->point(nx_-1,ii) = c_bxb * pmlArr_[kk].Bx_end_->point(nx_-1,ii) - c_bxe * (Ez_->point(nx_-1,ny_-1-ii+1)-Ez_->point(nx_-1,ny_-1-ii));
                                Hx_->point(nx_-1,ny_-1-ii) = c_hxh * Hx_->point(nx_-1,ny_-1-ii) + c_hxb1 * pmlArr_[kk].Bx_end_->point(nx_-1,ii) - c_hxb0 * bxstore;
                            }
                        }
                        else
                        {
                            for (int jj = xPML_; jj < nx_-xPML_; jj++)
                            {
                                eps    = objArr_[pmlArr_[kk].phys_Hx_->point(jj,0)].dielectric(1.0);
                                sigyx = pmlArr_[kk].sigma(0.5,eps);
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                eps    = objArr_[pmlArr_[kk].phys_Hy_->point(jj,0)].dielectric(1.0);
                                sigyy = pmlArr_[kk].sigma(0.0,eps);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                //Update both on the bottom side and Hy on the top
                                bxstore = pmlArr_[kk].Bx_->point(jj,0);
                                bystore = pmlArr_[kk].By_->point(jj,0);
                                pmlArr_[kk].Bx_->point(jj,0) = c_bxb * pmlArr_[kk].Bx_->point(jj,0) - c_bxe * (Ez_->point(jj,0+1)-Ez_->point(jj,0));
                                pmlArr_[kk].By_->point(jj,0) = c_byb * pmlArr_[kk].By_->point(jj,0) + c_bye * (Ez_->point(jj+1,0)-Ez_->point(jj,0));
                                //cout <<"H1"<<endl;
                                Hx_->point(jj,0) = c_hxh * Hx_->point(jj,0) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,0) - c_hxb0 * bxstore;
                                Hy_->point(jj,0) = c_hyh * Hy_->point(jj,0) + c_hyb1 * pmlArr_[kk].By_->point(jj,0) - c_hyb0 * bystore;

                                eps    = objArr_[pmlArr_[kk].phys_Hy_end_->point(jj,0)].dielectric(1.0);
                                sigyy = pmlArr_[kk].sigma(0.0,eps);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                bystore = pmlArr_[kk].By_end_->point(jj,0);
                                pmlArr_[kk].By_end_->point(jj,0) = c_byb * pmlArr_[kk].By_end_->point(jj,0) + c_bye * (Ez_->point(jj+1,ny_-1)-Ez_->point(jj,ny_-1));
                                Hy_->point(jj,ny_-1) = c_hyh * Hy_->point(jj,ny_-1) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,0) - c_hyb0 * bystore;
                            }
                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                for(int jj = xPML_; jj < nx_-xPML_; jj ++)
                                {
                                    //Update everything
                                    eps    = objArr_[pmlArr_[kk].phys_Hx_->point(jj,ii)].dielectric(1.0);
                                    sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) + 0.5),eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[kk].phys_Hy_->point(jj,ii)].dielectric(1.0);
                                    sigyy = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    bxstore = pmlArr_[kk].Bx_->point(jj,ii);
                                    bystore = pmlArr_[kk].By_->point(jj,ii);
                                    pmlArr_[kk].Bx_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_->point(jj,ii) - c_bxe * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                                    pmlArr_[kk].By_->point(jj,ii) = c_byb * pmlArr_[kk].By_->point(jj,ii) + c_bye * (Ez_->point(jj+1,ii)-Ez_->point(jj,ii));
                                    //cout <<"H1"<<endl;
                                    Hx_->point(jj,ii) = c_hxh * Hx_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,ii) - c_hxb0 * bxstore;
                                    Hy_->point(jj,ii) = c_hyh * Hy_->point(jj,ii) + c_hyb1 * pmlArr_[kk].By_->point(jj,ii) - c_hyb0 * bystore;

                                    eps    = objArr_[pmlArr_[kk].phys_Hx_end_->point(jj,ii)].dielectric(1.0);
                                    sigyx = pmlArr_[kk].sigma(static_cast<double>((ii)-0.5),eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[kk].phys_Hy_end_->point(jj,ii)].dielectric(1.0);
                                    sigyy = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    bxstore = pmlArr_[kk].Bx_end_->point(jj,ii);
                                    bystore = pmlArr_[kk].By_end_->point(jj,ii);
                                    pmlArr_[kk].Bx_end_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_end_->point(jj,ii) - c_bxe * (Ez_->point(jj,ny_-1-ii+1)-Ez_->point(jj,ny_-1-ii));
                                    pmlArr_[kk].By_end_->point(jj,ii) = c_byb * pmlArr_[kk].By_end_->point(jj,ii) + c_bye * (Ez_->point(jj+1,ny_-1-ii)-Ez_->point(jj,ny_-1-ii));
                                    //cout <<"H2"<<endl;
                                    Hx_->point(jj,ny_-1-ii) = c_hxh * Hx_->point(jj,ny_-1-ii) + c_hxb1 * pmlArr_[kk].Bx_end_->point(jj,ii) - c_hxb0 * bxstore;
                                    Hy_->point(jj,ny_-1-ii) = c_hyh * Hy_->point(jj,ny_-1-ii) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,ii) - c_hyb0 * bystore;
                                }
                            }
                            if(kk==0)
                            {
                                kapz = 1.0; sigz = 0.0; eps = 1.0;
                                kapx = 1.0; kapy = 1.0;

                                eps    = objArr_[pmlArr_[1].phys_Hx_->point(0,0)].dielectric(1.0);
                                sigxx = pmlArr_[1].sigma(0.0,eps);
                                sigyx = pmlArr_[0].sigma(0.5,eps);
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                eps    = objArr_[pmlArr_[1].phys_Hy_->point(0,0)].dielectric(1.0);
                                sigxy = pmlArr_[1].sigma(0.5,eps);
                                sigyy = pmlArr_[0].sigma(0.0,eps);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                // Bot Left
                                bxstore = pmlArr_[1].Bx_->point(0,0);
                                bystore = pmlArr_[1].By_->point(0,0);
                                pmlArr_[1].Bx_->point(0,0) = c_bxb * pmlArr_[1].Bx_->point(0,0) - c_bxe * (Ez_->point(0,0+1)-Ez_->point(0,0));
                                pmlArr_[1].By_->point(0,0) = c_byb * pmlArr_[1].By_->point(0,0) + c_bye * (Ez_->point(0+1,0)-Ez_->point(0,0));
                                Hx_->point(0,0) = c_hxh * Hx_->point(0,0) + c_hxb1 * pmlArr_[1].Bx_->point(0,0) - c_hxb0 * bxstore;
                                Hy_->point(0,0) = c_hyh * Hy_->point(0,0) + c_hyb1 * pmlArr_[1].By_->point(0,0) - c_hyb0 * bystore;
                                //Bot Right
                                eps    = objArr_[pmlArr_[1].phys_Hx_end_->point(0,0)].dielectric(1.0);
                                sigxx = pmlArr_[1].sigma(0.0,eps);
                                sigyx = pmlArr_[0].sigma(0.5,eps);
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                bxstore = pmlArr_[1].Bx_end_->point(0,0);
                                pmlArr_[1].Bx_end_->point(0,0) = c_bxb * pmlArr_[1].Bx_end_->point(0,0) - c_bxe * (Ez_->point(nx_-1,0+1)-Ez_->point(nx_-1,0));
                                Hx_->point(nx_-1,0) = c_hxh * Hx_->point(nx_-1,0) + c_hxb1 * pmlArr_[1].Bx_end_->point(0,0) - c_hxb0 * bxstore;
                                //Top Left

                                eps    = objArr_[pmlArr_[1].phys_Hy_->point(0,ny_-1)].dielectric(1.0);
                                sigxy = pmlArr_[1].sigma(0.5,eps);
                                sigyy = pmlArr_[0].sigma(0.0,eps);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                bystore = pmlArr_[1].By_->point(0,ny_-1);
                                pmlArr_[1].By_->point(0,ny_-1) = c_byb * pmlArr_[1].By_->point(0,ny_-1) + c_bye * (Ez_->point(0+1,ny_-1)-Ez_->point(0,ny_-1));
                                Hy_->point(0,ny_-1) = c_hyh * Hy_->point(0,ny_-1) + c_hyb1 * pmlArr_[1].By_->point(0,ny_-1) - c_hyb0 * bystore;
                                for(int ii = 1; ii < pmlArr_[1].thickness(); ii++)
                                {
                                    kapz = 1.0; sigz = 0.0; eps = 1.0;
                                    kapx = 1.0; kapy = 1.0;

                                    eps    = objArr_[pmlArr_[1].phys_Hx_->point(ii,0)].dielectric(1.0);
                                    sigxx  = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                    sigyx  = pmlArr_[1].sigma(0.5,eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[1].phys_Hy_->point(ii,0)].dielectric(1.0);
                                    sigxy  = pmlArr_[1].sigma(static_cast<double>(ii) + 0.5,eps);
                                    sigyy  = pmlArr_[1].sigma(0.0,eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    // Bot Left
                                    bxstore = pmlArr_[1].Bx_->point(ii,0);
                                    bystore = pmlArr_[1].By_->point(ii,0);
                                    pmlArr_[1].Bx_->point(ii,0) = c_bxb * pmlArr_[1].Bx_->point(ii,0) - c_bxe * (Ez_->point(ii,0+1)-Ez_->point(ii,0));
                                    pmlArr_[1].By_->point(ii,0) = c_byb * pmlArr_[1].By_->point(ii,0) + c_bye * (Ez_->point(ii+1,0)-Ez_->point(ii,0));
                                    Hx_->point(ii,0) = c_hxh * Hx_->point(ii,0) + c_hxb1 * pmlArr_[1].Bx_->point(ii,0) - c_hxb0 * bxstore;
                                    Hy_->point(ii,0) = c_hyh * Hy_->point(ii,0) + c_hyb1 * pmlArr_[1].By_->point(ii,0) - c_hyb0 * bystore;
                                    //Bot Right
                                    eps    = objArr_[pmlArr_[1].phys_Hx_end_->point(ii,0)].dielectric(1.0);
                                    sigxx  = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                    sigyx  = pmlArr_[1].sigma(0.5,eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[1].phys_Hy_end_->point(ii,0)].dielectric(1.0);
                                    sigxy = pmlArr_[1].sigma(static_cast<double>(ii) - 0.5,eps);
                                    sigyy = pmlArr_[0].sigma(0.0,eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    bxstore = pmlArr_[1].Bx_end_->point(ii,0);
                                    bystore = pmlArr_[1].By_end_->point(ii,0);
                                    pmlArr_[1].Bx_end_->point(ii,0) = c_bxb * pmlArr_[1].Bx_end_->point(ii,0) - c_bxe * (Ez_->point(nx_-1-ii,0+1)-Ez_->point(nx_-1-ii,0));
                                    pmlArr_[1].By_end_->point(ii,0) = c_byb * pmlArr_[1].By_end_->point(ii,0) + c_bye * (Ez_->point(nx_-1-ii+1,0)-Ez_->point(nx_-1-ii,0));
                                    Hx_->point(nx_-1-ii,0) = c_hxh * Hx_->point(nx_-1-ii,0) + c_hxb1 * pmlArr_[1].Bx_end_->point(ii,0) - c_hxb0 * bxstore;
                                    Hy_->point(nx_-1-ii,0) = c_hyh * Hy_->point(nx_-1-ii,0) + c_hyb1 * pmlArr_[1].By_end_->point(ii,0) - c_hyb0 * bystore;
                                    //Top Right
                                    eps    = objArr_[pmlArr_[1].phys_Hy_end_->point(ii,ny_-1)].dielectric(1.0);
                                    sigxy = pmlArr_[1].sigma(static_cast<double>(ii) - 0.5,eps);
                                    sigyy = pmlArr_[0].sigma(0.0,eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    bystore = pmlArr_[1].By_end_->point(ii,ny_-1);
                                    pmlArr_[1].By_end_->point(ii,ny_-1) = c_byb * pmlArr_[1].By_end_->point(ii,ny_-1) + c_bye * (Ez_->point(nx_-1-ii+1,ny_-1)-Ez_->point(nx_-1-ii,ny_-1));
                                    Hy_->point(nx_-1-ii,ny_-1) = c_hyh * Hy_->point(nx_-1-ii,ny_-1) + c_hyb1 * pmlArr_[1].By_end_->point(ii,ny_-1) - c_hyb0 * bystore;
                                    //Top Left
                                    eps    = objArr_[pmlArr_[1].phys_Hy_->point(ii,ny_-1)].dielectric(1.0);
                                    sigxy = pmlArr_[1].sigma(static_cast<double>(ii) + 0.5,eps);
                                    sigyy = pmlArr_[0].sigma(0.0,eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    bystore = pmlArr_[1].By_->point(ii,ny_-1);
                                    pmlArr_[1].By_->point(ii,ny_-1) = c_byb * pmlArr_[1].By_->point(ii,ny_-1) + c_bye * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
                                    Hy_->point(ii,ny_-1) = c_hyh * Hy_->point(ii,ny_-1) + c_hyb1 * pmlArr_[1].By_->point(ii,ny_ - 1) - c_hyb0 * bystore;
                                }
                                for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                                {
                                    kapz = 1.0; sigz = 0.0; eps = 1.0;
                                    kapx = 1.0; kapy = 1.0;

                                    eps    = objArr_[pmlArr_[1].phys_Hx_->point(0,jj)].dielectric(1.0);
                                    sigxx = pmlArr_[1].sigma(0.0,eps);
                                    sigyx = pmlArr_[0].sigma(static_cast<double>(jj) + 0.5,eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[1].phys_Hy_->point(0,jj)].dielectric(1.0);
                                    sigxy = pmlArr_[1].sigma(0.5,eps);
                                    sigyy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                    // Bot Left
                                    bxstore = pmlArr_[1].Bx_->point(0,jj);
                                    bystore = pmlArr_[1].By_->point(0,jj);
                                    pmlArr_[1].Bx_->point(0,jj) = c_bxb * pmlArr_[1].Bx_->point(0,jj) - c_bxe * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                                    pmlArr_[1].By_->point(0,jj) = c_byb * pmlArr_[1].By_->point(0,jj) + c_bye * (Ez_->point(0+1,jj)-Ez_->point(0,jj));
                                    Hx_->point(0,jj) = c_hxh * Hx_->point(0,jj) + c_hxb1 * pmlArr_[1].Bx_->point(0,jj) - c_hxb0 * bxstore;
                                    Hy_->point(0,jj) = c_hyh * Hy_->point(0,jj) + c_hyb1 * pmlArr_[1].By_->point(0,jj) - c_hyb0 * bystore;
                                    //Bot Right
                                    eps    = objArr_[pmlArr_[1].phys_Hx_end_->point(0,jj)].dielectric(1.0);
                                    sigxx = pmlArr_[1].sigma(0.0,eps);
                                    sigyx = pmlArr_[0].sigma(static_cast<double>(jj) + 0.5,eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    bxstore = pmlArr_[1].Bx_end_->point(0,jj);
                                    pmlArr_[1].Bx_end_->point(0,jj) = c_bxb * pmlArr_[1].Bx_end_->point(0,jj) - c_bxe * (Ez_->point(nx_-1,jj+1)-Ez_->point(nx_-1,jj));
                                    Hx_->point(nx_-1,jj) = c_hxh * Hx_->point(nx_-1,jj) + c_hxb1 * pmlArr_[1].Bx_end_->point(0,jj) - c_hxb0 * bxstore;
                                    //Top Rigt
                                    eps    = objArr_[pmlArr_[1].phys_Hx_end_->point(0,ny_-1-jj)].dielectric(1.0);
                                    sigxx = pmlArr_[1].sigma(0.0,eps);
                                    sigyx = pmlArr_[0].sigma(static_cast<double>(jj) - 0.5,eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    bxstore = pmlArr_[1].Bx_end_->point(0,ny_-1-jj);
                                    pmlArr_[1].Bx_end_->point(0,ny_-1-jj) = c_bxb * pmlArr_[1].Bx_end_->point(0,ny_-1-jj) - c_bxe * (Ez_->point(nx_-1,ny_-1-jj+1)-Ez_->point(nx_-1,ny_-1-jj));
                                    Hx_->point(nx_-1,ny_-1-jj) = c_hxh * Hx_->point(nx_-1,ny_-1-jj) + c_hxb1 * pmlArr_[1].Bx_end_->point(0,ny_-1-jj) - c_hxb0 * bxstore;
                                    //Top Left
                                    eps    = objArr_[pmlArr_[1].phys_Hx_->point(0,ny_-1-jj)].dielectric(1.0);
                                    sigxx = pmlArr_[1].sigma(0.0,eps);
                                    sigyx = pmlArr_[0].sigma(static_cast<double>(jj) - 0.5,eps);
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                    eps    = objArr_[pmlArr_[1].phys_Hy_->point(0,ny_-1-jj)].dielectric(1.0);
                                    sigxy = pmlArr_[1].sigma(0.5,eps);
                                    sigyy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    bxstore = pmlArr_[1].Bx_->point(0,ny_-1-jj);
                                    bystore = pmlArr_[1].By_->point(0,ny_-1-jj);
                                    pmlArr_[1].Bx_->point(0,ny_-1-jj) = c_bxb * pmlArr_[1].Bx_->point(0,ny_-1-jj) - c_bxe * (Ez_->point(0,ny_-1-jj+1)-Ez_->point(0,ny_-1-jj));
                                    pmlArr_[1].By_->point(0,ny_-1-jj) = c_byb * pmlArr_[1].By_->point(0,ny_-1-jj) + c_bye * (Ez_->point(0+1,ny_-1-jj)-Ez_->point(0,ny_-1-jj));
                                    Hx_->point(0,ny_-1-jj) = c_hxh * Hx_->point(0,ny_-1-jj) + c_hxb1 * pmlArr_[1].Bx_->point(0,ny_-1-jj) - c_hxb0 * bxstore;
                                    Hy_->point(0,ny_-1-jj) = c_hyh * Hy_->point(0,ny_-1-jj) + c_hyb1 * pmlArr_[1].By_->point(0,ny_-1-jj) - c_hyb0 * bystore;
                                }
                                for(int ii = 1; ii < pmlArr_[1].thickness(); ii++)
                                {
                                    for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                                    {
                                        kapz = 1.0; sigz = 0.0; eps = 1.0;
                                        kapx = 1.0; kapy = 1.0;

                                        eps    = objArr_[pmlArr_[1].phys_Hx_->point(ii,jj)].dielectric(1.0);
                                        sigxx = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                        sigyx = pmlArr_[0].sigma(static_cast<double>(jj) + 0.5,eps);
                                        c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                        c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                        c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                        eps    = objArr_[pmlArr_[1].phys_Hy_->point(ii,jj)].dielectric(1.0);
                                        sigxy = pmlArr_[1].sigma(static_cast<double>(ii) + 0.5,eps);
                                        sigyy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                        c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                        c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                        // Bot Left
                                        bxstore = pmlArr_[1].Bx_->point(ii,jj);
                                        bystore = pmlArr_[1].By_->point(ii,jj);
                                        pmlArr_[1].Bx_->point(ii,jj) = c_bxb * pmlArr_[1].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                        pmlArr_[1].By_->point(ii,jj) = c_byb * pmlArr_[1].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                                        Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[1].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                                        Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[1].By_->point(ii,jj) - c_hyb0 * bystore;
                                        //Bot Right
                                        eps    = objArr_[pmlArr_[1].phys_Hx_end_->point(ii,jj)].dielectric(1.0);
                                        sigxx = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                        sigyx = pmlArr_[0].sigma(static_cast<double>(jj) + 0.5,eps);
                                        c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                        c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                        c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                        eps    = objArr_[pmlArr_[1].phys_Hy_end_->point(ii,jj)].dielectric(1.0);
                                        sigxy = pmlArr_[1].sigma(static_cast<double>(ii) - 0.5,eps);
                                        sigyy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                        c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                        c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);

                                        bxstore = pmlArr_[1].Bx_end_->point(ii,jj);
                                        bystore = pmlArr_[1].By_end_->point(ii,jj);
                                        pmlArr_[1].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[1].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point(nx_-1-ii,jj+1)-Ez_->point(nx_-1-ii,jj));
                                        pmlArr_[1].By_end_->point(ii,jj) = c_byb * pmlArr_[1].By_end_->point(ii,jj) + c_bye * (Ez_->point(nx_-1-ii+1,jj)-Ez_->point(nx_-1-ii,jj));
                                        Hx_->point(nx_-1-ii,jj) = c_hxh * Hx_->point(nx_-1-ii,jj) + c_hxb1 * pmlArr_[1].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                                        Hy_->point(nx_-1-ii,jj) = c_hyh * Hy_->point(nx_-1-ii,jj) + c_hyb1 * pmlArr_[1].By_end_->point(ii,jj) - c_hyb0 * bystore;
                                        //Top Right
                                        eps    = objArr_[pmlArr_[1].phys_Hx_end_->point(ii,ny_-1-jj)].dielectric(1.0);
                                        sigyx = pmlArr_[0].sigma(static_cast<double>(jj) - 0.5,eps);
                                        sigxx = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                        c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                        c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                        c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                        eps    = objArr_[pmlArr_[1].phys_Hy_end_->point(ii,ny_-1-jj)].dielectric(1.0);
                                        sigxy = pmlArr_[1].sigma(static_cast<double>(ii) - 0.5,eps);
                                        sigyy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                        c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                        c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        bxstore = pmlArr_[1].Bx_end_->point(ii,ny_-1-jj);
                                        bystore = pmlArr_[1].By_end_->point(ii,ny_-1-jj);
                                        pmlArr_[1].Bx_end_->point(ii,ny_-1-jj) = c_bxb * pmlArr_[1].Bx_end_->point(ii,ny_-1-jj) - c_bxe * (Ez_->point(nx_-1-ii,ny_-1-jj+1)-Ez_->point(nx_-1-ii,ny_-1-jj));
                                        pmlArr_[1].By_end_->point(ii,ny_-1-jj) = c_byb * pmlArr_[1].By_end_->point(ii,ny_-1-jj) + c_bye * (Ez_->point(nx_-1-ii+1,ny_-1-jj)-Ez_->point(nx_-1-ii,ny_-1-jj));
                                        Hx_->point(nx_-1-ii,ny_-1-jj) = c_hxh * Hx_->point(nx_-1-ii,ny_-1-jj) + c_hxb1 * pmlArr_[1].Bx_end_->point(ii,ny_-1-jj) - c_hxb0 * bxstore;
                                        Hy_->point(nx_-1-ii,ny_-1-jj) = c_hyh * Hy_->point(nx_-1-ii,ny_-1-jj) + c_hyb1 * pmlArr_[1].By_end_->point(ii,ny_-1-jj) - c_hyb0 * bystore;
                                        //Top Left
                                        eps    = objArr_[pmlArr_[1].phys_Hx_->point(ii,ny_-1-jj)].dielectric(1.0);
                                        sigyx = pmlArr_[0].sigma(static_cast<double>(jj) - 0.5,eps);
                                        sigxx = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                        c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                        c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                        c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_);

                                        eps    = objArr_[pmlArr_[1].phys_Hy_->point(ii,ny_-1-jj)].dielectric(1.0);
                                        sigxy = pmlArr_[1].sigma(static_cast<double>(ii) + 0.5,eps);
                                        sigyy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                        c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                        c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                        c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_);
                                        bxstore = pmlArr_[1].Bx_->point(ii,ny_-1-jj);
                                        bystore = pmlArr_[1].By_->point(ii,ny_-1-jj);
                                        pmlArr_[1].Bx_->point(ii,ny_-1-jj) = c_bxb * pmlArr_[1].Bx_->point(ii,ny_-1-jj) - c_bxe * (Ez_->point(ii,ny_-1-jj+1)-Ez_->point(ii,ny_-1-jj));
                                        pmlArr_[1].By_->point(ii,ny_-1-jj) = c_byb * pmlArr_[1].By_->point(ii,ny_-1-jj) + c_bye * (Ez_->point(ii+1,ny_-1-jj)-Ez_->point(ii,ny_-1-jj));
                                        Hx_->point(ii,ny_-1-jj) = c_hxh * Hx_->point(ii,ny_-1-jj) + c_hxb1 * pmlArr_[1].Bx_->point(ii,ny_-1-jj) - c_hxb0 * bxstore;
                                        Hy_->point(ii,ny_-1-jj) = c_hyh * Hy_->point(ii,ny_-1-jj) + c_hyb1 * pmlArr_[1].By_->point(ii,ny_-1-jj) - c_hyb0 * bystore;
                                    }
                                }
                            }
                        }
                        break;
                    }
                    case Z:
                    {
                        throw logic_error("Z dir not implimented");
                        break;
                    }
                    default:
                        throw logic_error("how did you get here, it is a switch default try again");
                        break;
                }
            }
        }
    }*/
    /*else
    {
        if(Ez_)
        {
            for(int kk = 0; kk < pmlArr_.size(); kk++)
            {
                double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                double sigz = 0.0;
                switch(pmlArr_[kk].d())
                {
                    case X:
                    {
                        complex<double> bxstore(0.0,0.0); complex<double> bystore(0.0,0.0);
                        if (yPML_== 0)
                        {
                            for (int jj = 0; jj < ny_-1; jj++)
                            {
                                //Update both on the left side and Hx on the right
                                bxstore = pmlArr_[kk].Bx_->point(0,jj);
                                bystore = pmlArr_[kk].By_->point(0,jj);

                                pmlArr_[kk].Bx_->point(0,jj) = pmlArr_[kk].c_bxb_0_ -> point(0,jj) * pmlArr_[kk].Bx_->point(0,jj) - pmlArr_[kk].c_bxe_0_ -> point(0,jj) * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                                pmlArr_[kk].By_->point(0,jj) = pmlArr_[kk].c_byb_0_ -> point(0,jj) * pmlArr_[kk].By_->point(0,jj) + pmlArr_[kk].c_bye_0_ -> point(0,jj) * (Ez_->point(0+1,jj)-Ez_->point(0,jj));

                                Hx_->point(0,jj) = pmlArr_[kk].c_hxh_0_ -> point(0,jj) * Hx_->point(0,jj) + pmlArr_[kk].c_hxb1_0_ -> point(0,jj) * pmlArr_[kk].Bx_->point(0,jj) - pmlArr_[kk].c_hxb0_0_ -> point(0,jj) * bxstore;
                                Hy_->point(0,jj) = pmlArr_[kk].c_hyh_0_ -> point(0,jj) * Hy_->point(0,jj) + pmlArr_[kk].c_hyb1_0_ -> point(0,jj) * pmlArr_[kk].By_->point(0,jj) - pmlArr_[kk].c_hyb0_0_ -> point(0,jj) * bystore;

                                bxstore = pmlArr_[kk].Bx_end_->point(0,jj);
                                pmlArr_[kk].Bx_end_->point(0,jj) = pmlArr_[kk].c_bxb_n_ -> point(0,jj) * pmlArr_[kk].Bx_end_->point(0,jj) - pmlArr_[kk].c_bxe_n_ -> point(0,jj) * (Ez_->point(nx_-1, jj+1)-Ez_->point(nx_-1,jj));
                                Hx_->point(nx_-1,jj) = pmlArr_[kk].c_hxh_n_ -> point(0,jj) * Hx_->point(nx_-1,jj) + pmlArr_[kk].c_hxb1_n_ -> point(0,jj) * pmlArr_[kk].Bx_end_->point(0,jj) - pmlArr_[kk].c_hxb0_n_ -> point(0,jj) * bxstore;
                            }
                            //Top left corner
                            bystore = pmlArr_[kk].By_->point(0,ny_-1);
                            pmlArr_[kk].By_->point(0,ny_-1) = pmlArr_[kk].c_byb_0_ -> point(0,ny_-1) * pmlArr_[kk].By_->point(0,ny_-1) + pmlArr_[kk].c_bye_0_ -> point(0,ny_-1) * (Ez_->point(0+1,ny_-1)-Ez_->point(0,ny_-1));
                            Hy_->point(0,ny_-1) = pmlArr_[kk].c_hyh_0_ -> point(0,ny_-1) * Hy_->point(0,ny_-1) + pmlArr_[kk].c_hyb1_0_ -> point(0,ny_-1) * pmlArr_[kk].By_->point(0,ny_-1) - pmlArr_[kk].c_hyb0_0_ -> point(0,ny_-1) * bystore;

                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                //top row
                                bystore = pmlArr_[kk].By_->point(ii,ny_-1);
                                pmlArr_[kk].By_->point(ii,ny_-1) = pmlArr_[kk].c_byb_0_ -> point(ii,ny_-1) * pmlArr_[kk].By_->point(ii,ny_-1) + pmlArr_[kk].c_bye_0_ -> point(ii,ny_-1) * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
                                Hy_->point(ii,ny_-1) = pmlArr_[kk].c_hyh_0_ -> point(ii,ny_-1) * Hy_->point(ii,ny_-1) + pmlArr_[kk].c_hyb1_0_ -> point(ii,ny_-1) * pmlArr_[kk].By_->point(ii,ny_-1) - pmlArr_[kk].c_hyb0_0_ -> point(ii,ny_-1) * bystore;

                                for(int jj = 0; jj < ny_-1; jj ++)
                                {
                                    //Update everything
                                    bxstore = pmlArr_[kk].Bx_->point(ii,jj);
                                    bystore = pmlArr_[kk].By_->point(ii,jj);

                                    pmlArr_[kk].Bx_->point(ii,jj) = pmlArr_[kk].c_bxb_0_ -> point(ii,jj) * pmlArr_[kk].Bx_->point(ii,jj) - pmlArr_[kk].c_bxe_0_ -> point(ii,jj) * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                    pmlArr_[kk].By_->point(ii,jj) = pmlArr_[kk].c_byb_0_ -> point(ii,jj) * pmlArr_[kk].By_->point(ii,jj) + pmlArr_[kk].c_bye_0_ -> point(ii,jj) * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));

                                    Hx_->point(ii,jj) = pmlArr_[kk].c_hxh_0_ -> point(ii,jj) * Hx_->point(ii,jj) + pmlArr_[kk].c_hxb1_0_ -> point(ii,jj) * pmlArr_[kk].Bx_->point(ii,jj) - pmlArr_[kk].c_hxb0_0_ -> point(ii,jj) * bxstore;
                                    Hy_->point(ii,jj) = pmlArr_[kk].c_hyh_0_ -> point(ii,jj) * Hy_->point(ii,jj) + pmlArr_[kk].c_hyb1_0_ -> point(ii,jj) * pmlArr_[kk].By_->point(ii,jj) - pmlArr_[kk].c_hyb0_0_ -> point(ii,jj) * bystore;

                                    bxstore = pmlArr_[kk].Bx_end_->point(ii,jj);
                                    bystore = pmlArr_[kk].By_end_->point(ii,jj);

                                    pmlArr_[kk].Bx_end_->point(ii,jj) = pmlArr_[kk].c_bxb_n_ -> point(ii,jj) * pmlArr_[kk].Bx_end_->point(ii,jj) - pmlArr_[kk].c_bxe_n_ -> point(ii,jj) * (Ez_->point(nx_-1-ii, jj+1)-Ez_->point(nx_-1-ii,jj));
                                    pmlArr_[kk].By_end_->point(ii,jj) = pmlArr_[kk].c_byb_n_ -> point(ii,jj) * pmlArr_[kk].By_end_->point(ii,jj) + pmlArr_[kk].c_bye_n_ -> point(ii,jj) * (Ez_->point(nx_-1-ii+1, jj)-Ez_->point(nx_-1-ii,jj));

                                    Hx_->point(nx_-1-ii,jj) = pmlArr_[kk].c_hxh_n_ -> point(ii,jj) * Hx_->point((nx_-1) - ii,jj) + pmlArr_[kk].c_hxb1_n_ -> point(ii,jj) * pmlArr_[kk].Bx_end_->point(ii,jj) - pmlArr_[kk].c_hxb0_n_ -> point(ii,jj) * bxstore;
                                    Hy_->point(nx_-1-ii,jj) = pmlArr_[kk].c_hyh_n_ -> point(ii,jj) * Hy_->point((nx_-1) - ii,jj) + pmlArr_[kk].c_hyb1_n_ -> point(ii,jj) * pmlArr_[kk].By_end_->point(ii,jj) - pmlArr_[kk].c_hyb0_n_ -> point(ii,jj) * bystore;
                                }
                                // for the top row only upadate Hy

                                bystore = pmlArr_[kk].By_end_->point(ii,ny_-1);
                                pmlArr_[kk].By_end_->point(ii,ny_-1) = pmlArr_[kk].c_byb_n_ -> point(ii,ny_-1) * pmlArr_[kk].By_end_->point(ii,ny_-1) + pmlArr_[kk].c_bye_n_ -> point(ii,ny_-1) * (Ez_->point(nx_-1-ii+1, ny_-1)-Ez_->point(nx_-1-ii,ny_-1));
                                Hy_->point(nx_-1-ii,ny_-1) = pmlArr_[kk].c_hyh_n_ -> point(ii,ny_-1) * Hy_->point(nx_-1-ii,ny_-1) + pmlArr_[kk].c_hyb1_n_ -> point(ii,ny_-1) * pmlArr_[kk].By_end_->point(ii,ny_-1) - pmlArr_[kk].c_hyb0_n_ -> point(ii,ny_-1) * bystore;
                            }
                        }
                        else
                        {
                            for (int jj = yPML_; jj < ny_- yPML_; jj++)
                            {
                                //Update both on the left side and Hx on the right
                                bxstore = pmlArr_[kk].Bx_->point(0,jj);
                                bystore = pmlArr_[kk].By_->point(0,jj);

                                pmlArr_[kk].Bx_->point(0,jj) = pmlArr_[kk].c_bxb_0_ -> point(0,jj) * pmlArr_[kk].Bx_->point(0,jj) - pmlArr_[kk].c_bxe_0_ -> point(0,jj) * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                                pmlArr_[kk].By_->point(0,jj) = pmlArr_[kk].c_byb_0_ -> point(0,jj) * pmlArr_[kk].By_->point(0,jj) + pmlArr_[kk].c_bye_0_ -> point(0,jj) * (Ez_->point(0+1,jj)-Ez_->point(0,jj));

                                Hx_->point(0,jj) = pmlArr_[kk].c_hxh_0_ -> point(0,jj) * Hx_->point(0,jj) + pmlArr_[kk].c_hxb1_0_ -> point(0,jj) * pmlArr_[kk].Bx_->point(0,jj) - pmlArr_[kk].c_hxb0_0_ -> point(0,jj) * bxstore;
                                Hy_->point(0,jj) = pmlArr_[kk].c_hyh_0_ -> point(0,jj) * Hy_->point(0,jj) + pmlArr_[kk].c_hyb1_0_ -> point(0,jj) * pmlArr_[kk].By_->point(0,jj) - pmlArr_[kk].c_hyb0_0_ -> point(0,jj) * bystore;

                                bxstore = pmlArr_[kk].Bx_end_->point(0,jj);
                                pmlArr_[kk].Bx_end_->point(0,jj) = pmlArr_[kk].c_bxb_n_ -> point(0,jj) * pmlArr_[kk].Bx_end_->point(0,jj) - pmlArr_[kk].c_bxe_n_ -> point(0,jj) * (Ez_->point(nx_-1, jj+1)-Ez_->point(nx_-1,jj));
                                Hx_->point(nx_-1,jj) = pmlArr_[kk].c_hxh_n_ -> point(0,jj) * Hx_->point(nx_-1,jj) + pmlArr_[kk].c_hxb1_n_ -> point(0,jj) * pmlArr_[kk].Bx_end_->point(0,jj) - pmlArr_[kk].c_hxb0_n_ -> point(0,jj) * bxstore;
                            }
                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                for(int jj = yPML_; jj < ny_- yPML_; jj ++)
                                {
                                    //Update everything
                                    bxstore = pmlArr_[kk].Bx_->point(ii,jj);
                                    bystore = pmlArr_[kk].By_->point(ii,jj);

                                    pmlArr_[kk].Bx_->point(ii,jj) = pmlArr_[kk].c_bxb_0_ -> point(ii,jj) * pmlArr_[kk].Bx_->point(ii,jj) - pmlArr_[kk].c_bxe_0_ -> point(ii,jj) * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                    pmlArr_[kk].By_->point(ii,jj) = pmlArr_[kk].c_byb_0_ -> point(ii,jj) * pmlArr_[kk].By_->point(ii,jj) + pmlArr_[kk].c_bye_0_ -> point(ii,jj) * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));

                                    Hx_->point(ii,jj) = pmlArr_[kk].c_hxh_0_ -> point(ii,jj) * Hx_->point(ii,jj) + pmlArr_[kk].c_hxb1_0_ -> point(ii,jj) * pmlArr_[kk].Bx_->point(ii,jj) - pmlArr_[kk].c_hxb0_0_ -> point(ii,jj) * bxstore;
                                    Hy_->point(ii,jj) = pmlArr_[kk].c_hyh_0_ -> point(ii,jj) * Hy_->point(ii,jj) + pmlArr_[kk].c_hyb1_0_ -> point(ii,jj) * pmlArr_[kk].By_->point(ii,jj) - pmlArr_[kk].c_hyb0_0_ -> point(ii,jj) * bystore;

                                    bxstore = pmlArr_[kk].Bx_end_->point(ii,jj);
                                    bystore = pmlArr_[kk].By_end_->point(ii,jj);

                                    pmlArr_[kk].Bx_end_->point(ii,jj) = pmlArr_[kk].c_bxb_n_ -> point(ii,jj) * pmlArr_[kk].Bx_end_->point(ii,jj) - pmlArr_[kk].c_bxe_n_ -> point(ii,jj) * (Ez_->point(nx_-1-ii, jj+1)-Ez_->point(nx_-1-ii,jj));
                                    pmlArr_[kk].By_end_->point(ii,jj) = pmlArr_[kk].c_byb_n_ -> point(ii,jj) * pmlArr_[kk].By_end_->point(ii,jj) + pmlArr_[kk].c_bye_n_ -> point(ii,jj) * (Ez_->point(nx_-1-ii+1, jj)-Ez_->point(nx_-1-ii,jj));

                                    Hx_->point(nx_-1-ii,jj) = pmlArr_[kk].c_hxh_n_ -> point(ii,jj) * Hx_->point(nx_-1-ii,jj) + pmlArr_[kk].c_hxb1_n_ -> point(ii,jj) * pmlArr_[kk].Bx_end_->point(ii,jj) - pmlArr_[kk].c_hxb0_n_ -> point(ii,jj) * bxstore;
                                    Hy_->point(nx_-1-ii,jj) = pmlArr_[kk].c_hyh_n_ -> point(ii,jj) * Hy_->point(nx_-1-ii,jj) + pmlArr_[kk].c_hyb1_n_ -> point(ii,jj) * pmlArr_[kk].By_end_->point(ii,jj) - pmlArr_[kk].c_hyb0_n_ -> point(ii,jj) * bystore;
                                }
                            }
                            if(kk==0)
                            {
                                // Bot Left
                                bxstore = pmlArr_[0].Bx_->point(0,0);
                                bystore = pmlArr_[0].By_->point(0,0);
                                pmlArr_[0].Bx_->point(0,0) = pmlArr_[0].c_bxb_0_ -> point(0,0) * pmlArr_[0].Bx_->point(0,0) - pmlArr_[0].c_bxe_0_ -> point(0,0) * (Ez_->point(0,0+1)-Ez_->point(0,0));
                                pmlArr_[0].By_->point(0,0) = pmlArr_[0].c_byb_0_ -> point(0,0) * pmlArr_[0].By_->point(0,0) + pmlArr_[0].c_bye_0_ -> point(0,0) * (Ez_->point(0+1,0)-Ez_->point(0,0));
                                Hx_->point(0,0) = pmlArr_[0].c_hxh_0_ -> point(0,0) * Hx_->point(0,0) + pmlArr_[0].c_hxb1_0_ -> point(0,0) * pmlArr_[0].Bx_->point(0,0) - pmlArr_[0].c_hxb0_0_ -> point(0,0) * bxstore;
                                Hy_->point(0,0) = pmlArr_[0].c_hyh_0_ -> point(0,0) * Hy_->point(0,0) + pmlArr_[0].c_hyb1_0_ -> point(0,0) * pmlArr_[0].By_->point(0,0) - pmlArr_[0].c_hyb0_0_ -> point(0,0) * bystore;
                                //Bot Right
                                bxstore = pmlArr_[0].Bx_end_->point(0,0);
                                pmlArr_[0].Bx_end_->point(0,0) = pmlArr_[0].c_bxb_n_ -> point(0,0) * pmlArr_[0].Bx_end_->point(0,0) - pmlArr_[0].c_bxe_n_ -> point(0,0) * (Ez_->point(nx_-1,0+1)-Ez_->point(nx_-1,0));
                                Hx_->point(nx_-1,0) = pmlArr_[0].c_hxh_n_ -> point(0,0) * Hx_->point(nx_-1,0) + pmlArr_[0].c_hxb1_n_ -> point(0,0) * pmlArr_[0].Bx_end_->point(0,0) - pmlArr_[0].c_hxb0_n_ -> point(0,0) * bxstore;
                                //Top Left
                                bystore = pmlArr_[0].By_->point(0,ny_-1);
                                pmlArr_[0].By_->point(0,ny_-1) = pmlArr_[0].c_byb_0_ -> point(0,ny_-1) * pmlArr_[0].By_->point(0,ny_-1) + pmlArr_[0].c_bye_0_ -> point(0,ny_-1) * (Ez_->point(0+1,ny_-1)-Ez_->point(0,ny_-1));
                                Hy_->point(0,ny_-1) = pmlArr_[0].c_hyh_0_ -> point(0,ny_-1) * Hy_->point(0,ny_-1) + pmlArr_[0].c_hyb1_0_ -> point(0,ny_-1) * pmlArr_[0].By_->point(0,ny_-1) - pmlArr_[0].c_hyb0_0_ -> point(0,ny_-1) * bystore;

                                for(int ii = 1; ii < pmlArr_[0].thickness(); ii++)
                                {
                                    // Bot Left
                                    bxstore = pmlArr_[0].Bx_->point(ii,0);
                                    bystore = pmlArr_[0].By_->point(ii,0);
                                    pmlArr_[0].Bx_->point(ii,0) = pmlArr_[0].c_bxb_0_ -> point(ii,0) * pmlArr_[0].Bx_->point(ii,0) - pmlArr_[0].c_bxe_0_ -> point(ii,0) * (Ez_->point(ii,0+1)-Ez_->point(ii,0));
                                    pmlArr_[0].By_->point(ii,0) = pmlArr_[0].c_byb_0_ -> point(ii,0) * pmlArr_[0].By_->point(ii,0) + pmlArr_[0].c_bye_0_ -> point(ii,0) * (Ez_->point(ii+1,0)-Ez_->point(ii,0));
                                    Hx_->point(ii,0) = pmlArr_[0].c_hxh_0_ -> point(ii,0) * Hx_->point(ii,0) + pmlArr_[0].c_hxb1_0_ -> point(ii,0) * pmlArr_[0].Bx_->point(ii,0) - pmlArr_[0].c_hxb0_0_ -> point(ii,0) * bxstore;
                                    Hy_->point(ii,0) = pmlArr_[0].c_hyh_0_ -> point(ii,0) * Hy_->point(ii,0) + pmlArr_[0].c_hyb1_0_ -> point(ii,0) * pmlArr_[0].By_->point(ii,0) - pmlArr_[0].c_hyb0_0_ -> point(ii,0) * bystore;
                                    //Bot Right
                                    bxstore = pmlArr_[0].Bx_end_->point(ii,0);
                                    bystore = pmlArr_[0].By_end_->point(ii,0);
                                    pmlArr_[0].Bx_end_->point(ii,0) = pmlArr_[0].c_bxb_n_ -> point(ii,0) * pmlArr_[0].Bx_end_->point(ii,0) - pmlArr_[0].c_bxe_n_ -> point(ii,0) * (Ez_->point(nx_-1-ii,0+1)-Ez_->point(nx_-1-ii,0));
                                    pmlArr_[0].By_end_->point(ii,0) = pmlArr_[0].c_byb_n_ -> point(ii,0) * pmlArr_[0].By_end_->point(ii,0) + pmlArr_[0].c_bye_n_ -> point(ii,0) * (Ez_->point(nx_-1-ii+1,0)-Ez_->point(nx_-1-ii,0));
                                    Hx_->point(nx_-1-ii,0) = pmlArr_[0].c_hxh_n_ -> point(ii,0) * Hx_->point(nx_-1-ii,0) + pmlArr_[0].c_hxb1_n_ -> point(ii,0) * pmlArr_[0].Bx_end_->point(ii,0) - pmlArr_[0].c_hxb0_n_ -> point(ii,0) * bxstore;
                                    Hy_->point(nx_-1-ii,0) = pmlArr_[0].c_hyh_n_ -> point(ii,0) * Hy_->point(nx_-1-ii,0) + pmlArr_[0].c_hyb1_n_ -> point(ii,0) * pmlArr_[0].By_end_->point(ii,0) - pmlArr_[0].c_hyb0_n_ -> point(ii,0) * bystore;
                                    //Top Right
                                    bystore = pmlArr_[0].By_end_->point(ii,ny_-1);
                                    pmlArr_[0].By_end_->point(ii,ny_-1) = pmlArr_[0].c_byb_n_ -> point(ii,ny_-1) * pmlArr_[0].By_end_->point(ii,ny_-1) + pmlArr_[0].c_bye_n_ -> point(ii,ny_-1) * (Ez_->point(nx_-1-ii+1,ny_-1)-Ez_->point(nx_-1-ii,ny_-1));
                                    Hy_->point(nx_-1-ii,ny_-1) = pmlArr_[0].c_hyh_n_ -> point(ii,ny_-1) * Hy_->point(nx_-1-ii,ny_-1) + pmlArr_[0].c_hyb1_n_ -> point(ii,ny_-1) * pmlArr_[0].By_end_->point(ii,ny_-1) - pmlArr_[0].c_hyb0_n_ -> point(ii,ny_-1) * bystore;
                                    //Top Left
                                    bystore = pmlArr_[0].By_->point(ii,ny_-1);
                                    pmlArr_[0].By_->point(ii,ny_-1) = pmlArr_[0].c_byb_0_ -> point(ii,ny_-1) * pmlArr_[0].By_->point(ii,ny_-1) + pmlArr_[0].c_bye_0_ -> point(ii,ny_-1) * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
                                    Hy_->point(ii,ny_-1) = pmlArr_[0].c_hyh_0_ -> point(ii,ny_-1) * Hy_->point(ii,ny_-1) + pmlArr_[0].c_hyb1_0_ -> point(ii,ny_-1) * pmlArr_[0].By_->point(ii,ny_ - 1) - pmlArr_[0].c_hyb0_0_ -> point(ii,ny_-1) * bystore;
                                }
                                for(int jj = 1; jj < pmlArr_[0].thickness(); jj++)
                                {
                                    // Bot Left
                                    bxstore = pmlArr_[0].Bx_->point(0,jj);
                                    bystore = pmlArr_[0].By_->point(0,jj);
                                    pmlArr_[0].Bx_->point(0,jj) = pmlArr_[0].c_bxb_0_ -> point(0,jj) * pmlArr_[0].Bx_->point(0,jj) - pmlArr_[0].c_bxe_0_ -> point(0,jj) * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                                    pmlArr_[0].By_->point(0,jj) = pmlArr_[0].c_byb_0_ -> point(0,jj) * pmlArr_[0].By_->point(0,jj) + pmlArr_[0].c_bye_0_ -> point(0,jj) * (Ez_->point(0+1,jj)-Ez_->point(0,jj));
                                    Hx_->point(0,jj) = pmlArr_[0].c_hxh_0_ -> point(0,jj) * Hx_->point(0,jj) + pmlArr_[0].c_hxb1_0_ -> point(0,jj) * pmlArr_[0].Bx_->point(0,jj) - pmlArr_[0].c_hxb0_0_ -> point(0,jj) * bxstore;
                                    Hy_->point(0,jj) = pmlArr_[0].c_hyh_0_ -> point(0,jj) * Hy_->point(0,jj) + pmlArr_[0].c_hyb1_0_ -> point(0,jj) * pmlArr_[0].By_->point(0,jj) - pmlArr_[0].c_hyb0_0_ -> point(0,jj) * bystore;
                                    //Bot Right
                                    bxstore = pmlArr_[0].Bx_end_->point(0,jj);
                                    pmlArr_[0].Bx_end_->point(0,jj) = pmlArr_[0].c_bxb_n_ -> point(0,jj) * pmlArr_[0].Bx_end_->point(0,jj) - pmlArr_[0].c_bxe_n_ -> point(0,jj) * (Ez_->point(nx_-1,jj+1)-Ez_->point(nx_-1,jj));
                                    Hx_->point(nx_-1,jj) = pmlArr_[0].c_hxh_n_ -> point(0,jj) * Hx_->point(nx_-1,jj) + pmlArr_[0].c_hxb1_n_ -> point(0,jj) * pmlArr_[0].Bx_end_->point(0,jj) - pmlArr_[0].c_hxb0_n_ -> point(0,jj) * bxstore;
                                    //Top Rigt
                                    bxstore = pmlArr_[0].Bx_end_->point(0,ny_-1-jj);
                                    pmlArr_[0].Bx_end_->point(0,ny_-1-jj) = pmlArr_[0].c_bxb_n_ -> point(0,ny_-1-jj) * pmlArr_[0].Bx_end_->point(0,ny_-1-jj) - pmlArr_[0].c_bxe_n_ -> point(0,ny_-1-jj) * (Ez_->point(nx_-1,ny_-1-jj+1)-Ez_->point(nx_-1,ny_-1-jj));
                                    Hx_->point(nx_-1,ny_-1-jj) = pmlArr_[0].c_hxh_n_ -> point(0,ny_-1-jj) * Hx_->point(nx_-1,ny_-1-jj) + pmlArr_[0].c_hxb1_n_ -> point(0,ny_-1-jj) * pmlArr_[0].Bx_end_->point(0,ny_-1-jj) - pmlArr_[0].c_hxb0_n_ -> point(0,ny_-1-jj) * bxstore;
                                    //Top Left
                                    bxstore = pmlArr_[0].Bx_->point(0,ny_-1-jj);
                                    bystore = pmlArr_[0].By_->point(0,ny_-1-jj);
                                    pmlArr_[0].Bx_->point(0,ny_-1-jj) = pmlArr_[0].c_bxb_0_ -> point(0,ny_-1-jj) * pmlArr_[0].Bx_->point(0,ny_-1-jj) - pmlArr_[0].c_bxe_0_ -> point(0,ny_-1-jj) * (Ez_->point(0,ny_-1-jj+1)-Ez_->point(0,ny_-1-jj));
                                    pmlArr_[0].By_->point(0,ny_-1-jj) = pmlArr_[0].c_byb_0_ -> point(0,ny_-1-jj) * pmlArr_[0].By_->point(0,ny_-1-jj) + pmlArr_[0].c_bye_0_ -> point(0,ny_-1-jj) * (Ez_->point(0+1,ny_-1-jj)-Ez_->point(0,ny_-1-jj));
                                    Hx_->point(0,ny_-1-jj) = pmlArr_[0].c_hxh_0_ -> point(0,ny_-1-jj) * Hx_->point(0,ny_-1-jj) + pmlArr_[0].c_hxb1_0_ -> point(0,ny_-1-jj) * pmlArr_[0].Bx_->point(0,ny_-1-jj) - pmlArr_[0].c_hxb0_0_ -> point(0,ny_-1-jj) * bxstore;
                                    Hy_->point(0,ny_-1-jj) = pmlArr_[0].c_hyh_0_ -> point(0,ny_-1-jj) * Hy_->point(0,ny_-1-jj) + pmlArr_[0].c_hyb1_0_ -> point(0,ny_-1-jj) * pmlArr_[0].By_->point(0,ny_-1-jj) - pmlArr_[0].c_hyb0_0_ -> point(0,ny_-1-jj) * bystore;
                                }
                                for(int ii = 1; ii < pmlArr_[0].thickness(); ii++)
                                {
                                    for(int jj = 1; jj < pmlArr_[0].thickness(); jj++)
                                    {
                                        // Bot Left
                                        bxstore = pmlArr_[0].Bx_->point(ii,jj);
                                        bystore = pmlArr_[0].By_->point(ii,jj);
                                        pmlArr_[0].Bx_->point(ii,jj) = pmlArr_[0].c_bxb_0_ -> point(ii,jj) * pmlArr_[0].Bx_->point(ii,jj) - pmlArr_[0].c_bxe_0_ -> point(ii,jj) * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                        pmlArr_[0].By_->point(ii,jj) = pmlArr_[0].c_byb_0_ -> point(ii,jj) * pmlArr_[0].By_->point(ii,jj) + pmlArr_[0].c_bye_0_ -> point(ii,jj) * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                                        Hx_->point(ii,jj) = pmlArr_[0].c_hxh_0_ -> point(ii,jj) * Hx_->point(ii,jj) + pmlArr_[0].c_hxb1_0_ -> point(ii,jj) * pmlArr_[0].Bx_->point(ii,jj) - pmlArr_[0].c_hxb0_0_ -> point(ii,jj) * bxstore;
                                        Hy_->point(ii,jj) = pmlArr_[0].c_hyh_0_ -> point(ii,jj) * Hy_->point(ii,jj) + pmlArr_[0].c_hyb1_0_ -> point(ii,jj) * pmlArr_[0].By_->point(ii,jj) - pmlArr_[0].c_hyb0_0_ -> point(ii,jj) * bystore;
                                        //Bot Right
                                        bxstore = pmlArr_[0].Bx_end_->point(ii,jj);
                                        bystore = pmlArr_[0].By_end_->point(ii,jj);
                                        pmlArr_[0].Bx_end_->point(ii,jj) = pmlArr_[0].c_bxb_n_ -> point(ii,jj) * pmlArr_[0].Bx_end_->point(ii,jj) - pmlArr_[0].c_bxe_n_ -> point(ii,jj) * (Ez_->point(nx_-1-ii,jj+1)-Ez_->point(nx_-1-ii,jj));
                                        pmlArr_[0].By_end_->point(ii,jj) = pmlArr_[0].c_byb_n_ -> point(ii,jj) * pmlArr_[0].By_end_->point(ii,jj) + pmlArr_[0].c_bye_n_ -> point(ii,jj) * (Ez_->point(nx_-1-ii+1,jj)-Ez_->point(nx_-1-ii,jj));
                                        Hx_->point(nx_-1-ii,jj) = pmlArr_[0].c_hxh_n_ -> point(ii,jj) * Hx_->point(nx_-1-ii,jj) + pmlArr_[0].c_hxb1_n_ -> point(ii,jj) * pmlArr_[0].Bx_end_->point(ii,jj) - pmlArr_[0].c_hxb0_n_ -> point(ii,jj) * bxstore;
                                        Hy_->point(nx_-1-ii,jj) = pmlArr_[0].c_hyh_n_ -> point(ii,jj) * Hy_->point(nx_-1-ii,jj) + pmlArr_[0].c_hyb1_n_ -> point(ii,jj) * pmlArr_[0].By_end_->point(ii,jj) - pmlArr_[0].c_hyb0_n_ -> point(ii,jj) * bystore;
                                        //Top Right
                                        bxstore = pmlArr_[0].Bx_end_->point(ii,ny_-1-jj);
                                        bystore = pmlArr_[0].By_end_->point(ii,ny_-1-jj);
                                        pmlArr_[0].Bx_end_->point(ii,ny_-1-jj) = pmlArr_[0].c_bxb_n_ -> point(ii,ny_-1-jj) * pmlArr_[0].Bx_end_->point(ii,ny_-1-jj) - pmlArr_[0].c_bxe_n_ -> point(ii,ny_-1-jj) * (Ez_->point(nx_-1-ii,ny_-1-jj+1)-Ez_->point(nx_-1-ii,ny_-1-jj));
                                        pmlArr_[0].By_end_->point(ii,ny_-1-jj) = pmlArr_[0].c_byb_n_ -> point(ii,ny_-1-jj) * pmlArr_[0].By_end_->point(ii,ny_-1-jj) + pmlArr_[0].c_bye_n_ -> point(ii,ny_-1-jj) * (Ez_->point(nx_-1-ii+1,ny_-1-jj)-Ez_->point(nx_-1-ii,ny_-1-jj));
                                        Hx_->point(nx_-1-ii,ny_-1-jj) = pmlArr_[0].c_hxh_n_ -> point(ii,ny_-1-jj) * Hx_->point(nx_-1-ii,ny_-1-jj) + pmlArr_[0].c_hxb1_n_ -> point(ii,ny_-1-jj) * pmlArr_[0].Bx_end_->point(ii,ny_-1-jj) - pmlArr_[0].c_hxb0_n_ -> point(ii,ny_-1-jj) * bxstore;
                                        Hy_->point(nx_-1-ii,ny_-1-jj) = pmlArr_[0].c_hyh_n_ -> point(ii,ny_-1-jj) * Hy_->point(nx_-1-ii,ny_-1-jj) + pmlArr_[0].c_hyb1_n_ -> point(ii,ny_-1-jj) * pmlArr_[0].By_end_->point(ii,ny_-1-jj) - pmlArr_[0].c_hyb0_n_ -> point(ii,ny_-1-jj) * bystore;
                                        //Top Left
                                        bxstore = pmlArr_[0].Bx_->point(ii,ny_-1-jj);
                                        bystore = pmlArr_[0].By_->point(ii,ny_-1-jj);
                                        pmlArr_[0].Bx_->point(ii,ny_-1-jj) = pmlArr_[0].c_bxb_0_ -> point(ii,ny_-1-jj) * pmlArr_[0].Bx_->point(ii,ny_-1-jj) - pmlArr_[0].c_bxe_0_ -> point(ii,ny_-1-jj) * (Ez_->point(ii,ny_-1-jj+1)-Ez_->point(ii,ny_-1-jj));
                                        pmlArr_[0].By_->point(ii,ny_-1-jj) = pmlArr_[0].c_byb_0_ -> point(ii,ny_-1-jj) * pmlArr_[0].By_->point(ii,ny_-1-jj) + pmlArr_[0].c_bye_0_ -> point(ii,ny_-1-jj) * (Ez_->point(ii+1,ny_-1-jj)-Ez_->point(ii,ny_-1-jj));
                                        Hx_->point(ii,ny_-1-jj) = pmlArr_[0].c_hxh_0_ -> point(ii,ny_-1-jj) * Hx_->point(ii,ny_-1-jj) + pmlArr_[0].c_hxb1_0_ -> point(ii,ny_-1-jj) * pmlArr_[0].Bx_->point(ii,ny_-1-jj) - pmlArr_[0].c_hxb0_0_ -> point(ii,ny_-1-jj) * bxstore;
                                        Hy_->point(ii,ny_-1-jj) = pmlArr_[0].c_hyh_0_ -> point(ii,ny_-1-jj) * Hy_->point(ii,ny_-1-jj) + pmlArr_[0].c_hyb1_0_ -> point(ii,ny_-1-jj) * pmlArr_[0].By_->point(ii,ny_-1-jj) - pmlArr_[0].c_hyb0_0_ -> point(ii,ny_-1-jj) * bystore;
                                    }
                                }
                            }
                        }
                        break;
                    }
                    case Y:
                    {
                        complex<double>bxstore(0.0,0.0); complex<double>bystore(0.0,0.0);
                        if (xPML_== 0)
                        {
                            for (int jj = 0; jj < nx_-1; jj++)
                            {
                                //Update both on the bottom side and Hy on the top
                                bxstore = pmlArr_[kk].Bx_->point(jj,0);
                                bystore = pmlArr_[kk].By_->point(jj,0);
                                pmlArr_[kk].Bx_->point(jj,0) = pmlArr_[kk].c_bxb_0_ -> point(jj,0) * pmlArr_[kk].Bx_->point(jj,0) - pmlArr_[kk].c_bxe_0_ -> point(jj,0) * (Ez_->point(jj,0+1)-Ez_->point(jj,0));
                                pmlArr_[kk].By_->point(jj,0) = pmlArr_[kk].c_byb_0_ -> point(jj,0) * pmlArr_[kk].By_->point(jj,0) + pmlArr_[kk].c_bye_0_ -> point(jj,0) * (Ez_->point(jj+1,0)-Ez_->point(jj,0));
                                //cout <<"H1"<<endl;
                                Hx_->point(jj,0) = pmlArr_[kk].c_hxh_0_ -> point(jj,0) * Hx_->point(jj,0) + pmlArr_[kk].c_hxb1_0_ -> point(jj,0) * pmlArr_[kk].Bx_->point(jj,0) - pmlArr_[kk].c_hxb0_0_ -> point(jj,0) * bxstore;
                                Hy_->point(jj,0) = pmlArr_[kk].c_hyh_0_ -> point(jj,0) * Hy_->point(jj,0) + pmlArr_[kk].c_hyb1_0_ -> point(jj,0) * pmlArr_[kk].By_->point(jj,0) - pmlArr_[kk].c_hyb0_0_ -> point(jj,0) * bystore;

                                bystore = pmlArr_[kk].By_end_->point(jj,0);
                                pmlArr_[kk].By_end_->point(jj,0) = pmlArr_[kk].c_byb_n_ -> point(jj,0) * pmlArr_[kk].By_end_->point(jj,0) + pmlArr_[kk].c_bye_n_ -> point(jj,0) * (Ez_->point(jj+1, ny_-1)-Ez_->point(jj,ny_-1));
                                Hy_->point(jj,ny_-1) = pmlArr_[kk].c_hyh_n_ -> point(jj,0) * Hy_->point(jj,ny_-1) + pmlArr_[kk].c_hyb1_n_ -> point(jj,0) * pmlArr_[kk].By_end_->point(jj,0) - pmlArr_[kk].c_hyb0_n_ -> point(jj,0) * bystore;
                            }
                            //Top left corner
                            bxstore = pmlArr_[kk].Bx_->point(nx_-1,0);
                            pmlArr_[kk].Bx_->point(nx_-1,0) = pmlArr_[kk].c_bxb_0_ -> point(nx_-1,0) * pmlArr_[kk].Bx_->point(nx_-1,0) - pmlArr_[kk].c_bxe_0_ -> point(nx_-1,0) * (Ez_->point(nx_-1,0+1)-Ez_->point(nx_-1,0));
                            Hx_->point(nx_-1,0) = pmlArr_[kk].c_hxh_0_ -> point(nx_-1,0) * Hx_->point(nx_-1,0) + pmlArr_[kk].c_hxb1_0_ -> point(nx_-1,0) * pmlArr_[kk].Bx_->point(nx_-1,0) - pmlArr_[kk].c_hxb0_0_ -> point(nx_-1,0) * bxstore;
                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                //right col
                                bxstore = pmlArr_[kk].Bx_->point(nx_-1,ii);
                                pmlArr_[kk].Bx_->point(nx_-1,ii) = pmlArr_[kk].c_bxb_0_ -> point(nx_-1,ii) * pmlArr_[kk].Bx_->point(nx_-1,ii) - pmlArr_[kk].c_bxe_0_ -> point(nx_-1,ii) * (Ez_->point(nx_-1,ii+1)-Ez_->point(nx_-1,ii));
                                Hx_->point(nx_-1,ii) = pmlArr_[kk].c_hxh_0_ -> point(nx_-1,ii) * Hx_->point(nx_-1,ii) + pmlArr_[kk].c_hxb1_0_ -> point(nx_-1,ii) * pmlArr_[kk].Bx_->point(nx_-1,ii) - pmlArr_[kk].c_hxb0_0_ -> point(nx_-1,ii) * bxstore;
                                for(int jj = 0; jj < nx_-1; jj ++)
                                {
                                    //Update everything
                                    bxstore = pmlArr_[kk].Bx_->point(jj,ii);
                                    bystore = pmlArr_[kk].By_->point(jj,ii);
                                    pmlArr_[kk].Bx_->point(jj,ii) = pmlArr_[kk].c_bxb_0_ -> point(jj,ii) * pmlArr_[kk].Bx_->point(jj,ii) - pmlArr_[kk].c_bxe_0_ -> point(jj,ii) * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                                    pmlArr_[kk].By_->point(jj,ii) = pmlArr_[kk].c_byb_0_ -> point(jj,ii) * pmlArr_[kk].By_->point(jj,ii) + pmlArr_[kk].c_bye_0_ -> point(jj,ii) * (Ez_->point(jj+1,ii)-Ez_->point(jj,ii));
                                    //cout <<"H1"<<endl;
                                    Hx_->point(jj,ii) = pmlArr_[kk].c_hxh_0_ -> point(jj,ii) * Hx_->point(jj,ii) + pmlArr_[kk].c_hxb1_0_ -> point(jj,ii) * pmlArr_[kk].Bx_->point(jj,ii) - pmlArr_[kk].c_hxb0_0_ -> point(jj,ii) * bxstore;
                                    Hy_->point(jj,ii) = pmlArr_[kk].c_hyh_0_ -> point(jj,ii) * Hy_->point(jj,ii) + pmlArr_[kk].c_hyb1_0_ -> point(jj,ii) * pmlArr_[kk].By_->point(jj,ii) - pmlArr_[kk].c_hyb0_0_ -> point(jj,ii) * bystore;

                                    bxstore = pmlArr_[kk].Bx_end_->point(jj,ii);
                                    bystore = pmlArr_[kk].By_end_->point(jj,ii);
                                    pmlArr_[kk].Bx_end_->point(jj,ii) = pmlArr_[kk].c_bxb_n_ -> point(jj,ii) * pmlArr_[kk].Bx_end_->point(jj,ii) - pmlArr_[kk].c_bxe_n_ -> point(jj,ii) * (Ez_->point(jj,ny_-1-ii+1)-Ez_->point(jj,ny_-1-ii));
                                    pmlArr_[kk].By_end_->point(jj,ii) = pmlArr_[kk].c_byb_n_ -> point(jj,ii) * pmlArr_[kk].By_end_->point(jj,ii) + pmlArr_[kk].c_bye_n_ -> point(jj,ii) * (Ez_->point(jj+1,ny_-1-ii)-Ez_->point(jj,ny_-1-ii));
                                    //cout <<"H2"<<endl;
                                    Hx_->point(jj,ny_-1-ii) = pmlArr_[kk].c_hxh_n_ -> point(jj,ii) * Hx_->point(jj,ny_-1-ii) + pmlArr_[kk].c_hxb1_n_ -> point(jj,ii) * pmlArr_[kk].Bx_end_->point(jj,ii) - pmlArr_[kk].c_hxb0_n_ -> point(jj,ii) * bxstore;
                                    Hy_->point(jj,ny_-1-ii) = pmlArr_[kk].c_hyh_n_ -> point(jj,ii) * Hy_->point(jj,ny_-1-ii) + pmlArr_[kk].c_hyb1_n_ -> point(jj,ii) * pmlArr_[kk].By_end_->point(jj,ii) - pmlArr_[kk].c_hyb0_n_ -> point(jj,ii) * bystore;
                                }
                                // for the right only upadate Hx
                                bxstore = pmlArr_[kk].Bx_end_->point(nx_-1,ii);
                                pmlArr_[kk].Bx_end_->point(nx_-1,ii) = pmlArr_[kk].c_bxb_n_ -> point(nx_-1,ii) * pmlArr_[kk].Bx_end_->point(nx_-1,ii) - pmlArr_[kk].c_bxe_n_ -> point(nx_-1,ii) * (Ez_->point(nx_-1,ny_-1-ii+1)-Ez_->point(nx_-1,ny_-1-ii));
                                Hx_->point(nx_-1,ny_-1-ii) = pmlArr_[kk].c_hxh_n_ -> point(nx_-1,ii) * Hx_->point(nx_-1,ny_-1-ii) + pmlArr_[kk].c_hxb1_n_ -> point(nx_-1,ii) * pmlArr_[kk].Bx_end_->point(nx_-1,ii) - pmlArr_[kk].c_hxb0_n_ -> point(nx_-1,ii) * bxstore;
                            }
                        }
                        else
                        {
                            for (int jj = xPML_; jj < nx_-xPML_; jj++)
                            {
                                //Update both on the bottom side and Hy on the top
                                bxstore = pmlArr_[kk].Bx_->point(jj,0);
                                bystore = pmlArr_[kk].By_->point(jj,0);
                                pmlArr_[kk].Bx_->point(jj,0) = pmlArr_[kk].c_bxb_0_ -> point(jj,0) * pmlArr_[kk].Bx_->point(jj,0) - pmlArr_[kk].c_bxe_0_ -> point(jj,0) * (Ez_->point(jj,0+1)-Ez_->point(jj,0));
                                pmlArr_[kk].By_->point(jj,0) = pmlArr_[kk].c_byb_0_ -> point(jj,0) * pmlArr_[kk].By_->point(jj,0) + pmlArr_[kk].c_bye_0_ -> point(jj,0) * (Ez_->point(jj+1,0)-Ez_->point(jj,0));
                                //cout <<"H1"<<endl;
                                Hx_->point(jj,0) = pmlArr_[kk].c_hxh_0_ -> point(jj,0) * Hx_->point(jj,0) + pmlArr_[kk].c_hxb1_0_ -> point(jj,0) * pmlArr_[kk].Bx_->point(jj,0) - pmlArr_[kk].c_hxb0_0_ -> point(jj,0) * bxstore;
                                Hy_->point(jj,0) = pmlArr_[kk].c_hyh_0_ -> point(jj,0) * Hy_->point(jj,0) + pmlArr_[kk].c_hyb1_0_ -> point(jj,0) * pmlArr_[kk].By_->point(jj,0) - pmlArr_[kk].c_hyb0_0_ -> point(jj,0) * bystore;

                                bystore = pmlArr_[kk].By_end_->point(jj,0);
                                pmlArr_[kk].By_end_->point(jj,0) = pmlArr_[kk].c_byb_n_ -> point(jj,0) * pmlArr_[kk].By_end_->point(jj,0) + pmlArr_[kk].c_bye_n_ -> point(jj,0) * (Ez_->point(jj+1,ny_-1)-Ez_->point(jj,ny_-1));
                                Hy_->point(jj,ny_-1) = pmlArr_[kk].c_hyh_n_ -> point(jj,0) * Hy_->point(jj,ny_-1) + pmlArr_[kk].c_hyb1_n_ -> point(jj,0) * pmlArr_[kk].By_end_->point(jj,0) - pmlArr_[kk].c_hyb0_n_ -> point(jj,0) * bystore;
                            }
                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                for(int jj = xPML_; jj < nx_-xPML_; jj ++)
                                {
                                    //Update everything
                                    bxstore = pmlArr_[kk].Bx_->point(jj,ii);
                                    bystore = pmlArr_[kk].By_->point(jj,ii);
                                    pmlArr_[kk].Bx_->point(jj,ii) = pmlArr_[kk].c_bxb_0_ -> point(jj,ii) * pmlArr_[kk].Bx_->point(jj,ii) - pmlArr_[kk].c_bxe_0_ -> point(jj,ii) * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                                    pmlArr_[kk].By_->point(jj,ii) = pmlArr_[kk].c_byb_0_ -> point(jj,ii) * pmlArr_[kk].By_->point(jj,ii) + pmlArr_[kk].c_bye_0_ -> point(jj,ii) * (Ez_->point(jj+1,ii)-Ez_->point(jj,ii));
                                    //cout <<"H1"<<endl;
                                    Hx_->point(jj,ii) = pmlArr_[kk].c_hxh_0_ -> point(jj,ii) * Hx_->point(jj,ii) + pmlArr_[kk].c_hxb1_0_ -> point(jj,ii) * pmlArr_[kk].Bx_->point(jj,ii) - pmlArr_[kk].c_hxb0_0_ -> point(jj,ii) * bxstore;
                                    Hy_->point(jj,ii) = pmlArr_[kk].c_hyh_0_ -> point(jj,ii) * Hy_->point(jj,ii) + pmlArr_[kk].c_hyb1_0_ -> point(jj,ii) * pmlArr_[kk].By_->point(jj,ii) - pmlArr_[kk].c_hyb0_0_ -> point(jj,ii) * bystore;

                                    bxstore = pmlArr_[kk].Bx_end_->point(jj,ii);
                                    bystore = pmlArr_[kk].By_end_->point(jj,ii);
                                    pmlArr_[kk].Bx_end_->point(jj,ii) = pmlArr_[kk].c_bxb_n_ -> point(jj,ii) * pmlArr_[kk].Bx_end_->point(jj,ii) - pmlArr_[kk].c_bxe_n_ -> point(jj,ii) * (Ez_->point(jj,ny_-1-ii+1)-Ez_->point(jj,ny_-1-ii));
                                    pmlArr_[kk].By_end_->point(jj,ii) = pmlArr_[kk].c_byb_n_ -> point(jj,ii) * pmlArr_[kk].By_end_->point(jj,ii) + pmlArr_[kk].c_bye_n_ -> point(jj,ii) * (Ez_->point(jj+1,ny_-1-ii)-Ez_->point(jj,ny_-1-ii));
                                    //cout <<"H2"<<endl;
                                    Hx_->point(jj,ny_-1-ii) = pmlArr_[kk].c_hxh_n_ -> point(jj,ii) * Hx_->point(jj,ny_-1-ii) + pmlArr_[kk].c_hxb1_n_ -> point(jj,ii) * pmlArr_[kk].Bx_end_->point(jj,ii) - pmlArr_[kk].c_hxb0_n_ -> point(jj,ii) * bxstore;
                                    Hy_->point(jj,ny_-1-ii) = pmlArr_[kk].c_hyh_n_ -> point(jj,ii) * Hy_->point(jj,ny_-1-ii) + pmlArr_[kk].c_hyb1_n_ -> point(jj,ii) * pmlArr_[kk].By_end_->point(jj,ii) - pmlArr_[kk].c_hyb0_n_ -> point(jj,ii) * bystore;
                                }
                            }
                            if(kk==0)
                            {
                                // Bot Left
                                bxstore = pmlArr_[1].Bx_->point(0,0);
                                bystore = pmlArr_[1].By_->point(0,0);
                                pmlArr_[1].Bx_->point(0,0) = pmlArr_[1].c_bxb_0_ -> point(0,0) * pmlArr_[1].Bx_->point(0,0) - pmlArr_[1].c_bxe_0_ -> point(0,0) * (Ez_->point(0,0+1)-Ez_->point(0,0));
                                pmlArr_[1].By_->point(0,0) = pmlArr_[1].c_byb_0_ -> point(0,0) * pmlArr_[1].By_->point(0,0) + pmlArr_[1].c_bye_0_ -> point(0,0) * (Ez_->point(0+1,0)-Ez_->point(0,0));
                                Hx_->point(0,0) = pmlArr_[1].c_hxh_0_ -> point(0,0) * Hx_->point(0,0) + pmlArr_[1].c_hxb1_0_ -> point(0,0) * pmlArr_[1].Bx_->point(0,0) - pmlArr_[1].c_hxb0_0_ -> point(0,0) * bxstore;
                                Hy_->point(0,0) = pmlArr_[1].c_hyh_0_ -> point(0,0) * Hy_->point(0,0) + pmlArr_[1].c_hyb1_0_ -> point(0,0) * pmlArr_[1].By_->point(0,0) - pmlArr_[1].c_hyb0_0_ -> point(0,0) * bystore;
                                //Bot Right
                                bxstore = pmlArr_[1].Bx_end_->point(0,0);
                                pmlArr_[1].Bx_end_->point(0,0) = pmlArr_[1].c_bxb_n_ -> point(0,0) * pmlArr_[1].Bx_end_->point(0,0) - pmlArr_[1].c_bxe_n_ -> point(0,0) * (Ez_->point(nx_-1,0+1)-Ez_->point(nx_-1,0));
                                Hx_->point(nx_-1,0) = pmlArr_[1].c_hxh_n_ -> point(0,0) * Hx_->point(nx_-1,0) + pmlArr_[1].c_hxb1_n_ -> point(0,0) * pmlArr_[1].Bx_end_->point(0,0) - pmlArr_[1].c_hxb0_n_ -> point(0,0) * bxstore;
                                //Top Left
                                bystore = pmlArr_[1].By_->point(0,ny_-1);
                                pmlArr_[1].By_->point(0,ny_-1) = pmlArr_[1].c_byb_0_ -> point(0,ny_-1) * pmlArr_[1].By_->point(0,ny_-1) + pmlArr_[1].c_bye_0_ -> point(0,ny_-1) * (Ez_->point(0+1,ny_-1)-Ez_->point(0,ny_-1));
                                Hy_->point(0,ny_-1) = pmlArr_[1].c_hyh_0_ -> point(0,ny_-1) * Hy_->point(0,ny_-1) + pmlArr_[1].c_hyb1_0_ -> point(0,ny_-1) * pmlArr_[1].By_->point(0,ny_-1) - pmlArr_[1].c_hyb0_0_ -> point(0,ny_-1) * bystore;

                                for(int ii = 1; ii < pmlArr_[1].thickness(); ii++)
                                {
                                    // Bot Left
                                    bxstore = pmlArr_[1].Bx_->point(ii,0);
                                    bystore = pmlArr_[1].By_->point(ii,0);
                                    pmlArr_[1].Bx_->point(ii,0) = pmlArr_[1].c_bxb_0_ -> point(ii,0) * pmlArr_[1].Bx_->point(ii,0) - pmlArr_[1].c_bxe_0_ -> point(ii,0) * (Ez_->point(ii,0+1)-Ez_->point(ii,0));
                                    pmlArr_[1].By_->point(ii,0) = pmlArr_[1].c_byb_0_ -> point(ii,0) * pmlArr_[1].By_->point(ii,0) + pmlArr_[1].c_bye_0_ -> point(ii,0) * (Ez_->point(ii+1,0)-Ez_->point(ii,0));
                                    Hx_->point(ii,0) = pmlArr_[1].c_hxh_0_ -> point(ii,0) * Hx_->point(ii,0) + pmlArr_[1].c_hxb1_0_ -> point(ii,0) * pmlArr_[1].Bx_->point(ii,0) - pmlArr_[1].c_hxb0_0_ -> point(ii,0) * bxstore;
                                    Hy_->point(ii,0) = pmlArr_[1].c_hyh_0_ -> point(ii,0) * Hy_->point(ii,0) + pmlArr_[1].c_hyb1_0_ -> point(ii,0) * pmlArr_[1].By_->point(ii,0) - pmlArr_[1].c_hyb0_0_ -> point(ii,0) * bystore;
                                    //Bot Right
                                    bxstore = pmlArr_[1].Bx_end_->point(ii,0);
                                    bystore = pmlArr_[1].By_end_->point(ii,0);
                                    pmlArr_[1].Bx_end_->point(ii,0) = pmlArr_[1].c_bxb_n_ -> point(ii,0) * pmlArr_[1].Bx_end_->point(ii,0) - pmlArr_[1].c_bxe_n_ -> point(ii,0) * (Ez_->point(nx_-1-ii,0+1)-Ez_->point(nx_-1-ii,0));
                                    pmlArr_[1].By_end_->point(ii,0) = pmlArr_[1].c_byb_n_ -> point(ii,0) * pmlArr_[1].By_end_->point(ii,0) + pmlArr_[1].c_bye_n_ -> point(ii,0) * (Ez_->point(nx_-1-ii+1,0)-Ez_->point(nx_-1-ii,0));
                                    Hx_->point(nx_-1-ii,0) = pmlArr_[1].c_hxh_n_ -> point(ii,0) * Hx_->point(nx_-1-ii,0) + pmlArr_[1].c_hxb1_n_ -> point(ii,0) * pmlArr_[1].Bx_end_->point(ii,0) - pmlArr_[1].c_hxb0_n_ -> point(ii,0) * bxstore;
                                    Hy_->point(nx_-1-ii,0) = pmlArr_[1].c_hyh_n_ -> point(ii,0) * Hy_->point(nx_-1-ii,0) + pmlArr_[1].c_hyb1_n_ -> point(ii,0) * pmlArr_[1].By_end_->point(ii,0) - pmlArr_[1].c_hyb0_n_ -> point(ii,0) * bystore;
                                    //Top Right
                                    bystore = pmlArr_[1].By_end_->point(ii,ny_-1);
                                    pmlArr_[1].By_end_->point(ii,ny_-1) = pmlArr_[1].c_byb_n_ -> point(ii,ny_-1) * pmlArr_[1].By_end_->point(ii,ny_-1) + pmlArr_[1].c_bye_n_ -> point(ii,ny_-1) * (Ez_->point(nx_-1-ii+1,ny_-1)-Ez_->point(nx_-1-ii,ny_-1));
                                    Hy_->point(nx_-1-ii,ny_-1) = pmlArr_[1].c_hyh_n_ -> point(ii,ny_-1) * Hy_->point(nx_-1-ii,ny_-1) + pmlArr_[1].c_hyb1_n_ -> point(ii,ny_-1) * pmlArr_[1].By_end_->point(ii,ny_-1) - pmlArr_[1].c_hyb0_n_ -> point(ii,ny_-1) * bystore;
                                    //Top Left
                                    bystore = pmlArr_[1].By_->point(ii,ny_-1);
                                    pmlArr_[1].By_->point(ii,ny_-1) = pmlArr_[1].c_byb_0_ -> point(ii,ny_-1) * pmlArr_[1].By_->point(ii,ny_-1) + pmlArr_[1].c_bye_0_ -> point(ii,ny_-1) * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
                                    Hy_->point(ii,ny_-1) = pmlArr_[1].c_hyh_0_ -> point(ii,ny_-1) * Hy_->point(ii,ny_-1) + pmlArr_[1].c_hyb1_0_ -> point(ii,ny_-1) * pmlArr_[1].By_->point(ii,ny_ - 1) - pmlArr_[1].c_hyb0_0_ -> point(ii,ny_-1) * bystore;
                                }
                                for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                                {
                                    // Bot Left
                                    bxstore = pmlArr_[1].Bx_->point(0,jj);
                                    bystore = pmlArr_[1].By_->point(0,jj);
                                    pmlArr_[1].Bx_->point(0,jj) = pmlArr_[1].c_bxb_0_ -> point(0,jj) * pmlArr_[1].Bx_->point(0,jj) - pmlArr_[1].c_bxe_0_ -> point(0,jj) * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                                    pmlArr_[1].By_->point(0,jj) = pmlArr_[1].c_byb_0_ -> point(0,jj) * pmlArr_[1].By_->point(0,jj) + pmlArr_[1].c_bye_0_ -> point(0,jj) * (Ez_->point(0+1,jj)-Ez_->point(0,jj));
                                    Hx_->point(0,jj) = pmlArr_[1].c_hxh_0_ -> point(0,jj) * Hx_->point(0,jj) + pmlArr_[1].c_hxb1_0_ -> point(0,jj) * pmlArr_[1].Bx_->point(0,jj) - pmlArr_[1].c_hxb0_0_ -> point(0,jj) * bxstore;
                                    Hy_->point(0,jj) = pmlArr_[1].c_hyh_0_ -> point(0,jj) * Hy_->point(0,jj) + pmlArr_[1].c_hyb1_0_ -> point(0,jj) * pmlArr_[1].By_->point(0,jj) - pmlArr_[1].c_hyb0_0_ -> point(0,jj) * bystore;
                                    //Bot Right
                                    bxstore = pmlArr_[1].Bx_end_->point(0,jj);
                                    pmlArr_[1].Bx_end_->point(0,jj) = pmlArr_[1].c_bxb_n_ -> point(0,jj) * pmlArr_[1].Bx_end_->point(0,jj) - pmlArr_[1].c_bxe_n_ -> point(0,jj) * (Ez_->point(nx_-1,jj+1)-Ez_->point(nx_-1,jj));
                                    Hx_->point(nx_-1,jj) = pmlArr_[1].c_hxh_n_ -> point(0,jj) * Hx_->point(nx_-1,jj) + pmlArr_[1].c_hxb1_n_ -> point(0,jj) * pmlArr_[1].Bx_end_->point(0,jj) - pmlArr_[1].c_hxb0_n_ -> point(0,jj) * bxstore;
                                    //Top Rigt
                                    bxstore = pmlArr_[1].Bx_end_->point(0,ny_-1-jj);
                                    pmlArr_[1].Bx_end_->point(0,ny_-1-jj) = pmlArr_[1].c_bxb_n_ -> point(0,ny_-1-jj) * pmlArr_[1].Bx_end_->point(0,ny_-1-jj) - pmlArr_[1].c_bxe_n_ -> point(0,ny_-1-jj) * (Ez_->point(nx_-1,ny_-1-jj+1)-Ez_->point(nx_-1,ny_-1-jj));
                                    Hx_->point(nx_-1,ny_-1-jj) = pmlArr_[1].c_hxh_n_ -> point(0,ny_-1-jj) * Hx_->point(nx_-1,ny_-1-jj) + pmlArr_[1].c_hxb1_n_ -> point(0,ny_-1-jj) * pmlArr_[1].Bx_end_->point(0,ny_-1-jj) - pmlArr_[1].c_hxb0_n_ -> point(0,ny_-1-jj) * bxstore;
                                    //Top Left
                                    bxstore = pmlArr_[1].Bx_->point(0,ny_-1-jj);
                                    bystore = pmlArr_[1].By_->point(0,ny_-1-jj);
                                    pmlArr_[1].Bx_->point(0,ny_-1-jj) = pmlArr_[1].c_bxb_0_ -> point(0,ny_-1-jj) * pmlArr_[1].Bx_->point(0,ny_-1-jj) - pmlArr_[1].c_bxe_0_ -> point(0,ny_-1-jj) * (Ez_->point(0,ny_-1-jj+1)-Ez_->point(0,ny_-1-jj));
                                    pmlArr_[1].By_->point(0,ny_-1-jj) = pmlArr_[1].c_byb_0_ -> point(0,ny_-1-jj) * pmlArr_[1].By_->point(0,ny_-1-jj) + pmlArr_[1].c_bye_0_ -> point(0,ny_-1-jj) * (Ez_->point(0+1,ny_-1-jj)-Ez_->point(0,ny_-1-jj));
                                    Hx_->point(0,ny_-1-jj) = pmlArr_[1].c_hxh_0_ -> point(0,ny_-1-jj) * Hx_->point(0,ny_-1-jj) + pmlArr_[1].c_hxb1_0_ -> point(0,ny_-1-jj) * pmlArr_[1].Bx_->point(0,ny_-1-jj) - pmlArr_[1].c_hxb0_0_ -> point(0,ny_-1-jj) * bxstore;
                                    Hy_->point(0,ny_-1-jj) = pmlArr_[1].c_hyh_0_ -> point(0,ny_-1-jj) * Hy_->point(0,ny_-1-jj) + pmlArr_[1].c_hyb1_0_ -> point(0,ny_-1-jj) * pmlArr_[1].By_->point(0,ny_-1-jj) - pmlArr_[1].c_hyb0_0_ -> point(0,ny_-1-jj) * bystore;
                                }
                                for(int ii = 1; ii < pmlArr_[1].thickness(); ii++)
                                {
                                    for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                                    {
                                        // Bot Left
                                        bxstore = pmlArr_[1].Bx_->point(ii,jj);
                                        bystore = pmlArr_[1].By_->point(ii,jj);
                                        pmlArr_[1].Bx_->point(ii,jj) = pmlArr_[1].c_bxb_0_ -> point(ii,jj) * pmlArr_[1].Bx_->point(ii,jj) - pmlArr_[1].c_bxe_0_ -> point(ii,jj) * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                        pmlArr_[1].By_->point(ii,jj) = pmlArr_[1].c_byb_0_ -> point(ii,jj) * pmlArr_[1].By_->point(ii,jj) + pmlArr_[1].c_bye_0_ -> point(ii,jj) * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                                        Hx_->point(ii,jj) = pmlArr_[1].c_hxh_0_ -> point(ii,jj) * Hx_->point(ii,jj) + pmlArr_[1].c_hxb1_0_ -> point(ii,jj) * pmlArr_[1].Bx_->point(ii,jj) - pmlArr_[1].c_hxb0_0_ -> point(ii,jj) * bxstore;
                                        Hy_->point(ii,jj) = pmlArr_[1].c_hyh_0_ -> point(ii,jj) * Hy_->point(ii,jj) + pmlArr_[1].c_hyb1_0_ -> point(ii,jj) * pmlArr_[1].By_->point(ii,jj) - pmlArr_[1].c_hyb0_0_ -> point(ii,jj) * bystore;
                                        //Bot Right
                                        bxstore = pmlArr_[1].Bx_end_->point(ii,jj);
                                        bystore = pmlArr_[1].By_end_->point(ii,jj);
                                        pmlArr_[1].Bx_end_->point(ii,jj) = pmlArr_[1].c_bxb_n_ -> point(ii,jj) * pmlArr_[1].Bx_end_->point(ii,jj) - pmlArr_[1].c_bxe_n_ -> point(ii,jj) * (Ez_->point(nx_-1-ii,jj+1)-Ez_->point(nx_-1-ii,jj));
                                        pmlArr_[1].By_end_->point(ii,jj) = pmlArr_[1].c_byb_n_ -> point(ii,jj) * pmlArr_[1].By_end_->point(ii,jj) + pmlArr_[1].c_bye_n_ -> point(ii,jj) * (Ez_->point(nx_-1-ii+1,jj)-Ez_->point(nx_-1-ii,jj));
                                        Hx_->point(nx_-1-ii,jj) = pmlArr_[1].c_hxh_n_ -> point(ii,jj) * Hx_->point(nx_-1-ii,jj) + pmlArr_[1].c_hxb1_n_ -> point(ii,jj) * pmlArr_[1].Bx_end_->point(ii,jj) - pmlArr_[1].c_hxb0_n_ -> point(ii,jj) * bxstore;
                                        Hy_->point(nx_-1-ii,jj) = pmlArr_[1].c_hyh_n_ -> point(ii,jj) * Hy_->point(nx_-1-ii,jj) + pmlArr_[1].c_hyb1_n_ -> point(ii,jj) * pmlArr_[1].By_end_->point(ii,jj) - pmlArr_[1].c_hyb0_n_ -> point(ii,jj) * bystore;
                                        //Top Right
                                        bxstore = pmlArr_[1].Bx_end_->point(ii,ny_-1-jj);
                                        bystore = pmlArr_[1].By_end_->point(ii,ny_-1-jj);
                                        pmlArr_[1].Bx_end_->point(ii,ny_-1-jj) = pmlArr_[1].c_bxb_n_ -> point(ii,ny_-1-jj) * pmlArr_[1].Bx_end_->point(ii,ny_-1-jj) - pmlArr_[1].c_bxe_n_ -> point(ii,ny_-1-jj) * (Ez_->point(nx_-1-ii,ny_-1-jj+1)-Ez_->point(nx_-1-ii,ny_-1-jj));
                                        pmlArr_[1].By_end_->point(ii,ny_-1-jj) = pmlArr_[1].c_byb_n_ -> point(ii,ny_-1-jj) * pmlArr_[1].By_end_->point(ii,ny_-1-jj) + pmlArr_[1].c_bye_n_ -> point(ii,ny_-1-jj) * (Ez_->point(nx_-1-ii+1,ny_-1-jj)-Ez_->point(nx_-1-ii,ny_-1-jj));
                                        Hx_->point(nx_-1-ii,ny_-1-jj) = pmlArr_[1].c_hxh_n_ -> point(ii,ny_-1-jj) * Hx_->point(nx_-1-ii,ny_-1-jj) + pmlArr_[1].c_hxb1_n_ -> point(ii,ny_-1-jj) * pmlArr_[1].Bx_end_->point(ii,ny_-1-jj) - pmlArr_[1].c_hxb0_n_ -> point(ii,ny_-1-jj) * bxstore;
                                        Hy_->point(nx_-1-ii,ny_-1-jj) = pmlArr_[1].c_hyh_n_ -> point(ii,ny_-1-jj) * Hy_->point(nx_-1-ii,ny_-1-jj) + pmlArr_[1].c_hyb1_n_ -> point(ii,ny_-1-jj) * pmlArr_[1].By_end_->point(ii,ny_-1-jj) - pmlArr_[1].c_hyb0_n_ -> point(ii,ny_-1-jj) * bystore;
                                        //Top Left
                                        bxstore = pmlArr_[1].Bx_->point(ii,ny_-1-jj);
                                        bystore = pmlArr_[1].By_->point(ii,ny_-1-jj);
                                        pmlArr_[1].Bx_->point(ii,ny_-1-jj) = pmlArr_[1].c_bxb_0_ -> point(ii,ny_-1-jj) * pmlArr_[1].Bx_->point(ii,ny_-1-jj) - pmlArr_[1].c_bxe_0_ -> point(ii,ny_-1-jj) * (Ez_->point(ii,ny_-1-jj+1)-Ez_->point(ii,ny_-1-jj));
                                        pmlArr_[1].By_->point(ii,ny_-1-jj) = pmlArr_[1].c_byb_0_ -> point(ii,ny_-1-jj) * pmlArr_[1].By_->point(ii,ny_-1-jj) + pmlArr_[1].c_bye_0_ -> point(ii,ny_-1-jj) * (Ez_->point(ii+1,ny_-1-jj)-Ez_->point(ii,ny_-1-jj));
                                        Hx_->point(ii,ny_-1-jj) = pmlArr_[1].c_hxh_0_ -> point(ii,ny_-1-jj) * Hx_->point(ii,ny_-1-jj) + pmlArr_[1].c_hxb1_0_ -> point(ii,ny_-1-jj) * pmlArr_[1].Bx_->point(ii,ny_-1-jj) - pmlArr_[1].c_hxb0_0_ -> point(ii,ny_-1-jj) * bxstore;
                                        Hy_->point(ii,ny_-1-jj) = pmlArr_[1].c_hyh_0_ -> point(ii,ny_-1-jj) * Hy_->point(ii,ny_-1-jj) + pmlArr_[1].c_hyb1_0_ -> point(ii,ny_-1-jj) * pmlArr_[1].By_->point(ii,ny_-1-jj) - pmlArr_[1].c_hyb0_0_ -> point(ii,ny_-1-jj) * bystore;
                                    }
                                }
                            }
                        }
                        break;
                    }
                    case Z:
                    {
                        throw logic_error("Z dir not implimented");
                        break;
                    }
                    default:
                        throw logic_error("how did you get here, it is a switch default try again");
                        break;
                }
            }
        }
    }*/
    if(Ez_)
    {
        for(int kk =0; kk < pmlArr_.size(); kk++)
        {
            int stride = 1; int stride_rel = 1;
            int ni = 0; int nj = 0; int d = 0;
            if(pmlArr_[kk].d() == X)
            {
                stride_rel = nx_;
                ni         = nx_ - 1;
                stride     = pmlArr_[kk].thickness();
            }
            else
            {
                d          = 1;
                nj         = ny_ - 1;
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxHx_.size();zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxHx_[zz];
                int xx = static_cast<int>(zaxArr[0]);
                int yy = static_cast<int>(zaxArr[1]);
                int nZax = static_cast<int>(zaxArr[2]);

                vector<complex<double>> bxstore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].Bx_ -> point(xx,yy), stride, bxstore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].Bx_ -> point(xx,yy), stride);
                zscal_(nZax, zaxArr[6],             &Hx_ -> point(xx,yy), stride_rel);

                zaxpy_(nZax,      zaxArr[5], &Ez_ -> point(xx,yy  ), stride_rel, &pmlArr_[kk].Bx_ -> point(xx,yy), stride);
                zaxpy_(nZax, -1.0*zaxArr[5], &Ez_ -> point(xx,yy+1), stride_rel, &pmlArr_[kk].Bx_ -> point(xx,yy), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].Bx_ -> point(xx,yy), stride, &Hx_ -> point(xx,yy), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], bxstore.data()                  , 1     , &Hx_ -> point(xx,yy), stride_rel);
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxHx_end_.size();zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxHx_end_[zz];
                int x_rel = ni + pow(-1, 1-d) * static_cast<int>(zaxArr[0]);
                int y_rel = nj + pow(-1, d)   * static_cast<int>(zaxArr[1]); //the -1^(d()) is to account for changing signs
                int xx = static_cast<int>(zaxArr[0]);
                int yy = static_cast<int>(zaxArr[1]);
                int nZax = static_cast<int>(zaxArr[2]);

                vector<complex<double>> bxstore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].Bx_end_ -> point(xx,yy), stride, bxstore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].Bx_end_ -> point(xx   ,yy   ), stride);
                zscal_(nZax, zaxArr[6],                 &Hx_ -> point(x_rel,y_rel), stride_rel);

                zaxpy_(nZax,      zaxArr[5], &Ez_ -> point(x_rel,y_rel  ), stride_rel, &pmlArr_[kk].Bx_end_ -> point(xx,yy), stride);
                zaxpy_(nZax, -1.0*zaxArr[5], &Ez_ -> point(x_rel,y_rel+1), stride_rel, &pmlArr_[kk].Bx_end_ -> point(xx,yy), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].Bx_end_ -> point(xx,yy), stride, &Hx_ -> point(x_rel,y_rel), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], bxstore.data()                      , 1     , &Hx_ -> point(x_rel,y_rel), stride_rel);
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxHy_.size();zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxHy_[zz];
                int xx = static_cast<int>(zaxArr[0]);
                int yy = static_cast<int>(zaxArr[1]);
                int nZax = static_cast<int>(zaxArr[2]);
                vector<complex<double>> bystore(nZax, 0.0);

                zcopy_(nZax, &pmlArr_[kk].By_ -> point(xx,yy), stride, bystore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].By_ -> point(xx,yy), stride);
                zscal_(nZax, zaxArr[6],             &Hy_ -> point(xx,yy), stride_rel);

                zaxpy_(nZax,      zaxArr[5], &Ez_ -> point(xx+1,yy), stride_rel, &pmlArr_[kk].By_ -> point(xx,yy), stride);
                zaxpy_(nZax, -1.0*zaxArr[5], &Ez_ -> point(xx,yy  ), stride_rel, &pmlArr_[kk].By_ -> point(xx,yy), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].By_ -> point(xx,yy), stride, &Hy_ -> point(xx,yy), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], bystore.data()                  , 1     , &Hy_ -> point(xx,yy), stride_rel);
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxHy_end_.size();zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxHy_end_[zz];
                int x_rel = ni + pow(-1, 1-d) * static_cast<int>(zaxArr[0]);
                int y_rel = nj + pow(-1, d)   * static_cast<int>(zaxArr[1]); //the -1^(d()) is to account for changing signs
                int xx = static_cast<int>(zaxArr[0]);
                int yy = static_cast<int>(zaxArr[1]);
                int nZax = static_cast<int>(zaxArr[2]);

                vector<complex<double>> bystore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].By_end_ -> point(xx,yy), stride, bystore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].By_end_ -> point(xx   ,yy   ), stride);
                zscal_(nZax, zaxArr[6],                 &Hy_ -> point(x_rel,y_rel), stride_rel);

                zaxpy_(nZax,      zaxArr[5], &Ez_ -> point(x_rel+1,y_rel), stride_rel, &pmlArr_[kk].By_end_ -> point(xx,yy), stride);
                zaxpy_(nZax, -1.0*zaxArr[5], &Ez_ -> point(x_rel  ,y_rel), stride_rel, &pmlArr_[kk].By_end_ -> point(xx,yy), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].By_end_ -> point(xx,yy), stride, &Hy_ -> point(x_rel,y_rel), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], bystore.data()                      , 1     , &Hy_ -> point(x_rel,y_rel), stride_rel);
            }
        }
    }
}

/**
 * @brief Updates the E field components for both free space and inside the PML
 * @details Te following equations use teh FDTD update equations to update the E-field to the next time step.
 * For Normal Space TM modes \\
 * \f$E_z^{q+1}\left[i,j\right] = \frac{1 - \frac{\sigma\left[i,j\right]*d_t}{2*\epsilon\left[i,j\right]}}{1+\frac{\sigma\left[i,j\right]*d_t}{2\epsilon\left[i,j\right]}} E_z^q\left[i,j]\right] +\frac{\frac{dt}{\epsilon\left[i,j\right]}}{1+\frac{\sigma\left[i,j\right]*d_t}{2*\epsilon\left[i,j\right]}}\left\{\frac{1}{dx} \left(H_y^{q+\frac{1}{2}}\left[i+\frac{1}{2},j+\frac{1}{2}\right] - H_y^{q+\frac{1}{2}}\left[i-\frac{1}{2},j+\frac{1}{2}\right]\right) - \frac{1}{dy}\left(H_x^{q+\frac{1}{2}}\left[i+\frac{1}{2},j+\frac{1}{2}\right] - H_x^{q+\frac{1}{2}}\left[i+\frac{1}{2},j-\frac{1}{2}\right]\right) \right\}  \f$
 *
 * For UPML space TM mode \\
 * \f$D_z^{q+1}\left[i,j\right] = \frac{2 \epsilon \kappa_x- \sigma_x dt}{2 \epsilon  \kappa_x + \sigma_x dt} D_z^{q}\left[i,j\right] +  \frac{2  \epsilon  dt}{2  \epsilon  \kappa_x + \sigma_x  dt} \left\{\frac{1}{dx} \left(H_y^{q+\frac{1}{2}}\left[i+\frac{1}{2},j+\frac{1}{2}\right] - H_y^{q+\frac{1}{2}}\left[i-\frac{1}{2},j+\frac{1}{2}\right]\right) - \frac{1}{dy}\left(H_x^{q+\frac{1}{2}}\left[i+\frac{1}{2},j+\frac{1}{2}\right] - H_x^{q+\frac{1}{2}}\left[i+\frac{1}{2},j-\frac{1}{2}\right]\right) \right\}\f$
 *
 *\f$E_z^{q+1}\left[i,j\right] =  \frac{2 \epsilon \kappa_y- \sigma_y dt}{2 \epsilon  \kappa_y + \sigma_y dt} E_z^{q}\left[i,j\right] + \frac{1}{\left(2\epsilon \kappa_y +\sigma_y dt\right)\epsilon} \left\{ \left(2\epsilon\kappa_z + \sigma_z dt\right) D_z^{q+1}\left[i,j\right]  -\left(2\epsilon\kappa_z - \sigma_z dt\right) D_z^{q}\left[i,j\right] \right\} \f$
 *
 */
void FDTDField::updateE()
{
    int srcX = srcArr_[0].loc()[0];
    int srcY = srcArr_[0].loc()[1];
    if(Ez_)
    {
        double eps =1.0;
        complex<double> c_eze(1.0,0.0);
        double c_ezh = dt_/(eps*dx_);
        if(xPML_ != 0 && yPML_ !=0)
        {
            for(int kk = 0; kk < zaxEzList_.size(); kk++)
            {
                eps = objArr_[zaxEzList_[kk][3]].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                zscal_(zaxEzList_[kk][2], c_eze, &Ez_ ->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]-1), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hy_->point(zaxEzList_[kk][0]-1,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hy_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
            }
            /*for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
                {
                    eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0);
                    c_ezh = dt_/(eps*dx_);
                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                }
            }*/
        }
        else if(xPML_ != 0)
        {
            /*for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                for(int jj = 1; jj < ny_ - 1; jj ++)
                {
                    eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0);
                    c_ezh = dt_/(eps*dx_);
                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                }
            }*/
            for(int kk = 0; kk < y0EdgeInd_; kk++)
            {
                eps = objArr_[zaxEzList_[kk][3]].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                zscal_(zaxEzList_[kk][2], c_eze, &Ez_ ->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]-1), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hy_->point(zaxEzList_[kk][0]-1,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hy_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
            }
            for(int kk = y0EdgeInd_; kk < ynEdgeInd_; kk++)
            {
                eps = objArr_[zaxEzList_[kk][3]].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                zscal_(zaxEzList_[kk][2], c_eze, &Ez_ ->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hy_->point(zaxEzList_[kk][0]-1,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hy_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
            }
            for(int kk = ynEdgeInd_; kk < zaxEzList_.size(); kk++)
            {
                eps = objArr_[zaxEzList_[kk][3]].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                zscal_(zaxEzList_[kk][2], c_eze, &Ez_ ->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]-1), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hy_->point(zaxEzList_[kk][0]-1,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hy_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
            }
            /*for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                eps = objArr_[phys_Ez_->point(ii,0)].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                Ez_->point(ii,0) = c_eze * Ez_->point(ii,0) + c_ezh * ((Hy_->point(ii,0)-Hy_->point(ii-1,0)) - (Hx_->point(ii,0)));
                eps = objArr_[phys_Ez_->point(ii,ny_-1)].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                Ez_->point(ii,ny_-1) = c_eze * Ez_->point(ii,ny_-1) + c_ezh * ((Hy_->point(ii,ny_-1)-Hy_->point(ii-1,ny_-1)) - (-1.0*Hx_->point(ii,ny_-1-1)));
            }*/
        }
        else if(yPML_ != 0)
        {
            /*for(int ii = 1; ii < nx_ - 1; ii ++)
            {
                for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
                {
                    eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0);
                    c_ezh = dt_/(eps*dx_);
                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                }
            }*/
            for(int kk = 0; kk < x0EdgeInd_; kk++)
            {
                eps = objArr_[zaxEzList_[kk][3]].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                zscal_(zaxEzList_[kk][2], c_eze, &Ez_ ->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]-1), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hy_->point(zaxEzList_[kk][0]-1,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hy_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
            }
            for(int kk = x0EdgeInd_; kk < xnEdgeInd_; kk++)
            {
                eps = objArr_[zaxEzList_[kk][3]].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                zscal_(zaxEzList_[kk][2], c_eze, &Ez_ ->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), nx_,   &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),nx_);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]-1), nx_,   &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),nx_);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hy_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), nx_-1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),nx_);
            }
            for(int kk = xnEdgeInd_; kk < zaxEzList_.size(); kk++)
            {
                eps = objArr_[zaxEzList_[kk][3]].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                zscal_(zaxEzList_[kk][2], c_eze, &Ez_ ->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]-1), nx_,   &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),nx_);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), nx_,   &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),nx_);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hy_->point(zaxEzList_[kk][0]-1,zaxEzList_[kk][1]  ), nx_-1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),nx_);
            }
            /*for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
            {
                eps = objArr_[phys_Ez_->point(0,jj)].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                Ez_->point(0,jj) = c_eze * Ez_->point(0,jj) + c_ezh * ((Hy_->point(0,jj)) - (Hx_->point(0,jj)-Hx_->point(0,jj-1)));
                eps = objArr_[phys_Ez_->point(nx_-1,jj)].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                Ez_->point(nx_-1,jj) = c_eze * Ez_->point(nx_-1,jj) + c_ezh * ((-1.0*Hy_->point(nx_-1-1,jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1)));
            }*/
        }
        else
        {
            /*for(int ii = 1; ii < nx_ - 1; ii ++)
            {
                for(int jj = 1; jj < ny_ - 1; jj ++)
                {
                    eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0);
                    c_ezh = dt_/(eps*dx_);
                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                }
            }*/
            for(int kk = 0; kk < y0EdgeInd_; kk++)
            {
                eps = objArr_[zaxEzList_[kk][3]].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                zscal_(zaxEzList_[kk][2], c_eze, &Ez_ ->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]-1), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hy_->point(zaxEzList_[kk][0]-1,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hy_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
            }
            /*for(int jj = 1; jj < ny_ - 1; jj ++)
            {
                eps = objArr_[phys_Ez_->point(0,jj)].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                Ez_->point(0,jj) = c_eze * Ez_->point(0,jj) + c_ezh * ((Hy_->point(0,jj)) - (Hx_->point(0,jj)-Hx_->point(0,jj-1)));
                eps = objArr_[phys_Ez_->point(nx_-1,jj)].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                Ez_->point(nx_-1,jj) = c_eze * Ez_->point(nx_-1,jj) + c_ezh * ((-1.0*Hy_->point(nx_-1-1,jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1)));
            }*/
            for(int kk = y0EdgeInd_; kk < ynEdgeInd_; kk++)
            {
                eps = objArr_[zaxEzList_[kk][3]].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                zscal_(zaxEzList_[kk][2], c_eze, &Ez_ ->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hy_->point(zaxEzList_[kk][0]-1,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hy_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
            }
            for(int kk = ynEdgeInd_; kk < x0EdgeInd_; kk++)
            {
                eps = objArr_[zaxEzList_[kk][3]].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                zscal_(zaxEzList_[kk][2], c_eze, &Ez_ ->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hy_->point(zaxEzList_[kk][0]-1,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hy_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), 1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
            }
            for(int kk = x0EdgeInd_; kk < xnEdgeInd_; kk++)
            {
                eps = objArr_[zaxEzList_[kk][3]].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                zscal_(zaxEzList_[kk][2], c_eze, &Ez_ ->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), nx_,   &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),nx_);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]-1), nx_,   &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),nx_);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hy_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), nx_-1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),nx_);
            }
            for(int kk = xnEdgeInd_; kk < zaxEzList_.size(); kk++)
            {
                eps = objArr_[zaxEzList_[kk][3]].dielectric(1.0);
                c_ezh = dt_/(eps*dx_);
                zscal_(zaxEzList_[kk][2], c_eze, &Ez_ ->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),1);
                zaxpy_(zaxEzList_[kk][2],     c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]-1), nx_,   &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),nx_);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hx_->point(zaxEzList_[kk][0]  ,zaxEzList_[kk][1]  ), nx_,   &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),nx_);
                zaxpy_(zaxEzList_[kk][2],-1.0*c_ezh, &Hy_->point(zaxEzList_[kk][0]-1,zaxEzList_[kk][1]  ), nx_-1, &Ez_->point(zaxEzList_[kk][0],zaxEzList_[kk][1]),nx_);
            }
            eps = objArr_[phys_Ez_->point(0,0)].dielectric(1.0);
            c_ezh = dt_/(eps*dx_);
            Ez_->point(0,0) = c_eze * Ez_->point(0,0) + c_ezh * ((Hy_->point(0,0)) - (Hx_->point(0,0)));
            eps = objArr_[phys_Ez_->point(0,ny_-1)].dielectric(1.0);
            c_ezh = dt_/(eps*dx_);
            Ez_->point(0,ny_-1) = c_eze * Ez_->point(0,ny_-1) + c_ezh * ((Hy_->point(0,ny_-1)) - (-1.0*Hx_->point(0,ny_-1-1)));
            eps = objArr_[phys_Ez_->point(nx_-1,0)].dielectric(1.0);
            c_ezh = dt_/(eps*dx_);
            Ez_->point(nx_-1,0) = c_eze * Ez_->point(nx_-1,0) + c_ezh * ((-1.0*Hy_->point(nx_-1-1,0)) - (Hx_->point(nx_-1,0)));
            eps = objArr_[phys_Ez_->point(nx_-1,ny_-1)].dielectric(1.0);
            c_ezh = dt_/(eps*dx_);
            Ez_->point(nx_-1,ny_-1) = c_eze * Ez_->point(nx_-1,ny_-1) + c_ezh * ((-1.0*Hy_->point(nx_-1-1,ny_-1)) - (-1.0*Hx_->point(nx_-1,ny_-1-1)));
        }
    }
    if(precalcPML_ == false)
    {
        if(Ez_)
        {
            for(int kk = 0; kk < pmlArr_.size(); kk++)
            {
                switch(pmlArr_[kk].d())
                {
                    case X:
                    {
                        double eps = 1.0;
                        double kapx = 1.0;
                        double kapy = 1.0;
                        double kapz = 1.0;
                        double sigx = 0.0;
                        double sigy = 0.0;
                        double sigz = 0.0;
                        double c_dzd = 0.0; double c_dzh = 0.0; double c_eze = 0.0; double c_ezd1 = 0.0; double c_ezd0 = 0.0;
                        complex<double> dzstore(0.0,0.0);
                        if(yPML_ == 0)
                        {
                            for(int jj = 1; jj < ny_ - 1; jj++)
                            {
                                eps = objArr_[phys_Ez_->point(0,jj)].dielectric(1.0);
                                sigx = pmlArr_[kk].sigma(0.0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;

                                dzstore = pmlArr_[kk].Dz_->point(0,jj);
                                pmlArr_[kk].Dz_->point(0,jj) = c_dzd * pmlArr_[kk].Dz_->point(0,jj) + c_dzh * ((Hy_->point(0,jj)) - (Hx_->point(0,jj)-Hx_->point(0,jj-1))); // 0 is boundary condtion
                                Ez_->point(0,jj) = c_eze * Ez_->point(0,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(0,jj) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(nx_-1,jj)].dielectric(1.0);
                                sigx = pmlArr_[kk].sigma(0.0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;

                                dzstore = pmlArr_[kk].Dz_end_->point(0,jj);
                                pmlArr_[kk].Dz_end_->point(0,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(0,jj) + c_dzh * ((-1.0 *Hy_->point((nx_-2),jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1))); // 0 is boundary condition
                                Ez_->point(nx_-1,jj) = c_eze * Ez_->point(nx_-1,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(0,jj) - c_ezd0 * dzstore;
                            }
                            eps = objArr_[phys_Ez_->point(0,ny_-1)].dielectric(1.0);
                            sigx = pmlArr_[kk].sigma(0.0,eps);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;

                            dzstore = pmlArr_[kk].Dz_->point(0,ny_-1);
                            pmlArr_[kk].Dz_->point(0,ny_-1) = c_dzd * pmlArr_[kk].Dz_->point(0,ny_-1) + c_dzh * ((Hy_->point(0,ny_-1)) - (-1.0*Hx_->point(0,ny_-2))); // 0 is boundary condtion
                            Ez_->point(0,ny_-1) = c_eze * Ez_->point(0,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_->point(0,ny_-1) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(nx_-1,ny_-1)].dielectric(1.0);
                            sigx = pmlArr_[kk].sigma(0.0,eps);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_end_->point(0,ny_-1);
                            pmlArr_[kk].Dz_end_->point(0,ny_-1) = c_dzd * pmlArr_[kk].Dz_end_->point(0,ny_-1) + c_dzh * ((-1.0 *Hy_->point(nx_-2,ny_-1)) - (-1.0*Hx_->point(nx_-1,ny_-2))); // 0 is boundary condition
                            Ez_->point(nx_-1,ny_-1) = c_eze * Ez_->point(nx_-1,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_end_->point(0,ny_-1) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(0,0)].dielectric(1.0);
                            sigx = pmlArr_[kk].sigma(0.0,eps);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_->point(0,0);
                            pmlArr_[kk].Dz_->point(0,0) = c_dzd * pmlArr_[kk].Dz_->point(0,0) + c_dzh * ((Hy_->point(0,0)) - (Hx_->point(0,0))); // 0 is boundary condtion
                            Ez_->point(0,0) = c_eze * Ez_->point(0,0) + c_ezd1 * pmlArr_[kk].Dz_->point(0,0) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(nx_-1,0)].dielectric(1.0);
                            sigx = pmlArr_[kk].sigma(0.0,eps);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_end_->point(0,0);
                            pmlArr_[kk].Dz_end_->point(0,0) = c_dzd * pmlArr_[kk].Dz_end_->point(0,0) + c_dzh * ((-1.0 *Hy_->point(nx_-2,0)) - (Hx_->point(nx_-1,0))); // 0 is boundary condition
                            Ez_->point(nx_-1,0) = c_eze * Ez_->point(nx_-1,0) + c_ezd1 * pmlArr_[kk].Dz_end_->point(0,0) - c_ezd0 * dzstore;

                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                sigy = 0.0; sigz = 0.0;
                                //double kapx = pmlArr_[kk].kappa(ii);
                                //Kappas change throughout
                                for(int jj =  1; jj < ny_ - 1; jj ++)
                                {
                                    eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0);
                                    sigx = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[kk].Dz_->point(ii,jj);
                                    pmlArr_[kk].Dz_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,jj) - c_ezd0 * dzstore;

                                    eps = objArr_[phys_Ez_->point(nx_-1-ii,jj)].dielectric(1.0);
                                    sigx = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[kk].Dz_end_->point(ii,jj);
                                    pmlArr_[kk].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point(nx_-1-ii,jj)-Hy_->point(nx_-1-ii-1,jj)) - (Hx_->point(nx_-1-ii,jj)-Hx_->point(nx_-1-ii,jj-1)));
                                    Ez_->point(nx_-1-ii,jj) = c_eze * Ez_->point(nx_-1-ii,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                                eps = objArr_[phys_Ez_->point(ii,ny_-1)].dielectric(1.0);
                                sigx = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_->point(ii,ny_-1);
                                pmlArr_[kk].Dz_->point(ii,ny_-1) = c_dzd * pmlArr_[kk].Dz_->point(ii,ny_-1) + c_dzh * ((Hy_->point(ii,ny_-1)-Hy_->point(ii-1,ny_-1)) - (-1.0*Hx_->point(ii,ny_-1-1)));
                                Ez_->point(ii,ny_-1) = c_eze * Ez_->point(ii,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,ny_-1) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(nx_-1-ii,ny_-1)].dielectric(1.0);
                                sigx = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_end_->point(ii,ny_-1);
                                pmlArr_[kk].Dz_end_->point(ii,ny_-1) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,ny_-1) + c_dzh * ((Hy_->point(nx_-1-ii,ny_-1)-Hy_->point(nx_-1-ii-1,ny_-1)) - (-1.0*Hx_->point(nx_-1-ii,ny_-1-1)));
                                Ez_->point(nx_-1-ii,ny_-1) = c_eze * Ez_->point(nx_-1-ii,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,ny_-1) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(ii,0)].dielectric(1.0);
                                sigx = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_->point(ii,0);
                                pmlArr_[kk].Dz_->point(ii,0) = c_dzd * pmlArr_[kk].Dz_->point(ii,0) + c_dzh * ((Hy_->point(ii,0)-Hy_->point(ii-1,0)) - (Hx_->point(ii,0)));
                                Ez_->point(ii,0) = c_eze * Ez_->point(ii,0) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,0) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(nx_-1-ii,0)].dielectric(1.0);
                                sigx = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_end_->point(ii,0);
                                pmlArr_[kk].Dz_end_->point(ii,0) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,0) + c_dzh * ((Hy_->point(nx_-1-ii,0)-Hy_->point(nx_-1-ii-1,0)) - (Hx_->point(nx_-1-ii,0)));
                                Ez_->point(nx_-1-ii,0) = c_eze * Ez_->point(nx_-1-ii,0) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,0) - c_ezd0 * dzstore;
                            }
                        }
                        else
                        {
                            for(int jj = yPML_; jj < ny_ - yPML_; jj++)
                            {
                                eps = objArr_[phys_Ez_->point(0,jj)].dielectric(1.0);
                                sigx = pmlArr_[kk].sigma(0.0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_->point(0,jj);
                                pmlArr_[kk].Dz_->point(0,jj) = c_dzd * pmlArr_[kk].Dz_->point(0,jj) + c_dzh * ((Hy_->point(0,jj)) - (Hx_->point(0,jj)-Hx_->point(0,jj-1))); // 0 is boundary condtion
                                Ez_->point(0,jj) = c_eze * Ez_->point(0,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(0,jj) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(nx_-1,jj)].dielectric(1.0);
                                sigx = pmlArr_[kk].sigma(0.0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_end_->point(0,jj);
                                pmlArr_[kk].Dz_end_->point(0,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(0,jj) + c_dzh * ((-1.0 *Hy_->point((nx_-2),jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1))); // 0 is boundary condition
                                Ez_->point(nx_-1,jj) = c_eze * Ez_->point(nx_-1,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(0,jj) - c_ezd0 * dzstore;
                            }
                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                sigy = 0.0; sigz = 0.0;
                                for(int jj =  yPML_; jj < ny_ - yPML_; jj ++)
                                {
                                    eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0);
                                    sigx = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[kk].Dz_->point(ii,jj);
                                    pmlArr_[kk].Dz_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,jj) - c_ezd0 * dzstore;

                                    eps = objArr_[phys_Ez_->point(nx_-1-ii,jj)].dielectric(1.0);
                                    sigx = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[kk].Dz_end_->point(ii,jj);
                                    pmlArr_[kk].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point(nx_-1-ii,jj)-Hy_->point(nx_-1-ii-1,jj)) - (Hx_->point(nx_-1-ii,jj)-Hx_->point(nx_-1-ii,jj-1)));
                                    Ez_->point(nx_-1-ii,jj) = c_eze * Ez_->point(nx_-1-ii,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                }
                            }
                            if(kk == 0)
                            {
                                kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                sigz = 0.0;

                                eps = objArr_[phys_Ez_->point(0,0)].dielectric(1.0);
                                sigx = pmlArr_[0].sigma(0,eps);
                                sigy = pmlArr_[1].sigma(0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                //Bot Left
                                dzstore = pmlArr_[0].Dz_->point(0,0);
                                pmlArr_[0].Dz_->point(0,0) = c_dzd * pmlArr_[0].Dz_->point(0,0) + c_dzh * ((Hy_->point(0,0)) - (Hx_->point(0,0)));
                                Ez_->point(0,0) = c_eze * Ez_->point(0,0) + c_ezd1 * pmlArr_[0].Dz_->point(0,0) - c_ezd0 * dzstore;
                                //Bot Right
                                eps = objArr_[phys_Ez_->point(nx_-1,0)].dielectric(1.0);
                                sigx = pmlArr_[0].sigma(0,eps);
                                sigy = pmlArr_[1].sigma(0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[0].Dz_end_->point(0,0);
                                pmlArr_[0].Dz_end_->point(0,0) = c_dzd * pmlArr_[0].Dz_end_->point(0,0) + c_dzh * (-1.0*Hy_->point(nx_-1-1,0)) - (Hx_->point(nx_-1,0));
                                Ez_->point(nx_-1,0) = c_eze * Ez_->point(nx_-1,0) + c_ezd1 * pmlArr_[0].Dz_end_->point(0,0) - c_ezd0 * dzstore;
                                //Top Right
                                eps = objArr_[phys_Ez_->point(nx_-1,ny_-1)].dielectric(1.0);
                                sigx = pmlArr_[0].sigma(0,eps);
                                sigy = pmlArr_[1].sigma(0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[0].Dz_end_->point(0,ny_-1);
                                pmlArr_[0].Dz_end_->point(0,ny_-1) = c_dzd * pmlArr_[0].Dz_end_->point(0,ny_-1) + c_dzh * ((-1.0*Hy_->point(nx_-1-1,ny_-1)) - (-1.0*Hx_->point(nx_-1,ny_-1-1)));
                                Ez_->point(nx_-1,ny_-1) = c_eze * Ez_->point(nx_-1,ny_-1) + c_ezd1 * pmlArr_[0].Dz_end_->point(0,ny_-1) - c_ezd0 * dzstore;
                                //Top Left
                                eps = objArr_[phys_Ez_->point(0,ny_-1)].dielectric(1.0);
                                sigx = pmlArr_[0].sigma(0,eps);
                                sigy = pmlArr_[1].sigma(0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[0].Dz_->point(0,ny_-1);
                                pmlArr_[0].Dz_->point(0,ny_-1) = c_dzd * pmlArr_[0].Dz_->point(0,ny_-1) + c_dzh * ((Hy_->point(0,ny_-1)) - (-1.0*Hx_->point(0,ny_-1-1)));
                                Ez_->point(0,ny_-1) = c_eze * Ez_->point(0,ny_-1) + c_ezd1 * pmlArr_[0].Dz_->point(0,ny_-1) - c_ezd0 * dzstore;

                                for(int ii = 1; ii < pmlArr_[0].thickness(); ii++)
                                {
                                    kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                    sigz = 0.0;

                                    eps = objArr_[phys_Ez_->point(ii,0)].dielectric(1.0);
                                    sigx = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                    sigy = pmlArr_[1].sigma(0.0,eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    //Bot Left
                                    dzstore = pmlArr_[0].Dz_->point(ii,0);
                                    pmlArr_[0].Dz_->point(ii,0) = c_dzd * pmlArr_[0].Dz_->point(ii,0) + c_dzh * ((Hy_->point(ii,0)-Hy_->point(ii-1,0)) - (Hx_->point(ii,0)));
                                    Ez_->point(ii,0) = c_eze * Ez_->point(ii,0) + c_ezd1 * pmlArr_[0].Dz_->point(ii,0) - c_ezd0 * dzstore;
                                    //Bot Right
                                    eps = objArr_[phys_Ez_->point(nx_-1-ii,0)].dielectric(1.0);
                                    sigx = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                    sigy = pmlArr_[1].sigma(0.0,eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[0].Dz_end_->point(ii,0);
                                    pmlArr_[0].Dz_end_->point(ii,0) = c_dzd * pmlArr_[0].Dz_end_->point(ii,0) + c_dzh * ((Hy_->point(nx_-1-ii,0)-Hy_->point(nx_-1-ii-1,0)) - (Hx_->point(nx_-1-ii,0)));
                                    Ez_->point(nx_-1-ii,0) = c_eze * Ez_->point(nx_-1-ii,0) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,0) - c_ezd0 * dzstore;
                                    //Top Right
                                    eps = objArr_[phys_Ez_->point(nx_-1-ii,ny_-1)].dielectric(1.0);
                                    sigx = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                    sigy = pmlArr_[1].sigma(0.0,eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[0].Dz_end_->point(ii,ny_-1);
                                    pmlArr_[0].Dz_end_->point(ii,ny_-1) = c_dzd * pmlArr_[0].Dz_end_->point(ii,ny_-1) + c_dzh * ((Hy_->point(nx_-1-ii,ny_-1)-Hy_->point(nx_-1-ii-1,ny_-1)) - (-1.0*Hx_->point(nx_-1-ii,ny_-1-1)));
                                    Ez_->point(nx_-1-ii,ny_-1) = c_eze * Ez_->point(nx_-1-ii,ny_-1) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,ny_-1) - c_ezd0 * dzstore;
                                    //Top Left
                                    eps = objArr_[phys_Ez_->point(ii,ny_-1)].dielectric(1.0);
                                    sigx = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                    sigy = pmlArr_[1].sigma(0.0,eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[0].Dz_->point(ii,ny_-1);
                                    pmlArr_[0].Dz_->point(ii,ny_-1) = c_dzd * pmlArr_[0].Dz_->point(ii,ny_-1) + c_dzh * ((Hy_->point(ii,ny_-1)-Hy_->point(ii-1,ny_-1)) - (-1.0*Hx_->point(ii,ny_-1-1)));
                                    Ez_->point(ii,ny_-1) = c_eze * Ez_->point(ii,ny_-1) + c_ezd1 * pmlArr_[0].Dz_->point(ii,ny_-1) - c_ezd0 * dzstore;
                                }
                                for(int jj = 1; jj < pmlArr_[0].thickness(); jj++)
                                {
                                    kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                    sigz = 0.0;

                                    eps = objArr_[phys_Ez_->point(0,jj)].dielectric(1.0);
                                    sigx = pmlArr_[0].sigma(0.0,eps);
                                    sigy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    //Bot Left
                                    dzstore = pmlArr_[0].Dz_->point(0,jj);
                                    pmlArr_[0].Dz_->point(0,jj) = c_dzd * pmlArr_[0].Dz_->point(0,jj) + c_dzh * ((Hy_->point(0,jj)) - (Hx_->point(0,jj)-Hx_->point(0,jj-1)));
                                    Ez_->point(0,jj) = c_eze * Ez_->point(0,jj) + c_ezd1 * pmlArr_[0].Dz_->point(0,jj) - c_ezd0 * dzstore;
                                    //Bot Right
                                    eps = objArr_[phys_Ez_->point(nx_-1,jj)].dielectric(1.0);
                                    sigx = pmlArr_[0].sigma(0.0,eps);
                                    sigy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[0].Dz_end_->point(0,jj);
                                    pmlArr_[0].Dz_end_->point(0,jj) = c_dzd * pmlArr_[0].Dz_end_->point(0,jj) + c_dzh * ((-1.0*Hy_->point(nx_-1-1,jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1)));
                                    Ez_->point(nx_-1,jj) = c_eze * Ez_->point(nx_-1,jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(0,jj) - c_ezd0 * dzstore;
                                    //Top Right
                                    eps = objArr_[phys_Ez_->point(nx_-1,ny_-1-jj)].dielectric(1.0);
                                    sigx = pmlArr_[0].sigma(0.0,eps);
                                    sigy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[0].Dz_end_->point(0,ny_-1-jj);
                                    pmlArr_[0].Dz_end_->point(0,ny_-1-jj) = c_dzd * pmlArr_[0].Dz_end_->point(0,ny_-1-jj) + c_dzh * ((-1.0*Hy_->point(nx_-1-1,ny_-1-jj)) - (Hx_->point(nx_-1,ny_-1-jj)-Hx_->point(nx_-1,ny_-1-jj-1)));
                                    Ez_->point(nx_-1,ny_-1-jj) = c_eze * Ez_->point(nx_-1,ny_-1-jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(0,ny_-1-jj) - c_ezd0 * dzstore;
                                    //Top Left
                                    eps = objArr_[phys_Ez_->point(0,ny_-1-jj)].dielectric(1.0);
                                    sigx = pmlArr_[0].sigma(0.0,eps);
                                    sigy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[0].Dz_->point(0,ny_-1-jj);
                                    pmlArr_[0].Dz_->point(0,ny_-1-jj) = c_dzd * pmlArr_[0].Dz_->point(0,ny_-1-jj) + c_dzh * ((Hy_->point(0,ny_-1-jj)) - (Hx_->point(0,ny_-1-jj)-Hx_->point(0,ny_-1-jj-1)));
                                    Ez_->point(0,ny_-1-jj) = c_eze * Ez_->point(0,ny_-1-jj) + c_ezd1 * pmlArr_[0].Dz_->point(0,ny_-1-jj) - c_ezd0 * dzstore;
                                }
                                for(int ii = 1; ii < pmlArr_[0].thickness(); ii++)
                                {
                                    for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                                    {
                                        kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                        sigz = 0.0;

                                        eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0);
                                        sigx = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                        sigy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        //Bot Left
                                        dzstore = pmlArr_[0].Dz_->point(ii,jj);
                                        pmlArr_[0].Dz_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                        Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                        //Bot Right
                                        eps = objArr_[phys_Ez_->point(nx_-1-ii,jj)].dielectric(1.0);
                                        sigx = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                        sigy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        dzstore = pmlArr_[0].Dz_end_->point(ii,jj);
                                        pmlArr_[0].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point(nx_-1-ii,jj)-Hy_->point(nx_-1-ii-1,jj)) - (Hx_->point(nx_-1-ii,jj)-Hx_->point(nx_-1-ii,jj-1)));
                                        Ez_->point(nx_-1-ii,jj) = c_eze * Ez_->point(nx_-1-ii,jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                        //Top Right
                                        eps = objArr_[phys_Ez_->point(nx_-1-ii,ny_-1-jj)].dielectric(1.0);
                                        sigx = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                        sigy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        dzstore = pmlArr_[0].Dz_end_->point(ii,ny_-1-jj);
                                        pmlArr_[0].Dz_end_->point(ii,ny_-1-jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,ny_-1-jj) + c_dzh * ((Hy_->point(nx_-1-ii,ny_-1-jj)-Hy_->point(nx_-1-ii-1,ny_-1-jj)) - (Hx_->point(nx_-1-ii,ny_-1-jj)-Hx_->point(nx_-1-ii,ny_-1-jj-1)));
                                        Ez_->point(nx_-1-ii,ny_-1-jj) = c_eze * Ez_->point(nx_-1-ii,ny_-1-jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,ny_-1-jj) - c_ezd0 * dzstore;
                                        //Top Left
                                        eps = objArr_[phys_Ez_->point(ii,ny_-1-jj)].dielectric(1.0);
                                        sigx = pmlArr_[0].sigma(static_cast<double>(ii),eps);
                                        sigy = pmlArr_[1].sigma(static_cast<double>(jj),eps);
                                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        dzstore = pmlArr_[0].Dz_->point(ii,ny_-1-jj);
                                        pmlArr_[0].Dz_->point(ii,ny_-1-jj) = c_dzd * pmlArr_[0].Dz_->point(ii,ny_-1-jj) + c_dzh * ((Hy_->point(ii,ny_-1-jj)-Hy_->point(ii-1,ny_-1-jj)) - (Hx_->point(ii,ny_-1-jj)-Hx_->point(ii,ny_-1-jj-1)));
                                        Ez_->point(ii,ny_-1-jj) = c_eze * Ez_->point(ii,ny_-1-jj) + c_ezd1 * pmlArr_[0].Dz_->point(ii,ny_-1-jj) - c_ezd0 * dzstore;
                                    }
                                }
                            }
                        }
                        break;
                    }
                    case Y:
                    {
                        double eps = 1.0;
                        double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                        double sigx = 0.0;
                        double sigy = 0.0; double sigz = 0.0;

                        double c_dzd = 0.0; double c_dzh = 0.0; double c_eze = 0.0; double c_ezd1 = 0.0; double c_ezd0 = 0.0;
                        complex<double> dzstore(0.0,0.0);
                        if(xPML_ == 0)
                        {
                            for(int jj = 1; jj < nx_ - 1; jj++)
                            {
                                eps = objArr_[phys_Ez_->point(jj,0)].dielectric(1.0);
                                sigy = pmlArr_[kk].sigma(0.0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;

                                dzstore = pmlArr_[kk].Dz_->point(jj,0);
                                pmlArr_[kk].Dz_->point(jj,0) = c_dzd * pmlArr_[kk].Dz_->point(jj,0) + c_dzh * ((Hy_->point(jj,0)-Hy_->point(jj-1,0)) - (Hx_->point(jj,0))); // 0 is boundary condition
                                Ez_->point(jj,0) = c_eze * Ez_->point(jj,0) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,0) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(jj,ny_-1)].dielectric(1.0);
                                sigy = pmlArr_[kk].sigma(0.0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_end_->point(jj,0);
                                pmlArr_[kk].Dz_end_->point(jj,0) = c_dzd * pmlArr_[kk].Dz_end_->point(jj,0) + c_dzh * ((Hy_->point(jj,ny_-1) - Hy_->point(jj - 1,ny_-1)) - (-1.0 * Hx_->point(jj,ny_-1-1))); // 0 is boundary condtion
                                Ez_->point(jj,ny_-1) = c_eze * Ez_->point(jj,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_end_->point(jj,0) - c_ezd0 * dzstore;
                            }
                            eps = objArr_[phys_Ez_->point(nx_-1,0)].dielectric(1.0);
                            sigy = pmlArr_[kk].sigma(0.0,eps);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_->point(nx_-1,0);
                            pmlArr_[kk].Dz_->point(nx_-1,0) = c_dzd * pmlArr_[kk].Dz_->point(nx_-1,0) + c_dzh * ((-1.0*Hy_->point(nx_-1-1,0)) - (Hx_->point(nx_-1,0))); // 0 is boundary condition
                            Ez_->point(nx_-1,0) = c_eze * Ez_->point(nx_-1,0) + c_ezd1 * pmlArr_[kk].Dz_->point(nx_-1,0) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(nx_-1,ny_-1)].dielectric(1.0);
                            sigy = pmlArr_[kk].sigma(0.0,eps);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_end_->point(nx_-1,0);
                            pmlArr_[kk].Dz_end_->point(nx_-1,0) = c_dzd * pmlArr_[kk].Dz_end_->point(nx_-1,0) + c_dzh * ((-1.0 * Hy_->point(nx_-1 - 1,ny_-1)) - (-1.0 * Hx_->point(nx_-1,ny_-1-1))); // 0 is boundary condtion
                            Ez_->point(nx_-1,ny_-1) = c_eze * Ez_->point(nx_-1,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_end_->point(nx_-1,0) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(0,0)].dielectric(1.0);
                            sigy = pmlArr_[kk].sigma(0.0,eps);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_->point(0,0);
                            pmlArr_[kk].Dz_->point(0,0) = c_dzd * pmlArr_[kk].Dz_->point(0,0) + c_dzh * ((Hy_->point(0,0)) - (Hx_->point(0,0))); // 0 is boundary condition
                            Ez_->point(0,0) = c_eze * Ez_->point(0,0) + c_ezd1 * pmlArr_[kk].Dz_->point(0,0) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(0,ny_-1)].dielectric(1.0);
                            sigy = pmlArr_[kk].sigma(0.0,eps);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_end_->point(0,0);
                            pmlArr_[kk].Dz_end_->point(0,0) = c_dzd * pmlArr_[kk].Dz_end_->point(0,0) + c_dzh * ((Hy_->point(0,ny_-1)) - (-1.0 * Hx_->point(0,ny_-1-1))); // 0 is boundary condtion
                            Ez_->point(0,ny_-1) = c_eze * Ez_->point(0,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_end_->point(0,0) - c_ezd0 * dzstore;

                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                eps = 1.0;
                                kapx = 1.0; kapy = 1.0;
                                kapz = 1.0;
                                sigx = 0.0;
                                sigz = 0.0;

                                for(int jj =  1; jj < nx_ - 1; jj ++)
                                {
                                    eps = objArr_[phys_Ez_->point(jj,ii)].dielectric(1.0);
                                    sigy = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[kk].Dz_->point(jj,ii);
                                    pmlArr_[kk].Dz_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_->point(jj,ii) + c_dzh * ((Hy_->point(jj,ii)-Hy_->point(jj-1,ii)) - (Hx_->point(jj,ii)-Hx_->point(jj,ii-1)));
                                    Ez_->point(jj,ii) = c_eze * Ez_->point(jj,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,ii) - c_ezd0 * dzstore;

                                    eps = objArr_[phys_Ez_->point(jj,ny_-1-ii)].dielectric(1.0);
                                    sigy = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[kk].Dz_end_->point(jj,ii);
                                    pmlArr_[kk].Dz_end_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_end_->point(jj,ii) + c_dzh * ((Hy_->point(jj,ny_-1-ii) - Hy_->point(jj - 1,ny_-1-ii)) - (Hx_->point(jj,ny_-1-ii) - Hx_->point(jj,ny_-1-ii-1)));
                                    Ez_->point(jj,ny_-1-ii) = c_eze * Ez_->point(jj,ny_-1-ii) + c_ezd1 * pmlArr_[kk].Dz_end_->point(jj,ii) - c_ezd0 * dzstore;
                                }
                                eps = objArr_[phys_Ez_->point(nx_-1,ii)].dielectric(1.0);
                                sigy = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_->point(nx_-1,ii);
                                pmlArr_[kk].Dz_->point(nx_-1,ii) = c_dzd * pmlArr_[kk].Dz_->point(nx_-1,ii) + c_dzh * ((-1.0*Hy_->point(nx_-1-1,ii)) - (Hx_->point(nx_-1,ii)-Hx_->point(nx_-1,ii-1)));
                                Ez_->point(nx_-1,ii) = c_eze * Ez_->point(nx_-1,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(nx_-1,ii) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(nx_-1, ny_-1-ii)].dielectric(1.0);
                                sigy = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_end_->point(nx_-1,ii);
                                pmlArr_[kk].Dz_end_->point(nx_-1,ii) = c_dzd * pmlArr_[kk].Dz_end_->point(nx_-1,ii) + c_dzh * ((-1.0 * Hy_->point(nx_-1-1,ny_-1-ii)) - (Hx_->point(nx_-1,ny_-1-ii) - Hx_->point(nx_-1,ny_-1-ii-1)));
                                Ez_->point(nx_-1,ny_-1-ii) = c_eze * Ez_->point(nx_-1,ny_-1-ii) + c_ezd1 * pmlArr_[kk].Dz_end_->point(nx_-1,ii) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(0, ii)].dielectric(1.0);
                                sigy = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_->point(0,ii);
                                pmlArr_[kk].Dz_->point(0,ii) = c_dzd * pmlArr_[kk].Dz_->point(0,ii) + c_dzh * ((Hy_->point(0,ii)) - (Hx_->point(0,ii)-Hx_->point(0,ii-1)));
                                Ez_->point(0,ii) = c_eze * Ez_->point(0,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(0,ii) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(0, ny_-1-ii)].dielectric(1.0);
                                sigy = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_end_->point(0,ii);
                                pmlArr_[kk].Dz_end_->point(0,ii) = c_dzd * pmlArr_[kk].Dz_end_->point(0,ii) + c_dzh * ((Hy_->point(0,ny_-1-ii)) - (Hx_->point(0,ny_-1-ii) - Hx_->point(0,ny_-1-ii-1)));
                                Ez_->point(0,ny_-1-ii) = c_eze * Ez_->point(0,ny_-1-ii) + c_ezd1 * pmlArr_[kk].Dz_end_->point(0,ii) - c_ezd0 * dzstore;
                            }
                        }
                        else
                        {
                            for(int jj = xPML_; jj < nx_ - xPML_; jj++)
                            {
                                eps = objArr_[phys_Ez_->point(jj,0)].dielectric(1.0);
                                sigy = pmlArr_[kk].sigma(0.0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_->point(jj,0);
                                pmlArr_[kk].Dz_->point(jj,0) = c_dzd * pmlArr_[kk].Dz_->point(jj,0) + c_dzh * ((Hy_->point(jj,0)-Hy_->point(jj-1,0)) - (Hx_->point(jj,0))); // 0 is boundary condition
                                Ez_->point(jj,0) = c_eze * Ez_->point(jj,0) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,0) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(jj, ny_-1)].dielectric(1.0);
                                sigy = pmlArr_[kk].sigma(0.0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_end_->point(jj,0);
                                pmlArr_[kk].Dz_end_->point(jj,0) = c_dzd * pmlArr_[kk].Dz_end_->point(jj,0) + c_dzh * ((Hy_->point(jj,ny_-1) - Hy_->point(jj - 1,ny_-1)) - (-1.0 * Hx_->point(jj,ny_-1-1))); // 0 is boundary condtion
                                Ez_->point(jj,ny_-1) = c_eze * Ez_->point(jj,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_end_->point(jj,0) - c_ezd0 * dzstore;
                            }
                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                kapx = 1.0; kapy = 1.0;
                                kapz = 1.0;
                                sigx = 0.0;
                                sigz = 0.0;

                                for(int jj =  xPML_; jj < nx_ - xPML_; jj ++)
                                {
                                    eps = objArr_[phys_Ez_->point(jj, ii)].dielectric(1.0);
                                    sigy = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[kk].Dz_->point(jj,ii);
                                    pmlArr_[kk].Dz_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_->point(jj,ii) + c_dzh * ((Hy_->point(jj,ii)-Hy_->point(jj-1,ii)) - (Hx_->point(jj,ii)-Hx_->point(jj,ii-1)));
                                    Ez_->point(jj,ii) = c_eze * Ez_->point(jj,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,ii) - c_ezd0 * dzstore;

                                    eps = objArr_[phys_Ez_->point(jj, ny_-1-ii)].dielectric(1.0);
                                    sigy = pmlArr_[kk].sigma(static_cast<double>(ii),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[kk].Dz_end_->point(jj,ii);
                                    pmlArr_[kk].Dz_end_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_end_->point(jj,ii) + c_dzh * ((Hy_->point(jj,ny_-1-ii) - Hy_->point(jj - 1,ny_-1-ii)) - (Hx_->point(jj,ny_-1-ii) - Hx_->point(jj,ny_-1-ii-1)));
                                    Ez_->point(jj,ny_-1-ii) = c_eze * Ez_->point(jj,ny_-1-ii) + c_ezd1 * pmlArr_[kk].Dz_end_->point(jj,ii) - c_ezd0 * dzstore;
                                }
                            }
                            if(kk == 0)
                            {
                                kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                sigz = 0.0;

                                eps = objArr_[phys_Ez_->point(0,0)].dielectric(1.0);
                                sigx = pmlArr_[1].sigma(0,eps);
                                sigy = pmlArr_[0].sigma(0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                //Bot Left
                                dzstore = pmlArr_[1].Dz_->point(0,0);
                                pmlArr_[1].Dz_->point(0,0) = c_dzd * pmlArr_[1].Dz_->point(0,0) + c_dzh * ((Hy_->point(0,0)) - (Hx_->point(0,0)));
                                Ez_->point(0,0) = c_eze * Ez_->point(0,0) + c_ezd1 * pmlArr_[1].Dz_->point(0,0) - c_ezd0 * dzstore;
                                //Bot Right
                                eps = objArr_[phys_Ez_->point(nx_-1,0)].dielectric(1.0);
                                sigx = pmlArr_[1].sigma(0,eps);
                                sigy = pmlArr_[0].sigma(0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[1].Dz_end_->point(0,0);
                                pmlArr_[1].Dz_end_->point(0,0) = c_dzd * pmlArr_[1].Dz_end_->point(0,0) + c_dzh * (-1.0*Hy_->point(nx_-1-1,0)) - (Hx_->point(nx_-1,0));
                                Ez_->point(nx_-1,0) = c_eze * Ez_->point(nx_-1,0) + c_ezd1 * pmlArr_[1].Dz_end_->point(0,0) - c_ezd0 * dzstore;
                                //Top Right
                                eps = objArr_[phys_Ez_->point(nx_-1,ny_-1)].dielectric(1.0);
                                sigx = pmlArr_[1].sigma(0,eps);
                                sigy = pmlArr_[0].sigma(0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[1].Dz_end_->point(0,ny_-1);
                                pmlArr_[1].Dz_end_->point(0,ny_-1) = c_dzd * pmlArr_[1].Dz_end_->point(0,ny_-1) + c_dzh * ((-1.0*Hy_->point(nx_-1-1,ny_-1)) - (-1.0*Hx_->point(nx_-1,ny_-1-1)));
                                Ez_->point(nx_-1,ny_-1) = c_eze * Ez_->point(nx_-1,ny_-1) + c_ezd1 * pmlArr_[1].Dz_end_->point(0,ny_-1) - c_ezd0 * dzstore;
                                //Top Left
                                eps = objArr_[phys_Ez_->point(0,ny_-1)].dielectric(1.0);
                                sigx = pmlArr_[1].sigma(0,eps);
                                sigy = pmlArr_[0].sigma(0,eps);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[1].Dz_->point(0,ny_-1);
                                pmlArr_[1].Dz_->point(0,ny_-1) = c_dzd * pmlArr_[1].Dz_->point(0,ny_-1) + c_dzh * ((Hy_->point(0,ny_-1)) - (-1.0*Hx_->point(0,ny_-1-1)));
                                Ez_->point(0,ny_-1) = c_eze * Ez_->point(0,ny_-1) + c_ezd1 * pmlArr_[1].Dz_->point(0,ny_-1) - c_ezd0 * dzstore;

                                for(int ii = 1; ii < pmlArr_[1].thickness(); ii++)
                                {
                                    kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                    sigz = 0.0;

                                    eps = objArr_[phys_Ez_->point(ii,0)].dielectric(1.0);
                                    sigx = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                    sigy = pmlArr_[0].sigma(0.0,eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    //Bot Left
                                    dzstore = pmlArr_[1].Dz_->point(ii,0);
                                    pmlArr_[1].Dz_->point(ii,0) = c_dzd * pmlArr_[1].Dz_->point(ii,0) + c_dzh * ((Hy_->point(ii,0)-Hy_->point(ii-1,0)) - (Hx_->point(ii,0)));
                                    Ez_->point(ii,0) = c_eze * Ez_->point(ii,0) + c_ezd1 * pmlArr_[1].Dz_->point(ii,0) - c_ezd0 * dzstore;
                                    //Bot Right
                                    eps = objArr_[phys_Ez_->point(nx_-1-ii,0)].dielectric(1.0);
                                    sigx = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                    sigy = pmlArr_[0].sigma(0.0,eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[1].Dz_end_->point(ii,0);
                                    pmlArr_[1].Dz_end_->point(ii,0) = c_dzd * pmlArr_[1].Dz_end_->point(ii,0) + c_dzh * ((Hy_->point(nx_-1-ii,0)-Hy_->point(nx_-1-ii-1,0)) - (Hx_->point(nx_-1-ii,0)));
                                    Ez_->point(nx_-1-ii,0) = c_eze * Ez_->point(nx_-1-ii,0) + c_ezd1 * pmlArr_[1].Dz_end_->point(ii,0) - c_ezd0 * dzstore;
                                    //Top Right
                                    eps = objArr_[phys_Ez_->point(nx_-1-ii,ny_-1)].dielectric(1.0);
                                    sigx = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                    sigy = pmlArr_[0].sigma(0.0,eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[1].Dz_end_->point(ii,ny_-1);
                                    pmlArr_[1].Dz_end_->point(ii,ny_-1) = c_dzd * pmlArr_[1].Dz_end_->point(ii,ny_-1) + c_dzh * ((Hy_->point(nx_-1-ii,ny_-1)-Hy_->point(nx_-1-ii-1,ny_-1)) - (-1.0*Hx_->point(nx_-1-ii,ny_-1-1)));
                                    Ez_->point(nx_-1-ii,ny_-1) = c_eze * Ez_->point(nx_-1-ii,ny_-1) + c_ezd1 * pmlArr_[1].Dz_end_->point(ii,ny_-1) - c_ezd0 * dzstore;
                                    //Top Left
                                    eps = objArr_[phys_Ez_->point(ii,ny_-1)].dielectric(1.0);
                                    sigx = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                    sigy = pmlArr_[0].sigma(0.0,eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[1].Dz_->point(ii,ny_-1);
                                    pmlArr_[1].Dz_->point(ii,ny_-1) = c_dzd * pmlArr_[1].Dz_->point(ii,ny_-1) + c_dzh * ((Hy_->point(ii,ny_-1)-Hy_->point(ii-1,ny_-1)) - (-1.0*Hx_->point(ii,ny_-1-1)));
                                    Ez_->point(ii,ny_-1) = c_eze * Ez_->point(ii,ny_-1) + c_ezd1 * pmlArr_[1].Dz_->point(ii,ny_-1) - c_ezd0 * dzstore;
                                }
                                for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                                {
                                    kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                    sigz = 0.0;

                                    eps = objArr_[phys_Ez_->point(0,jj)].dielectric(1.0);
                                    sigx = pmlArr_[1].sigma(0.0,eps);
                                    sigy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    //Bot Left
                                    dzstore = pmlArr_[1].Dz_->point(0,jj);
                                    pmlArr_[1].Dz_->point(0,jj) = c_dzd * pmlArr_[1].Dz_->point(0,jj) + c_dzh * ((Hy_->point(0,jj)) - (Hx_->point(0,jj)-Hx_->point(0,jj-1)));
                                    Ez_->point(0,jj) = c_eze * Ez_->point(0,jj) + c_ezd1 * pmlArr_[1].Dz_->point(0,jj) - c_ezd0 * dzstore;
                                    //Bot Right
                                    eps = objArr_[phys_Ez_->point(nx_-1,jj)].dielectric(1.0);
                                    sigx = pmlArr_[1].sigma(0.0,eps);
                                    sigy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[1].Dz_end_->point(0,jj);
                                    pmlArr_[1].Dz_end_->point(0,jj) = c_dzd * pmlArr_[1].Dz_end_->point(0,jj) + c_dzh * ((-1.0*Hy_->point(nx_-1-1,jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1)));
                                    Ez_->point(nx_-1,jj) = c_eze * Ez_->point(nx_-1,jj) + c_ezd1 * pmlArr_[1].Dz_end_->point(0,jj) - c_ezd0 * dzstore;
                                    //Top Right
                                    eps = objArr_[phys_Ez_->point(nx_-1,ny_-1-jj)].dielectric(1.0);
                                    sigx = pmlArr_[1].sigma(0.0,eps);
                                    sigy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[1].Dz_end_->point(0,ny_-1-jj);
                                    pmlArr_[1].Dz_end_->point(0,ny_-1-jj) = c_dzd * pmlArr_[1].Dz_end_->point(0,ny_-1-jj) + c_dzh * ((-1.0*Hy_->point(nx_-1-1,ny_-1-jj)) - (Hx_->point(nx_-1,ny_-1-jj)-Hx_->point(nx_-1,ny_-1-jj-1)));
                                    Ez_->point(nx_-1,ny_-1-jj) = c_eze * Ez_->point(nx_-1,ny_-1-jj) + c_ezd1 * pmlArr_[1].Dz_end_->point(0,ny_-1-jj) - c_ezd0 * dzstore;
                                    //Top Left
                                    eps = objArr_[phys_Ez_->point(0,ny_-1-jj)].dielectric(1.0);
                                    sigx = pmlArr_[1].sigma(0.0,eps);
                                    sigy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[1].Dz_->point(0,ny_-1-jj);
                                    pmlArr_[1].Dz_->point(0,ny_-1-jj) = c_dzd * pmlArr_[1].Dz_->point(0,ny_-1-jj) + c_dzh * ((Hy_->point(0,ny_-1-jj)) - (Hx_->point(0,ny_-1-jj)-Hx_->point(0,ny_-1-jj-1)));
                                    Ez_->point(0,ny_-1-jj) = c_eze * Ez_->point(0,ny_-1-jj) + c_ezd1 * pmlArr_[1].Dz_->point(0,ny_-1-jj) - c_ezd0 * dzstore;
                                }
                                for(int ii = 1; ii < pmlArr_[1].thickness(); ii++)
                                {
                                    for(int jj = 1; jj < pmlArr_[0].thickness(); jj++)
                                    {
                                        kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                        sigz = 0.0;

                                        eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0);
                                        sigx = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                        sigy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        //Bot Left
                                        dzstore = pmlArr_[1].Dz_->point(ii,jj);
                                        pmlArr_[1].Dz_->point(ii,jj) = c_dzd * pmlArr_[1].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                        Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[1].Dz_->point(ii,jj) - c_ezd0 * dzstore;
                                        //Bot Right
                                        eps = objArr_[phys_Ez_->point(nx_-1-ii,jj)].dielectric(1.0);
                                        sigx = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                        sigy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        dzstore = pmlArr_[1].Dz_end_->point(ii,jj);
                                        pmlArr_[1].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[1].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point(nx_-1-ii,jj)-Hy_->point(nx_-1-ii-1,jj)) - (Hx_->point(nx_-1-ii,jj)-Hx_->point(nx_-1-ii,jj-1)));
                                        Ez_->point(nx_-1-ii,jj) = c_eze * Ez_->point(nx_-1-ii,jj) + c_ezd1 * pmlArr_[1].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                        //Top Right
                                        eps = objArr_[phys_Ez_->point(nx_-1-ii,ny_-1-jj)].dielectric(1.0);
                                        sigx = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                        sigy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        dzstore = pmlArr_[1].Dz_end_->point(ii,ny_-1-jj);
                                        pmlArr_[1].Dz_end_->point(ii,ny_-1-jj) = c_dzd * pmlArr_[1].Dz_end_->point(ii,ny_-1-jj) + c_dzh * ((Hy_->point(nx_-1-ii,ny_-1-jj)-Hy_->point(nx_-1-ii-1,ny_-1-jj)) - (Hx_->point(nx_-1-ii,ny_-1-jj)-Hx_->point(nx_-1-ii,ny_-1-jj-1)));
                                        Ez_->point(nx_-1-ii,ny_-1-jj) = c_eze * Ez_->point(nx_-1-ii,ny_-1-jj) + c_ezd1 * pmlArr_[1].Dz_end_->point(ii,ny_-1-jj) - c_ezd0 * dzstore;
                                        //Top Left
                                        eps = objArr_[phys_Ez_->point(ii,ny_-1-jj)].dielectric(1.0);
                                        sigx = pmlArr_[1].sigma(static_cast<double>(ii),eps);
                                        sigy = pmlArr_[0].sigma(static_cast<double>(jj),eps);
                                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                        dzstore = pmlArr_[1].Dz_->point(ii,ny_-1-jj);
                                        pmlArr_[1].Dz_->point(ii,ny_-1-jj) = c_dzd * pmlArr_[1].Dz_->point(ii,ny_-1-jj) + c_dzh * ((Hy_->point(ii,ny_-1-jj)-Hy_->point(ii-1,ny_-1-jj)) - (Hx_->point(ii,ny_-1-jj)-Hx_->point(ii,ny_-1-jj-1)));
                                        Ez_->point(ii,ny_-1-jj) = c_eze * Ez_->point(ii,ny_-1-jj) + c_ezd1 * pmlArr_[1].Dz_->point(ii,ny_-1-jj) - c_ezd0 * dzstore;
                                    }
                                }
                            }
                        }
                        break;
                    }
                    case Z:
                        throw logic_error("z not implimented");
                        break;
                    default:
                        throw logic_error("hit default");
                        break;
                }
            }
        }
    }
    /*else
    {
        if(Ez_)
        {
            for(int kk = 0; kk < pmlArr_.size(); kk++)
            {
                switch(pmlArr_[kk].d())
                {
                    case X:
                    {
                        complex<double> dzstore(0.0,0.0);
                        if(yPML_ == 0)
                        {
                            for(int jj = 1; jj < ny_ - 1; jj++)
                            {
                                dzstore = pmlArr_[kk].Dz_->point(0,jj);
                                pmlArr_[kk].Dz_->point(0,jj) = pmlArr_[kk].c_dzd_0_ -> point(0,jj) * pmlArr_[kk].Dz_->point(0,jj) + pmlArr_[kk].c_dzh_0_ -> point(0,jj) * ((Hy_->point(0,jj)) - (Hx_->point(0,jj)-Hx_->point(0,jj-1))); // 0 is boundary condtion
                                Ez_->point(0,jj) = pmlArr_[kk].c_eze_0_ -> point(0,jj) * Ez_->point(0,jj) + pmlArr_[kk].c_ezd1_0_ -> point(0,jj) * pmlArr_[kk].Dz_->point(0,jj) - pmlArr_[kk].c_ezd0_0_ -> point(0,jj) * dzstore;

                                dzstore = pmlArr_[kk].Dz_end_->point(0,jj);
                                pmlArr_[kk].Dz_end_->point(0,jj) = pmlArr_[kk].c_dzd_n_ -> point(0,jj) * pmlArr_[kk].Dz_end_->point(0,jj) + pmlArr_[kk].c_dzh_n_ -> point(0,jj) * ((-1.0 *Hy_->point((nx_-2),jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1))); // 0 is boundary condition
                                Ez_->point(nx_-1,jj) = pmlArr_[kk].c_eze_n_ -> point(0,jj) * Ez_->point(nx_-1,jj) + pmlArr_[kk].c_ezd1_n_ -> point(0,jj) * pmlArr_[kk].Dz_end_->point(0,jj) - pmlArr_[kk].c_ezd0_n_ -> point(0,jj) * dzstore;
                            }

                            dzstore = pmlArr_[kk].Dz_->point(0,ny_-1);
                            pmlArr_[kk].Dz_->point(0,ny_-1) = pmlArr_[kk].c_dzd_0_ -> point(0,ny_-1) * pmlArr_[kk].Dz_->point(0,ny_-1) + pmlArr_[kk].c_dzh_0_ -> point(0,ny_-1) * ((Hy_->point(0,ny_-1)) - (-1.0*Hx_->point(0,ny_-2))); // 0 is boundary condtion
                            Ez_->point(0,ny_-1) = pmlArr_[kk].c_eze_0_ -> point(0,ny_-1) * Ez_->point(0,ny_-1) + pmlArr_[kk].c_ezd1_0_ -> point(0,ny_-1) * pmlArr_[kk].Dz_->point(0,ny_-1) - pmlArr_[kk].c_ezd0_0_ -> point(0,ny_-1) * dzstore;

                            dzstore = pmlArr_[kk].Dz_end_->point(0,ny_-1);
                            pmlArr_[kk].Dz_end_->point(0,ny_-1) = pmlArr_[kk].c_dzd_n_ -> point(0,ny_-1) * pmlArr_[kk].Dz_end_->point(0,ny_-1) + pmlArr_[kk].c_dzh_n_ -> point(0,ny_-1) * ((-1.0 *Hy_->point(nx_-2,ny_-1)) - (-1.0*Hx_->point(nx_-1,ny_-2))); // 0 is boundary condition
                            Ez_->point(nx_-1,ny_-1) = pmlArr_[kk].c_eze_n_ -> point(0,ny_-1) * Ez_->point(nx_-1,ny_-1) + pmlArr_[kk].c_ezd1_n_ -> point(0,ny_-1) * pmlArr_[kk].Dz_end_->point(0,ny_-1) - pmlArr_[kk].c_ezd0_n_ -> point(0,ny_-1) * dzstore;

                            dzstore = pmlArr_[kk].Dz_->point(0,0);
                            pmlArr_[kk].Dz_->point(0,0) = pmlArr_[kk].c_dzd_0_ -> point(0,0) * pmlArr_[kk].Dz_->point(0,0) + pmlArr_[kk].c_dzh_0_ -> point(0,0) * ((Hy_->point(0,0)) - (Hx_->point(0,0))); // 0 is boundary condtion
                            Ez_->point(0,0) = pmlArr_[kk].c_eze_0_ -> point(0,0) * Ez_->point(0,0) + pmlArr_[kk].c_ezd1_0_ -> point(0,0) * pmlArr_[kk].Dz_->point(0,0) - pmlArr_[kk].c_ezd0_0_ -> point(0,0) * dzstore;

                            dzstore = pmlArr_[kk].Dz_end_->point(0,0);
                            pmlArr_[kk].Dz_end_->point(0,0) = pmlArr_[kk].c_dzd_n_ -> point(0,0) * pmlArr_[kk].Dz_end_->point(0,0) + pmlArr_[kk].c_dzh_n_ -> point(0,0) * ((-1.0 *Hy_->point(nx_-2,0)) - (Hx_->point(nx_-1,0))); // 0 is boundary condition
                            Ez_->point(nx_-1,0) = pmlArr_[kk].c_eze_n_ -> point(0,0) * Ez_->point(nx_-1,0) + pmlArr_[kk].c_ezd1_n_ -> point(0,0) * pmlArr_[kk].Dz_end_->point(0,0) - pmlArr_[kk].c_ezd0_n_ -> point(0,0) * dzstore;

                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                //double kapx = pmlArr_[kk].kappa(ii);
                                //Kappas change throughout
                                for(int jj =  1; jj < ny_ - 1; jj ++)
                                {
                                    dzstore = pmlArr_[kk].Dz_->point(ii,jj);
                                    pmlArr_[kk].Dz_->point(ii,jj) = pmlArr_[kk].c_dzd_0_ -> point(ii,jj) * pmlArr_[kk].Dz_->point(ii,jj) + pmlArr_[kk].c_dzh_0_ -> point(ii,jj) * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                    Ez_->point(ii,jj) = pmlArr_[kk].c_eze_0_ -> point(ii,jj) * Ez_->point(ii,jj) + pmlArr_[kk].c_ezd1_0_ -> point(ii,jj) * pmlArr_[kk].Dz_->point(ii,jj) - pmlArr_[kk].c_ezd0_0_ -> point(ii,jj) * dzstore;

                                    dzstore = pmlArr_[kk].Dz_end_->point(ii,jj);
                                    pmlArr_[kk].Dz_end_->point(ii,jj) = pmlArr_[kk].c_dzd_n_ -> point(ii,jj) * pmlArr_[kk].Dz_end_->point(ii,jj) + pmlArr_[kk].c_dzh_n_ -> point(ii,jj) * ((Hy_->point(nx_-1-ii,jj)-Hy_->point(nx_-1-ii-1,jj)) - (Hx_->point(nx_-1-ii,jj)-Hx_->point(nx_-1-ii,jj-1)));
                                    Ez_->point(nx_-1-ii,jj) = pmlArr_[kk].c_eze_n_ -> point(ii,jj) * Ez_->point(nx_-1-ii,jj) + pmlArr_[kk].c_ezd1_n_ -> point(ii,jj) * pmlArr_[kk].Dz_end_->point(ii,jj) - pmlArr_[kk].c_ezd0_n_ -> point(ii,jj) * dzstore;
                                }
                                dzstore = pmlArr_[kk].Dz_->point(ii,ny_-1);
                                pmlArr_[kk].Dz_->point(ii,ny_-1) = pmlArr_[kk].c_dzd_0_ -> point(ii,ny_-1) * pmlArr_[kk].Dz_->point(ii,ny_-1) + pmlArr_[kk].c_dzh_0_ -> point(ii,ny_-1) * ((Hy_->point(ii,ny_-1)-Hy_->point(ii-1,ny_-1)) - (-1.0*Hx_->point(ii,ny_-1-1)));
                                Ez_->point(ii,ny_-1) = pmlArr_[kk].c_eze_0_ -> point(ii,ny_-1) * Ez_->point(ii,ny_-1) + pmlArr_[kk].c_ezd1_0_ -> point(ii,ny_-1) * pmlArr_[kk].Dz_->point(ii,ny_-1) - pmlArr_[kk].c_ezd0_0_ -> point(ii,ny_-1) * dzstore;

                                dzstore = pmlArr_[kk].Dz_end_->point(ii,ny_-1);
                                pmlArr_[kk].Dz_end_->point(ii,ny_-1) = pmlArr_[kk].c_dzd_n_ -> point(ii,ny_-1) * pmlArr_[kk].Dz_end_->point(ii,ny_-1) + pmlArr_[kk].c_dzh_n_ -> point(ii,ny_-1) * ((Hy_->point(nx_-1-ii,ny_-1)-Hy_->point(nx_-1-ii-1,ny_-1)) - (-1.0*Hx_->point(nx_-1-ii,ny_-1-1)));
                                Ez_->point(nx_-1-ii,ny_-1) = pmlArr_[kk].c_eze_n_ -> point(ii,ny_-1) * Ez_->point(nx_-1-ii,ny_-1) + pmlArr_[kk].c_ezd1_n_ -> point(ii,ny_-1) * pmlArr_[kk].Dz_end_->point(ii,ny_-1) - pmlArr_[kk].c_ezd0_n_ -> point(ii,ny_-1) * dzstore;

                                dzstore = pmlArr_[kk].Dz_->point(ii,0);
                                pmlArr_[kk].Dz_->point(ii,0) = pmlArr_[kk].c_dzd_0_ -> point(ii,0) * pmlArr_[kk].Dz_->point(ii,0) + pmlArr_[kk].c_dzh_0_ -> point(ii,0) * ((Hy_->point(ii,0)-Hy_->point(ii-1,0)) - (Hx_->point(ii,0)));
                                Ez_->point(ii,0) = pmlArr_[kk].c_eze_0_ -> point(ii,0) * Ez_->point(ii,0) + pmlArr_[kk].c_ezd1_0_ -> point(ii,0) * pmlArr_[kk].Dz_->point(ii,0) - pmlArr_[kk].c_ezd0_0_ -> point(ii,0) * dzstore;

                                dzstore = pmlArr_[kk].Dz_end_->point(ii,0);
                                pmlArr_[kk].Dz_end_->point(ii,0) = pmlArr_[kk].c_dzd_n_ -> point(ii,0) * pmlArr_[kk].Dz_end_->point(ii,0) + pmlArr_[kk].c_dzh_n_ -> point(ii,0) * ((Hy_->point(nx_-1-ii,0)-Hy_->point(nx_-1-ii-1,0)) - (Hx_->point(nx_-1-ii,0)));
                                Ez_->point(nx_-1-ii,0) = pmlArr_[kk].c_eze_n_ -> point(ii,0) * Ez_->point(nx_-1-ii,0) + pmlArr_[kk].c_ezd1_n_ -> point(ii,0) * pmlArr_[kk].Dz_end_->point(ii,0) - pmlArr_[kk].c_ezd0_n_ -> point(ii,0) * dzstore;
                            }
                        }
                        else
                        {
                            for(int jj = yPML_; jj < ny_ - yPML_; jj++)
                            {
                                dzstore = pmlArr_[kk].Dz_->point(0,jj);
                                pmlArr_[kk].Dz_->point(0,jj) = pmlArr_[kk].c_dzd_0_ -> point(0,jj) * pmlArr_[kk].Dz_->point(0,jj) + pmlArr_[kk].c_dzh_0_ -> point(0,jj) * ((Hy_->point(0,jj)) - (Hx_->point(0,jj)-Hx_->point(0,jj-1))); // 0 is boundary condtion
                                Ez_->point(0,jj) = pmlArr_[kk].c_eze_0_ -> point(0,jj) * Ez_->point(0,jj) + pmlArr_[kk].c_ezd1_0_ -> point(0,jj) * pmlArr_[kk].Dz_->point(0,jj) - pmlArr_[kk].c_ezd0_0_ -> point(0,jj) * dzstore;

                                dzstore = pmlArr_[kk].Dz_end_->point(0,jj);
                                pmlArr_[kk].Dz_end_->point(0,jj) = pmlArr_[kk].c_dzd_n_ -> point(0,jj) * pmlArr_[kk].Dz_end_->point(0,jj) + pmlArr_[kk].c_dzh_n_ -> point(0,jj) * ((-1.0 *Hy_->point((nx_-2),jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1))); // 0 is boundary condition
                                Ez_->point(nx_-1,jj) = pmlArr_[kk].c_eze_n_ -> point(0,jj) * Ez_->point(nx_-1,jj) + pmlArr_[kk].c_ezd1_n_ -> point(0,jj) * pmlArr_[kk].Dz_end_->point(0,jj) - pmlArr_[kk].c_ezd0_n_ -> point(0,jj) * dzstore;
                            }
                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                for(int jj =  yPML_; jj < ny_ - yPML_; jj ++)
                                {
                                    dzstore = pmlArr_[kk].Dz_->point(ii,jj);
                                    pmlArr_[kk].Dz_->point(ii,jj) = pmlArr_[kk].c_dzd_0_ -> point(ii,jj) * pmlArr_[kk].Dz_->point(ii,jj) + pmlArr_[kk].c_dzh_0_ -> point(ii,jj) * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                    Ez_->point(ii,jj) = pmlArr_[kk].c_eze_0_ -> point(ii,jj) * Ez_->point(ii,jj) + pmlArr_[kk].c_ezd1_0_ -> point(ii,jj) * pmlArr_[kk].Dz_->point(ii,jj) - pmlArr_[kk].c_ezd0_0_ -> point(ii,jj) * dzstore;

                                    dzstore = pmlArr_[kk].Dz_end_->point(ii,jj);
                                    pmlArr_[kk].Dz_end_->point(ii,jj) = pmlArr_[kk].c_dzd_n_ -> point(ii,jj) * pmlArr_[kk].Dz_end_->point(ii,jj) + pmlArr_[kk].c_dzh_n_ -> point(ii,jj) * ((Hy_->point(nx_-1-ii,jj)-Hy_->point(nx_-1-ii-1,jj)) - (Hx_->point(nx_-1-ii,jj)-Hx_->point(nx_-1-ii,jj-1)));
                                    Ez_->point(nx_-1-ii,jj) = pmlArr_[kk].c_eze_n_ -> point(ii,jj) * Ez_->point(nx_-1-ii,jj) + pmlArr_[kk].c_ezd1_n_ -> point(ii,jj) * pmlArr_[kk].Dz_end_->point(ii,jj) - pmlArr_[kk].c_ezd0_n_ -> point(ii,jj) * dzstore;
                                }
                            }
                            if(kk == 0)
                            {
                                //Bot Left
                                dzstore = pmlArr_[0].Dz_->point(0,0);
                                pmlArr_[0].Dz_->point(0,0) = pmlArr_[0].c_dzd_0_ -> point(0,0) * pmlArr_[0].Dz_->point(0,0) + pmlArr_[0].c_dzh_0_ -> point(0,0) * ((Hy_->point(0,0)) - (Hx_->point(0,0)));
                                Ez_->point(0,0) = pmlArr_[0].c_eze_0_ -> point(0,0) * Ez_->point(0,0) + pmlArr_[0].c_ezd1_0_ -> point(0,0) * pmlArr_[0].Dz_->point(0,0) - pmlArr_[0].c_ezd0_0_ -> point(0,0) * dzstore;
                                //Bot Right
                                dzstore = pmlArr_[0].Dz_end_->point(0,0);
                                pmlArr_[0].Dz_end_->point(0,0) = pmlArr_[0].c_dzd_n_ -> point(0,0) * pmlArr_[0].Dz_end_->point(0,0) + pmlArr_[0].c_dzh_n_ -> point(0,0) * (-1.0*Hy_->point(nx_-1-1,0)) - (Hx_->point(nx_-1,0));
                                Ez_->point(nx_-1,0) = pmlArr_[0].c_eze_n_ -> point(0,0) * Ez_->point(nx_-1,0) + pmlArr_[0].c_ezd1_n_ -> point(0,0) * pmlArr_[0].Dz_end_->point(0,0) - pmlArr_[0].c_ezd0_n_ -> point(0,0) * dzstore;
                                //Top Right
                                dzstore = pmlArr_[0].Dz_end_->point(0,ny_-1);
                                pmlArr_[0].Dz_end_->point(0,ny_-1) = pmlArr_[0].c_dzd_n_ -> point(0,ny_-1) * pmlArr_[0].Dz_end_->point(0,ny_-1) + pmlArr_[0].c_dzh_n_ -> point(0,ny_-1) * ((-1.0*Hy_->point(nx_-1-1,ny_-1)) - (-1.0*Hx_->point(nx_-1,ny_-1-1)));
                                Ez_->point(nx_-1,ny_-1) = pmlArr_[0].c_eze_n_ -> point(0,ny_-1) * Ez_->point(nx_-1,ny_-1) + pmlArr_[0].c_ezd1_n_ -> point(0,ny_-1) * pmlArr_[0].Dz_end_->point(0,ny_-1) - pmlArr_[0].c_ezd0_n_ -> point(0,ny_-1) * dzstore;
                                //Top Left
                                dzstore = pmlArr_[0].Dz_->point(0,ny_-1);
                                pmlArr_[0].Dz_->point(0,ny_-1) = pmlArr_[0].c_dzd_0_ -> point(0,ny_-1) * pmlArr_[0].Dz_->point(0,ny_-1) + pmlArr_[0].c_dzh_0_ -> point(0,ny_-1) * ((Hy_->point(0,ny_-1)) - (-1.0*Hx_->point(0,ny_-1-1)));
                                Ez_->point(0,ny_-1) = pmlArr_[0].c_eze_0_ -> point(0,ny_-1) * Ez_->point(0,ny_-1) + pmlArr_[0].c_ezd1_0_ -> point(0,ny_-1) * pmlArr_[0].Dz_->point(0,ny_-1) - pmlArr_[0].c_ezd0_0_ -> point(0,ny_-1) * dzstore;

                                for(int ii = 1; ii < pmlArr_[0].thickness(); ii++)
                                {
                                    //Bot Left
                                    dzstore = pmlArr_[0].Dz_->point(ii,0);
                                    pmlArr_[0].Dz_->point(ii,0) = pmlArr_[0].c_dzd_0_ -> point(ii,0) * pmlArr_[0].Dz_->point(ii,0) + pmlArr_[0].c_dzh_0_ -> point(ii,0) * ((Hy_->point(ii,0)-Hy_->point(ii-1,0)) - (Hx_->point(ii,0)));
                                    Ez_->point(ii,0) = pmlArr_[0].c_eze_0_ -> point(ii,0) * Ez_->point(ii,0) + pmlArr_[0].c_ezd1_0_ -> point(ii,0) * pmlArr_[0].Dz_->point(ii,0) - pmlArr_[0].c_ezd0_0_ -> point(ii,0) * dzstore;
                                    //Bot Right
                                    dzstore = pmlArr_[0].Dz_end_->point(ii,0);
                                    pmlArr_[0].Dz_end_->point(ii,0) = pmlArr_[0].c_dzd_n_ -> point(ii,0) * pmlArr_[0].Dz_end_->point(ii,0) + pmlArr_[0].c_dzh_n_ -> point(ii,0) * ((Hy_->point(nx_-1-ii,0)-Hy_->point(nx_-1-ii-1,0)) - (Hx_->point(nx_-1-ii,0)));
                                    Ez_->point(nx_-1-ii,0) = pmlArr_[0].c_eze_n_ -> point(ii,0) * Ez_->point(nx_-1-ii,0) + pmlArr_[0].c_ezd1_n_ -> point(ii,0) * pmlArr_[0].Dz_end_->point(ii,0) - pmlArr_[0].c_ezd0_n_ -> point(ii,0) * dzstore;
                                    //Top Right
                                    dzstore = pmlArr_[0].Dz_end_->point(ii,ny_-1);
                                    pmlArr_[0].Dz_end_->point(ii,ny_-1) = pmlArr_[0].c_dzd_n_ -> point(ii,ny_-1) * pmlArr_[0].Dz_end_->point(ii,ny_-1) + pmlArr_[0].c_dzh_n_ -> point(ii,ny_-1) * ((Hy_->point(nx_-1-ii,ny_-1)-Hy_->point(nx_-1-ii-1,ny_-1)) - (-1.0*Hx_->point(nx_-1-ii,ny_-1-1)));
                                    Ez_->point(nx_-1-ii,ny_-1) = pmlArr_[0].c_eze_n_ -> point(ii,ny_-1) * Ez_->point(nx_-1-ii,ny_-1) + pmlArr_[0].c_ezd1_n_ -> point(ii,ny_-1) * pmlArr_[0].Dz_end_->point(ii,ny_-1) - pmlArr_[0].c_ezd0_n_ -> point(ii,ny_-1) * dzstore;
                                    //Top Left
                                    dzstore = pmlArr_[0].Dz_->point(ii,ny_-1);
                                    pmlArr_[0].Dz_->point(ii,ny_-1) = pmlArr_[0].c_dzd_0_ -> point(ii,ny_-1) * pmlArr_[0].Dz_->point(ii,ny_-1) + pmlArr_[0].c_dzh_0_ -> point(ii,ny_-1) * ((Hy_->point(ii,ny_-1)-Hy_->point(ii-1,ny_-1)) - (-1.0*Hx_->point(ii,ny_-1-1)));
                                    Ez_->point(ii,ny_-1) = pmlArr_[0].c_eze_0_ -> point(ii,ny_-1) * Ez_->point(ii,ny_-1) + pmlArr_[0].c_ezd1_0_ -> point(ii,ny_-1) * pmlArr_[0].Dz_->point(ii,ny_-1) - pmlArr_[0].c_ezd0_0_ -> point(ii,ny_-1) * dzstore;
                                }
                                for(int jj = 1; jj < pmlArr_[0].thickness(); jj++)
                                {
                                    //Bot Left
                                    dzstore = pmlArr_[0].Dz_->point(0,jj);
                                    pmlArr_[0].Dz_->point(0,jj) = pmlArr_[0].c_dzd_0_ -> point(0,jj) * pmlArr_[0].Dz_->point(0,jj) + pmlArr_[0].c_dzh_0_ -> point(0,jj) * ((Hy_->point(0,jj)) - (Hx_->point(0,jj)-Hx_->point(0,jj-1)));
                                    Ez_->point(0,jj) = pmlArr_[0].c_eze_0_ -> point(0,jj) * Ez_->point(0,jj) + pmlArr_[0].c_ezd1_0_ -> point(0,jj) * pmlArr_[0].Dz_->point(0,jj) - pmlArr_[0].c_ezd0_0_ -> point(0,jj) * dzstore;
                                    //Bot Right
                                    dzstore = pmlArr_[0].Dz_end_->point(0,jj);
                                    pmlArr_[0].Dz_end_->point(0,jj) = pmlArr_[0].c_dzd_n_ -> point(0,jj) * pmlArr_[0].Dz_end_->point(0,jj) + pmlArr_[0].c_dzh_n_ -> point(0,jj) * ((-1.0*Hy_->point(nx_-1-1,jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1)));
                                    Ez_->point(nx_-1,jj) = pmlArr_[0].c_eze_n_ -> point(0,jj) * Ez_->point(nx_-1,jj) + pmlArr_[0].c_ezd1_n_ -> point(0,jj) * pmlArr_[0].Dz_end_->point(0,jj) - pmlArr_[0].c_ezd0_n_ -> point(0,jj) * dzstore;
                                    //Top Right
                                    dzstore = pmlArr_[0].Dz_end_->point(0,ny_-1-jj);
                                    pmlArr_[0].Dz_end_->point(0,ny_-1-jj) = pmlArr_[0].c_dzd_n_ -> point(0,ny_-1-jj) * pmlArr_[0].Dz_end_->point(0,ny_-1-jj) + pmlArr_[0].c_dzh_n_ -> point(0,ny_-1-jj) * ((-1.0*Hy_->point(nx_-1-1,ny_-1-jj)) - (Hx_->point(nx_-1,ny_-1-jj)-Hx_->point(nx_-1,ny_-1-jj-1)));
                                    Ez_->point(nx_-1,ny_-1-jj) = pmlArr_[0].c_eze_n_ -> point(0,ny_-1-jj) * Ez_->point(nx_-1,ny_-1-jj) + pmlArr_[0].c_ezd1_n_ -> point(0,ny_-1-jj) * pmlArr_[0].Dz_end_->point(0,ny_-1-jj) - pmlArr_[0].c_ezd0_n_ -> point(0,ny_-1-jj) * dzstore;
                                    //Top Left
                                    dzstore = pmlArr_[0].Dz_->point(0,ny_-1-jj);
                                    pmlArr_[0].Dz_->point(0,ny_-1-jj) = pmlArr_[0].c_dzd_0_ -> point(0,ny_-1-jj) * pmlArr_[0].Dz_->point(0,ny_-1-jj) + pmlArr_[0].c_dzh_0_ -> point(0,ny_-1-jj) * ((Hy_->point(0,ny_-1-jj)) - (Hx_->point(0,ny_-1-jj)-Hx_->point(0,ny_-1-jj-1)));
                                    Ez_->point(0,ny_-1-jj) = pmlArr_[0].c_eze_0_ -> point(0,ny_-1-jj) * Ez_->point(0,ny_-1-jj) + pmlArr_[0].c_ezd1_0_ -> point(0,ny_-1-jj) * pmlArr_[0].Dz_->point(0,ny_-1-jj) - pmlArr_[0].c_ezd0_0_ -> point(0,ny_-1-jj) * dzstore;
                                }
                                for(int ii = 1; ii < pmlArr_[0].thickness(); ii++)
                                {
                                    for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                                    {
                                        //Bot Left
                                        dzstore = pmlArr_[0].Dz_->point(ii,jj);
                                        pmlArr_[0].Dz_->point(ii,jj) = pmlArr_[0].c_dzd_0_ -> point(ii,jj) * pmlArr_[0].Dz_->point(ii,jj) + pmlArr_[0].c_dzh_0_ -> point(ii,jj) * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                        Ez_->point(ii,jj) = pmlArr_[0].c_eze_0_ -> point(ii,jj) * Ez_->point(ii,jj) + pmlArr_[0].c_ezd1_0_ -> point(ii,jj) * pmlArr_[0].Dz_->point(ii,jj) - pmlArr_[0].c_ezd0_0_ -> point(ii,jj) * dzstore;
                                        //Bot Right
                                        dzstore = pmlArr_[0].Dz_end_->point(ii,jj);
                                        pmlArr_[0].Dz_end_->point(ii,jj) = pmlArr_[0].c_dzd_n_ -> point(ii,jj) * pmlArr_[0].Dz_end_->point(ii,jj) + pmlArr_[0].c_dzh_n_ -> point(ii,jj) * ((Hy_->point(nx_-1-ii,jj)-Hy_->point(nx_-1-ii-1,jj)) - (Hx_->point(nx_-1-ii,jj)-Hx_->point(nx_-1-ii,jj-1)));
                                        Ez_->point(nx_-1-ii,jj) = pmlArr_[0].c_eze_n_ -> point(ii,jj) * Ez_->point(nx_-1-ii,jj) + pmlArr_[0].c_ezd1_n_ -> point(ii,jj) * pmlArr_[0].Dz_end_->point(ii,jj) - pmlArr_[0].c_ezd0_n_ -> point(ii,jj) * dzstore;
                                        //Top Right
                                        dzstore = pmlArr_[0].Dz_end_->point(ii,ny_-1-jj);
                                        pmlArr_[0].Dz_end_->point(ii,ny_-1-jj) = pmlArr_[0].c_dzd_n_ -> point(ii,ny_-1-jj) * pmlArr_[0].Dz_end_->point(ii,ny_-1-jj) + pmlArr_[0].c_dzh_n_ -> point(ii,ny_-1-jj) * ((Hy_->point(nx_-1-ii,ny_-1-jj)-Hy_->point(nx_-1-ii-1,ny_-1-jj)) - (Hx_->point(nx_-1-ii,ny_-1-jj)-Hx_->point(nx_-1-ii,ny_-1-jj-1)));
                                        Ez_->point(nx_-1-ii,ny_-1-jj) = pmlArr_[0].c_eze_n_ -> point(ii,ny_-1-jj) * Ez_->point(nx_-1-ii,ny_-1-jj) + pmlArr_[0].c_ezd1_n_ -> point(ii,ny_-1-jj) * pmlArr_[0].Dz_end_->point(ii,ny_-1-jj) - pmlArr_[0].c_ezd0_n_ -> point(ii,ny_-1-jj) * dzstore;
                                        //Top Left
                                        dzstore = pmlArr_[0].Dz_->point(ii,ny_-1-jj);
                                        pmlArr_[0].Dz_->point(ii,ny_-1-jj) = pmlArr_[0].c_dzd_0_ -> point(ii,ny_-1-jj) * pmlArr_[0].Dz_->point(ii,ny_-1-jj) + pmlArr_[0].c_dzh_0_ -> point(ii,ny_-1-jj) * ((Hy_->point(ii,ny_-1-jj)-Hy_->point(ii-1,ny_-1-jj)) - (Hx_->point(ii,ny_-1-jj)-Hx_->point(ii,ny_-1-jj-1)));
                                        Ez_->point(ii,ny_-1-jj) = pmlArr_[0].c_eze_0_ -> point(ii,ny_-1-jj) * Ez_->point(ii,ny_-1-jj) + pmlArr_[0].c_ezd1_0_ -> point(ii,ny_-1-jj) * pmlArr_[0].Dz_->point(ii,ny_-1-jj) - pmlArr_[0].c_ezd0_0_ -> point(ii,ny_-1-jj) * dzstore;
                                    }
                                }
                            }
                        }
                        break;
                    }
                    case Y:
                    {
                        complex<double> dzstore(0.0,0.0);
                        if(xPML_ == 0)
                        {
                            for(int jj = 1; jj < nx_ - 1; jj++)
                            {

                                dzstore = pmlArr_[kk].Dz_->point(jj,0);
                                pmlArr_[kk].Dz_->point(jj,0) = pmlArr_[kk].c_dzd_0_ -> point(jj,0) * pmlArr_[kk].Dz_->point(jj,0) + pmlArr_[kk].c_dzh_0_ -> point(jj,0) * ((Hy_->point(jj,0)-Hy_->point(jj-1,0)) - (Hx_->point(jj,0))); // 0 is boundary condition
                                Ez_->point(jj,0) = pmlArr_[kk].c_eze_0_ -> point(jj,0) * Ez_->point(jj,0) + pmlArr_[kk].c_ezd1_0_ -> point(jj,0) * pmlArr_[kk].Dz_->point(jj,0) - pmlArr_[kk].c_ezd0_0_ -> point(jj,0) * dzstore;

                                dzstore = pmlArr_[kk].Dz_end_->point(jj,0);
                                pmlArr_[kk].Dz_end_->point(jj,0) = pmlArr_[kk].c_dzd_n_ -> point(jj,0) * pmlArr_[kk].Dz_end_->point(jj,0) + pmlArr_[kk].c_dzh_n_ -> point(jj,0) * ((Hy_->point(jj,ny_-1) - Hy_->point(jj - 1,ny_-1)) - (-1.0 * Hx_->point(jj,ny_-1-1))); // 0 is boundary condtion
                                Ez_->point(jj,ny_-1) = pmlArr_[kk].c_eze_n_ -> point(jj,0) * Ez_->point(jj,ny_-1) + pmlArr_[kk].c_ezd1_n_ -> point(jj,0) * pmlArr_[kk].Dz_end_->point(jj,0) - pmlArr_[kk].c_ezd0_n_ -> point(jj,0) * dzstore;
                            }
                            dzstore = pmlArr_[kk].Dz_->point(nx_-1,0);
                            pmlArr_[kk].Dz_->point(nx_-1,0) = pmlArr_[kk].c_dzd_0_ -> point(nx_-1,0) * pmlArr_[kk].Dz_->point(nx_-1,0) + pmlArr_[kk].c_dzh_0_ -> point(nx_-1,0) * ((-1.0*Hy_->point(nx_-1-1,0)) - (Hx_->point(nx_-1,0))); // 0 is boundary condition
                            Ez_->point(nx_-1,0) = pmlArr_[kk].c_eze_0_ -> point(nx_-1,0) * Ez_->point(nx_-1,0) + pmlArr_[kk].c_ezd1_0_ -> point(nx_-1,0) * pmlArr_[kk].Dz_->point(nx_-1,0) - pmlArr_[kk].c_ezd0_0_ -> point(nx_-1,0) * dzstore;

                            dzstore = pmlArr_[kk].Dz_end_->point(nx_-1,0);
                            pmlArr_[kk].Dz_end_->point(nx_-1,0) = pmlArr_[kk].c_dzd_n_ -> point(nx_-1,0) * pmlArr_[kk].Dz_end_->point(nx_-1,0) + pmlArr_[kk].c_dzh_n_ -> point(nx_-1,0) * ((-1.0 * Hy_->point(nx_-1 - 1,ny_-1)) - (-1.0 * Hx_->point(nx_-1,ny_-1-1))); // 0 is boundary condtion
                            Ez_->point(nx_-1,ny_-1) = pmlArr_[kk].c_eze_n_ -> point(nx_-1,0) * Ez_->point(nx_-1,ny_-1) + pmlArr_[kk].c_ezd1_n_ -> point(nx_-1,0) * pmlArr_[kk].Dz_end_->point(nx_-1,0) - pmlArr_[kk].c_ezd0_n_ -> point(nx_-1,0) * dzstore;

                            dzstore = pmlArr_[kk].Dz_->point(0,0);
                            pmlArr_[kk].Dz_->point(0,0) = pmlArr_[kk].c_dzd_0_ -> point(0,0) * pmlArr_[kk].Dz_->point(0,0) + pmlArr_[kk].c_dzh_0_ -> point(0,0) * ((Hy_->point(0,0)) - (Hx_->point(0,0))); // 0 is boundary condition
                            Ez_->point(0,0) = pmlArr_[kk].c_eze_0_ -> point(0,0) * Ez_->point(0,0) + pmlArr_[kk].c_ezd1_0_ -> point(0,0) * pmlArr_[kk].Dz_->point(0,0) - pmlArr_[kk].c_ezd0_0_ -> point(0,0) * dzstore;

                            dzstore = pmlArr_[kk].Dz_end_->point(0,0);
                            pmlArr_[kk].Dz_end_->point(0,0) = pmlArr_[kk].c_dzd_n_ -> point(0,0) * pmlArr_[kk].Dz_end_->point(0,0) + pmlArr_[kk].c_dzh_n_ -> point(0,0) * ((Hy_->point(0,ny_-1)) - (-1.0 * Hx_->point(0,ny_-1-1))); // 0 is boundary condtion
                            Ez_->point(0,ny_-1) = pmlArr_[kk].c_eze_n_ -> point(0,0) * Ez_->point(0,ny_-1) + pmlArr_[kk].c_ezd1_n_ -> point(0,0) * pmlArr_[kk].Dz_end_->point(0,0) - pmlArr_[kk].c_ezd0_n_ -> point(0,0) * dzstore;

                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                for(int jj =  1; jj < nx_ - 1; jj ++)
                                {
                                    dzstore = pmlArr_[kk].Dz_->point(jj,ii);
                                    pmlArr_[kk].Dz_->point(jj,ii) = pmlArr_[kk].c_dzd_0_ -> point(jj,ii) * pmlArr_[kk].Dz_->point(jj,ii) + pmlArr_[kk].c_dzh_0_ -> point(jj,ii) * ((Hy_->point(jj,ii)-Hy_->point(jj-1,ii)) - (Hx_->point(jj,ii)-Hx_->point(jj,ii-1)));
                                    Ez_->point(jj,ii) = pmlArr_[kk].c_eze_0_ -> point(jj,ii) * Ez_->point(jj,ii) + pmlArr_[kk].c_ezd1_0_ -> point(jj,ii) * pmlArr_[kk].Dz_->point(jj,ii) - pmlArr_[kk].c_ezd0_0_ -> point(jj,ii) * dzstore;

                                    dzstore = pmlArr_[kk].Dz_end_->point(jj,ii);
                                    pmlArr_[kk].Dz_end_->point(jj,ii) = pmlArr_[kk].c_dzd_n_ -> point(jj,ii) * pmlArr_[kk].Dz_end_->point(jj,ii) + pmlArr_[kk].c_dzh_n_ -> point(jj,ii) * ((Hy_->point(jj,ny_-1-ii) - Hy_->point(jj - 1,ny_-1-ii)) - (Hx_->point(jj,ny_-1-ii) - Hx_->point(jj,ny_-1-ii-1)));
                                    Ez_->point(jj,ny_-1-ii) = pmlArr_[kk].c_eze_n_ -> point(jj,ii) * Ez_->point(jj,ny_-1-ii) + pmlArr_[kk].c_ezd1_n_ -> point(jj,ii) * pmlArr_[kk].Dz_end_->point(jj,ii) - pmlArr_[kk].c_ezd0_n_ -> point(jj,ii) * dzstore;
                                }
                                dzstore = pmlArr_[kk].Dz_->point(nx_-1,ii);
                                pmlArr_[kk].Dz_->point(nx_-1,ii) = pmlArr_[kk].c_dzd_0_ -> point(nx_-1,ii) * pmlArr_[kk].Dz_->point(nx_-1,ii) + pmlArr_[kk].c_dzh_0_ -> point(nx_-1,ii) * ((-1.0*Hy_->point(nx_-1-1,ii)) - (Hx_->point(nx_-1,ii)-Hx_->point(nx_-1,ii-1)));
                                Ez_->point(nx_-1,ii) = pmlArr_[kk].c_eze_0_ -> point(nx_-1,ii) * Ez_->point(nx_-1,ii) + pmlArr_[kk].c_ezd1_0_ -> point(nx_-1,ii) * pmlArr_[kk].Dz_->point(nx_-1,ii) - pmlArr_[kk].c_ezd0_0_ -> point(nx_-1,ii) * dzstore;

                                dzstore = pmlArr_[kk].Dz_end_->point(nx_-1,ii);
                                pmlArr_[kk].Dz_end_->point(nx_-1,ii) = pmlArr_[kk].c_dzd_n_ -> point(nx_-1,ii) * pmlArr_[kk].Dz_end_->point(nx_-1,ii) + pmlArr_[kk].c_dzh_n_ -> point(nx_-1,ii) * ((-1.0 * Hy_->point(nx_-1-1,ny_-1-ii)) - (Hx_->point(nx_-1,ny_-1-ii) - Hx_->point(nx_-1,ny_-1-ii-1)));
                                Ez_->point(nx_-1,ny_-1-ii) = pmlArr_[kk].c_eze_n_ -> point(nx_-1,ii) * Ez_->point(nx_-1,ny_-1-ii) + pmlArr_[kk].c_ezd1_n_ -> point(nx_-1,ii) * pmlArr_[kk].Dz_end_->point(nx_-1,ii) - pmlArr_[kk].c_ezd0_n_ -> point(nx_-1,ii) * dzstore;

                                dzstore = pmlArr_[kk].Dz_->point(0,ii);
                                pmlArr_[kk].Dz_->point(0,ii) = pmlArr_[kk].c_dzd_0_ -> point(0,ii) * pmlArr_[kk].Dz_->point(0,ii) + pmlArr_[kk].c_dzh_0_ -> point(0,ii) * ((Hy_->point(0,ii)) - (Hx_->point(0,ii)-Hx_->point(0,ii-1)));
                                Ez_->point(0,ii) = pmlArr_[kk].c_eze_0_ -> point(0,ii) * Ez_->point(0,ii) + pmlArr_[kk].c_ezd1_0_ -> point(0,ii) * pmlArr_[kk].Dz_->point(0,ii) - pmlArr_[kk].c_ezd0_0_ -> point(0,ii) * dzstore;

                                dzstore = pmlArr_[kk].Dz_end_->point(0,ii);
                                pmlArr_[kk].Dz_end_->point(0,ii) = pmlArr_[kk].c_dzd_n_ -> point(0,ii) * pmlArr_[kk].Dz_end_->point(0,ii) + pmlArr_[kk].c_dzh_n_ -> point(0,ii) * ((Hy_->point(0,ny_-1-ii)) - (Hx_->point(0,ny_-1-ii) - Hx_->point(0,ny_-1-ii-1)));
                                Ez_->point(0,ny_-1-ii) = pmlArr_[kk].c_eze_n_ -> point(0,ii) * Ez_->point(0,ny_-1-ii) + pmlArr_[kk].c_ezd1_n_ -> point(0,ii) * pmlArr_[kk].Dz_end_->point(0,ii) - pmlArr_[kk].c_ezd0_n_ -> point(0,ii) * dzstore;
                            }
                        }
                        else
                        {
                            for(int jj = xPML_; jj < nx_ - xPML_; jj++)
                            {
                                dzstore = pmlArr_[kk].Dz_->point(jj,0);
                                pmlArr_[kk].Dz_->point(jj,0) = pmlArr_[kk].c_dzd_0_ -> point(jj,0) * pmlArr_[kk].Dz_->point(jj,0) + pmlArr_[kk].c_dzh_0_ -> point(jj,0) * ((Hy_->point(jj,0)-Hy_->point(jj-1,0)) - (Hx_->point(jj,0))); // 0 is boundary condition
                                Ez_->point(jj,0) = pmlArr_[kk].c_eze_0_ -> point(jj,0) * Ez_->point(jj,0) + pmlArr_[kk].c_ezd1_0_ -> point(jj,0) * pmlArr_[kk].Dz_->point(jj,0) - pmlArr_[kk].c_ezd0_0_ -> point(jj,0) * dzstore;

                                dzstore = pmlArr_[kk].Dz_end_->point(jj,0);
                                pmlArr_[kk].Dz_end_->point(jj,0) = pmlArr_[kk].c_dzd_n_ -> point(jj,0) * pmlArr_[kk].Dz_end_->point(jj,0) + pmlArr_[kk].c_dzh_n_ -> point(jj,0) * ((Hy_->point(jj,ny_-1) - Hy_->point(jj - 1,ny_-1)) - (-1.0 * Hx_->point(jj,ny_-1-1))); // 0 is boundary condtion
                                Ez_->point(jj,ny_-1) = pmlArr_[kk].c_eze_n_ -> point(jj,0) * Ez_->point(jj,ny_-1) + pmlArr_[kk].c_ezd1_n_ -> point(jj,0) * pmlArr_[kk].Dz_end_->point(jj,0) - pmlArr_[kk].c_ezd0_n_ -> point(jj,0) * dzstore;
                            }
                            for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                            {
                                for(int jj =  xPML_; jj < nx_ - xPML_; jj ++)
                                {
                                    dzstore = pmlArr_[kk].Dz_->point(jj,ii);
                                    pmlArr_[kk].Dz_->point(jj,ii) = pmlArr_[kk].c_dzd_0_ -> point(jj,ii) * pmlArr_[kk].Dz_->point(jj,ii) + pmlArr_[kk].c_dzh_0_ -> point(jj,ii) * ((Hy_->point(jj,ii)-Hy_->point(jj-1,ii)) - (Hx_->point(jj,ii)-Hx_->point(jj,ii-1)));
                                    Ez_->point(jj,ii) = pmlArr_[kk].c_eze_0_ -> point(jj,ii) * Ez_->point(jj,ii) + pmlArr_[kk].c_ezd1_0_ -> point(jj,ii) * pmlArr_[kk].Dz_->point(jj,ii) - pmlArr_[kk].c_ezd0_0_ -> point(jj,ii) * dzstore;

                                    dzstore = pmlArr_[kk].Dz_end_->point(jj,ii);
                                    pmlArr_[kk].Dz_end_->point(jj,ii) = pmlArr_[kk].c_dzd_n_ -> point(jj,ii) * pmlArr_[kk].Dz_end_->point(jj,ii) + pmlArr_[kk].c_dzh_n_ -> point(jj,ii) * ((Hy_->point(jj,ny_-1-ii) - Hy_->point(jj - 1,ny_-1-ii)) - (Hx_->point(jj,ny_-1-ii) - Hx_->point(jj,ny_-1-ii-1)));
                                    Ez_->point(jj,ny_-1-ii) = pmlArr_[kk].c_eze_n_ -> point(jj,ii) * Ez_->point(jj,ny_-1-ii) + pmlArr_[kk].c_ezd1_n_ -> point(jj,ii) * pmlArr_[kk].Dz_end_->point(jj,ii) - pmlArr_[kk].c_ezd0_n_ -> point(jj,ii) * dzstore;
                                }
                            }
                            if(kk == 0)
                            {
                                //Bot Left
                                dzstore = pmlArr_[1].Dz_->point(0,0);
                                pmlArr_[1].Dz_->point(0,0) = pmlArr_[1].c_dzd_0_ -> point(0,0) * pmlArr_[1].Dz_->point(0,0) + pmlArr_[1].c_dzh_0_ -> point(0,0) * ((Hy_->point(0,0)) - (Hx_->point(0,0)));
                                Ez_->point(0,0) = pmlArr_[1].c_eze_0_ -> point(0,0) * Ez_->point(0,0) + pmlArr_[1].c_ezd1_0_ -> point(0,0) * pmlArr_[1].Dz_->point(0,0) - pmlArr_[1].c_ezd0_0_ -> point(0,0) * dzstore;
                                //Bot Right
                                dzstore = pmlArr_[1].Dz_end_->point(0,0);
                                pmlArr_[1].Dz_end_->point(0,0) = pmlArr_[1].c_dzd_n_ -> point(0,0) * pmlArr_[1].Dz_end_->point(0,0) + pmlArr_[1].c_dzh_n_ -> point(0,0) * (-1.0*Hy_->point(nx_-1-1,0)) - (Hx_->point(nx_-1,0));
                                Ez_->point(nx_-1,0) = pmlArr_[1].c_eze_n_ -> point(0,0) * Ez_->point(nx_-1,0) + pmlArr_[1].c_ezd1_n_ -> point(0,0) * pmlArr_[1].Dz_end_->point(0,0) - pmlArr_[1].c_ezd0_n_ -> point(0,0) * dzstore;
                                //Top Right
                                dzstore = pmlArr_[1].Dz_end_->point(0,ny_-1);
                                pmlArr_[1].Dz_end_->point(0,ny_-1) = pmlArr_[1].c_dzd_n_ -> point(0,ny_-1) * pmlArr_[1].Dz_end_->point(0,ny_-1) + pmlArr_[1].c_dzh_n_ -> point(0,ny_-1) * ((-1.0*Hy_->point(nx_-1-1,ny_-1)) - (-1.0*Hx_->point(nx_-1,ny_-1-1)));
                                Ez_->point(nx_-1,ny_-1) = pmlArr_[1].c_eze_n_ -> point(0,ny_-1) * Ez_->point(nx_-1,ny_-1) + pmlArr_[1].c_ezd1_n_ -> point(0,ny_-1) * pmlArr_[1].Dz_end_->point(0,ny_-1) - pmlArr_[1].c_ezd0_n_ -> point(0,ny_-1) * dzstore;
                                //Top Left
                                dzstore = pmlArr_[1].Dz_->point(0,ny_-1);
                                pmlArr_[1].Dz_->point(0,ny_-1) = pmlArr_[1].c_dzd_0_ -> point(0,ny_-1) * pmlArr_[1].Dz_->point(0,ny_-1) + pmlArr_[1].c_dzh_0_ -> point(0,ny_-1) * ((Hy_->point(0,ny_-1)) - (-1.0*Hx_->point(0,ny_-1-1)));
                                Ez_->point(0,ny_-1) = pmlArr_[1].c_eze_0_ -> point(0,ny_-1) * Ez_->point(0,ny_-1) + pmlArr_[1].c_ezd1_0_ -> point(0,ny_-1) * pmlArr_[1].Dz_->point(0,ny_-1) - pmlArr_[1].c_ezd0_0_ -> point(0,ny_-1) * dzstore;

                                for(int ii = 1; ii < pmlArr_[1].thickness(); ii++)
                                {
                                    //Bot Left
                                    dzstore = pmlArr_[1].Dz_->point(ii,0);
                                    pmlArr_[1].Dz_->point(ii,0) = pmlArr_[1].c_dzd_0_ -> point(ii,0) * pmlArr_[1].Dz_->point(ii,0) + pmlArr_[1].c_dzh_0_ -> point(ii,0) * ((Hy_->point(ii,0)-Hy_->point(ii-1,0)) - (Hx_->point(ii,0)));
                                    Ez_->point(ii,0) = pmlArr_[1].c_eze_0_ -> point(ii,0) * Ez_->point(ii,0) + pmlArr_[1].c_ezd1_0_ -> point(ii,0) * pmlArr_[1].Dz_->point(ii,0) - pmlArr_[1].c_ezd0_0_ -> point(ii,0) * dzstore;
                                    //Bot Right
                                    dzstore = pmlArr_[1].Dz_end_->point(ii,0);
                                    pmlArr_[1].Dz_end_->point(ii,0) = pmlArr_[1].c_dzd_n_ -> point(ii,0) * pmlArr_[1].Dz_end_->point(ii,0) + pmlArr_[1].c_dzh_n_ -> point(ii,0) * ((Hy_->point(nx_-1-ii,0)-Hy_->point(nx_-1-ii-1,0)) - (Hx_->point(nx_-1-ii,0)));
                                    Ez_->point(nx_-1-ii,0) = pmlArr_[1].c_eze_n_ -> point(ii,0) * Ez_->point(nx_-1-ii,0) + pmlArr_[1].c_ezd1_n_ -> point(ii,0) * pmlArr_[1].Dz_end_->point(ii,0) - pmlArr_[1].c_ezd0_n_ -> point(ii,0) * dzstore;
                                    //Top Right
                                    dzstore = pmlArr_[1].Dz_end_->point(ii,ny_-1);
                                    pmlArr_[1].Dz_end_->point(ii,ny_-1) = pmlArr_[1].c_dzd_n_ -> point(ii,ny_-1) * pmlArr_[1].Dz_end_->point(ii,ny_-1) + pmlArr_[1].c_dzh_n_ -> point(ii,ny_-1) * ((Hy_->point(nx_-1-ii,ny_-1)-Hy_->point(nx_-1-ii-1,ny_-1)) - (-1.0*Hx_->point(nx_-1-ii,ny_-1-1)));
                                    Ez_->point(nx_-1-ii,ny_-1) = pmlArr_[1].c_eze_n_ -> point(ii,ny_-1) * Ez_->point(nx_-1-ii,ny_-1) + pmlArr_[1].c_ezd1_n_ -> point(ii,ny_-1) * pmlArr_[1].Dz_end_->point(ii,ny_-1) - pmlArr_[1].c_ezd0_n_ -> point(ii,ny_-1) * dzstore;
                                    //Top Left
                                    dzstore = pmlArr_[1].Dz_->point(ii,ny_-1);
                                    pmlArr_[1].Dz_->point(ii,ny_-1) = pmlArr_[1].c_dzd_0_ -> point(ii,ny_-1) * pmlArr_[1].Dz_->point(ii,ny_-1) + pmlArr_[1].c_dzh_0_ -> point(ii,ny_-1) * ((Hy_->point(ii,ny_-1)-Hy_->point(ii-1,ny_-1)) - (-1.0*Hx_->point(ii,ny_-1-1)));
                                    Ez_->point(ii,ny_-1) = pmlArr_[1].c_eze_0_ -> point(ii,ny_-1) * Ez_->point(ii,ny_-1) + pmlArr_[1].c_ezd1_0_ -> point(ii,ny_-1) * pmlArr_[1].Dz_->point(ii,ny_-1) - pmlArr_[1].c_ezd0_0_ -> point(ii,ny_-1) * dzstore;
                                }
                                for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                                {
                                    //Bot Left
                                    dzstore = pmlArr_[1].Dz_->point(0,jj);
                                    pmlArr_[1].Dz_->point(0,jj) = pmlArr_[1].c_dzd_0_ -> point(0,jj) * pmlArr_[1].Dz_->point(0,jj) + pmlArr_[1].c_dzh_0_ -> point(0,jj) * ((Hy_->point(0,jj)) - (Hx_->point(0,jj)-Hx_->point(0,jj-1)));
                                    Ez_->point(0,jj) = pmlArr_[1].c_eze_0_ -> point(0,jj) * Ez_->point(0,jj) + pmlArr_[1].c_ezd1_0_ -> point(0,jj) * pmlArr_[1].Dz_->point(0,jj) - pmlArr_[1].c_ezd0_0_ -> point(0,jj) * dzstore;
                                    //Bot Right
                                    dzstore = pmlArr_[1].Dz_end_->point(0,jj);
                                    pmlArr_[1].Dz_end_->point(0,jj) = pmlArr_[1].c_dzd_n_ -> point(0,jj) * pmlArr_[1].Dz_end_->point(0,jj) + pmlArr_[1].c_dzh_n_ -> point(0,jj) * ((-1.0*Hy_->point(nx_-1-1,jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1)));
                                    Ez_->point(nx_-1,jj) = pmlArr_[1].c_eze_n_ -> point(0,jj) * Ez_->point(nx_-1,jj) + pmlArr_[1].c_ezd1_n_ -> point(0,jj) * pmlArr_[1].Dz_end_->point(0,jj) - pmlArr_[1].c_ezd0_n_ -> point(0,jj) * dzstore;
                                    //Top Right
                                    dzstore = pmlArr_[1].Dz_end_->point(0,ny_-1-jj);
                                    pmlArr_[1].Dz_end_->point(0,ny_-1-jj) = pmlArr_[1].c_dzd_n_ -> point(0,ny_-1-jj) * pmlArr_[1].Dz_end_->point(0,ny_-1-jj) + pmlArr_[1].c_dzh_n_ -> point(0,ny_-1-jj) * ((-1.0*Hy_->point(nx_-1-1,ny_-1-jj)) - (Hx_->point(nx_-1,ny_-1-jj)-Hx_->point(nx_-1,ny_-1-jj-1)));
                                    Ez_->point(nx_-1,ny_-1-jj) = pmlArr_[1].c_eze_n_ -> point(0,ny_-1-jj) * Ez_->point(nx_-1,ny_-1-jj) + pmlArr_[1].c_ezd1_n_ -> point(0,ny_-1-jj) * pmlArr_[1].Dz_end_->point(0,ny_-1-jj) - pmlArr_[1].c_ezd0_n_ -> point(0,ny_-1-jj) * dzstore;
                                    //Top Left
                                    dzstore = pmlArr_[1].Dz_->point(0,ny_-1-jj);
                                    pmlArr_[1].Dz_->point(0,ny_-1-jj) = pmlArr_[1].c_dzd_0_ -> point(0,ny_-1-jj) * pmlArr_[1].Dz_->point(0,ny_-1-jj) + pmlArr_[1].c_dzh_0_ -> point(0,ny_-1-jj) * ((Hy_->point(0,ny_-1-jj)) - (Hx_->point(0,ny_-1-jj)-Hx_->point(0,ny_-1-jj-1)));
                                    Ez_->point(0,ny_-1-jj) = pmlArr_[1].c_eze_0_ -> point(0,ny_-1-jj) * Ez_->point(0,ny_-1-jj) + pmlArr_[1].c_ezd1_0_ -> point(0,ny_-1-jj) * pmlArr_[1].Dz_->point(0,ny_-1-jj) - pmlArr_[1].c_ezd0_0_ -> point(0,ny_-1-jj) * dzstore;
                                }
                                for(int ii = 1; ii < pmlArr_[1].thickness(); ii++)
                                {
                                    for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                                    {
                                        //Bot Left
                                        dzstore = pmlArr_[1].Dz_->point(ii,jj);
                                        pmlArr_[1].Dz_->point(ii,jj) = pmlArr_[1].c_dzd_0_ -> point(ii,jj) * pmlArr_[1].Dz_->point(ii,jj) + pmlArr_[1].c_dzh_0_ -> point(ii,jj) * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                        Ez_->point(ii,jj) = pmlArr_[1].c_eze_0_ -> point(ii,jj) * Ez_->point(ii,jj) + pmlArr_[1].c_ezd1_0_ -> point(ii,jj) * pmlArr_[1].Dz_->point(ii,jj) - pmlArr_[1].c_ezd0_0_ -> point(ii,jj) * dzstore;
                                        //Bot Right
                                        dzstore = pmlArr_[1].Dz_end_->point(ii,jj);
                                        pmlArr_[1].Dz_end_->point(ii,jj) = pmlArr_[1].c_dzd_n_ -> point(ii,jj) * pmlArr_[1].Dz_end_->point(ii,jj) + pmlArr_[1].c_dzh_n_ -> point(ii,jj) * ((Hy_->point(nx_-1-ii,jj)-Hy_->point(nx_-1-ii-1,jj)) - (Hx_->point(nx_-1-ii,jj)-Hx_->point(nx_-1-ii,jj-1)));
                                        Ez_->point(nx_-1-ii,jj) = pmlArr_[1].c_eze_n_ -> point(ii,jj) * Ez_->point(nx_-1-ii,jj) + pmlArr_[1].c_ezd1_n_ -> point(ii,jj) * pmlArr_[1].Dz_end_->point(ii,jj) - pmlArr_[1].c_ezd0_n_ -> point(ii,jj) * dzstore;
                                        //Top Right
                                        dzstore = pmlArr_[1].Dz_end_->point(ii,ny_-1-jj);
                                        pmlArr_[1].Dz_end_->point(ii,ny_-1-jj) = pmlArr_[1].c_dzd_n_ -> point(ii,ny_-1-jj) * pmlArr_[1].Dz_end_->point(ii,ny_-1-jj) + pmlArr_[1].c_dzh_n_ -> point(ii,ny_-1-jj) * ((Hy_->point(nx_-1-ii,ny_-1-jj)-Hy_->point(nx_-1-ii-1,ny_-1-jj)) - (Hx_->point(nx_-1-ii,ny_-1-jj)-Hx_->point(nx_-1-ii,ny_-1-jj-1)));
                                        Ez_->point(nx_-1-ii,ny_-1-jj) = pmlArr_[1].c_eze_n_ -> point(ii,ny_-1-jj) * Ez_->point(nx_-1-ii,ny_-1-jj) + pmlArr_[1].c_ezd1_n_ -> point(ii,ny_-1-jj) * pmlArr_[1].Dz_end_->point(ii,ny_-1-jj) - pmlArr_[1].c_ezd0_n_ -> point(ii,ny_-1-jj) * dzstore;
                                        //Top Left
                                        dzstore = pmlArr_[1].Dz_->point(ii,ny_-1-jj);
                                        pmlArr_[1].Dz_->point(ii,ny_-1-jj) = pmlArr_[1].c_dzd_0_ -> point(ii,ny_-1-jj) * pmlArr_[1].Dz_->point(ii,ny_-1-jj) + pmlArr_[1].c_dzh_0_ -> point(ii,ny_-1-jj) * ((Hy_->point(ii,ny_-1-jj)-Hy_->point(ii-1,ny_-1-jj)) - (Hx_->point(ii,ny_-1-jj)-Hx_->point(ii,ny_-1-jj-1)));
                                        Ez_->point(ii,ny_-1-jj) = pmlArr_[1].c_eze_0_ -> point(ii,ny_-1-jj) * Ez_->point(ii,ny_-1-jj) + pmlArr_[1].c_ezd1_0_ -> point(ii,ny_-1-jj) * pmlArr_[1].Dz_->point(ii,ny_-1-jj) - pmlArr_[1].c_ezd0_0_ -> point(ii,ny_-1-jj) * dzstore;
                                    }
                                }
                            }
                        }
                        break;
                    }
                    case Z:
                        throw logic_error("z not implimented");
                        break;
                    default:
                        throw logic_error("hit default");
                        break;
                }
            }
        }
    }*/
}

/**
 * @brief Updates all fields to the next time step
 * @details updates the FDTD cell to the next time step
 *
 */
void FDTDField::step()
{
    // Source
    for(int kk = 0; kk < srcArr_.size(); kk ++)
    {
        if(abs(real(srcArr_[kk].prof().pulse(static_cast<double>(t_step_)))) > 1.0e-70)
        {
            int ii = srcArr_[kk].loc()[0];
            int jj = srcArr_[kk].loc()[1];
            double eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0);
            double c_ezj = dt_/(eps);
            Ez_ -> point(ii,jj) += c_ezj * srcArr_[kk].prof().pulse(static_cast<double>(t_step_));
        }
        /*switch ( srcArr_[kk].pol() )
        {
            case EZ: //if(srcArr[kk].pol() == EZ)
                if(abs(srcArr_[kk].prof().pulse(static_cast<double>(t_step_)).real()) > 1.0e-70)
                    Ez_ -> point(ii,jj) = srcArr_[kk].prof().pulse(static_cast<double>(t_step_));
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
        }*/
    }
    updateH();

    updateE();
    for(int ii = 0; ii < dtcArr_.size(); ii ++)
        ouputField(dtcArr_[ii]);
    //if(abs(tcur_-floor(tcur_+0.5)) < 1e-7)
    if(false)
    {
        string fname("fout/Hx/HxField_t" + to_string(static_cast<int>(t_step_))+".dat");
        Hx_->gridOut(fname);
        fname = "fout/Hy/HyField_t" + to_string(static_cast<int>(t_step_))+".dat";
        Hy_->gridOut(fname);
        fname = "fout/Ez/EzField_t" + to_string(static_cast<int>(t_step_))+".dat";
        Ez_->gridOut(fname);
    }

    tcur_ += dt_;
    t_step_ ++;
}
