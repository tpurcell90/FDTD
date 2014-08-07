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
/**
 * @brief Constructs a FDTD field from a programInputs object
 * @details Constructs the FDTD field manger using the information from the inputs parameter. Sets dx, and dy to 1/res, dt_ is set to S/res, The unused fields are set to null, nx and ny are calculated by rouding the product of the physical cell size with the
 *
 * @param IP Input object created from an input file
 */
FDTDField::FDTDField(programInputs &IP)
{
    //Define cell parameters
    tcur_     = 0;
    res_      = IP.res_;
    dx_       = 1.0/res_;
    dy_       = 1.0/res_;
    dt_       = IP.courant_ * dx_;
    nx_       = floor(static_cast<double>(res_) * static_cast<double>(IP.x_size_)+ 0.5) + 1; //Better way here; + 1 to include the 0 point
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

/**
 * @brief Initializes the physical grid for materials look up
 * @details Initializes the physical grid for materials look up, sets the object to the number in the object array will overwrite to the last one if multiple objects exist at the same point
 *
 */
void FDTDField::initializeGrid()
{
    for(int kk = 0; kk < objArr_.size(); kk ++)
    {
        if(Hz_)
        {
            if(objArr_[kk].s() == sphere)
            {
                for(int ii = 0; ii < nx_-1;ii ++)
                {
                    for(int jj = 0; jj < ny_-1; jj ++)
                    {
                        vector<double> pt({(ii+0.5)*dx_,jj*dy_});
                        if(objArr_[kk].isObj(pt))
                            phys_Ex_->point(ii,jj) = kk;
                        pt[1] += 0.5*dy_;
                        if(objArr_[kk].isObj(pt))
                            phys_Hz_->point(ii,jj) = kk;
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
                        vector<double> pt({(ii+0.5)*dx_,jj*dy_});
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Ex_->point(ii,jj) = kk;
                        pt[1] += 0.5*dy_;
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Hz_->point(ii,jj) = kk;
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
                        vector<double> pt({(ii+0.5-(nx_-1)/2.0)*dx_,(jj-static_cast<double>(ny_-1)/2.0)*dy_});
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Hy_->point(ii,jj) = kk;
                        pt[0] -= 0.5*dx_;
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Ez_->point(ii,jj) = kk;
                        pt[1] += 0.5*dy_;
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Hx_->point(ii,jj) = kk;
                    }
                }
                for(int ii = 0; ii < nx_-1;ii ++)
                {
                    vector<double> pt({(ii+0.5-(nx_-1)/2.0)*dx_,(ny_-1-(ny_-1)/2.0)*dy_});
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Hy_->point(ii,ny_-1) = kk;
                    pt[0] -= 0.5*dx_;
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ez_->point(ii,ny_-1) = kk;
                }
                for(int jj = 0; jj < ny_-1; jj ++)
                {
                    vector<double> pt({(nx_-1-(nx_-1)/2.0)*dx_,(jj-(ny_-1)/2.0)*dy_});
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ez_->point(nx_-1,jj) = kk;
                    pt[1] += 0.5*dy_;
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Hx_->point(nx_-1,jj) = kk;
                }
                vector<double> pt({(nx_-1-(nx_-1)/2.0)*dx_,(ny_-1-(ny_-1)/2.0)*dy_});
                if(objArr_[kk].isObj(pt)==true)
                    phys_Ez_->point(nx_-1,ny_-1) = kk;
            }
            else if(objArr_[kk].s() == block)
            {
                for(int ii = 0; ii < nx_-1;ii ++)
                {
                    for(int jj = 0; jj < ny_-1; jj ++)
                    {
                        vector<double> pt({(ii+0.5-(nx_-1)/2.0)*dx_,(jj-static_cast<double>(ny_-1)/2.0)*dy_});
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Hy_->point(ii,jj) = kk;
                        pt[0] -= 0.5*dx_;
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Ez_->point(ii,jj) = kk;
                        pt[1] += 0.5*dy_;
                        if(objArr_[kk].isObj(pt)==true)
                            phys_Hx_->point(ii,jj) = kk;
                    }
                }
                for(int ii = 0; ii < nx_-1;ii ++)
                {
                    vector<double> pt({(ii+0.5-(nx_-1)/2.0)*dx_,(ny_-1-(ny_-1)/2.0)*dy_});
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Hy_->point(ii,ny_-1) = kk;
                    pt[0] -= 0.5*dx_;
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ez_->point(ii,ny_-1) = kk;
                }
                for(int jj = 0; jj < ny_-1; jj ++)
                {
                    vector<double> pt({(nx_-1-(nx_-1)/2.0)*dx_,(jj-(ny_-1)/2.0)*dy_});
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ez_->point(nx_-1,jj) = kk;
                    pt[1] += 0.5*dy_;
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Hx_->point(nx_-1,jj) = kk;
                }
                vector<double> pt({(nx_-1-(nx_-1)/2.0)*dx_,(ny_-1-(ny_-1)/2.0)*dy_});
                if(objArr_[kk].isObj(pt)==true)
                    phys_Ez_->point(nx_-1,ny_-1) = kk;
            }
        }
    }
    /*string fname("fout/Hx/HxField_t" + to_string(static_cast<int>(tcur_/dt_))+".dat");
    phys_Hx_->gridOut(fname);
    fname = "fout/Hy/HyField_t" + to_string(static_cast<int>(tcur_/dt_))+".dat";
    phys_Hy_->gridOut(fname);
    fname = "fout/Ez/EzField_t" + to_string(static_cast<int>(tcur_/dt_))+".dat";
    phys_Ez_->gridOut(fname);*/
}
/**
 * @brief Outputs the relevant field information to an output file specified by the detector
 * @details Outputs the relevant field information to an output file specified by the detector
 *
 * @param d The Detector used for the output
 */
void FDTDField::ouputField(Detector<double> d) //iostream as input parameter?
{
    ofstream outFile;
    outFile.open(d.outfile(), ios_base::app);
    double eps = 1.00;
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
/**
 * @brief This performs one time update using the FDTD update equations and a UPML
 * @details The function updates each point according to what it is. Because The units are scalable eps_0, mu_0 and c are all set to 1 leading to the following update equations
 * For Normal Space TM modes
 * \(E_z^{q+1}\left[i,j\right] = \frac{1 - \frac{\sigma\left[i,j\right]*d_t}{2*\epsilon\left[i,j\right]}}{1+\frac{\sigma\left[i,j\right]*d_t}{2\epsilon\left[i,j\right]}} E_z^q\left[i,j]\right] +\frac{\frac{dt}{\epsilon\left[i,j\right]}}{1+\frac{\sigma\left[i,j\right]*d_t}{2*\epsilon\left[i,j\right]}}\left\{\frac{1}{dx} \left(H_y^{q+\frac{1}{2}}\left[i+\frac{1}{2},j+\frac{1}{2}\right] - H_y^{q+\frac{1}{2}}\left[i-\frac{1}{2},j+\frac{1}{2}\right]\right) - \frac{1}{dy}\left(H_x^{q+\frac{1}{2}}\left[i+\frac{1}{2},j+\frac{1}{2}\right] - H_x^{q+\frac{1}{2}}\left[i+\frac{1}{2},j-\frac{1}{2}\right]\right) \right\}  \)
 * \(H_x^{q+\frac{3}{2}}\left[i,j+\frac{1}{2}\right] =  \frac{1 - \frac{\sigma^*\left[i,j+\frac{1}{2}\right] d_t}{2*mu\left[i,j+\frac{1}{2}\right]d_t}{2\mu\left[i,j+\frac{1}{2}\right]}} H_x^{q+\frac{3}{2}}\left[i,j+\frac{1}{2}\right] + \frac{\frac{dt}{mu\left[i,j+\frac{1}{2}\right]}}{\sigma^*\left[i,j+\frac{1}{2}\right] d_t}{2*mu\left[i,j+\frac{1}{2}\right]d_t}{2\mu\left[i,j+\frac{1}{2}\right]}} \left\{ \frac{1}{dy} \left( E_z^{q+1}\left[i,j+1\right] - E_z^{q+1}\left[i,j\right]\right)  \right\}\)
 * \(H_y^{q+\frac{3}{2}}\left[i+\frac{1}{2},j\right] =  \frac{1 - \frac{\sigma^*\left[i+\frac{1}{2},j\right] d_t}{2*mu\left[i+\frac{1}{2},j\right]d_t}{2\mu\left[i+\frac{1}{2},j\right]}} H_y^{q+\frac{3}{2}}\left[i+\frac{1}{2},j\right] + \frac{\frac{dt}{mu\left[i+\frac{1}{2},j\right]}}{\sigma^*\left[i+\frac{1}{2},j\right] d_t}{2*mu\left[i+\frac{1}{2},j\right]d_t}{2\mu\left[i+\frac{1}{2},j\right]}} \left\{ \frac{1}{dy} \left( E_z^{q+1}\left[i+1,j\right] - E_z^{q+1}\left[i,j\right]\right)  \right\}\)
 *
 * For UPML space TM mode
 * \(D_z^{q+1}\left[i,j\right] = \frac{2 \epsilon \kappa_x - \sigma_x dt}{2 \epsilon  \kappa_x + \sigma_x dt} D_z^{q}\left[i,j\right] +  \frac{2  \epsilon  dt}{2  \epsilon  \kappa_x + \sigma_x  dt} \left\{\frac{1}{dx} \left(H_y^{q+\frac{1}{2}}\left[i+\frac{1}{2},j+\frac{1}{2}\right] - H_y^{q+\frac{1}{2}}\left[i-\frac{1}{2},j+\frac{1}{2}\right]\right) - \frac{1}{dy}\left(H_x^{q+\frac{1}{2}}\left[i+\frac{1}{2},j+\frac{1}{2}\right] - H_x^{q+\frac{1}{2}}\left[i+\frac{1}{2},j-\frac{1}{2}\right]\right) \right\}\)
 *
 *\(E_z^{q+1}\left[i,j\right] =  \frac{2 \epsilon \kappa_y - \sigma_y dt}{2 \epsilon  \kappa_y + \sigma_y dt} E_z^{q}\left[i,j\right] + \frac{1}{\left(2\epsilon \kappa_y +\sigma_y dt\right)\epsilon} \left\{ \left(2\epsilon\kappa_z + \sigma_z dt\right) D_z^{q+1}\left[i,j\right]  -\left(2\epsilon\kappa_z - \sigma_z dt\right) D_z^{q}\left[i,j\right] \right\} \)
 *
 *\(B_x^{q+\frac{3}{2}}\left[i,j+\frac{1}{2}\right] = \frac{2 \epsilon \kappa_y - \sigma_y dt}{2 \epsilon  \kappa_y + \sigma_y dt} B_x^{q+\frac{3}{2}}\left[i,j+\frac{1}{2}\right] - \frac{2  \epsilon  dt}{2  \epsilon  \kappa_y + \sigma_y  dt} \left\{ \frac{1}{dy} \left( E_z^{q+1}\left[i,j+1\right] - E_z^{q+1}\left[i,j\right]\right)  \right\} \)
 *
 *\(H_x^{q+\frac{3}{2}}\left[i,j+\frac{1}{2}\right] =  \frac{2 \epsilon \kappa_z - \sigma_z dt}{2 \epsilon  \kappa_z + \sigma_z dt} H_x^{q+\frac{1}{2}}\left[i,j+\frac{1}{2}\right] + \frac{1}{\left(2\epsilon \kappa_z +\sigma_z dt\right)\epsilon} \left\{ \left(2\epsilon\kappa_x + \sigma_x dt\right) B_x^{q+\frac{3}{2}}\left[i,j+\frac{1}{2}\right] -\left(2\epsilon\kappa_x - \sigma_x dt\right) B_x^{q+\frac{1}{2}}\left[i,j+\frac{1}{2}\right]\right\} \)
 *
 *\(B_y^{q+\frac{3}{2}}\left[i+\frac{1}{2},j\right] = \frac{2 \epsilon \kappa_z - \sigma_z dt}{2 \epsilon  \kappa_z + \sigma_z dt} B_y^{q+\frac{3}{2}}\left[i+\frac{1}{2},j\right] +  \frac{2  \epsilon  dt}{2  \epsilon  \kappa_z + \sigma_z  dt} \left\{ \frac{1}{dz} \frac{1}{dx} \left(E_z^{q+1}\left[i+1,j\right] - E_z^{q+1}\left[i,j\right)\right]     \right\} \)
 *
 *\(H_y^{q+\frac{3}{2}}\left[i+\frac{1}{2},j,\right]  = \frac{2 \epsilon \kappa_x - \sigma_x dt}{2 \epsilon  \kappa_x + \sigma_x dt} H_y^{q+\frac{1}{2}}\left[i+\frac{1}{2},j\right] + \frac{1}{\left(2\epsilon \kappa_x +\sigma_x dt\right)\epsilon} \left\{ \left(2\epsilon\kappa_y + \sigma_y dt\right) B_y^{q+\frac{3}{2}}\left[i+\frac{1}{2},j\right]  -\left(2\epsilon\kappa_y - \sigma_y dt\right) B_y^{q+\frac{1}{2}}\left[i+\frac{1}{2},j\right] \right\}  \)
 *
 */
void FDTDField::step()
{
    // Source
    for(int kk = 0; kk < srcArr_.size(); kk ++)
    {
        int ii = srcArr_[kk].loc()[0];
        int jj = srcArr_[kk].loc()[1];
        switch ( srcArr_[kk].pol() )
        {
            case EZ: //if(srcArr[kk].pol() == EZ)
                if(srcArr_[kk].prof().pulse(tcur_) != 0.0)
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
    if(Ez_)
    {
        if(xPML_ != 0 && yPML_ !=0)
        {
            for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
                {
                    //Always true in a magnetic lossless enviroment, with mu0 = 1, sigma_m = 0;
                    double c_hxh = 1.0;
                    double c_hxe = 1.0 * dt_/dx_;
                    double c_hyh = 1.0;
                    double c_hye = 1.0 * dt_/dy_;
                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) - c_hxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                }
            }
        }
        else if(xPML_ != 0)
        {
            for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                for(int jj = 0; jj < ny_ - 1; jj ++)
                {
                    double c_hxh = 1.0;
                    double c_hxe = 1.0 * dt_/dx_;
                    double c_hyh = 1.0;
                    double c_hye = 1.0 * dt_/dy_;
                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) - c_hxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                }
            }
            for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                double c_hyh = 1.0;
                double c_hye = 1.0 * dt_/dy_;
                Hy_->point(ii,ny_-1) = c_hyh * Hy_->point(ii,ny_-1) + c_hye * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
            }
        }
        else if(yPML_ != 0)
        {
            for(int ii = 0; ii < nx_ - 1; ii ++)
            {
                for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
                {
                    double c_hxh = 1.0;
                    double c_hxe = 1.0 * dt_/dx_;
                    double c_hyh = 1.0;
                    double c_hye = 1.0 * dt_/dy_;
                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) - c_hxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                }
            }
            for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
            {
                double c_hxh = 1.0;
                double c_hxe = 1.0 * dt_/dx_;
                Hx_->point(nx_-1,jj) = c_hxh * Hx_->point(nx_-1,jj) - c_hxe * (Ez_->point(nx_-1,jj+1)-Ez_->point(nx_-1,jj));
            }
        }
        else
        {
            for(int ii = 0; ii < nx_ - 1; ii ++)
            {
                for(int jj = 0; jj < ny_ - 1; jj ++)
                {
                    double c_hxh = 1.0;
                    double c_hxe = 1.0 * dt_/dx_;
                    double c_hyh = 1.0;
                    double c_hye = 1.0 * dt_/dy_;
                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) - c_hxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                }
            }
            for(int jj = 0; jj < nx_ - 1; jj ++)
            {
                double c_hxh = 1.0;
                double c_hxe = 1.0 * dt_/dx_;
                Hx_->point(nx_-1,jj) = c_hxh * Hx_->point(nx_-1,jj) - c_hxe * (Ez_->point(nx_-1,jj+1)-Ez_->point(nx_-1,jj));
            }
            for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                double c_hyh = 1.0;
                double c_hye = 1.0 * dt_/dy_;
                Hy_->point(ii,ny_-1) = c_hyh * Hy_->point(ii,ny_-1) + c_hye * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
            }
        }
    }
    else
    {
        for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
        {
            for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
            {
                //Always true in a magnetic lossless enviroment, with mu0 = 1, sigma_m = 0;
                double c_hzh = 1.0;
                double c_hze = 1.0 * dt_/dx_;
                Hz_->point(ii,jj) = c_hzh * Hz_->point(ii,jj) + c_hze * ((Ex_->point(ii,jj+1) - Ex_->point(ii,jj)) - (Ey_->point(ii+1,jj)-Ey_->point(ii,jj)));
            }
        }
    }
    if(Ez_)
    {
        for(int kk = 0; kk < pmlArr_.size(); kk++)
        {
            double eps = 1.0;
            double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
            double sigz = 0.0;
            switch(pmlArr_[kk].d())
            {
                case X:
                {
                    double sigxx = pmlArr_[kk].sigma(0.0);
                    double sigxy = pmlArr_[kk].sigma(0.5);
                    double sigyx = 0.0;
                    double sigyy = 0.0;

                    double c_bxb; double c_bxe; double c_hxh; double c_hxb0; double c_hxb1;
                    double c_byb; double c_bye; double c_hyh;double c_hyb0; double c_hyb1;

                    if (yPML_== 0)
                    {
                        for (int jj = 0; jj < ny_-1; jj++)
                        {
                            eps    = objArr_[phys_Hx_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            eps    = objArr_[phys_Hy_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            //Update both on the left side and Hx on the right
                            double bxstore = pmlArr_[kk].Bx_->point(0,jj);
                            double bystore = pmlArr_[kk].By_->point(0,jj);

                            pmlArr_[kk].Bx_->point(0,jj) = c_bxb * pmlArr_[kk].Bx_->point(0,jj) - c_bxe * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                            pmlArr_[kk].By_->point(0,jj) = c_byb * pmlArr_[kk].By_->point(0,jj) + c_bye * (Ez_->point(0+1,jj)-Ez_->point(0,jj));

                            Hx_->point(0,jj) = c_hxh * Hx_->point(0,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(0,jj) - c_hxb0 * bxstore;
                            Hy_->point(0,jj) = c_hyh * Hy_->point(0,jj) + c_hyb1 * pmlArr_[kk].By_->point(0,jj) - c_hyb0 * bystore;

                            eps    = objArr_[phys_Hx_->point(nx_-1,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            bxstore = pmlArr_[kk].Bx_end_->point(0,jj);
                            pmlArr_[kk].Bx_end_->point(0,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(0,jj) - c_bxe * (Ez_->point(nx_-1, jj+1)-Ez_->point(nx_-1,jj));
                            Hx_->point(nx_-1,jj) = c_hxh * Hx_->point(nx_-1,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(0,jj) - c_hxb0 * bxstore;
                        }
                        //Top left corner
                        eps    = objArr_[phys_Hy_->point(0,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                        c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                        c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                        c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                        c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                        c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                        double bystore = pmlArr_[kk].By_->point(0,ny_-1);
                        pmlArr_[kk].By_->point(0,ny_-1) = c_byb * pmlArr_[kk].By_->point(0,ny_-1) + c_bye * (Ez_->point(0+1,ny_-1)-Ez_->point(0,ny_-1));
                        Hy_->point(0,ny_-1) = c_hyh * Hy_->point(0,ny_-1) + c_hyb1 * pmlArr_[kk].By_->point(0,ny_-1) - c_hyb0 * bystore;
                        for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                        {
                            sigxx = pmlArr_[kk].sigma(static_cast<double>(ii));
                            sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5);

                            sigyx = 0.0;
                            sigyy = 0.0;

                            eps    = objArr_[phys_Hy_->point(ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                            //top row
                            bystore = pmlArr_[kk].By_->point(ii,ny_-1);
                            pmlArr_[kk].By_->point(ii,ny_-1) = c_byb * pmlArr_[kk].By_->point(ii,ny_-1) + c_bye * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
                            Hy_->point(ii,ny_-1) = c_hyh * Hy_->point(ii,ny_-1) + c_hyb1 * pmlArr_[kk].By_->point(ii,ny_-1) - c_hyb0 * bystore;

                            for(int jj = 0; jj < ny_-1; jj ++)
                            {
                                sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5);

                                eps    = objArr_[phys_Hx_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                //Update everything

                                double bxstore = pmlArr_[kk].Bx_->point(ii,jj);
                                bystore = pmlArr_[kk].By_->point(ii,jj);

                                pmlArr_[kk].Bx_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                pmlArr_[kk].By_->point(ii,jj) = c_byb * pmlArr_[kk].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));

                                Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                                Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[kk].By_->point(ii,jj) - c_hyb0 * bystore;

                                sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) - 0.5);
                                eps    = objArr_[phys_Hx_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                bxstore = pmlArr_[kk].Bx_end_->point(ii,jj);
                                bystore = pmlArr_[kk].By_end_->point(ii,jj);

                                pmlArr_[kk].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point(nx_-1-ii, jj+1)-Ez_->point(nx_-1-ii,jj));
                                pmlArr_[kk].By_end_->point(ii,jj) = c_byb * pmlArr_[kk].By_end_->point(ii,jj) + c_bye * (Ez_->point(nx_-1-ii+1, jj)-Ez_->point(nx_-1-ii,jj));

                                Hx_->point(nx_-1-ii,jj) = c_hxh * Hx_->point((nx_-1) - ii,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                                Hy_->point(nx_-1-ii,jj) = c_hyh * Hy_->point((nx_-1) - ii,jj) + c_hyb1 * pmlArr_[kk].By_end_->point(ii,jj) - c_hyb0 * bystore;
                            }
                            // for the top row only upadate Hy

                            sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) - 0.5);
                            eps    = objArr_[phys_Hy_->point(nx_-1-ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                            bystore = pmlArr_[kk].By_end_->point(ii,ny_-1);
                            pmlArr_[kk].By_end_->point(ii,ny_-1) = c_byb * pmlArr_[kk].By_end_->point(ii,ny_-1) + c_bye * (Ez_->point(nx_-1-ii+1, ny_-1)-Ez_->point(nx_-1-ii,ny_-1));
                            Hy_->point(nx_-1-ii,ny_-1) = c_hyh * Hy_->point(nx_-1-ii,ny_-1) + c_hyb1 * pmlArr_[kk].By_end_->point(ii,ny_-1) - c_hyb0 * bystore;
                        }
                    }
                    else
                    {
                        for (int jj = yPML_; jj < ny_- yPML_; jj++)
                        {
                            eps    = objArr_[phys_Hx_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            eps    = objArr_[phys_Hy_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            //Update both on the left side and Hx on the right
                            double bxstore = pmlArr_[kk].Bx_->point(0,jj);
                            double bystore = pmlArr_[kk].By_->point(0,jj);

                            pmlArr_[kk].Bx_->point(0,jj) = c_bxb * pmlArr_[kk].Bx_->point(0,jj) - c_bxe * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                            pmlArr_[kk].By_->point(0,jj) = c_byb * pmlArr_[kk].By_->point(0,jj) + c_bye * (Ez_->point(0+1,jj)-Ez_->point(0,jj));

                            Hx_->point(0,jj) = c_hxh * Hx_->point(0,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(0,jj) - c_hxb0 * bxstore;
                            Hy_->point(0,jj) = c_hyh * Hy_->point(0,jj) + c_hyb1 * pmlArr_[kk].By_->point(0,jj) - c_hyb0 * bystore;

                            eps    = objArr_[phys_Hx_->point(nx_-1,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            bxstore = pmlArr_[kk].Bx_end_->point(0,jj);
                            pmlArr_[kk].Bx_end_->point(0,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(0,jj) - c_bxe * (Ez_->point(nx_-1, jj+1)-Ez_->point(nx_-1,jj));
                            Hx_->point(nx_-1,jj) = c_hxh * Hx_->point(nx_-1,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(0,jj) - c_hxb0 * bxstore;
                        }
                        for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                        {
                            for(int jj = yPML_; jj < ny_- yPML_; jj ++)
                            {
                                //Update everything
                                sigxx = pmlArr_[kk].sigma(static_cast<double>(ii));
                                sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5);
                                sigyx = 0.0;
                                sigyy = 0.0;

                                eps    = objArr_[phys_Hx_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                double bxstore = pmlArr_[kk].Bx_->point(ii,jj);
                                double bystore = pmlArr_[kk].By_->point(ii,jj);

                                pmlArr_[kk].Bx_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                pmlArr_[kk].By_->point(ii,jj) = c_byb * pmlArr_[kk].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));

                                Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                                Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[kk].By_->point(ii,jj) - c_hyb0 * bystore;

                                eps    = objArr_[phys_Hx_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

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
                            sigxx = pmlArr_[0].sigma(0.0);
                            sigxy = pmlArr_[0].sigma(0.5);
                            sigyx = pmlArr_[1].sigma(0.5);
                            sigyy = pmlArr_[1].sigma(0.0);

                            eps    = objArr_[phys_Hx_->point(0,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            eps    = objArr_[phys_Hy_->point(0,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                            // Bot Left
                            double bxstore = pmlArr_[0].Bx_->point(0,0);
                            double bystore = pmlArr_[0].By_->point(0,0);
                            pmlArr_[0].Bx_->point(0,0) = c_bxb * pmlArr_[0].Bx_->point(0,0) - c_bxe * (Ez_->point(0,0+1)-Ez_->point(0,0));
                            pmlArr_[0].By_->point(0,0) = c_byb * pmlArr_[0].By_->point(0,0) + c_bye * (Ez_->point(0+1,0)-Ez_->point(0,0));
                            Hx_->point(0,0) = c_hxh * Hx_->point(0,0) + c_hxb1 * pmlArr_[0].Bx_->point(0,0) - c_hxb0 * bxstore;
                            Hy_->point(0,0) = c_hyh * Hy_->point(0,0) + c_hyb1 * pmlArr_[0].By_->point(0,0) - c_hyb0 * bystore;
                            //Bot Right
                            eps    = objArr_[phys_Hx_->point(nx_-1,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            bxstore = pmlArr_[0].Bx_end_->point(0,0);
                            pmlArr_[0].Bx_end_->point(0,0) = c_bxb * pmlArr_[0].Bx_end_->point(0,0) - c_bxe * (Ez_->point(nx_-1,0+1)-Ez_->point(nx_-1,0));
                            Hx_->point(nx_-1,0) = c_hxh * Hx_->point(nx_-1,0) + c_hxb1 * pmlArr_[0].Bx_end_->point(0,0) - c_hxb0 * bxstore;
                            //Top Left

                            eps    = objArr_[phys_Hy_->point(0,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                            bystore = pmlArr_[0].By_->point(0,ny_-1);
                            pmlArr_[0].By_->point(0,ny_-1) = c_byb * pmlArr_[0].By_->point(0,ny_-1) + c_bye * (Ez_->point(0+1,ny_-1)-Ez_->point(0,ny_-1));
                            Hy_->point(0,ny_-1) = c_hyh * Hy_->point(0,ny_-1) + c_hyb1 * pmlArr_[0].By_->point(0,ny_-1) - c_hyb0 * bystore;
                            for(int ii = 1; ii < pmlArr_[0].thickness(); ii++)
                            {
                                kapz = 1.0; sigz = 0.0; eps = 1.0;
                                kapx = 1.0; kapy = 1.0;
                                sigxx = pmlArr_[0].sigma(static_cast<double>(ii));
                                sigxy = pmlArr_[0].sigma(static_cast<double>(ii) + 0.5);
                                sigyx = pmlArr_[1].sigma(0.5);
                                sigyy = pmlArr_[1].sigma(0.0);

                                eps    = objArr_[phys_Hx_->point(ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                // Bot Left
                                double bxstore = pmlArr_[0].Bx_->point(ii,0);
                                double bystore = pmlArr_[0].By_->point(ii,0);
                                pmlArr_[0].Bx_->point(ii,0) = c_bxb * pmlArr_[0].Bx_->point(ii,0) - c_bxe * (Ez_->point(ii,0+1)-Ez_->point(ii,0));
                                pmlArr_[0].By_->point(ii,0) = c_byb * pmlArr_[0].By_->point(ii,0) + c_bye * (Ez_->point(ii+1,0)-Ez_->point(ii,0));
                                Hx_->point(ii,0) = c_hxh * Hx_->point(ii,0) + c_hxb1 * pmlArr_[0].Bx_->point(ii,0) - c_hxb0 * bxstore;
                                Hy_->point(ii,0) = c_hyh * Hy_->point(ii,0) + c_hyb1 * pmlArr_[0].By_->point(ii,0) - c_hyb0 * bystore;
                                //Bot Right
                                sigxy = pmlArr_[0].sigma(static_cast<double>(ii) - 0.5);
                                eps    = objArr_[phys_Hx_->point(nx_-1-ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(nx_-1-ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                bxstore = pmlArr_[0].Bx_end_->point(ii,0);
                                bystore = pmlArr_[0].By_end_->point(ii,0);
                                pmlArr_[0].Bx_end_->point(ii,0) = c_bxb * pmlArr_[0].Bx_end_->point(ii,0) - c_bxe * (Ez_->point(nx_-1-ii,0+1)-Ez_->point(nx_-1-ii,0));
                                pmlArr_[0].By_end_->point(ii,0) = c_byb * pmlArr_[0].By_end_->point(ii,0) + c_bye * (Ez_->point(nx_-1-ii+1,0)-Ez_->point(nx_-1-ii,0));
                                Hx_->point(nx_-1-ii,0) = c_hxh * Hx_->point(nx_-1-ii,0) + c_hxb1 * pmlArr_[0].Bx_end_->point(ii,0) - c_hxb0 * bxstore;
                                Hy_->point(nx_-1-ii,0) = c_hyh * Hy_->point(nx_-1-ii,0) + c_hyb1 * pmlArr_[0].By_end_->point(ii,0) - c_hyb0 * bystore;
                                //Top Right
                                eps    = objArr_[phys_Hy_->point(nx_-1-ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                bystore = pmlArr_[0].By_end_->point(ii,ny_-1);
                                pmlArr_[0].By_end_->point(ii,ny_-1) = c_byb * pmlArr_[0].By_end_->point(ii,ny_-1) + c_bye * (Ez_->point(nx_-1-ii+1,ny_-1)-Ez_->point(nx_-1-ii,ny_-1));
                                Hy_->point(nx_-1-ii,ny_-1) = c_hyh * Hy_->point(nx_-1-ii,ny_-1) + c_hyb1 * pmlArr_[0].By_end_->point(ii,ny_-1) - c_hyb0 * bystore;
                                //Top Left
                                sigxy = pmlArr_[0].sigma(static_cast<double>(ii) + 0.5);
                                eps    = objArr_[phys_Hy_->point(ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                bystore = pmlArr_[0].By_->point(ii,ny_-1);
                                pmlArr_[0].By_->point(ii,ny_-1) = c_byb * pmlArr_[0].By_->point(ii,ny_-1) + c_bye * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
                                Hy_->point(ii,ny_-1) = c_hyh * Hy_->point(ii,ny_-1) + c_hyb1 * pmlArr_[0].By_->point(ii,ny_ - 1 - 0) - c_hyb0 * bystore;
                            }
                            for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                            {
                                kapz = 1.0; sigz = 0.0; eps = 1.0;
                                kapx = 1.0; kapy = 1.0;
                                sigxx = pmlArr_[0].sigma(0.0);
                                sigxy = pmlArr_[0].sigma(0.5);
                                sigyx = pmlArr_[1].sigma(static_cast<double>(jj) + 0.5);
                                sigyy = pmlArr_[1].sigma(static_cast<double>(jj));

                                eps    = objArr_[phys_Hx_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                // Bot Left
                                double bxstore = pmlArr_[0].Bx_->point(0,jj);
                                double bystore = pmlArr_[0].By_->point(0,jj);
                                pmlArr_[0].Bx_->point(0,jj) = c_bxb * pmlArr_[0].Bx_->point(0,jj) - c_bxe * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                                pmlArr_[0].By_->point(0,jj) = c_byb * pmlArr_[0].By_->point(0,jj) + c_bye * (Ez_->point(0+1,jj)-Ez_->point(0,jj));
                                Hx_->point(0,jj) = c_hxh * Hx_->point(0,jj) + c_hxb1 * pmlArr_[0].Bx_->point(0,jj) - c_hxb0 * bxstore;
                                Hy_->point(0,jj) = c_hyh * Hy_->point(0,jj) + c_hyb1 * pmlArr_[0].By_->point(0,jj) - c_hyb0 * bystore;
                                //Bot Right
                                eps    = objArr_[phys_Hx_->point(nx_-1,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                bxstore = pmlArr_[0].Bx_end_->point(0,jj);
                                pmlArr_[0].Bx_end_->point(0,jj) = c_bxb * pmlArr_[0].Bx_end_->point(0,jj) - c_bxe * (Ez_->point(nx_-1,jj+1)-Ez_->point(nx_-1,jj));
                                Hx_->point(nx_-1,jj) = c_hxh * Hx_->point(nx_-1,jj) + c_hxb1 * pmlArr_[0].Bx_end_->point(0,jj) - c_hxb0 * bxstore;
                                //Top Rigt
                                sigyx = pmlArr_[1].sigma(static_cast<double>(jj) - 0.5);
                                eps    = objArr_[phys_Hx_->point(nx_-1,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                bxstore = pmlArr_[0].Bx_end_->point(0,ny_-1-jj);
                                pmlArr_[0].Bx_end_->point(0,ny_-1-jj) = c_bxb * pmlArr_[0].Bx_end_->point(0,ny_-1-jj) - c_bxe * (Ez_->point(nx_-1,ny_-1-jj+1)-Ez_->point(nx_-1,ny_-1-jj));
                                Hx_->point(nx_-1,ny_-1-jj) = c_hxh * Hx_->point(nx_-1,ny_-1-jj) + c_hxb1 * pmlArr_[0].Bx_end_->point(0,ny_-1-jj) - c_hxb0 * bxstore;
                                //Top Left
                                eps    = objArr_[phys_Hx_->point(0,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(0,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                bxstore = pmlArr_[0].Bx_->point(0,ny_-1-jj);
                                bystore = pmlArr_[0].By_->point(0,ny_-1-jj);
                                pmlArr_[0].Bx_->point(0,ny_-1-jj) = c_bxb * pmlArr_[0].Bx_->point(0,ny_-1-jj) - c_bxe * (Ez_->point(0,ny_-1-jj+1)-Ez_->point(0,ny_-1-jj));
                                pmlArr_[0].By_->point(0,ny_-1-jj) = c_byb * pmlArr_[0].By_->point(0,ny_-1-jj) + c_bye * (Ez_->point(0+1,ny_-1-jj)-Ez_->point(0,ny_-1-jj));
                                Hx_->point(0,ny_-1-jj) = c_hxh * Hx_->point(0,ny_-1-jj) + c_hxb1 * pmlArr_[0].Bx_->point(0,ny_-1-jj) - c_hxb0 * bxstore;
                                Hy_->point(0,ny_-1-jj) = c_hyh * Hy_->point(0,ny_-1-jj) + c_hyb1 * pmlArr_[0].By_->point(0,ny_-1-jj) - c_hyb0 * bystore;
                            }
                            for(int ii = 1; ii < pmlArr_[0].thickness(); ii++)
                            {
                                for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                                {
                                    kapz = 1.0; sigz = 0.0; eps = 1.0;
                                    kapx = 1.0; kapy = 1.0;
                                    sigxx = pmlArr_[0].sigma(static_cast<double>(ii));
                                    sigxy = pmlArr_[0].sigma(static_cast<double>(ii) + 0.5);
                                    sigyx = pmlArr_[1].sigma(static_cast<double>(jj) + 0.5);
                                    sigyy = pmlArr_[1].sigma(static_cast<double>(jj));

                                    eps    = objArr_[phys_Hx_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                    eps    = objArr_[phys_Hy_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                    // Bot Left
                                    double bxstore = pmlArr_[0].Bx_->point(ii,jj);
                                    double bystore = pmlArr_[0].By_->point(ii,jj);
                                    pmlArr_[0].Bx_->point(ii,jj) = c_bxb * pmlArr_[0].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                    pmlArr_[0].By_->point(ii,jj) = c_byb * pmlArr_[0].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[0].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[0].By_->point(ii,jj) - c_hyb0 * bystore;
                                    //Bot Right
                                    sigxy = pmlArr_[0].sigma(static_cast<double>(ii) - 0.5);
                                    eps    = objArr_[phys_Hx_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                    eps    = objArr_[phys_Hy_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                    bxstore = pmlArr_[0].Bx_end_->point(ii,jj);
                                    bystore = pmlArr_[0].By_end_->point(ii,jj);
                                    pmlArr_[0].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[0].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point(nx_-1-ii,jj+1)-Ez_->point(nx_-1-ii,jj));
                                    pmlArr_[0].By_end_->point(ii,jj) = c_byb * pmlArr_[0].By_end_->point(ii,jj) + c_bye * (Ez_->point(nx_-1-ii+1,jj)-Ez_->point(nx_-1-ii,jj));
                                    Hx_->point(nx_-1-ii,jj) = c_hxh * Hx_->point(nx_-1-ii,jj) + c_hxb1 * pmlArr_[0].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                                    Hy_->point(nx_-1-ii,jj) = c_hyh * Hy_->point(nx_-1-ii,jj) + c_hyb1 * pmlArr_[0].By_end_->point(ii,jj) - c_hyb0 * bystore;
                                    //Top Right
                                    sigyx = pmlArr_[1].sigma(static_cast<double>(jj) - 0.5);
                                    eps    = objArr_[phys_Hx_->point(nx_-1-ii,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                    eps    = objArr_[phys_Hy_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                    bxstore = pmlArr_[0].Bx_end_->point(ii,ny_-1-jj);
                                    bystore = pmlArr_[0].By_end_->point(ii,ny_-1-jj);
                                    pmlArr_[0].Bx_end_->point(ii,ny_-1-jj) = c_bxb * pmlArr_[0].Bx_end_->point(ii,ny_-1-jj) - c_bxe * (Ez_->point(nx_-1-ii,ny_-1-jj+1)-Ez_->point(nx_-1-ii,ny_-1-jj));
                                    pmlArr_[0].By_end_->point(ii,ny_-1-jj) = c_byb * pmlArr_[0].By_end_->point(ii,ny_-1-jj) + c_bye * (Ez_->point(nx_-1-ii+1,ny_-1-jj)-Ez_->point(nx_-1-ii,ny_-1-jj));
                                    Hx_->point(nx_-1-ii,ny_-1-jj) = c_hxh * Hx_->point(nx_-1-ii,ny_-1-jj) + c_hxb1 * pmlArr_[0].Bx_end_->point(ii,ny_-1-jj) - c_hxb0 * bxstore;
                                    Hy_->point(nx_-1-ii,ny_-1-jj) = c_hyh * Hy_->point(nx_-1-ii,ny_-1-jj) + c_hyb1 * pmlArr_[0].By_end_->point(ii,ny_-1-jj) - c_hyb0 * bystore;
                                    //Top Left
                                    sigxy = pmlArr_[0].sigma(static_cast<double>(ii) + 0.5);
                                    eps    = objArr_[phys_Hx_->point(ii,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                    eps    = objArr_[phys_Hy_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
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
                    double sigyx = pmlArr_[kk].sigma(0.5);
                    double sigyy = pmlArr_[kk].sigma(0.0);

                    double c_bxb; double c_bxe; double c_hxh; double c_hxb0; double c_hxb1;
                    double c_byb; double c_bye; double c_hyh; double c_hyb0; double c_hyb1;

                    if (xPML_== 0)
                    {
                        for (int jj = 0; jj < nx_-1; jj++)
                        {
                            eps    = objArr_[phys_Hx_->point(jj,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            eps    = objArr_[phys_Hy_->point(jj,0)].dielectric(1.0);
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            //Update both on the bottom side and Hy on the top
                            double bxstore = pmlArr_[kk].Bx_->point(jj,0);
                            double bystore = pmlArr_[kk].By_->point(jj,0);
                            pmlArr_[kk].Bx_->point(jj,0) = c_bxb * pmlArr_[kk].Bx_->point(jj,0) - c_bxe * (Ez_->point(jj,0+1)-Ez_->point(jj,0));
                            pmlArr_[kk].By_->point(jj,0) = c_byb * pmlArr_[kk].By_->point(jj,0) + c_bye * (Ez_->point(jj+1,0)-Ez_->point(jj,0));
                            //cout <<"H1"<<endl;
                            Hx_->point(jj,0) = c_hxh * Hx_->point(jj,0) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,0) - c_hxb0 * bxstore;
                            Hy_->point(jj,0) = c_hyh * Hy_->point(jj,0) + c_hyb1 * pmlArr_[kk].By_->point(jj,0) - c_hyb0 * bystore;

                            eps    = objArr_[phys_Hy_->point(jj,ny_-1)].dielectric(1.0);
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            bystore = pmlArr_[kk].By_end_->point(jj,0);
                            pmlArr_[kk].By_end_->point(jj,0) = c_byb * pmlArr_[kk].By_end_->point(jj,0) + c_bye * (Ez_->point(jj+1, ny_-1)-Ez_->point(jj,ny_-1));
                            Hy_->point(jj,ny_-1) = c_hyh * Hy_->point(jj,ny_-1) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,0) - c_hyb0 * bystore;
                        }
                        //Top left corner
                        eps    = objArr_[phys_Hx_->point(nx_-1,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                        c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                        c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                        c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                        c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                        c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                        double bxstore = pmlArr_[kk].Bx_->point(nx_-1,0);
                        pmlArr_[kk].Bx_->point(nx_-1,0) = c_bxb * pmlArr_[kk].Bx_->point(nx_-1,0) - c_bxe * (Ez_->point(nx_-1,0+1)-Ez_->point(nx_-1,0));
                        Hx_->point(nx_-1,0) = c_hxh * Hx_->point(nx_-1,0) + c_hxb1 * pmlArr_[kk].Bx_->point(nx_-1,0) - c_hxb0 * bxstore;
                        for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                        {
                            sigyy = pmlArr_[kk].sigma(static_cast<double>(ii));
                            sigyx = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5);
                            eps    = objArr_[phys_Hx_->point(nx_-1,ii)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            //right col
                            bxstore = pmlArr_[kk].Bx_->point(nx_-1,ii);
                            pmlArr_[kk].Bx_->point(nx_-1,ii) = c_bxb * pmlArr_[kk].Bx_->point(nx_-1,ii) - c_bxe * (Ez_->point(nx_-1,ii+1)-Ez_->point(nx_-1,ii));
                            Hx_->point(nx_-1,ii) = c_hxh * Hx_->point(nx_-1,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(nx_-1,ii) - c_hxb0 * bxstore;
                            for(int jj = 0; jj < nx_-1; jj ++)
                            {
                                //Update everything
                                sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) + 0.5));
                                eps    = objArr_[phys_Hx_->point(jj,ii)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(jj,ii)].dielectric(1.0);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                bxstore = pmlArr_[kk].Bx_->point(jj,ii);
                                double bystore = pmlArr_[kk].By_->point(jj,ii);
                                pmlArr_[kk].Bx_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_->point(jj,ii) - c_bxe * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                                pmlArr_[kk].By_->point(jj,ii) = c_byb * pmlArr_[kk].By_->point(jj,ii) + c_bye * (Ez_->point(jj+1,ii)-Ez_->point(jj,ii));
                                //cout <<"H1"<<endl;
                                Hx_->point(jj,ii) = c_hxh * Hx_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,ii) - c_hxb0 * bxstore;
                                Hy_->point(jj,ii) = c_hyh * Hy_->point(jj,ii) + c_hyb1 * pmlArr_[kk].By_->point(jj,ii) - c_hyb0 * bystore;

                                sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) - 0.5));
                                eps    = objArr_[phys_Hx_->point(jj,ny_-1-ii)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(jj,ny_-1-ii)].dielectric(1.0);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                bxstore = pmlArr_[kk].Bx_end_->point(jj,ii);
                                bystore = pmlArr_[kk].By_end_->point(jj,ii);
                                pmlArr_[kk].Bx_end_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_end_->point(jj,ii) - c_bxe * (Ez_->point(jj,ny_-1-ii+1)-Ez_->point(jj,ny_-1-ii));
                                pmlArr_[kk].By_end_->point(jj,ii) = c_byb * pmlArr_[kk].By_end_->point(jj,ii) + c_bye * (Ez_->point(jj+1,ny_-1-ii)-Ez_->point(jj,ny_-1-ii));
                                //cout <<"H2"<<endl;
                                Hx_->point(jj,ny_-1-ii) = c_hxh * Hx_->point(jj,ny_-1-ii) + c_hxb1 * pmlArr_[kk].Bx_end_->point(jj,ii) - c_hxb0 * bxstore;
                                Hy_->point(jj,ny_-1-ii) = c_hyh * Hy_->point(jj,ny_-1-ii) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,ii) - c_hyb0 * bystore;
                            }
                            // for the right only upadate Hx
                            eps    = objArr_[phys_Hx_->point(nx_-1,ny_-1-ii)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            bxstore = pmlArr_[kk].Bx_end_->point(nx_-1,ii);
                            pmlArr_[kk].Bx_end_->point(nx_-1,ii) = c_bxb * pmlArr_[kk].Bx_end_->point(nx_-1,ii) - c_bxe * (Ez_->point(nx_-1,ny_-1-ii+1)-Ez_->point(nx_-1,ny_-1-ii));
                            Hx_->point(nx_-1,ny_-1-ii) = c_hxh * Hx_->point(nx_-1,ny_-1-ii) + c_hxb1 * pmlArr_[kk].Bx_end_->point(nx_-1,ii) - c_hxb0 * bxstore;
                        }
                    }
                    else
                    {
                        for (int jj = xPML_; jj < nx_-xPML_; jj++)
                        {
                            eps    = objArr_[phys_Hx_->point(jj,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            eps    = objArr_[phys_Hy_->point(jj,0)].dielectric(1.0);
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            //Update both on the bottom side and Hy on the top
                            double bxstore = pmlArr_[kk].Bx_->point(jj,0);
                            double bystore = pmlArr_[kk].By_->point(jj,0);
                            pmlArr_[kk].Bx_->point(jj,0) = c_bxb * pmlArr_[kk].Bx_->point(jj,0) - c_bxe * (Ez_->point(jj,0+1)-Ez_->point(jj,0));
                            pmlArr_[kk].By_->point(jj,0) = c_byb * pmlArr_[kk].By_->point(jj,0) + c_bye * (Ez_->point(jj+1,0)-Ez_->point(jj,0));
                            //cout <<"H1"<<endl;
                            Hx_->point(jj,0) = c_hxh * Hx_->point(jj,0) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,0) - c_hxb0 * bxstore;
                            Hy_->point(jj,0) = c_hyh * Hy_->point(jj,0) + c_hyb1 * pmlArr_[kk].By_->point(jj,0) - c_hyb0 * bystore;

                            eps    = objArr_[phys_Hy_->point(jj,ny_-1)].dielectric(1.0);
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            bystore = pmlArr_[kk].By_end_->point(jj,0);
                            pmlArr_[kk].By_end_->point(jj,0) = c_byb * pmlArr_[kk].By_end_->point(jj,0) + c_bye * (Ez_->point(jj+1,ny_-1)-Ez_->point(jj,ny_-1));
                            Hy_->point(jj,ny_-1) = c_hyh * Hy_->point(jj,ny_-1) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,0) - c_hyb0 * bystore;
                        }
                        for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                        {
                            sigyy = pmlArr_[kk].sigma(static_cast<double>(ii));
                            sigyx = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5);

                            for(int jj = xPML_; jj < nx_-xPML_; jj ++)
                            {
                                //Update everything
                                sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) + 0.5));
                                eps    = objArr_[phys_Hx_->point(jj,ii)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(jj,ii)].dielectric(1.0);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                double bxstore = pmlArr_[kk].Bx_->point(jj,ii);
                                double bystore = pmlArr_[kk].By_->point(jj,ii);
                                pmlArr_[kk].Bx_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_->point(jj,ii) - c_bxe * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                                pmlArr_[kk].By_->point(jj,ii) = c_byb * pmlArr_[kk].By_->point(jj,ii) + c_bye * (Ez_->point(jj+1,ii)-Ez_->point(jj,ii));
                                //cout <<"H1"<<endl;
                                Hx_->point(jj,ii) = c_hxh * Hx_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,ii) - c_hxb0 * bxstore;
                                Hy_->point(jj,ii) = c_hyh * Hy_->point(jj,ii) + c_hyb1 * pmlArr_[kk].By_->point(jj,ii) - c_hyb0 * bystore;

                                sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) - 0.5));
                                eps    = objArr_[phys_Hx_->point(jj,ny_-1-ii)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(jj,ny_-1-ii)].dielectric(1.0);
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

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
                            sigxx = pmlArr_[1].sigma(0.0);
                            sigxy = pmlArr_[1].sigma(0.5);
                            sigyx = pmlArr_[0].sigma(0.5);
                            sigyy = pmlArr_[0].sigma(0.0);

                            eps    = objArr_[phys_Hx_->point(0,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            eps    = objArr_[phys_Hy_->point(0,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                            // Bot Left
                            double bxstore = pmlArr_[1].Bx_->point(0,0);
                            double bystore = pmlArr_[1].By_->point(0,0);
                            pmlArr_[1].Bx_->point(0,0) = c_bxb * pmlArr_[1].Bx_->point(0,0) - c_bxe * (Ez_->point(0,0+1)-Ez_->point(0,0));
                            pmlArr_[1].By_->point(0,0) = c_byb * pmlArr_[1].By_->point(0,0) + c_bye * (Ez_->point(0+1,0)-Ez_->point(0,0));
                            Hx_->point(0,0) = c_hxh * Hx_->point(0,0) + c_hxb1 * pmlArr_[1].Bx_->point(0,0) - c_hxb0 * bxstore;
                            Hy_->point(0,0) = c_hyh * Hy_->point(0,0) + c_hyb1 * pmlArr_[1].By_->point(0,0) - c_hyb0 * bystore;
                            //Bot Right
                            eps    = objArr_[phys_Hx_->point(nx_-1,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                            c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                            c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                            c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                            bxstore = pmlArr_[1].Bx_end_->point(0,0);
                            pmlArr_[1].Bx_end_->point(0,0) = c_bxb * pmlArr_[1].Bx_end_->point(0,0) - c_bxe * (Ez_->point(nx_-1,0+1)-Ez_->point(nx_-1,0));
                            Hx_->point(nx_-1,0) = c_hxh * Hx_->point(nx_-1,0) + c_hxb1 * pmlArr_[1].Bx_end_->point(0,0) - c_hxb0 * bxstore;
                            //Top Left

                            eps    = objArr_[phys_Hy_->point(0,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                            c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                            c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                            c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                            c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                            bystore = pmlArr_[1].By_->point(0,ny_-1);
                            pmlArr_[1].By_->point(0,ny_-1) = c_byb * pmlArr_[1].By_->point(0,ny_-1) + c_bye * (Ez_->point(0+1,ny_-1)-Ez_->point(0,ny_-1));
                            Hy_->point(0,ny_-1) = c_hyh * Hy_->point(0,ny_-1) + c_hyb1 * pmlArr_[1].By_->point(0,ny_-1) - c_hyb0 * bystore;
                            for(int ii = 1; ii < pmlArr_[1].thickness(); ii++)
                            {
                                kapz = 1.0; sigz = 0.0; eps = 1.0;
                                kapx = 1.0; kapy = 1.0;
                                sigxx = pmlArr_[1].sigma(static_cast<double>(ii));
                                sigxy = pmlArr_[1].sigma(static_cast<double>(ii) + 0.5);
                                sigyx = pmlArr_[0].sigma(0.5);
                                sigyy = pmlArr_[0].sigma(0.0);

                                eps    = objArr_[phys_Hx_->point(ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                // Bot Left
                                double bxstore = pmlArr_[1].Bx_->point(ii,0);
                                double bystore = pmlArr_[1].By_->point(ii,0);
                                pmlArr_[1].Bx_->point(ii,0) = c_bxb * pmlArr_[1].Bx_->point(ii,0) - c_bxe * (Ez_->point(ii,0+1)-Ez_->point(ii,0));
                                pmlArr_[1].By_->point(ii,0) = c_byb * pmlArr_[1].By_->point(ii,0) + c_bye * (Ez_->point(ii+1,0)-Ez_->point(ii,0));
                                Hx_->point(ii,0) = c_hxh * Hx_->point(ii,0) + c_hxb1 * pmlArr_[1].Bx_->point(ii,0) - c_hxb0 * bxstore;
                                Hy_->point(ii,0) = c_hyh * Hy_->point(ii,0) + c_hyb1 * pmlArr_[1].By_->point(ii,0) - c_hyb0 * bystore;
                                //Bot Right
                                sigxy = pmlArr_[1].sigma(static_cast<double>(ii) - 0.5);
                                eps    = objArr_[phys_Hx_->point(ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                bxstore = pmlArr_[1].Bx_end_->point(ii,0);
                                bystore = pmlArr_[1].By_end_->point(ii,0);
                                pmlArr_[1].Bx_end_->point(ii,0) = c_bxb * pmlArr_[1].Bx_end_->point(ii,0) - c_bxe * (Ez_->point(nx_-1-ii,0+1)-Ez_->point(nx_-1-ii,0));
                                pmlArr_[1].By_end_->point(ii,0) = c_byb * pmlArr_[1].By_end_->point(ii,0) + c_bye * (Ez_->point(nx_-1-ii+1,0)-Ez_->point(nx_-1-ii,0));
                                Hx_->point(nx_-1-ii,0) = c_hxh * Hx_->point(nx_-1-ii,0) + c_hxb1 * pmlArr_[1].Bx_end_->point(ii,0) - c_hxb0 * bxstore;
                                Hy_->point(nx_-1-ii,0) = c_hyh * Hy_->point(nx_-1-ii,0) + c_hyb1 * pmlArr_[1].By_end_->point(ii,0) - c_hyb0 * bystore;
                                //Top Right
                                eps    = objArr_[phys_Hy_->point(nx_-ii-1,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                bystore = pmlArr_[1].By_end_->point(ii,ny_-1);
                                pmlArr_[1].By_end_->point(ii,ny_-1) = c_byb * pmlArr_[1].By_end_->point(ii,ny_-1) + c_bye * (Ez_->point(nx_-1-ii+1,ny_-1)-Ez_->point(nx_-1-ii,ny_-1));
                                Hy_->point(nx_-1-ii,ny_-1) = c_hyh * Hy_->point(nx_-1-ii,ny_-1) + c_hyb1 * pmlArr_[1].By_end_->point(ii,ny_-1) - c_hyb0 * bystore;
                                //Top Left
                                sigxy = pmlArr_[1].sigma(static_cast<double>(ii) + 0.5);
                                eps    = objArr_[phys_Hy_->point(ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                bystore = pmlArr_[1].By_->point(ii,ny_-1);
                                pmlArr_[1].By_->point(ii,ny_-1) = c_byb * pmlArr_[1].By_->point(ii,ny_-1) + c_bye * (Ez_->point(ii+1,ny_-1)-Ez_->point(ii,ny_-1));
                                Hy_->point(ii,ny_-1) = c_hyh * Hy_->point(ii,ny_-1) + c_hyb1 * pmlArr_[1].By_->point(ii,ny_ - 1 - 0) - c_hyb0 * bystore;
                            }
                            for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                            {
                                kapz = 1.0; sigz = 0.0; eps = 1.0;
                                kapx = 1.0; kapy = 1.0;
                                sigxx = pmlArr_[1].sigma(0.0);
                                sigxy = pmlArr_[1].sigma(0.5);
                                sigyx = pmlArr_[0].sigma(static_cast<double>(jj) + 0.5);
                                sigyy = pmlArr_[0].sigma(static_cast<double>(jj));

                                eps    = objArr_[phys_Hx_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                // Bot Left
                                double bxstore = pmlArr_[1].Bx_->point(0,jj);
                                double bystore = pmlArr_[1].By_->point(0,jj);
                                pmlArr_[1].Bx_->point(0,jj) = c_bxb * pmlArr_[1].Bx_->point(0,jj) - c_bxe * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                                pmlArr_[1].By_->point(0,jj) = c_byb * pmlArr_[1].By_->point(0,jj) + c_bye * (Ez_->point(0+1,jj)-Ez_->point(0,jj));
                                Hx_->point(0,jj) = c_hxh * Hx_->point(0,jj) + c_hxb1 * pmlArr_[1].Bx_->point(0,jj) - c_hxb0 * bxstore;
                                Hy_->point(0,jj) = c_hyh * Hy_->point(0,jj) + c_hyb1 * pmlArr_[1].By_->point(0,jj) - c_hyb0 * bystore;
                                //Bot Right
                                eps    = objArr_[phys_Hx_->point(nx_-1,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                bxstore = pmlArr_[1].Bx_end_->point(0,jj);
                                pmlArr_[1].Bx_end_->point(0,jj) = c_bxb * pmlArr_[1].Bx_end_->point(0,jj) - c_bxe * (Ez_->point(nx_-1,jj+1)-Ez_->point(nx_-1,jj));
                                Hx_->point(nx_-1,jj) = c_hxh * Hx_->point(nx_-1,jj) + c_hxb1 * pmlArr_[1].Bx_end_->point(0,jj) - c_hxb0 * bxstore;
                                //Top Rigt
                                sigyx = pmlArr_[0].sigma(static_cast<double>(jj) - 0.5);
                                eps    = objArr_[phys_Hx_->point(nx_-1,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                bxstore = pmlArr_[1].Bx_end_->point(0,ny_-1-jj);
                                pmlArr_[1].Bx_end_->point(0,ny_-1-jj) = c_bxb * pmlArr_[1].Bx_end_->point(0,ny_-1-jj) - c_bxe * (Ez_->point(nx_-1,ny_-1-jj+1)-Ez_->point(nx_-1,ny_-1-jj));
                                Hx_->point(nx_-1,ny_-1-jj) = c_hxh * Hx_->point(nx_-1,ny_-1-jj) + c_hxb1 * pmlArr_[1].Bx_end_->point(0,ny_-1-jj) - c_hxb0 * bxstore;
                                //Top Left
                                eps    = objArr_[phys_Hx_->point(0,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                eps    = objArr_[phys_Hy_->point(0,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
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
                                    sigxx = pmlArr_[1].sigma(static_cast<double>(ii));
                                    sigxy = pmlArr_[1].sigma(static_cast<double>(ii) + 0.5);
                                    sigyx = pmlArr_[0].sigma(static_cast<double>(jj) + 0.5);
                                    sigyy = pmlArr_[0].sigma(static_cast<double>(jj));

                                    eps    = objArr_[phys_Hx_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                    eps    = objArr_[phys_Hy_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                    // Bot Left
                                    double bxstore = pmlArr_[1].Bx_->point(ii,jj);
                                    double bystore = pmlArr_[1].By_->point(ii,jj);
                                    pmlArr_[1].Bx_->point(ii,jj) = c_bxb * pmlArr_[1].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                    pmlArr_[1].By_->point(ii,jj) = c_byb * pmlArr_[1].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                                    Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[1].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                                    Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[1].By_->point(ii,jj) - c_hyb0 * bystore;
                                    //Bot Right
                                    sigxy = pmlArr_[1].sigma(static_cast<double>(ii) - 0.5);
                                    eps    = objArr_[phys_Hx_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                    eps    = objArr_[phys_Hy_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;

                                    bxstore = pmlArr_[1].Bx_end_->point(ii,jj);
                                    bystore = pmlArr_[1].By_end_->point(ii,jj);
                                    pmlArr_[1].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[1].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point(nx_-1-ii,jj+1)-Ez_->point(nx_-1-ii,jj));
                                    pmlArr_[1].By_end_->point(ii,jj) = c_byb * pmlArr_[1].By_end_->point(ii,jj) + c_bye * (Ez_->point(nx_-1-ii+1,jj)-Ez_->point(nx_-1-ii,jj));
                                    Hx_->point(nx_-1-ii,jj) = c_hxh * Hx_->point(nx_-1-ii,jj) + c_hxb1 * pmlArr_[1].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                                    Hy_->point(nx_-1-ii,jj) = c_hyh * Hy_->point(nx_-1-ii,jj) + c_hyb1 * pmlArr_[1].By_end_->point(ii,jj) - c_hyb0 * bystore;
                                    //Top Right
                                    sigyx = pmlArr_[0].sigma(static_cast<double>(jj) - 0.5);
                                    eps    = objArr_[phys_Hx_->point(nx_-1-ii,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                    eps    = objArr_[phys_Hy_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                    bxstore = pmlArr_[1].Bx_end_->point(ii,ny_-1-jj);
                                    bystore = pmlArr_[1].By_end_->point(ii,ny_-1-jj);
                                    pmlArr_[1].Bx_end_->point(ii,ny_-1-jj) = c_bxb * pmlArr_[1].Bx_end_->point(ii,ny_-1-jj) - c_bxe * (Ez_->point(nx_-1-ii,ny_-1-jj+1)-Ez_->point(nx_-1-ii,ny_-1-jj));
                                    pmlArr_[1].By_end_->point(ii,ny_-1-jj) = c_byb * pmlArr_[1].By_end_->point(ii,ny_-1-jj) + c_bye * (Ez_->point(nx_-1-ii+1,ny_-1-jj)-Ez_->point(nx_-1-ii,ny_-1-jj));
                                    Hx_->point(nx_-1-ii,ny_-1-jj) = c_hxh * Hx_->point(nx_-1-ii,ny_-1-jj) + c_hxb1 * pmlArr_[1].Bx_end_->point(ii,ny_-1-jj) - c_hxb0 * bxstore;
                                    Hy_->point(nx_-1-ii,ny_-1-jj) = c_hyh * Hy_->point(nx_-1-ii,ny_-1-jj) + c_hyb1 * pmlArr_[1].By_end_->point(ii,ny_-1-jj) - c_hyb0 * bystore;
                                    //Top Left
                                    sigxy = pmlArr_[1].sigma(static_cast<double>(ii) + 0.5);
                                    eps    = objArr_[phys_Hx_->point(ii,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_bxb  = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                    c_bxe  = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_));
                                    c_hxh  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                                    c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                                    eps    = objArr_[phys_Hy_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_byb  = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                                    c_bye  = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_));
                                    c_hyh  = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                    c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                    c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
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
    else
    {

    }
    if(Ez_)
    {
        double eps;
        double c_eze = 1.0;
        double c_ezh;
        if(xPML_ != 0 && yPML_ !=0)
        {
            for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
                {
                    eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                    c_ezh = dt_/(eps*dx_);
                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                }
            }
        }
        else if(xPML_ != 0)
        {
            for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                for(int jj = 1; jj < ny_ - 1; jj ++)
                {
                    eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                    c_ezh = dt_/(eps*dx_);
                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                }
            }
            for(int ii = xPML_; ii < nx_ - xPML_; ii ++)
            {
                eps = objArr_[phys_Ez_->point(ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                c_ezh = dt_/(eps*dx_);
                Ez_->point(ii,0) = c_eze * Ez_->point(ii,0) + c_ezh * ((Hy_->point(ii,0)-Hy_->point(ii-1,0)) - (Hx_->point(ii,0)-0));
                eps = objArr_[phys_Ez_->point(ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                c_ezh = dt_/(eps*dx_);
                Ez_->point(ii,ny_-1) = c_eze * Ez_->point(ii,ny_-1) + c_ezh * ((Hy_->point(ii,ny_-1)-Hy_->point(ii-1,ny_-1)) - (0-Hx_->point(ii,ny_-1-1)));
            }
        }
        else if(yPML_ != 0)
        {
            for(int ii = 1; ii < nx_ - 1; ii ++)
            {
                for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
                {
                    eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                    c_ezh = dt_/(eps*dx_);
                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                }
            }
            for(int jj = yPML_; jj < ny_ - yPML_; jj ++)
            {
                eps = objArr_[phys_Ez_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                c_ezh = dt_/(eps*dx_);
                Ez_->point(0,jj) = c_eze * Ez_->point(0,jj) + c_ezh * ((Hy_->point(0,jj)-0) - (Hx_->point(0,jj)-Hx_->point(0,jj-1)));
                eps = objArr_[phys_Ez_->point(nx_-1,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                c_ezh = dt_/(eps*dx_);
                Ez_->point(nx_-1,jj) = c_eze * Ez_->point(nx_-1,jj) + c_ezh * ((0-Hy_->point(nx_-1-1,jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1)));
            }
        }
        else
        {
            for(int ii = 1; ii < nx_ - 1; ii ++)
            {
                for(int jj = 1; jj < ny_ - 1; jj ++)
                {
                    eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                    c_ezh = dt_/(eps*dx_);
                    Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                }
            }
            for(int jj = 1; jj < nx_ - 1; jj ++)
            {
                eps = objArr_[phys_Ez_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                c_ezh = dt_/(eps*dx_);
                Ez_->point(0,jj) = c_eze * Ez_->point(0,jj) + c_ezh * ((Hy_->point(0,jj)-0) - (Hx_->point(0,jj)-Hx_->point(0,jj-1)));
                eps = objArr_[phys_Ez_->point(nx_-1,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                c_ezh = dt_/(eps*dx_);
                Ez_->point(nx_-1,jj) = c_eze * Ez_->point(nx_-1,jj) + c_ezh * ((0-Hy_->point(nx_-1-1,jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1)));
            }
            for(int ii = 1; ii < nx_ - 1; ii ++)
            {
                eps = objArr_[phys_Ez_->point(ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                c_ezh = dt_/(eps*dx_);
                Ez_->point(ii,0) = c_eze * Ez_->point(ii,0) + c_ezh * ((Hy_->point(ii,0)-Hy_->point(ii-1,0)) - (Hx_->point(ii,0)-0));
                eps = objArr_[phys_Ez_->point(ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                c_ezh = dt_/(eps*dx_);
                Ez_->point(ii,ny_-1) = c_eze * Ez_->point(ii,ny_-1) + c_ezh * ((Hy_->point(ii,ny_-1)-Hy_->point(ii-1,ny_-1)) - (0-Hx_->point(ii,ny_-1-1)));
            }
            eps = objArr_[phys_Ez_->point(0,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
            c_ezh = dt_/(eps*dx_);
            Ez_->point(0,0) = c_eze * Ez_->point(0,0) + c_ezh * ((Hy_->point(0,0)-0) - (Hx_->point(0,0)-0));
            eps = objArr_[phys_Ez_->point(0,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
            c_ezh = dt_/(eps*dx_);
            Ez_->point(0,ny_-1) = c_eze * Ez_->point(0,ny_-1) + c_ezh * ((Hy_->point(0,ny_-1)-0) - (0-Hx_->point(0,ny_-1-1)));
            eps = objArr_[phys_Ez_->point(nx_-1,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
            c_ezh = dt_/(eps*dx_);
            Ez_->point(nx_-1,0) = c_eze * Ez_->point(nx_-1,0) + c_ezh * ((0-Hy_->point(nx_-1-1,0)) - (Hx_->point(nx_-1,0)-0));
            eps = objArr_[phys_Ez_->point(nx_-1,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
            c_ezh = dt_/(eps*dx_);
            Ez_->point(nx_-1,ny_-1) = c_eze * Ez_->point(nx_-1,ny_-1) + c_ezh * ((0-Hy_->point(nx_-1-1,ny_-1)) - (0-Hx_->point(nx_-1,ny_-1-1)));
        }
    }
    else
    {

    }
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
                    double sigx = pmlArr_[kk].sigma(0.0);
                    double sigy = 0.0;
                    double sigz = 0.0;
                    double c_dzd; double c_dzh; double c_eze; double c_ezd1; double c_ezd0;
                    if(yPML_ == 0)
                    {
                        for(int jj = 1; jj < ny_ - 1; jj++)
                        {
                            eps = objArr_[phys_Ez_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;

                            double dzstore = pmlArr_[kk].Dz_->point(0,jj);
                            pmlArr_[kk].Dz_->point(0,jj) = c_dzd * pmlArr_[kk].Dz_->point(0,jj) + c_dzh * ((Hy_->point(0,jj) - 0) - (Hx_->point(0,jj)-Hx_->point(0,jj-1))); // 0 is boundary condtion
                            Ez_->point(0,jj) = c_eze * Ez_->point(0,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(0,jj) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(nx_-1,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;

                            dzstore = pmlArr_[kk].Dz_end_->point(0,jj);
                            pmlArr_[kk].Dz_end_->point(0,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(0,jj) + c_dzh * ((0 -Hy_->point((nx_-2),jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1))); // 0 is boundary condition
                            Ez_->point(nx_-1,jj) = c_eze * Ez_->point(nx_-1,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(0,jj) - c_ezd0 * dzstore;
                        }
                        eps = objArr_[phys_Ez_->point(0,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;

                        double dzstore = pmlArr_[kk].Dz_->point(0,ny_-1);
                        pmlArr_[kk].Dz_->point(0,ny_-1) = c_dzd * pmlArr_[kk].Dz_->point(0,ny_-1) + c_dzh * ((Hy_->point(0,ny_-1) - 0) - (0-Hx_->point(0,ny_-2))); // 0 is boundary condtion
                        Ez_->point(0,ny_-1) = c_eze * Ez_->point(0,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_->point(0,ny_-1) - c_ezd0 * dzstore;

                        eps = objArr_[phys_Ez_->point(0,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        dzstore = pmlArr_[kk].Dz_end_->point(0,ny_-1);
                        pmlArr_[kk].Dz_end_->point(0,ny_-1) = c_dzd * pmlArr_[kk].Dz_end_->point(0,ny_-1) + c_dzh * ((0 -Hy_->point(nx_-2,ny_-1)) - (0-Hx_->point(nx_-1,ny_-2))); // 0 is boundary condition
                        Ez_->point(nx_-1,ny_-1) = c_eze * Ez_->point(nx_-1,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_end_->point(0,ny_-1) - c_ezd0 * dzstore;

                        eps = objArr_[phys_Ez_->point(0,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        dzstore = pmlArr_[kk].Dz_->point(0,0);
                        pmlArr_[kk].Dz_->point(0,0) = c_dzd * pmlArr_[kk].Dz_->point(0,0) + c_dzh * ((Hy_->point(0,0) - 0) - (Hx_->point(0,0)-0)); // 0 is boundary condtion
                        Ez_->point(0,0) = c_eze * Ez_->point(0,0) + c_ezd1 * pmlArr_[kk].Dz_->point(0,0) - c_ezd0 * dzstore;

                        eps = objArr_[phys_Ez_->point(nx_-1,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        dzstore = pmlArr_[kk].Dz_end_->point(0,0);
                        pmlArr_[kk].Dz_end_->point(0,0) = c_dzd * pmlArr_[kk].Dz_end_->point(0,0) + c_dzh * ((0 -Hy_->point(nx_-2,0)) - (Hx_->point(nx_-1,0)-0)); // 0 is boundary condition
                        Ez_->point(nx_-1,0) = c_eze * Ez_->point(nx_-1,0) + c_ezd1 * pmlArr_[kk].Dz_end_->point(0,0) - c_ezd0 * dzstore;

                        for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                        {
                            kapx = 1.0; kapy = 1.0; kapz = 1.0;
                            sigx = pmlArr_[kk].sigma(static_cast<double>(ii));
                            sigy = 0.0; sigz = 0.0;
                            //double kapx = pmlArr_[kk].kappa(ii);
                            //Kappas change throughout
                            for(int jj =  1; jj < ny_ - 1; jj ++)
                            {
                                eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                double dzstore = pmlArr_[kk].Dz_->point(ii,jj);
                                pmlArr_[kk].Dz_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,jj) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[kk].Dz_end_->point(ii,jj);
                                pmlArr_[kk].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point(nx_-1-ii,jj)-Hy_->point(nx_-1-ii-1,jj)) - (Hx_->point(nx_-1-ii,jj)-Hx_->point(nx_-1-ii,jj-1)));
                                Ez_->point(nx_-1-ii,jj) = c_eze * Ez_->point(nx_-1-ii,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                            }
                            eps = objArr_[phys_Ez_->point(ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            double dzstore = pmlArr_[kk].Dz_->point(ii,ny_-1);
                            pmlArr_[kk].Dz_->point(ii,ny_-1) = c_dzd * pmlArr_[kk].Dz_->point(ii,ny_-1) + c_dzh * ((Hy_->point(ii,ny_-1)-Hy_->point(ii-1,ny_-1)) - (0-Hx_->point(ii,ny_-1-1)));
                            Ez_->point(ii,ny_-1) = c_eze * Ez_->point(ii,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,ny_-1) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(nx_-1-ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_end_->point(ii,ny_-1);
                            pmlArr_[kk].Dz_end_->point(ii,ny_-1) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,ny_-1) + c_dzh * ((Hy_->point(nx_-1-ii,ny_-1)-Hy_->point(nx_-1-ii-1,ny_-1)) - (0-Hx_->point(nx_-1-ii,ny_-1-1)));
                            Ez_->point(nx_-1-ii,ny_-1) = c_eze * Ez_->point(nx_-1-ii,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,ny_-1) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_->point(ii,0);
                            pmlArr_[kk].Dz_->point(ii,0) = c_dzd * pmlArr_[kk].Dz_->point(ii,0) + c_dzh * ((Hy_->point(ii,0)-Hy_->point(ii-1,0)) - (Hx_->point(ii,0)-0));
                            Ez_->point(ii,0) = c_eze * Ez_->point(ii,0) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,0) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(nx_-1-ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_end_->point(ii,0);
                            pmlArr_[kk].Dz_end_->point(ii,0) = c_dzd * pmlArr_[kk].Dz_end_->point(ii,0) + c_dzh * ((Hy_->point(nx_-1-ii,0)-Hy_->point(nx_-1-ii-1,0)) - (Hx_->point(nx_-1-ii,0)-0));
                            Ez_->point(nx_-1-ii,0) = c_eze * Ez_->point(nx_-1-ii,0) + c_ezd1 * pmlArr_[kk].Dz_end_->point(ii,0) - c_ezd0 * dzstore;
                        }
                    }
                    else
                    {
                        for(int jj = yPML_; jj < ny_ - yPML_; jj++)
                        {
                            eps = objArr_[phys_Ez_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            double dzstore = pmlArr_[kk].Dz_->point(0,jj);
                            pmlArr_[kk].Dz_->point(0,jj) = c_dzd * pmlArr_[kk].Dz_->point(0,jj) + c_dzh * ((Hy_->point(0,jj) - 0) - (Hx_->point(0,jj)-Hx_->point(0,jj-1))); // 0 is boundary condtion
                            Ez_->point(0,jj) = c_eze * Ez_->point(0,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(0,jj) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(nx_-1,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_end_->point(0,jj);
                            pmlArr_[kk].Dz_end_->point(0,jj) = c_dzd * pmlArr_[kk].Dz_end_->point(0,jj) + c_dzh * ((0 -Hy_->point((nx_-2),jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1))); // 0 is boundary condition
                            Ez_->point(nx_-1,jj) = c_eze * Ez_->point(nx_-1,jj) + c_ezd1 * pmlArr_[kk].Dz_end_->point(0,jj) - c_ezd0 * dzstore;
                        }
                        for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                        {
                            kapx = 1.0; kapy = 1.0; kapz = 1.0;
                            sigx = pmlArr_[kk].sigma(static_cast<double>(ii));
                            sigy = 0.0; sigz = 0.0;
                            for(int jj =  yPML_; jj < ny_ - yPML_; jj ++)
                            {
                                eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                double dzstore = pmlArr_[kk].Dz_->point(ii,jj);
                                pmlArr_[kk].Dz_->point(ii,jj) = c_dzd * pmlArr_[kk].Dz_->point(ii,jj) + c_dzh * ((Hy_->point(ii,jj)-Hy_->point(ii-1,jj)) - (Hx_->point(ii,jj)-Hx_->point(ii,jj-1)));
                                Ez_->point(ii,jj) = c_eze * Ez_->point(ii,jj) + c_ezd1 * pmlArr_[kk].Dz_->point(ii,jj) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
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
                            sigx = pmlArr_[0].sigma(0);
                            sigy = pmlArr_[1].sigma(0);
                            sigz = 0.0;

                            eps = objArr_[phys_Ez_->point(0,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            //Bot Left
                            double dzstore = pmlArr_[0].Dz_->point(0,0);
                            pmlArr_[0].Dz_->point(0,0) = c_dzd * pmlArr_[0].Dz_->point(0,0) + c_dzh * ((Hy_->point(0,0)-0) - (Hx_->point(0,0)-0));
                            Ez_->point(0,0) = c_eze * Ez_->point(0,0) + c_ezd1 * pmlArr_[0].Dz_->point(0,0) - c_ezd0 * dzstore;
                            //Bot Right
                            eps = objArr_[phys_Ez_->point(nx_-1,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[0].Dz_end_->point(0,0);
                            pmlArr_[0].Dz_end_->point(0,0) = c_dzd * pmlArr_[0].Dz_end_->point(0,0) + c_dzh * (0-Hy_->point(nx_-1-1,0)) - (Hx_->point(nx_-1,0)-0);
                            Ez_->point(nx_-1,0) = c_eze * Ez_->point(nx_-1,0) + c_ezd1 * pmlArr_[0].Dz_end_->point(0,0) - c_ezd0 * dzstore;
                            //Top Right
                            eps = objArr_[phys_Ez_->point(nx_-1,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[0].Dz_end_->point(0,ny_-1);
                            pmlArr_[0].Dz_end_->point(0,ny_-1) = c_dzd * pmlArr_[0].Dz_end_->point(0,ny_-1) + c_dzh * ((0-Hy_->point(nx_-1-1,ny_-1)) - (0-Hx_->point(nx_-1,ny_-1-1)));
                            Ez_->point(nx_-1,ny_-1) = c_eze * Ez_->point(nx_-1,ny_-1) + c_ezd1 * pmlArr_[0].Dz_end_->point(0,ny_-1) - c_ezd0 * dzstore;
                            //Top Left
                            eps = objArr_[phys_Ez_->point(0,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[0].Dz_->point(0,ny_-1);
                            pmlArr_[0].Dz_->point(0,ny_-1) = c_dzd * pmlArr_[0].Dz_->point(0,ny_-1) + c_dzh * ((Hy_->point(0,ny_-1)-0) - (0-Hx_->point(0,ny_-1-1)));
                            Ez_->point(0,ny_-1) = c_eze * Ez_->point(0,ny_-1) + c_ezd1 * pmlArr_[0].Dz_->point(0,ny_-1) - c_ezd0 * dzstore;

                            for(int ii = 1; ii < pmlArr_[0].thickness(); ii++)
                            {
                                kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                sigx = pmlArr_[0].sigma(static_cast<double>(ii));
                                sigy = pmlArr_[1].sigma(0.0);
                                sigz = 0.0;

                                eps = objArr_[phys_Ez_->point(ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                //Bot Left
                                dzstore = pmlArr_[0].Dz_->point(ii,0);
                                pmlArr_[0].Dz_->point(ii,0) = c_dzd * pmlArr_[0].Dz_->point(ii,0) + c_dzh * ((Hy_->point(ii,0)-Hy_->point(ii-1,0)) - (Hx_->point(ii,0)-0));
                                Ez_->point(ii,0) = c_eze * Ez_->point(ii,0) + c_ezd1 * pmlArr_[0].Dz_->point(ii,0) - c_ezd0 * dzstore;
                                //Bot Right
                                eps = objArr_[phys_Ez_->point(nx_-1-ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[0].Dz_end_->point(ii,0);
                                pmlArr_[0].Dz_end_->point(ii,0) = c_dzd * pmlArr_[0].Dz_end_->point(ii,0) + c_dzh * ((Hy_->point(nx_-1-ii,0)-Hy_->point(nx_-1-ii-1,0)) - (Hx_->point(nx_-1-ii,0)-0));
                                Ez_->point(nx_-1-ii,0) = c_eze * Ez_->point(nx_-1-ii,0) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,0) - c_ezd0 * dzstore;
                                //Top Right
                                eps = objArr_[phys_Ez_->point(ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[0].Dz_end_->point(ii,ny_-1);
                                pmlArr_[0].Dz_end_->point(ii,ny_-1) = c_dzd * pmlArr_[0].Dz_end_->point(ii,ny_-1) + c_dzh * ((Hy_->point(nx_-1-ii,ny_-1)-Hy_->point(nx_-1-ii-1,ny_-1)) - (0-Hx_->point(nx_-1-ii,ny_-1-1)));
                                Ez_->point(nx_-1-ii,ny_-1) = c_eze * Ez_->point(nx_-1-ii,ny_-1) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,ny_-1) - c_ezd0 * dzstore;
                                //Top Left
                                eps = objArr_[phys_Ez_->point(nx_-1-ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[0].Dz_->point(ii,ny_-1);
                                pmlArr_[0].Dz_->point(ii,ny_-1) = c_dzd * pmlArr_[0].Dz_->point(ii,ny_-1) + c_dzh * ((Hy_->point(ii,ny_-1)-Hy_->point(ii-1,ny_-1)) - (0-Hx_->point(ii,ny_-1-1)));
                                Ez_->point(ii,ny_-1) = c_eze * Ez_->point(ii,ny_-1) + c_ezd1 * pmlArr_[0].Dz_->point(ii,ny_-1) - c_ezd0 * dzstore;
                            }
                            for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                            {
                                kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                sigx = pmlArr_[0].sigma(static_cast<double>(0));
                                sigy = pmlArr_[1].sigma(static_cast<double>(jj));
                                sigz = 0.0;

                                eps = objArr_[phys_Ez_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                //Bot Left
                                dzstore = pmlArr_[0].Dz_->point(0,jj);
                                pmlArr_[0].Dz_->point(0,jj) = c_dzd * pmlArr_[0].Dz_->point(0,jj) + c_dzh * ((Hy_->point(0,jj)-0) - (Hx_->point(0,jj)-Hx_->point(0,jj-1)));
                                Ez_->point(0,jj) = c_eze * Ez_->point(0,jj) + c_ezd1 * pmlArr_[0].Dz_->point(0,jj) - c_ezd0 * dzstore;
                                //Bot Right
                                eps = objArr_[phys_Ez_->point(nx_-1,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[0].Dz_end_->point(0,jj);
                                pmlArr_[0].Dz_end_->point(0,jj) = c_dzd * pmlArr_[0].Dz_end_->point(0,jj) + c_dzh * ((0-Hy_->point(nx_-1-1,jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1)));
                                Ez_->point(nx_-1,jj) = c_eze * Ez_->point(nx_-1,jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(0,jj) - c_ezd0 * dzstore;
                                //Top Right
                                eps = objArr_[phys_Ez_->point(nx_-1,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[0].Dz_end_->point(0,ny_-1-jj);
                                pmlArr_[0].Dz_end_->point(0,ny_-1-jj) = c_dzd * pmlArr_[0].Dz_end_->point(0,ny_-1-jj) + c_dzh * ((0-Hy_->point((nx_-1) - 0-1,ny_-1-jj)) - (Hx_->point((nx_-1) - 0,ny_-1-jj)-Hx_->point((nx_-1) - 0,ny_-1-jj-1)));
                                Ez_->point((nx_-1) - 0,ny_-1-jj) = c_eze * Ez_->point((nx_-1) - 0,ny_-1-jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(0,ny_-1-jj) - c_ezd0 * dzstore;
                                //Top Left
                                eps = objArr_[phys_Ez_->point(0,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[0].Dz_->point(0,ny_-1-jj);
                                pmlArr_[0].Dz_->point(0,ny_-1-jj) = c_dzd * pmlArr_[0].Dz_->point(0,ny_-1-jj) + c_dzh * ((Hy_->point(0,ny_-1-jj)-0) - (Hx_->point(0,ny_-1-jj)-Hx_->point(0,ny_-1-jj-1)));
                                Ez_->point(0,ny_-1-jj) = c_eze * Ez_->point(0,ny_-1-jj) + c_ezd1 * pmlArr_[0].Dz_->point(0,ny_-1-jj) - c_ezd0 * dzstore;
                            }
                            for(int ii = 1; ii < pmlArr_[0].thickness(); ii++)
                            {
                                for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                                {
                                    kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                    sigx = pmlArr_[0].sigma(static_cast<double>(ii));
                                    sigy = pmlArr_[1].sigma(static_cast<double>(jj));
                                    sigz = 0.0;

                                    eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
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
                                    eps = objArr_[phys_Ez_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[0].Dz_end_->point(ii,jj);
                                    pmlArr_[0].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point(nx_-1-ii,jj)-Hy_->point(nx_-1-ii-1,jj)) - (Hx_->point(nx_-1-ii,jj)-Hx_->point(nx_-1-ii,jj-1)));
                                    Ez_->point(nx_-1-ii,jj) = c_eze * Ez_->point(nx_-1-ii,jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                    //Top Right
                                    eps = objArr_[phys_Ez_->point(nx_-1-ii,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[0].Dz_end_->point(ii,ny_-1-jj);
                                    pmlArr_[0].Dz_end_->point(ii,ny_-1-jj) = c_dzd * pmlArr_[0].Dz_end_->point(ii,ny_-1-jj) + c_dzh * ((Hy_->point(nx_-1-ii,ny_-1-jj)-Hy_->point(nx_-1-ii-1,ny_-1-jj)) - (Hx_->point(nx_-1-ii,ny_-1-jj)-Hx_->point(nx_-1-ii,ny_-1-jj-1)));
                                    Ez_->point(nx_-1-ii,ny_-1-jj) = c_eze * Ez_->point(nx_-1-ii,ny_-1-jj) + c_ezd1 * pmlArr_[0].Dz_end_->point(ii,ny_-1-jj) - c_ezd0 * dzstore;
                                    //Top Left
                                    eps = objArr_[phys_Ez_->point(ii,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
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
                    double sigy = pmlArr_[kk].sigma(0.0); double sigz = 0.0;

                    double c_dzd; double c_dzh; double c_eze; double c_ezd1; double c_ezd0;
                    if(xPML_ == 0)
                    {
                        for(int jj = 1; jj < nx_ - 1; jj++)
                        {
                            eps = objArr_[phys_Ez_->point(jj,0)].dielectric(1.0);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;

                            double dzstore = pmlArr_[kk].Dz_->point(jj,0);
                            pmlArr_[kk].Dz_->point(jj,0) = c_dzd * pmlArr_[kk].Dz_->point(jj,0) + c_dzh * ((Hy_->point(jj,0)-Hy_->point(jj-1,0)) - (Hx_->point(jj,0)- 0)); // 0 is boundary condition
                            Ez_->point(jj,0) = c_eze * Ez_->point(jj,0) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,0) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(jj,ny_-1)].dielectric(1.0);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_end_->point(jj,0);
                            pmlArr_[kk].Dz_end_->point(jj,0) = c_dzd * pmlArr_[kk].Dz_end_->point(jj,0) + c_dzh * ((Hy_->point(jj,ny_-1) - Hy_->point(jj - 1,ny_-1)) - (0 - Hx_->point(jj,ny_-1-1))); // 0 is boundary condtion
                            Ez_->point(jj,ny_-1) = c_eze * Ez_->point(jj,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_end_->point(jj,0) - c_ezd0 * dzstore;
                        }
                        eps = objArr_[phys_Ez_->point(nx_-1,0)].dielectric(1.0);
                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        double dzstore = pmlArr_[kk].Dz_->point(nx_-1,0);
                        pmlArr_[kk].Dz_->point(nx_-1,0) = c_dzd * pmlArr_[kk].Dz_->point(nx_-1,0) + c_dzh * ((0-Hy_->point(nx_-1-1,0)) - (Hx_->point(nx_-1,0)- 0)); // 0 is boundary condition
                        Ez_->point(nx_-1,0) = c_eze * Ez_->point(nx_-1,0) + c_ezd1 * pmlArr_[kk].Dz_->point(nx_-1,0) - c_ezd0 * dzstore;

                        eps = objArr_[phys_Ez_->point(nx_-1,ny_-1)].dielectric(1.0);
                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        dzstore = pmlArr_[kk].Dz_end_->point(nx_-1,0);
                        pmlArr_[kk].Dz_end_->point(nx_-1,0) = c_dzd * pmlArr_[kk].Dz_end_->point(nx_-1,0) + c_dzh * ((0 - Hy_->point(nx_-1 - 1,ny_-1)) - (0 - Hx_->point(nx_-1,ny_-1-1))); // 0 is boundary condtion
                        Ez_->point(nx_-1,ny_-1) = c_eze * Ez_->point(nx_-1,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_end_->point(nx_-1,0) - c_ezd0 * dzstore;

                        eps = objArr_[phys_Ez_->point(0,0)].dielectric(1.0);
                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        dzstore = pmlArr_[kk].Dz_->point(0,0);
                        pmlArr_[kk].Dz_->point(0,0) = c_dzd * pmlArr_[kk].Dz_->point(0,0) + c_dzh * ((Hy_->point(0,0)-0) - (Hx_->point(0,0)- 0)); // 0 is boundary condition
                        Ez_->point(0,0) = c_eze * Ez_->point(0,0) + c_ezd1 * pmlArr_[kk].Dz_->point(0,0) - c_ezd0 * dzstore;

                        eps = objArr_[phys_Ez_->point(0,ny_-1)].dielectric(1.0);
                        c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                        c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                        c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                        c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                        dzstore = pmlArr_[kk].Dz_end_->point(0,0);
                        pmlArr_[kk].Dz_end_->point(0,0) = c_dzd * pmlArr_[kk].Dz_end_->point(0,0) + c_dzh * ((Hy_->point(0,ny_-1) - 0) - (0 - Hx_->point(0,ny_-1-1))); // 0 is boundary condtion
                        Ez_->point(0,ny_-1) = c_eze * Ez_->point(0,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_end_->point(0,0) - c_ezd0 * dzstore;

                        for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                        {
                            eps = 1.0;
                            kapx = 1.0; kapy = 1.0;
                            kapz = 1.0;
                            sigx = 0.0;
                            sigy = pmlArr_[kk].sigma(static_cast<double>(ii));
                            sigz = 0.0;

                            for(int jj =  1; jj < nx_ - 1; jj ++)
                            {
                                eps = objArr_[phys_Ez_->point(jj,ii)].dielectric(1.0);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                double dzstore = pmlArr_[kk].Dz_->point(jj,ii);
                                pmlArr_[kk].Dz_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_->point(jj,ii) + c_dzh * ((Hy_->point(jj,ii)-Hy_->point(jj-1,ii)) - (Hx_->point(jj,ii)-Hx_->point(jj,ii-1)));
                                Ez_->point(jj,ii) = c_eze * Ez_->point(jj,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,ii) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(jj,ny_-1-ii)].dielectric(1.0);
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
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            double dzstore = pmlArr_[kk].Dz_->point(nx_-1,ii);
                            pmlArr_[kk].Dz_->point(nx_-1,ii) = c_dzd * pmlArr_[kk].Dz_->point(nx_-1,ii) + c_dzh * ((0-Hy_->point(nx_-1-1,ii)) - (Hx_->point(nx_-1,ii)-Hx_->point(nx_-1,ii-1)));
                            Ez_->point(nx_-1,ii) = c_eze * Ez_->point(nx_-1,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(nx_-1,ii) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(nx_-1, ny_-1-ii)].dielectric(1.0);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_end_->point(nx_-1,ii);
                            pmlArr_[kk].Dz_end_->point(nx_-1,ii) = c_dzd * pmlArr_[kk].Dz_end_->point(nx_-1,ii) + c_dzh * ((0 - Hy_->point(nx_-1-1,ny_-1-ii)) - (Hx_->point(nx_-1,ny_-1-ii) - Hx_->point(nx_-1,ny_-1-ii-1)));
                            Ez_->point(nx_-1,ny_-1-ii) = c_eze * Ez_->point(nx_-1,ny_-1-ii) + c_ezd1 * pmlArr_[kk].Dz_end_->point(nx_-1,ii) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(0, ii)].dielectric(1.0);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_->point(0,ii);
                            pmlArr_[kk].Dz_->point(0,ii) = c_dzd * pmlArr_[kk].Dz_->point(0,ii) + c_dzh * ((Hy_->point(0,ii)-0) - (Hx_->point(0,ii)-Hx_->point(0,ii-1)));
                            Ez_->point(0,ii) = c_eze * Ez_->point(0,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(0,ii) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(0, ny_-1-ii)].dielectric(1.0);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_end_->point(0,ii);
                            pmlArr_[kk].Dz_end_->point(0,ii) = c_dzd * pmlArr_[kk].Dz_end_->point(0,ii) + c_dzh * ((Hy_->point(0,ny_-1-ii) - 0) - (Hx_->point(0,ny_-1-ii) - Hx_->point(0,ny_-1-ii-1)));
                            Ez_->point(0,ny_-1-ii) = c_eze * Ez_->point(0,ny_-1-ii) + c_ezd1 * pmlArr_[kk].Dz_end_->point(0,ii) - c_ezd0 * dzstore;
                        }
                    }
                    else
                    {
                        for(int jj = xPML_; jj < nx_ - xPML_; jj++)
                        {
                            eps = objArr_[phys_Ez_->point(jj,0)].dielectric(1.0);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            double dzstore = pmlArr_[kk].Dz_->point(jj,0);
                            pmlArr_[kk].Dz_->point(jj,0) = c_dzd * pmlArr_[kk].Dz_->point(jj,0) + c_dzh * ((Hy_->point(jj,0)-Hy_->point(jj-1,0)) - (Hx_->point(jj,0)- 0)); // 0 is boundary condition
                            Ez_->point(jj,0) = c_eze * Ez_->point(jj,0) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,0) - c_ezd0 * dzstore;

                            eps = objArr_[phys_Ez_->point(jj, ny_-1)].dielectric(1.0);
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[kk].Dz_end_->point(jj,0);
                            pmlArr_[kk].Dz_end_->point(jj,0) = c_dzd * pmlArr_[kk].Dz_end_->point(jj,0) + c_dzh * ((Hy_->point(jj,ny_-1) - Hy_->point(jj - 1,ny_-1)) - (0 - Hx_->point(jj,ny_-1-1))); // 0 is boundary condtion
                            Ez_->point(jj,ny_-1) = c_eze * Ez_->point(jj,ny_-1) + c_ezd1 * pmlArr_[kk].Dz_end_->point(jj,0) - c_ezd0 * dzstore;
                        }
                        for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                        {
                            kapx = 1.0; kapy = 1.0;
                            kapz = 1.0;
                            sigx = 0.0;
                            sigy = pmlArr_[kk].sigma(static_cast<double>(ii));
                            sigz = 0.0;

                            for(int jj =  xPML_; jj < nx_ - xPML_; jj ++)
                            {
                                eps = objArr_[phys_Ez_->point(jj, ii)].dielectric(1.0);
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                double dzstore = pmlArr_[kk].Dz_->point(jj,ii);
                                pmlArr_[kk].Dz_->point(jj,ii) = c_dzd * pmlArr_[kk].Dz_->point(jj,ii) + c_dzh * ((Hy_->point(jj,ii)-Hy_->point(jj-1,ii)) - (Hx_->point(jj,ii)-Hx_->point(jj,ii-1)));
                                Ez_->point(jj,ii) = c_eze * Ez_->point(jj,ii) + c_ezd1 * pmlArr_[kk].Dz_->point(jj,ii) - c_ezd0 * dzstore;

                                eps = objArr_[phys_Ez_->point(jj, ny_-1-ii)].dielectric(1.0);
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
                            sigx = pmlArr_[1].sigma(0);
                            sigy = pmlArr_[0].sigma(0);
                            sigz = 0.0;

                            eps = objArr_[phys_Ez_->point(0,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            //Bot Left
                            double dzstore = pmlArr_[1].Dz_->point(0,0);
                            pmlArr_[1].Dz_->point(0,0) = c_dzd * pmlArr_[1].Dz_->point(0,0) + c_dzh * ((Hy_->point(0,0)-0) - (Hx_->point(0,0)-0));
                            Ez_->point(0,0) = c_eze * Ez_->point(0,0) + c_ezd1 * pmlArr_[1].Dz_->point(0,0) - c_ezd0 * dzstore;
                            //Bot Right
                            eps = objArr_[phys_Ez_->point(nx_-1,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[1].Dz_end_->point(0,0);
                            pmlArr_[1].Dz_end_->point(0,0) = c_dzd * pmlArr_[1].Dz_end_->point(0,0) + c_dzh * (0-Hy_->point(nx_-1-1,0)) - (Hx_->point(nx_-1,0)-0);
                            Ez_->point(nx_-1,0) = c_eze * Ez_->point(nx_-1,0) + c_ezd1 * pmlArr_[1].Dz_end_->point(0,0) - c_ezd0 * dzstore;
                            //Top Right
                            eps = objArr_[phys_Ez_->point(nx_-1,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[1].Dz_end_->point(0,ny_-1);
                            pmlArr_[1].Dz_end_->point(0,ny_-1) = c_dzd * pmlArr_[1].Dz_end_->point(0,ny_-1) + c_dzh * ((0-Hy_->point(nx_-1-1,ny_-1)) - (0-Hx_->point(nx_-1,ny_-1-1)));
                            Ez_->point(nx_-1,ny_-1) = c_eze * Ez_->point(nx_-1,ny_-1) + c_ezd1 * pmlArr_[1].Dz_end_->point(0,ny_-1) - c_ezd0 * dzstore;
                            //Top Left
                            eps = objArr_[phys_Ez_->point(0,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                            c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                            c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                            c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                            c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                            dzstore = pmlArr_[1].Dz_->point(0,ny_-1);
                            pmlArr_[1].Dz_->point(0,ny_-1) = c_dzd * pmlArr_[1].Dz_->point(0,ny_-1) + c_dzh * ((Hy_->point(0,ny_-1)-0) - (0-Hx_->point(0,ny_-1-1)));
                            Ez_->point(0,ny_-1) = c_eze * Ez_->point(0,ny_-1) + c_ezd1 * pmlArr_[1].Dz_->point(0,ny_-1) - c_ezd0 * dzstore;

                            for(int ii = 1; ii < pmlArr_[1].thickness(); ii++)
                            {
                                kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                sigx = pmlArr_[1].sigma(static_cast<double>(ii));
                                sigy = pmlArr_[0].sigma(0.0);
                                sigz = 0.0;

                                eps = objArr_[phys_Ez_->point(ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                //Bot Left
                                dzstore = pmlArr_[1].Dz_->point(ii,0);
                                pmlArr_[1].Dz_->point(ii,0) = c_dzd * pmlArr_[1].Dz_->point(ii,0) + c_dzh * ((Hy_->point(ii,0)-Hy_->point(ii-1,0)) - (Hx_->point(ii,0)-0));
                                Ez_->point(ii,0) = c_eze * Ez_->point(ii,0) + c_ezd1 * pmlArr_[1].Dz_->point(ii,0) - c_ezd0 * dzstore;
                                //Bot Right
                                eps = objArr_[phys_Ez_->point(nx_-1-ii,0)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[1].Dz_end_->point(ii,0);
                                pmlArr_[1].Dz_end_->point(ii,0) = c_dzd * pmlArr_[1].Dz_end_->point(ii,0) + c_dzh * ((Hy_->point(nx_-1-ii,0)-Hy_->point(nx_-1-ii-1,0)) - (Hx_->point(nx_-1-ii,0)-0));
                                Ez_->point(nx_-1-ii,0) = c_eze * Ez_->point(nx_-1-ii,0) + c_ezd1 * pmlArr_[1].Dz_end_->point(ii,0) - c_ezd0 * dzstore;
                                //Top Right
                                eps = objArr_[phys_Ez_->point(ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[1].Dz_end_->point(ii,ny_-1);
                                pmlArr_[1].Dz_end_->point(ii,ny_-1) = c_dzd * pmlArr_[1].Dz_end_->point(ii,ny_-1) + c_dzh * ((Hy_->point(nx_-1-ii,ny_-1)-Hy_->point(nx_-1-ii-1,ny_-1)) - (0-Hx_->point(nx_-1-ii,ny_-1-1)));
                                Ez_->point(nx_-1-ii,ny_-1) = c_eze * Ez_->point(nx_-1-ii,ny_-1) + c_ezd1 * pmlArr_[1].Dz_end_->point(ii,ny_-1) - c_ezd0 * dzstore;
                                //Top Left
                                eps = objArr_[phys_Ez_->point(nx_-1-ii,ny_-1)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[1].Dz_->point(ii,ny_-1);
                                pmlArr_[1].Dz_->point(ii,ny_-1) = c_dzd * pmlArr_[1].Dz_->point(ii,ny_-1) + c_dzh * ((Hy_->point(ii,ny_-1)-Hy_->point(ii-1,ny_-1)) - (0-Hx_->point(ii,ny_-1-1)));
                                Ez_->point(ii,ny_-1) = c_eze * Ez_->point(ii,ny_-1) + c_ezd1 * pmlArr_[1].Dz_->point(ii,ny_-1) - c_ezd0 * dzstore;
                            }
                            for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                            {
                                kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                sigx = pmlArr_[1].sigma(static_cast<double>(0));
                                sigy = pmlArr_[0].sigma(static_cast<double>(jj));
                                sigz = 0.0;

                                eps = objArr_[phys_Ez_->point(0,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                //Bot Left
                                dzstore = pmlArr_[1].Dz_->point(0,jj);
                                pmlArr_[1].Dz_->point(0,jj) = c_dzd * pmlArr_[1].Dz_->point(0,jj) + c_dzh * ((Hy_->point(0,jj)-0) - (Hx_->point(0,jj)-Hx_->point(0,jj-1)));
                                Ez_->point(0,jj) = c_eze * Ez_->point(0,jj) + c_ezd1 * pmlArr_[1].Dz_->point(0,jj) - c_ezd0 * dzstore;
                                //Bot Right
                                eps = objArr_[phys_Ez_->point(nx_-1,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[1].Dz_end_->point(0,jj);
                                pmlArr_[1].Dz_end_->point(0,jj) = c_dzd * pmlArr_[1].Dz_end_->point(0,jj) + c_dzh * ((0-Hy_->point(nx_-1-1,jj)) - (Hx_->point(nx_-1,jj)-Hx_->point(nx_-1,jj-1)));
                                Ez_->point(nx_-1,jj) = c_eze * Ez_->point(nx_-1,jj) + c_ezd1 * pmlArr_[1].Dz_end_->point(0,jj) - c_ezd0 * dzstore;
                                //Top Right
                                eps = objArr_[phys_Ez_->point(nx_-1,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[1].Dz_end_->point(0,ny_-1-jj);
                                pmlArr_[1].Dz_end_->point(0,ny_-1-jj) = c_dzd * pmlArr_[1].Dz_end_->point(0,ny_-1-jj) + c_dzh * ((0-Hy_->point((nx_-1) - 0-1,ny_-1-jj)) - (Hx_->point((nx_-1) - 0,ny_-1-jj)-Hx_->point((nx_-1) - 0,ny_-1-jj-1)));
                                Ez_->point((nx_-1) - 0,ny_-1-jj) = c_eze * Ez_->point((nx_-1) - 0,ny_-1-jj) + c_ezd1 * pmlArr_[1].Dz_end_->point(0,ny_-1-jj) - c_ezd0 * dzstore;
                                //Top Left
                                eps = objArr_[phys_Ez_->point(0,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                dzstore = pmlArr_[1].Dz_->point(0,ny_-1-jj);
                                pmlArr_[1].Dz_->point(0,ny_-1-jj) = c_dzd * pmlArr_[1].Dz_->point(0,ny_-1-jj) + c_dzh * ((Hy_->point(0,ny_-1-jj)-0) - (Hx_->point(0,ny_-1-jj)-Hx_->point(0,ny_-1-jj-1)));
                                Ez_->point(0,ny_-1-jj) = c_eze * Ez_->point(0,ny_-1-jj) + c_ezd1 * pmlArr_[1].Dz_->point(0,ny_-1-jj) - c_ezd0 * dzstore;
                            }
                            for(int ii = 1; ii < pmlArr_[1].thickness(); ii++)
                            {
                                for(int jj = 1; jj < pmlArr_[1].thickness(); jj++)
                                {
                                    kapx = 1.0; kapy = 1.0; kapz = 1.0;
                                    sigx = pmlArr_[1].sigma(static_cast<double>(ii));
                                    sigy = pmlArr_[0].sigma(static_cast<double>(jj));
                                    sigz = 0.0;

                                    eps = objArr_[phys_Ez_->point(ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
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
                                    eps = objArr_[phys_Ez_->point(nx_-1-ii,jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[1].Dz_end_->point(ii,jj);
                                    pmlArr_[1].Dz_end_->point(ii,jj) = c_dzd * pmlArr_[1].Dz_end_->point(ii,jj) + c_dzh * ((Hy_->point(nx_-1-ii,jj)-Hy_->point(nx_-1-ii-1,jj)) - (Hx_->point(nx_-1-ii,jj)-Hx_->point(nx_-1-ii,jj-1)));
                                    Ez_->point(nx_-1-ii,jj) = c_eze * Ez_->point(nx_-1-ii,jj) + c_ezd1 * pmlArr_[1].Dz_end_->point(ii,jj) - c_ezd0 * dzstore;
                                    //Top Right
                                    eps = objArr_[phys_Ez_->point(nx_-1-ii,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
                                    c_dzd = (2*eps*kapx - sigx*dt_) / (2*eps*kapx + sigx*dt_);
                                    c_dzh = 2 * eps * dt_ / (dy_ * (2*eps*kapx + sigx*dt_));
                                    c_eze = (2*eps*kapy - sigy*dt_) / (2*eps*kapy + sigy*dt_);
                                    c_ezd1 = (2*eps*kapz + sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    c_ezd0 = (2*eps*kapz - sigz*dt_) / (2*eps*kapy + sigy*dt_) / eps;
                                    dzstore = pmlArr_[1].Dz_end_->point(ii,ny_-1-jj);
                                    pmlArr_[1].Dz_end_->point(ii,ny_-1-jj) = c_dzd * pmlArr_[1].Dz_end_->point(ii,ny_-1-jj) + c_dzh * ((Hy_->point(nx_-1-ii,ny_-1-jj)-Hy_->point(nx_-1-ii-1,ny_-1-jj)) - (Hx_->point(nx_-1-ii,ny_-1-jj)-Hx_->point(nx_-1-ii,ny_-1-jj-1)));
                                    Ez_->point(nx_-1-ii,ny_-1-jj) = c_eze * Ez_->point(nx_-1-ii,ny_-1-jj) + c_ezd1 * pmlArr_[1].Dz_end_->point(ii,ny_-1-jj) - c_ezd0 * dzstore;
                                    //Top Left
                                    eps = objArr_[phys_Ez_->point(ii,ny_-1-jj)].dielectric(1.0); ///< Will become freq dep once that becomes added into the equations
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
    for(int ii = 0; ii < dtcArr_.size(); ii ++)
        ouputField(dtcArr_[ii]);

    /*string fname("fout/Hx/HxField_t" + to_string(static_cast<int>(tcur_/dt_))+".dat");
    Hx_->gridOut(fname);
    fname = "fout/Hy/HyField_t" + to_string(static_cast<int>(tcur_/dt_))+".dat";
    Hy_->gridOut(fname);
    fname = "fout/Ez/EzField_t" + to_string(static_cast<int>(tcur_/dt_))+".dat";
    Ez_->gridOut(fname);*/

    tcur_ += dt_;
}
