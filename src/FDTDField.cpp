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
 * @details Creates the FDTD cell and related fields. The step sizes are set to 1/res, and the number of grid cells are calculated by multiplying the step size with the step sizes.
 * All fields have a one cell border around them to model the correct boundary conditions
 *
 * @param IP Input object created from an input file
 * @param dx_ and dy_: step size in the x and y direction respectively
 * @param dt_ time step
 * @param nx_ and ny_ number of steps in the x and y direction respectively
 * @param srcArr_ list of all the sources in the cell
 * @param objArr_ list of all the objects in the cell
 * @param pmlArr_ list of all the PML's n the cell
 * @param xPML_ and yPML_ thickness of the PML regions in number of grid points
 * @param periodic_ switch for implementing PBC
 * @param k_point_ the k vector for the PBC to correct the phase
 *
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
    nx_         = floor(static_cast<double>(res_) * static_cast<double>(IP.x_size_)+ 0.5) + 1;
    ny_         = floor(static_cast<double>(res_) * static_cast<double>(IP.y_size_)+ 0.5) + 1;
    srcArr_     = IP.srcArr_;
    objArr_     = IP.objArr_;
    dtcArr_     = IP.dctArr_;
    pmlArr_     = IP.pmlArr_;
    xPML_       = IP.xPml_;
    yPML_       = IP.yPml_;
    precalcPML_ = IP.pmlCalc_;
    periodic_   = IP.periodic_;
    k_point_    = IP.k_point_;
    if(IP.invCell_)
    {
        xDTC_   = 1;
        yDTC_   = 0;
    }
    else
    {
        xDTC_   = 0;
        yDTC_   = 1;
    }

    if(IP.pol_.compare("Hz") == 0 || IP.pol_.compare("Ey") == 0 || IP.pol_.compare("Ex") == 0)
    {
        Ex_      = make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_);
        Ey_      = make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_);
        Dx_      = make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_);
        Dy_      = make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_);
        Hz_      = make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_);
        phys_Ex_ = make_shared<Grid2D<int>>(nx_+2,ny_+2,dx_,dy_);
        phys_Ey_ = make_shared<Grid2D<int>>(nx_+2,ny_+2,dx_,dy_);
        //These are never used in the TE mode
        Hx_      = nullptr;
        Hy_      = nullptr;
        Ez_      = nullptr;
        Dz_      = nullptr;
        phys_Ez_ = nullptr;
    }
    else
    {
        Hx_ = make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_);
        Hy_ = make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_);
        Ez_ = make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_);
        Dz_ = make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_);

        phys_Ez_ = make_shared<Grid2D<int>>(nx_+2,ny_+2,dx_,dy_);
        // These are never used in the TM mode
        Ex_      = nullptr;
        Dx_      = nullptr;
        Ey_      = nullptr;
        Dy_      = nullptr;
        Hz_      = nullptr;
        phys_Ex_ = nullptr;
        phys_Ey_ = nullptr;
    }
}

/**
 * @brief Initializes the physical grid for materials look up
 * @details Initializes the lists that will be used to run ?axpy from the MKL libraries. First maps out where all the materials are, then constructs lists to run them all.
 * Ends by initializing the PMLs
 *
 */

void FDTDField::initializeGrid()
{
    // determines the boundaries of what to initialize
    int xmin = 0; int xmax = 0; int ymin =0; int ymax = 0;
    if(yPML_ != 0 && xPML_ != 0)
    {
        xmin = xPML_+1;
        xmax = nx_-xPML_+1;
        ymin = yPML_+1;
        ymax = ny_-yPML_+1;
    }
    else if(yPML_ != 0)
    {
        xmin = 1;
        xmax = nx_+1;
        ymin = yPML_+1;
        ymax = ny_-yPML_+1;
    }
    else if(xPML_ != 0)
    {
        xmin = xPML_+1;
        xmax = nx_-xPML_+1;
        ymin = 1;
        ymax = ny_+1;
    }
    else
    {
        xmin = 1;
        xmax = nx_+1;
        ymin = 1;
        ymax = ny_+1;
    }
    for(int kk = 0; kk < objArr_.size(); kk++)
    {
        // maps out where the objects are in the grid; Since no magnetically active materials are being modeled H fields don't need to be updated
        // Objects overwrite each other, last object listed gets modeled
        vector<double> pt(2,0.0);
        if(Hz_)
        {

            for(int ii = xmin; ii < xmax;ii ++)
            {
                for(int jj = ymin; jj < ymax; jj ++)
                {
                    pt[0] = ((ii-1) +0.5-(nx_-1)/2.0)*dx_;
                    pt[1] = ((jj-1)-static_cast<double>(ny_-1)/2.0)*dy_;
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ex_->point(ii,jj) = kk;
                    pt[0] -= 0.5*dx_;
                    pt[1] += 0.5*dy_;
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ey_->point(ii,jj) = kk;
                }
            }
        }
        else
        {
            for(int ii = xmin; ii < xmax; ii ++)
            {
                for(int jj = ymin; jj < ymax; jj ++)
                {
                    pt[0] = ((ii-1)-(nx_-1)/2.0)*dx_;
                    pt[1] = ((jj-1)-static_cast<double>(ny_-1)/2.0)*dy_;
                    if(objArr_[kk].isObj(pt)==true)
                        phys_Ez_->point(ii,jj) = kk;
                }
            }
        }
        // Sets up all constants and necessary fields
        objArr_[kk].setUpConsts(dt_);

        if(Hz_)
        {
            while((objArr_[kk].mat().size() - 1) / 3.0 > lorPx_.size())
            {
                    lorPx_.push_back(make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_));
                prevLorPx_.push_back(make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_));
                    lorPy_.push_back(make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_));
                prevLorPy_.push_back(make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_));
            }
            cout << lorPx_.size() << endl;
        }
        else
        {
            while((objArr_[kk].mat().size() - 1) / 3.0 > lorPz_.size())
            {
                    lorPz_.push_back(make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_));
                prevLorPz_.push_back(make_shared<Grid2D<complex<double>>>(nx_+2,ny_+2,dx_,dy_));
            }
        }
    }
    // D fields are for dispersion only, if no dispersive materials then no D-fields
    if(lorPx_.size() == 0)
    {
        Dx_ = nullptr;
        Dy_ = nullptr;
    }
    else if(lorPz_.size() ==0)
        Dz_ = nullptr;

    //Set up the parameters for the MKL calls for the Electric fields, conditionals are because of edge cases.
    if(Ez_)
    {
        for(int jj = ymin; jj < ymax; jj++)
        {
            int ii = xmin;
            while(ii < xmax)
            {
                int iistore = ii;
                while(ii < xmax-1 && phys_Ez_ -> point(ii,jj) == phys_Ez_ -> point(ii+1,jj) )
                    ii ++;
                array<int,4> tempArr = { iistore,jj,ii-iistore+1,phys_Ez_->point(iistore,jj)};
                if(objArr_[tempArr[3]].mat().size() <= 1)
                    zaxEz_.push_back(tempArr);
                else
                    zaxDz_.push_back(tempArr);
                ii++;
            }
        }
        // No magnetic Materials means H-Fields don't care about materials, cells have to account for difference in materials
        if(xPML_ != 0)
        {
            for(int jj = ymin; jj < ymax; jj++)
            {
                array<int,4> tempArr = {xPML_+1, jj, static_cast<int>(nx_)-2*xPML_,0};
                zaxHx_.push_back(tempArr);
                zaxHy_.push_back(tempArr);
            }
        }
        else
        {
            for(int jj = ymin; jj < ymax; jj++)
            {
                array<int,4> tempArr = {1, jj, static_cast<int>(nx_),0};
                zaxHx_.push_back(tempArr);
                tempArr = {1, jj, static_cast<int>(nx_-1),0};
                zaxHy_.push_back(tempArr);
            }
        }
        if(yPML_ == 0)
            zaxHx_.pop_back();
    }
    else
    {
        // sub is to account for field size differences

        int x_sub = 0; int y_sub = 0;

        if(yPML_ == 0)
            y_sub = 1;
        if(xPML_ == 0)
            x_sub = 1;

        for(int jj = ymin; jj < ymax; jj++)
        {
            int ii = xmin;
            while(ii < xmax-x_sub)
            {
                int iistore = ii;
                while(ii < xmax-1-x_sub && phys_Ex_ -> point(ii,jj) == phys_Ex_ -> point(ii+1,jj) )
                    ii ++;
                array<int,4> tempArr = { iistore,jj,ii-iistore+1,phys_Ex_->point(iistore,jj)};
                if(objArr_[tempArr[3]].mat().size() <= 1)
                    zaxEx_.push_back(tempArr);
                else
                    zaxDx_.push_back(tempArr);
                ii++;
            }
        }
        for(int jj = ymin; jj < ymax-y_sub; jj++)
        {
            int ii = xmin;
            while(ii < xmax)
            {
                int iistore = ii;
                while(ii < xmax-1 && phys_Ey_ -> point(ii,jj) == phys_Ey_ -> point(ii+1,jj) )
                    ii ++;
                array<int,4> tempArr = { iistore,jj,ii-iistore+1,phys_Ey_->point(iistore,jj)};
                if(objArr_[tempArr[3]].mat().size() <= 1)
                    zaxEy_.push_back(tempArr);
                else
                    zaxDy_.push_back(tempArr);
                ii++;
            }
        }
        if(xPML_ != 0)
        {
            for(int jj = ymin; jj < ymax-y_sub; jj++)
            {
                array<int,4> tempArr = {xPML_+1, jj, static_cast<int>(nx_)-2*xPML_,0};
                zaxHz_.push_back(tempArr);
            }

        }
        else
        {
            for(int jj = ymin; jj < ymax-y_sub; jj++)
            {
                array<int,4> tempArr = {1, jj, static_cast<int>(nx_-1),0};
                zaxHz_.push_back(tempArr);
            }
        }
    }
    //initialize PML's with or without an opposition sigma.
    if(pmlArr_.size() > 1)
    {
        for(int kk = 0; kk < pmlArr_.size(); kk++)
        {
            UPML<complex<double>>* opp = &pmlArr_[abs(kk-1)];
            pmlArr_[kk].initializeUPML(objArr_, nx_,ny_,dx_,dy_,dt_, yPML_, xPML_, opp);
        }
    }
    else if (pmlArr_.size() ==1)
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
            outFile << setw(9) << tcur_ << "\t" << (d.loc()[xDTC_]-1)*dx_ << "\t" << (d.loc()[yDTC_]-1)*dy_ << "\t" << setw(10) << d.output(Ez_,eps).real() << "\t" << setw(10) << d.output(Ez_,eps).imag() << "\t" << setw(10) << srcArr_[0].prof().pulse(t_step_).real() << endl;
            cout << setw(9) << tcur_ << "\t" << (d.loc()[xDTC_]-1)*dx_ << "\t" << (d.loc()[yDTC_]-1)*dy_ << "\t" << setw(10) << d.output(Ez_,eps).real() << "\t" << setw(10) << d.output(Ez_,eps).imag() << "\t" << setw(10) << srcArr_[0].prof().pulse(t_step_).real() << endl;
            break;
        case HX:
            outFile << setw(9) << tcur_ << "\t" << (d.loc()[xDTC_]-1)*dx_ << "\t" << (d.loc()[yDTC_]-1)*dy_ << "\t" << setw(10) << d.output(Hx_,eps).real() << "\t" << setw(10) << d.output(Hx_,eps).imag() << endl;
            cout << setw(9) << tcur_ << "\t" << (d.loc()[xDTC_]-1)*dx_ << "\t" << (d.loc()[yDTC_]-1)*dy_ << "\t" << setw(10) << d.output(Hx_,eps).real() << "\t" << setw(10) << d.output(Hx_,eps).imag() << endl;
            break;
        case HY:
            outFile << setw(9) << tcur_ << "\t" << (d.loc()[xDTC_]-1)*dx_ << "\t" << (d.loc()[yDTC_]-1)*dy_ << "\t" << setw(10) << d.output(Hy_,eps).real() << "\t" << setw(10) << d.output(Hy_,eps).imag() << endl;
            cout << setw(9) << tcur_ << "\t" << (d.loc()[xDTC_]-1)*dx_ << "\t" << (d.loc()[yDTC_]-1)*dy_ << "\t" << setw(10) << d.output(Hy_,eps).real() << "\t" << setw(10) << d.output(Hy_,eps).imag() << endl;
            break;
        case HZ:
            outFile << setw(9) << tcur_ << "\t" << (d.loc()[xDTC_]-1)*dx_ << "\t" << (d.loc()[yDTC_]-1)*dy_ << "\t" << setw(10) << d.output(Hz_,eps).real() << "\t" << setw(10) << d.output(Hz_,eps).imag() << endl;
            cout << setw(9) << tcur_ << "\t" << (d.loc()[xDTC_]-1)*dx_ << "\t" << (d.loc()[yDTC_]-1)*dy_ << "\t" << setw(10) << d.output(Hz_,eps).real() << "\t" << setw(10) << d.output(Hz_,eps).imag() << endl;
            break;
        case EX:
            outFile << setw(9) << tcur_ << "\t" << (d.loc()[xDTC_]-1)*dx_ << "\t" << (d.loc()[yDTC_]-1)*dy_ << "\t" << setw(10) << d.output(Ex_,eps).real() << "\t" << setw(10) << d.output(Ex_,eps).imag() << "\t" << setw(10) << srcArr_[0].prof().pulse(t_step_).real() << endl;
            cout << setw(9) << tcur_ << "\t" << (d.loc()[xDTC_]-1)*dx_ << "\t" << (d.loc()[yDTC_]-1)*dy_ << "\t" << setw(10) << d.output(Ex_,eps).real() << "\t" << setw(10) << d.output(Ex_,eps).imag() << "\t" << setw(10) << srcArr_[0].prof().pulse(t_step_).real() << endl;
            break;
        case EY:
            outFile << setw(9) << tcur_ << "\t" << (d.loc()[xDTC_]-1)*dx_ << "\t" << (d.loc()[yDTC_]-1)*dy_ << "\t" << setw(10) << d.output(Ey_,eps).real() << "\t" << setw(10) << d.output(Ey_,eps).imag() << "\t" << setw(10) << srcArr_[0].prof().pulse(t_step_).real() << endl;
            cout << setw(9) << tcur_ << "\t" << (d.loc()[xDTC_]-1)*dx_ << "\t" << (d.loc()[yDTC_]-1)*dy_ << "\t" << setw(10) << d.output(Ey_,eps).real() << "\t" << setw(10) << d.output(Ey_,eps).imag() << "\t" << setw(10) << srcArr_[0].prof().pulse(t_step_).real() << endl;
            break;
        default:
            throw logic_error("reached a default case in a switch state that should never happen!");
            break;
    }
    outFile.close();
}

/**
 * @brief Periodic Boundry correction
 * @details Determines the prefactor to correct for lattice vector distances for PBC
 *
 * @param r Point on the lattice to find the prefactor in
 * @return Prefactor correnction
 */
complex<double> FDTDField::per_factor(std::vector<double> r)
{
    complex<double> i1(0.0,-1.0);
    double ip = inner_product(r.begin(),r.end(), k_point_.begin(),0);
    return exp(i1*ip);
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
        // c_hxh is set in case of Magnetic materials get used
        double c_hxh = 1.0;
        double c_hxe = 1.0 * dt_/dx_;
        double c_hyh = 1.0;
        double c_hye = 1.0 * dt_/dy_;
        for(int kk = 0; kk < zaxHx_.size(); kk++)
        {
            int xx = zaxHx_[kk][0]; int yy = zaxHx_[kk][1]; int nZax = zaxHx_[kk][2];
            zscal_(nZax,      c_hxh, &Hx_->point(xx  ,yy  ), 1);
            zaxpy_(nZax, -1.0*c_hxe, &Ez_->point(xx  ,yy+1), 1, &Hx_->point(xx,yy), 1);
            zaxpy_(nZax,      c_hxe, &Ez_->point(xx  ,yy  ), 1, &Hx_->point(xx,yy), 1);

        }
        for(int kk = 0; kk < zaxHy_.size(); kk++)
        {
            int xx = zaxHy_[kk][0]; int yy = zaxHy_[kk][1]; int nZax = zaxHy_[kk][2];
            zscal_(nZax,      c_hyh, &Hy_->point(xx  ,yy  ), 1);
            zaxpy_(nZax,      c_hye, &Ez_->point(xx+1,yy  ), 1, &Hy_->point(xx,yy), 1);
            zaxpy_(nZax, -1.0*c_hye, &Ez_->point(xx  ,yy  ), 1, &Hy_->point(xx,yy), 1);
        }
        // PML
        for(int kk =0; kk < pmlArr_.size(); kk++)
        {
            // Set parameters for switches and general axpy stuff
            // stride is because XPML are constant in the direction opposite of continuous memory
            // d converts the direction of the PML into an int
            // ni/nj are to calculate the relative direction for the opposite end.
            // PML fields don't have the border which is why there is a xx-1 and yy-1 for all of them
            int stride = 1; int stride_rel = 1;
            int ni = 0; int nj = 0; int d = 0;
            if(pmlArr_[kk].d() == X)
            {
                stride_rel = nx_+2;
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
                zcopy_(nZax, &pmlArr_[kk].Bx_->point(xx-1,yy-1), stride, bxstore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].Bx_ -> point(xx-1,yy-1), stride);
                zscal_(nZax, zaxArr[6],             &Hx_ -> point(xx,yy), stride_rel);

                zaxpy_(nZax,      zaxArr[5], &Ez_ -> point(xx,yy  ), stride_rel, &pmlArr_[kk].Bx_->point(xx-1,yy-1), stride);
                zaxpy_(nZax, -1.0*zaxArr[5], &Ez_ -> point(xx,yy+1), stride_rel, &pmlArr_[kk].Bx_->point(xx-1,yy-1), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].Bx_->point(xx-1,yy-1), stride, &Hx_ -> point(xx,yy), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], bxstore.data()                    , 1     , &Hx_ -> point(xx,yy), stride_rel);
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxHx_end_.size();zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxHx_end_[zz];
                int xx = static_cast<int>(zaxArr[0]);
                int yy = static_cast<int>(zaxArr[1]);
                int xx_rel = ni + pow(-1, 1-d) * xx+2*(1-d);  /// If the pml is is the x direction xx_rel = nx - xx, else xrel = xx
                int yy_rel = nj + pow(-1, d)   * yy+2*d; /// If the pml is is the y direction yy_rel = ny - yy, else yrel = yy
                int nZax = static_cast<int>(zaxArr[2]);

                vector<complex<double>> bxstore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].Bx_end_->point(xx-1,yy-1), stride, bxstore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].Bx_end_ -> point(xx-1,yy-1), stride);
                zscal_(nZax, zaxArr[6],                 &Hx_ -> point(xx_rel,yy_rel), stride_rel);

                zaxpy_(nZax,      zaxArr[5], &Ez_ -> point(xx_rel,yy_rel  ), stride_rel, &pmlArr_[kk].Bx_end_->point(xx-1,yy-1), stride);
                zaxpy_(nZax, -1.0*zaxArr[5], &Ez_ -> point(xx_rel,yy_rel+1), stride_rel, &pmlArr_[kk].Bx_end_->point(xx-1,yy-1), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].Bx_end_->point(xx-1,yy-1), stride, &Hx_ -> point(xx_rel,yy_rel), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], bxstore.data()                      , 1     , &Hx_ -> point(xx_rel,yy_rel), stride_rel);
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxHy_.size();zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxHy_[zz];
                int xx = static_cast<int>(zaxArr[0]);
                int yy = static_cast<int>(zaxArr[1]);
                int nZax = static_cast<int>(zaxArr[2]);

                vector<complex<double>> bystore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].By_->point(xx-1,yy-1), stride, bystore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].By_ -> point(xx-1,yy-1), stride);
                zscal_(nZax, zaxArr[6],                 &Hy_ -> point(xx,yy), stride_rel);

                zaxpy_(nZax,      zaxArr[5], &Ez_ -> point(xx+1,yy), stride_rel, &pmlArr_[kk].By_->point(xx-1,yy-1), stride);
                zaxpy_(nZax, -1.0*zaxArr[5], &Ez_ -> point(xx  ,yy), stride_rel, &pmlArr_[kk].By_->point(xx-1,yy-1), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].By_->point(xx-1,yy-1), stride, &Hy_ -> point(xx,yy), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], bystore.data()                      , 1     , &Hy_ -> point(xx,yy), stride_rel);
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxHy_end_.size();zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxHy_end_[zz];
                int xx = static_cast<int>(zaxArr[0]);
                int yy = static_cast<int>(zaxArr[1]);
                int xx_rel = ni + pow(-1, 1-d) * xx+2*(1-d);  /// If the pml is is the x direction xx_rel = nx - xx, else xrel = xx
                int yy_rel = nj + pow(-1, d)   * yy+2*d; /// If the pml is is the y direction yy_rel = ny - yy, else yrel = yy
                int nZax = static_cast<int>(zaxArr[2]);

                vector<complex<double>> bystore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].By_end_->point(xx-1,yy-1), stride, bystore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].By_end_ -> point(xx-1,yy-1), stride);
                zscal_(nZax, zaxArr[6],                 &Hy_ -> point(xx_rel,yy_rel), stride_rel);

                zaxpy_(nZax,      zaxArr[5], &Ez_ -> point(xx_rel+1,yy_rel), stride_rel, &pmlArr_[kk].By_end_->point(xx-1,yy-1), stride);
                zaxpy_(nZax, -1.0*zaxArr[5], &Ez_ -> point(xx_rel  ,yy_rel), stride_rel, &pmlArr_[kk].By_end_->point(xx-1,yy-1), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].By_end_->point(xx-1,yy-1), stride, &Hy_ -> point(xx_rel,yy_rel), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], bystore.data()                      , 1     , &Hy_ -> point(xx_rel,yy_rel), stride_rel);
            }
        }
        //corners have to be done piecewise
        complex<double> bxstore(0.0,0.0); complex<double> bystore(0.0,0.0);
        if(pmlArr_.size() > 1)
        {
            // Setting everything such that the X-PML stores all the information for the corners
            int kx = 1;
            // Since PML corners are identical in either PML only store the information in one.
            shared_ptr<vector<vector<array<double,5>>>> c_hx_0_n;
            shared_ptr<vector<vector<array<double,5>>>> c_hx_n_0;
            shared_ptr<vector<vector<array<double,5>>>> c_hy_0_n;
            shared_ptr<vector<vector<array<double,5>>>> c_hy_n_0;
            // 0_0 and n_n corners will be the same
            shared_ptr<vector<vector<array<double,5>>>> c_hx_0_0 = pmlArr_[1].c_hx_0_0_;
            shared_ptr<vector<vector<array<double,5>>>> c_hx_n_n = pmlArr_[1].c_hx_n_n_;
            shared_ptr<vector<vector<array<double,5>>>> c_hy_0_0 = pmlArr_[1].c_hy_0_0_;
            shared_ptr<vector<vector<array<double,5>>>> c_hy_n_n = pmlArr_[1].c_hy_n_n_;
            // Ensures corners are correct in all cases 0_n and n_0 would be reversed depending on which direction the storing PML is in
            if(pmlArr_[1].d() == X)
            {
                c_hx_0_n = pmlArr_[1].c_hx_0_n_;
                c_hx_n_0 = pmlArr_[1].c_hx_n_0_;
                c_hy_0_n = pmlArr_[1].c_hy_0_n_;
                c_hy_n_0 = pmlArr_[1].c_hy_n_0_;
            }
            else
            {
                kx = 0;
                c_hx_0_n = pmlArr_[1].c_hx_n_0_;
                c_hx_n_0 = pmlArr_[1].c_hx_0_n_;
                c_hy_0_n = pmlArr_[1].c_hy_n_0_;
                c_hy_n_0 = pmlArr_[1].c_hy_0_n_;
            }
            complex<double> bxstore(0.0,0.0); complex<double> bystore(0.0,0.0);
            // Storage is done in the x direction PML why it's ii,yy-1
            for(int ii = 0; ii < xPML_; ii++)
            {
                for(int jj = 0; jj < yPML_; jj++)
                {
                    //Bot Left
                    int xx = ii+1; int yy = jj+1;
                    bxstore = pmlArr_[kx].Bx_->point(ii,yy-1);
                    bystore = pmlArr_[kx].By_->point(ii,yy-1);
                    pmlArr_[kx].Bx_->point(ii,yy-1) = c_hx_0_0->at(ii).at(jj)[0] * pmlArr_[kx].Bx_->point(ii,yy-1) - c_hx_0_0->at(ii).at(jj)[1] * (Ez_->point(xx,yy+1)-Ez_->point(xx,yy));
                    pmlArr_[kx].By_->point(ii,yy-1) = c_hy_0_0->at(ii).at(jj)[0] * pmlArr_[kx].By_->point(ii,yy-1) + c_hy_0_0->at(ii).at(jj)[1] * (Ez_->point(xx+1,yy)-Ez_->point(xx,yy));
                    Hx_->point(xx,yy) = c_hx_0_0->at(ii).at(jj)[2] * Hx_->point(xx,yy) + c_hx_0_0->at(ii).at(jj)[3] * pmlArr_[kx].Bx_->point(ii,yy-1) - c_hx_0_0->at(ii).at(jj)[4] * bxstore;
                    Hy_->point(xx,yy) = c_hy_0_0->at(ii).at(jj)[2] * Hy_->point(xx,yy) + c_hy_0_0->at(ii).at(jj)[3] * pmlArr_[kx].By_->point(ii,yy-1) - c_hy_0_0->at(ii).at(jj)[4] * bystore;

                    //Top Left
                    xx = ii+1; yy = ny_-jj;
                    bxstore = pmlArr_[kx].Bx_->point(ii,yy-1);
                    bystore = pmlArr_[kx].By_->point(ii,yy-1);
                    pmlArr_[kx].Bx_->point(ii,yy-1) = c_hx_0_n->at(ii).at(jj)[0] * pmlArr_[kx].Bx_->point(ii,yy-1) - c_hx_0_n->at(ii).at(jj)[1] * (Ez_->point(xx,yy+1)-Ez_->point(xx,yy));
                    pmlArr_[kx].By_->point(ii,yy-1) = c_hy_0_n->at(ii).at(jj)[0] * pmlArr_[kx].By_->point(ii,yy-1) + c_hy_0_n->at(ii).at(jj)[1] * (Ez_->point(xx+1,yy)-Ez_->point(xx,yy));
                    Hx_->point(xx,yy) = c_hx_0_n->at(ii).at(jj)[2] * Hx_->point(xx,yy) + c_hx_0_n->at(ii).at(jj)[3] * pmlArr_[kx].Bx_->point(ii,yy-1) - c_hx_0_n->at(ii).at(jj)[4] * bxstore;
                    Hy_->point(xx,yy) = c_hy_0_n->at(ii).at(jj)[2] * Hy_->point(xx,yy) + c_hy_0_n->at(ii).at(jj)[3] * pmlArr_[kx].By_->point(ii,yy-1) - c_hy_0_n->at(ii).at(jj)[4] * bystore;

                    //Top Right
                    xx = nx_-ii; yy = ny_-jj;
                    bxstore = pmlArr_[kx].Bx_end_->point(ii,yy-1);
                    bystore = pmlArr_[kx].By_end_->point(ii,yy-1);
                    pmlArr_[kx].Bx_end_->point(ii,yy-1) = c_hx_n_n->at(ii).at(jj)[0] * pmlArr_[kx].Bx_end_->point(ii,yy-1) - c_hx_n_n->at(ii).at(jj)[1] * (Ez_->point(xx,yy+1)-Ez_->point(xx,yy));
                    pmlArr_[kx].By_end_->point(ii,yy-1) = c_hy_n_n->at(ii).at(jj)[0] * pmlArr_[kx].By_end_->point(ii,yy-1) + c_hy_n_n->at(ii).at(jj)[1] * (Ez_->point(xx+1,yy)-Ez_->point(xx,yy));
                    Hx_->point(xx,yy) = c_hx_n_n->at(ii).at(jj)[2] * Hx_->point(xx,yy) + c_hx_n_n->at(ii).at(jj)[3] * pmlArr_[kx].Bx_end_->point(ii,yy-1) - c_hx_n_n->at(ii).at(jj)[4] * bxstore;
                    Hy_->point(xx,yy) = c_hy_n_n->at(ii).at(jj)[2] * Hy_->point(xx,yy) + c_hy_n_n->at(ii).at(jj)[3] * pmlArr_[kx].By_end_->point(ii,yy-1) - c_hy_n_n->at(ii).at(jj)[4] * bystore;

                    //Bot Right
                    xx = nx_-ii; yy = jj+1;
                    bxstore = pmlArr_[kx].Bx_end_->point(ii,yy-1);
                    bystore = pmlArr_[kx].By_end_->point(ii,yy-1);
                    pmlArr_[kx].Bx_end_->point(ii,yy-1) = c_hx_n_0->at(ii).at(jj)[0] * pmlArr_[kx].Bx_end_->point(ii,yy-1) - c_hx_n_0->at(ii).at(jj)[1] * (Ez_->point(xx,yy+1)-Ez_->point(xx,yy));
                    pmlArr_[kx].By_end_->point(ii,yy-1) = c_hy_n_0->at(ii).at(jj)[0] * pmlArr_[kx].By_end_->point(ii,yy-1) + c_hy_n_0->at(ii).at(jj)[1] * (Ez_->point(xx+1,yy)-Ez_->point(xx,yy));
                    Hx_->point(xx,yy) = c_hx_n_0->at(ii).at(jj)[2] * Hx_->point(xx,yy) + c_hx_n_0->at(ii).at(jj)[3] * pmlArr_[kx].Bx_end_->point(ii,yy-1) - c_hx_n_0->at(ii).at(jj)[4] * bxstore;
                    Hy_->point(xx,yy) = c_hy_n_0->at(ii).at(jj)[2] * Hy_->point(xx,yy) + c_hy_n_0->at(ii).at(jj)[3] * pmlArr_[kx].By_end_->point(ii,yy-1) - c_hy_n_0->at(ii).at(jj)[4] * bystore;
                }
            }
        }
    }
    else
    {
        // c_hxh is set in case of Magnetic materials get used
        double c_hzh = 1.0;
        double c_hze = 1.0 * dt_/dx_;
        for(int kk = 0; kk < zaxHz_.size(); kk++)
        {
            int xx = zaxHz_[kk][0]; int yy = zaxHz_[kk][1]; int nZax = zaxHz_[kk][2];
            zscal_(nZax,       c_hzh, &Hz_->point(xx  ,yy  ), 1);
            zaxpy_(nZax,  -1.0*c_hze, &Ex_->point(xx  ,yy  ), 1, &Hz_ ->point(xx,yy), 1);
            zaxpy_(nZax,       c_hze, &Ex_->point(xx  ,yy+1), 1, &Hz_ ->point(xx,yy), 1);
            zaxpy_(nZax,       c_hze, &Ey_->point(xx  ,yy  ), 1, &Hz_ ->point(xx,yy), 1);
            zaxpy_(nZax,  -1.0*c_hze, &Ey_->point(xx+1,yy  ), 1, &Hz_ ->point(xx,yy), 1);
            // array<complex<double>,5> upConsts = {c_hzh, -1.0*c_hze, c_hze, c_hze, -1.0*c_hze};
            // array<int, 4> axConsts = {zaxHz_[kk][0], zaxHz_[kk][1],zaxHz_[kk][2],1};
            // zFieldUpdate(Hz_, Ex_, Ey_, axConsts, upConsts);
        }
        for(int kk = 0; kk < pmlArr_.size(); kk++)
        {
            // Set parameters for switches and general axpy stuff
            // stride is because XPML are constant in the direction opposite of continuous memory
            // d converts the direction of the PML into an int
            // ni/nj are to calculate the relative direction for the opposite end.
            // PML fields don't have the border which is why there is a xx-1 and yy-1 for all of themint stride = 1; int stride_rel = 1;
            int stride = 1; int stride_rel = 1;
            int ni = 0; int nj = 0; int d = 0;
            if(pmlArr_[kk].d() == X)
            {
                stride_rel = nx_+2;
                ni         = nx_+1;
                stride     = pmlArr_[kk].thickness();
            }
            else
            {
                d          = 1;
                nj         = ny_+1;
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxHz_.size(); zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxHz_[zz];
                int xx   = static_cast<int>(zaxArr[0]);
                int yy   = static_cast<int>(zaxArr[1]);
                int nZax = static_cast<int>(zaxArr[2]);
                vector<complex<double>> bzstore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].Bz_->point(xx-1,yy-1), stride, bzstore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].Bz_->point(xx-1,yy-1), stride);
                zscal_(nZax, zaxArr[6],             &Hz_ -> point(xx,yy), stride_rel);

                zaxpy_(nZax,     zaxArr[5], &Ex_->point(xx  ,yy+1), stride_rel, &pmlArr_[kk].Bz_->point(xx-1,yy-1), stride);
                zaxpy_(nZax,-1.0*zaxArr[5], &Ex_->point(xx  ,yy  ), stride_rel, &pmlArr_[kk].Bz_->point(xx-1,yy-1), stride);
                zaxpy_(nZax,     zaxArr[5], &Ey_->point(xx  ,yy  ), stride_rel, &pmlArr_[kk].Bz_->point(xx-1,yy-1), stride);
                zaxpy_(nZax,-1.0*zaxArr[5], &Ey_->point(xx+1,yy  ), stride_rel, &pmlArr_[kk].Bz_->point(xx-1,yy-1), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].Bz_->point(xx-1,yy-1), stride, &Hz_ -> point(xx,yy), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], bzstore.data()                    , 1     , &Hz_ -> point(xx,yy), stride_rel);
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxHz_end_.size(); zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxHz_end_[zz];
                int xx    = static_cast<int>(zaxArr[0]);
                int yy    = static_cast<int>(zaxArr[1]);
                int nZax  = static_cast<int>(zaxArr[2]);
                int xx_rel = ni + pow(-1, 1-d) * xx; // If the pml is is the x direction xx_rel = nx - xx, else xrel = xx
                int yy_rel = nj + pow(-1, d)   * yy; // If the pml is is the y direction yy_rel = ny - yy, else yrel = yy

                vector<complex<double>> bzstore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].Bz_end_->point(xx-1,yy-1), stride, bzstore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].Bz_end_->point(xx-1,yy-1), stride);
                zscal_(nZax, zaxArr[6],             &Hz_     -> point(xx_rel,yy_rel), stride_rel);

                zaxpy_(nZax,     zaxArr[5], &Ex_->point(xx_rel  ,yy_rel+1), stride_rel, &pmlArr_[kk].Bz_end_->point(xx-1,yy-1), stride);
                zaxpy_(nZax,-1.0*zaxArr[5], &Ex_->point(xx_rel  ,yy_rel  ), stride_rel, &pmlArr_[kk].Bz_end_->point(xx-1,yy-1), stride);
                zaxpy_(nZax,     zaxArr[5], &Ey_->point(xx_rel  ,yy_rel  ), stride_rel, &pmlArr_[kk].Bz_end_->point(xx-1,yy-1), stride);
                zaxpy_(nZax,-1.0*zaxArr[5], &Ey_->point(xx_rel+1,yy_rel  ), stride_rel, &pmlArr_[kk].Bz_end_->point(xx-1,yy-1), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].Bz_end_->point(xx-1,yy-1), stride, &Hz_ -> point(xx_rel,yy_rel), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], bzstore.data()                      , 1     , &Hz_ -> point(xx_rel,yy_rel), stride_rel);
            }
        }
        if(pmlArr_.size() > 1)
        {
            int kx = 1;
            // Setting everything such that the X-PML stores all the information for the corners
            // Since PML corners are identical in either PML only store the information in one.
            shared_ptr<vector<vector<array<double,5>>>> c_hz_0_n;
            shared_ptr<vector<vector<array<double,5>>>> c_hz_n_0;
            // 0_0 and n_n corners will be the same
            shared_ptr<vector<vector<array<double,5>>>> c_hz_0_0 = pmlArr_[1].c_hz_0_0_;
            shared_ptr<vector<vector<array<double,5>>>> c_hz_n_n = pmlArr_[1].c_hz_n_n_;
            // Ensures corners are correct in all cases 0_n and n_0 would be reversed depending on which direction the storing PML is in
            if(pmlArr_[1].d() == X)
            {
                c_hz_0_n = pmlArr_[1].c_hz_0_n_;
                c_hz_n_0 = pmlArr_[1].c_hz_n_0_;
            }
            else
            {
                kx = 0;
                c_hz_0_n = pmlArr_[1].c_hz_n_0_;
                c_hz_n_0 = pmlArr_[1].c_hz_0_n_;
            }
            complex<double> bzstore(0.0,0.0);
            int xx = 1; int yy = 1;
            for(int ii = 0; ii< xPML_; ii++)
            {
                for(int jj = 0; jj < yPML_; jj++)
                {
                    // Storage is done in the x direction PML why it's ii,yy-1
                    //Bot Left
                    xx = ii+1; yy = jj+1;
                    bzstore = pmlArr_[kx].Bz_->point(ii,yy-1);
                    pmlArr_[kx].Bz_->point(ii,yy-1) = c_hz_0_0->at(ii).at(jj)[0] * pmlArr_[kx].Bz_->point(ii,yy-1) - c_hz_0_0->at(ii).at(jj)[1] * ((Ey_->point(xx+1,yy)-Ey_->point(xx,yy)) - (Ex_->point(xx,yy+1)-Ex_->point(xx,yy)));
                    Hz_->point(xx,yy) = c_hz_0_0->at(ii).at(jj)[2] * Hz_->point(xx,yy) + c_hz_0_0->at(ii).at(jj)[3] * pmlArr_[kx].Bz_->point(ii,yy-1) - c_hz_0_0->at(ii).at(jj)[4] * bzstore;

                    //Top Left
                    xx = ii+1; yy = ny_ - jj;
                    bzstore = pmlArr_[kx].Bz_->point(ii,yy-1);
                    pmlArr_[kx].Bz_->point(ii,yy-1) = c_hz_0_n->at(ii).at(jj)[0] * pmlArr_[kx].Bz_->point(ii,yy-1) - c_hz_0_n->at(ii).at(jj)[1] * ((Ey_->point(xx+1,yy)-Ey_->point(xx,yy)) - (Ex_->point(xx,yy+1)-Ex_->point(xx,yy)));
                    Hz_->point(xx,yy) = c_hz_0_n->at(ii).at(jj)[2] * Hz_->point(xx,yy) + c_hz_0_n->at(ii).at(jj)[3] * pmlArr_[kx].Bz_->point(ii,yy-1) - c_hz_0_n->at(ii).at(jj)[4] * bzstore;

                    //Top Right
                    xx = nx_ - ii; yy = ny_ - jj;
                    bzstore = pmlArr_[kx].Bz_end_->point(ii,yy-1);
                    pmlArr_[kx].Bz_end_->point(ii,yy-1) = c_hz_n_n->at(ii).at(jj)[0] * pmlArr_[kx].Bz_end_->point(ii,yy-1) - c_hz_n_n->at(ii).at(jj)[1] * ((Ey_->point(xx+1,yy)-Ey_->point(xx,yy)) - (Ex_->point(xx,yy+1)-Ex_->point(xx,yy)));
                    Hz_->point(xx,yy) = c_hz_n_n->at(ii).at(jj)[2] * Hz_->point(xx,yy) + c_hz_n_n->at(ii).at(jj)[3] * pmlArr_[kx].Bz_end_->point(ii,yy-1) - c_hz_n_n->at(ii).at(jj)[4] * bzstore;

                    //Bot Right
                    xx = nx_ - ii; yy = jj+1;
                    bzstore = pmlArr_[kx].Bz_end_->point(ii,yy-1);
                    pmlArr_[kx].Bz_end_->point(ii,yy-1) = c_hz_n_0->at(ii).at(jj)[0] * pmlArr_[kx].Bz_end_->point(ii,yy-1) - c_hz_n_0->at(ii).at(jj)[1] * ((Ey_->point(xx+1,yy)-Ey_->point(xx,yy)) - (Ex_->point(xx,yy+1)-Ex_->point(xx,yy)));
                    Hz_->point(xx,yy) = c_hz_n_0->at(ii).at(jj)[2] * Hz_->point(xx,yy) + c_hz_n_0->at(ii).at(jj)[3] * pmlArr_[kx].Bz_end_->point(ii,yy-1) - c_hz_n_0->at(ii).at(jj)[4] * bzstore;
                }
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
    if(Ez_)
    {
        double eps =1.0;
        double c_eze =1.0;
        double c_ezh = dt_/(eps*dx_);
        // Seperated out by PML because edge cases require special treatment
        for(int kk = 0; kk < zaxEz_.size(); kk++)
        {
            eps = objArr_[zaxEz_[kk][3]].dielectric();
            c_ezh = dt_/(eps*dx_);
            int xx = zaxEz_[kk][0]; int yy = zaxEz_[kk][1]; int nZax = zaxEz_[kk][2];
            zscal_(nZax,      c_eze, &Ez_->point(xx  ,yy  ), 1);
            zaxpy_(nZax,      c_ezh, &Hx_->point(xx  ,yy-1), 1, &Ez_ ->point(xx,yy), 1);
            zaxpy_(nZax, -1.0*c_ezh, &Hx_->point(xx  ,yy  ), 1, &Ez_ ->point(xx,yy), 1);
            zaxpy_(nZax, -1.0*c_ezh, &Hy_->point(xx-1,yy  ), 1, &Ez_ ->point(xx,yy), 1);
            zaxpy_(nZax,      c_ezh, &Hy_->point(xx  ,yy  ), 1, &Ez_ ->point(xx,yy), 1);

        }
         //PML
        for(int kk = 0; kk < pmlArr_.size(); kk++)
        {
            int stride = 1; int stride_rel = 1;
            int ni = 0; int nj = 0; int d = 0;
            if(pmlArr_[kk].d() == X)
            {
                stride_rel = nx_+2;
                ni         = nx_ - 1;
                stride     = pmlArr_[kk].thickness();
            }
            else
            {
                d          = 1;
                nj         = ny_ - 1;
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxEz_.size(); zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxEz_[zz];
                int xx    = static_cast<int>(zaxArr[0]);
                int yy    = static_cast<int>(zaxArr[1]);
                int nZax  = static_cast<int>(zaxArr[2]);

                vector<complex<double>> dzstore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].Dz_ -> point(xx-1,yy-1), stride, dzstore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].Dz_ -> point(xx-1,yy-1), stride);
                zscal_(nZax, zaxArr[6],             &Ez_     -> point(xx,yy), stride_rel);

                zaxpy_(nZax,-1.0*zaxArr[5], &Hx_->point(xx  ,yy  ), stride_rel, &pmlArr_[kk].Dz_->point(xx-1,yy-1), stride);
                zaxpy_(nZax,     zaxArr[5], &Hx_->point(xx  ,yy-1), stride_rel, &pmlArr_[kk].Dz_->point(xx-1,yy-1), stride);
                zaxpy_(nZax,-1.0*zaxArr[5], &Hy_->point(xx-1,yy  ), stride_rel, &pmlArr_[kk].Dz_->point(xx-1,yy-1), stride);
                zaxpy_(nZax,     zaxArr[5], &Hy_->point(xx  ,yy  ), stride_rel, &pmlArr_[kk].Dz_->point(xx-1,yy-1), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].Dz_ -> point(xx-1,yy-1), stride, &Ez_ -> point(xx,yy), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], dzstore.data()                      , 1     , &Ez_ -> point(xx,yy), stride_rel);
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxEz_end_.size(); zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxEz_end_[zz];
                int xx    = static_cast<int>(zaxArr[0]);
                int yy    = static_cast<int>(zaxArr[1]);
                int nZax  = static_cast<int>(zaxArr[2]);
                int xx_rel = ni + pow(-1, 1-d) * xx+2*(1-d);  /// If the pml is is the x direction xx_rel = nx - xx, else xrel = xx
                int yy_rel = nj + pow(-1, d)   * yy+2*d; /// If the pml is is the y direction yy_rel = ny - yy, else yrel = yy

                vector<complex<double>> dzstore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].Dz_end_ -> point(xx-1,yy-1), stride, dzstore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].Dz_end_ -> point(xx-1,yy-1), stride);
                zscal_(nZax, zaxArr[6],             &Ez_     -> point(xx_rel,yy_rel), stride_rel);

                zaxpy_(nZax,-1.0*zaxArr[5], &Hx_->point(xx_rel  ,yy_rel  ), stride_rel, &pmlArr_[kk].Dz_end_->point(xx-1,yy-1), stride);
                zaxpy_(nZax,     zaxArr[5], &Hx_->point(xx_rel  ,yy_rel-1), stride_rel, &pmlArr_[kk].Dz_end_->point(xx-1,yy-1), stride);
                zaxpy_(nZax,-1.0*zaxArr[5], &Hy_->point(xx_rel-1,yy_rel  ), stride_rel, &pmlArr_[kk].Dz_end_->point(xx-1,yy-1), stride);
                zaxpy_(nZax,     zaxArr[5], &Hy_->point(xx_rel  ,yy_rel  ), stride_rel, &pmlArr_[kk].Dz_end_->point(xx-1,yy-1), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].Dz_end_ -> point(xx-1,yy-1), stride, &Ez_ -> point(xx_rel,yy_rel), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], dzstore.data()                      , 1     , &Ez_ -> point(xx_rel,yy_rel), stride_rel);
            }
        }
        if(pmlArr_.size() > 1)
        {
            int kx = 1;
            shared_ptr<vector<vector<array<double,5>>>> c_ez_0_n;
            shared_ptr<vector<vector<array<double,5>>>> c_ez_n_0;

            shared_ptr<vector<vector<array<double,5>>>> c_ez_0_0 = pmlArr_[1].c_ez_0_0_;
            shared_ptr<vector<vector<array<double,5>>>> c_ez_n_n = pmlArr_[1].c_ez_n_n_;
            if(pmlArr_[1].d() == X)
            {
                c_ez_0_n = pmlArr_[1].c_ez_0_n_;
                c_ez_n_0 = pmlArr_[1].c_ez_n_0_;
            }
            else
            {
                kx = 0;
                c_ez_0_n = pmlArr_[1].c_ez_n_0_;
                c_ez_n_0 = pmlArr_[1].c_ez_0_n_;
            }
            complex<double> dzstore(0.0,0.0);
            int xx = 0; int yy = 0;
            for(int ii = 0; ii< xPML_; ii++)
            {
                for(int jj = 0; jj < yPML_; jj++)
                {
                    //Bot Left
                    xx = ii+1; yy = jj+1;
                    dzstore = pmlArr_[kx].Dz_->point(ii,yy-1);
                    pmlArr_[kx].Dz_->point(ii,yy-1) = c_ez_0_0->at(ii).at(jj)[0] * pmlArr_[kx].Dz_->point(ii,yy-1) + c_ez_0_0->at(ii).at(jj)[1] * ((Hy_->point(xx,yy)-Hy_->point(xx-1,yy)) - (Hx_->point(xx,yy)-Hx_->point(xx,yy-1)));
                    Ez_->point(xx,yy) = c_ez_0_0->at(ii).at(jj)[2] * Ez_->point(xx,yy) + c_ez_0_0->at(ii).at(jj)[3] * pmlArr_[kx].Dz_->point(ii,yy-1) - c_ez_0_0->at(ii).at(jj)[4] * dzstore;

                    //Top Left
                    xx = ii+1; yy = ny_- jj;
                    dzstore = pmlArr_[kx].Dz_->point(ii,yy-1);
                    pmlArr_[kx].Dz_->point(ii,yy-1) = c_ez_0_n->at(ii).at(jj)[0] * pmlArr_[kx].Dz_->point(ii,yy-1) + c_ez_0_n->at(ii).at(jj)[1] * ((Hy_->point(xx,yy)-Hy_->point(xx-1,yy)) - (Hx_->point(xx,yy)-Hx_->point(xx,yy-1)));
                    Ez_->point(xx,yy) = c_ez_0_n->at(ii).at(jj)[2] * Ez_->point(xx,yy) + c_ez_0_n->at(ii).at(jj)[3] * pmlArr_[kx].Dz_->point(ii,yy-1) - c_ez_0_n->at(ii).at(jj)[4] * dzstore;

                    //Top Right
                    xx = nx_ - ii; yy = ny_- jj;
                    dzstore = pmlArr_[kx].Dz_end_->point(ii,yy-1);
                    pmlArr_[kx].Dz_end_->point(ii,yy-1) = c_ez_n_n->at(ii).at(jj)[0] * pmlArr_[kx].Dz_end_->point(ii,yy-1) + c_ez_n_n->at(ii).at(jj)[1] * ((Hy_->point(xx,yy)-Hy_->point(xx-1,yy)) - (Hx_->point(xx,yy)-Hx_->point(xx,yy-1)));
                    Ez_->point(xx,yy) = c_ez_n_n->at(ii).at(jj)[2] * Ez_->point(xx,yy) + c_ez_n_n->at(ii).at(jj)[3] * pmlArr_[kx].Dz_end_->point(ii,yy-1) - c_ez_n_n->at(ii).at(jj)[4] * dzstore;

                    //Bot Right
                    xx = nx_- ii; yy = jj+1;
                    dzstore = pmlArr_[kx].Dz_end_->point(ii,yy-1);
                    pmlArr_[kx].Dz_end_->point(ii,yy-1) = c_ez_n_0->at(ii).at(jj)[0] * pmlArr_[kx].Dz_end_->point(ii,yy-1) + c_ez_n_0->at(ii).at(jj)[1] * ((Hy_->point(xx,yy)-Hy_->point(xx-1,yy)) - (Hx_->point(xx,yy)-Hx_->point(xx,yy-1)));
                    Ez_->point(xx,yy) = c_ez_n_0->at(ii).at(jj)[2] * Ez_->point(xx,yy) + c_ez_n_0->at(ii).at(jj)[3] * pmlArr_[kx].Dz_end_->point(ii,yy-1) - c_ez_n_0->at(ii).at(jj)[4] * dzstore;
                }
            }
        }
    }
    else
    {
        double c_exe = 1.0;
        double c_eye = 1.0;
        double c_exh = 0.0;
        double c_eyh = 0.0;
        double eps = 1.0;
        for(int kk = 0; kk < zaxEx_.size(); kk++)
        {
            eps = objArr_[zaxEx_[kk][3]].dielectric();
            c_exh = dt_/(eps*dy_);
            int xx = zaxEx_[kk][0]; int yy = zaxEx_[kk][1]; int nZax = zaxEx_[kk][2];
            zscal_(nZax,      c_exe, &Ex_->point(xx  ,yy  ), 1);
            zaxpy_(nZax,      c_exh, &Hz_->point(xx  ,yy  ), 1, &Ex_->point(xx,yy), 1);
            zaxpy_(nZax, -1.0*c_exh, &Hz_->point(xx  ,yy-1), 1, &Ex_->point(xx,yy), 1);
            if(objArr_[zaxEx_[kk][3]].mat().size() > 1)
                cout << "Visit here for the update eqs" << endl;
        }
        for(int kk = 0; kk < zaxEy_.size(); kk++)
        {
            eps = objArr_[zaxEy_[kk][3]].dielectric();
            c_eyh = dt_/(eps*dy_);
            int xx = zaxEy_[kk][0]; int yy = zaxEy_[kk][1]; int nZax = zaxEy_[kk][2];
            zscal_(nZax,      c_eye, &Ey_->point(xx  ,yy  ), 1);
            zaxpy_(nZax, -1.0*c_eyh, &Hz_->point(xx  ,yy  ), 1, &Ey_->point(xx,yy), 1);
            zaxpy_(nZax,      c_eyh, &Hz_->point(xx-1,yy  ), 1, &Ey_->point(xx,yy), 1);
            if(objArr_[zaxEy_[kk][3]].mat().size() > 1)
                cout << "Visit here for the update eqs" << endl;
        }
        //PML
        for(int kk =0; kk < pmlArr_.size(); kk++)
        {
            int stride = 1; int stride_rel = 1;
            int ni = 0; int nj = 0; int d = 0;
            if(pmlArr_[kk].d() == X)
            {
                stride_rel = nx_+2;
                ni         = nx_ - 1;
                stride     = pmlArr_[kk].thickness();
            }
            else
            {
                d          = 1;
                nj         = ny_ - 1;
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxEx_.size();zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxEx_[zz];
                int xx = static_cast<int>(zaxArr[0]);
                int yy = static_cast<int>(zaxArr[1]);
                int nZax = static_cast<int>(zaxArr[2]);
                vector<complex<double>> dxstore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].Dx_ -> point(xx-1,yy-1), stride, dxstore.data(), 1);

                zscal_(nZax, zaxArr[4],     &pmlArr_[kk].Dx_ -> point(xx-1,yy-1), stride);
                zscal_(nZax, zaxArr[6],                 &Ex_ -> point(xx,yy)    , stride_rel);

                zaxpy_(nZax, -1.0*zaxArr[5], &Hz_ -> point(xx,yy-1), stride_rel, &pmlArr_[kk].Dx_ -> point(xx-1,yy-1), stride);
                zaxpy_(nZax,      zaxArr[5], &Hz_ -> point(xx,yy  ), stride_rel, &pmlArr_[kk].Dx_ -> point(xx-1,yy-1), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].Dx_ -> point(xx-1,yy-1), stride, &Ex_ -> point(xx,yy), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], dxstore.data()                      , 1     , &Ex_ -> point(xx,yy), stride_rel);
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxEx_end_.size();zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxEx_end_[zz];
                int xx = static_cast<int>(zaxArr[0]);
                int yy = static_cast<int>(zaxArr[1]);
                int xx_rel = ni + pow(-1, 1-d) * xx+2*(1-d);  /// If the pml is is the x direction xx_rel = nx - xx, else xrel = xx
                int yy_rel = nj + pow(-1, d)   * yy+2*d; /// If the pml is is the y direction yy_rel = ny - yy, else yrel = yy
                int nZax = static_cast<int>(zaxArr[2]);

                vector<complex<double>> dxstore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].Dx_end_ -> point(xx-1,yy-1), stride, dxstore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].Dx_end_ -> point(xx-1,yy-1), stride);
                zscal_(nZax, zaxArr[6],                 &Ex_ -> point(xx_rel,yy_rel), stride_rel);

                zaxpy_(nZax, -1.0*zaxArr[5], &Hz_ -> point(xx_rel,yy_rel-1), stride_rel, &pmlArr_[kk].Dx_end_ -> point(xx-1,yy-1), stride);
                zaxpy_(nZax,      zaxArr[5], &Hz_ -> point(xx_rel,yy_rel  ), stride_rel, &pmlArr_[kk].Dx_end_ -> point(xx-1,yy-1), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].Dx_end_ -> point(xx-1,yy-1), stride, &Ex_ -> point(xx_rel,yy_rel), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], dxstore.data()                      , 1     , &Ex_ -> point(xx_rel,yy_rel), stride_rel);
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxEy_.size();zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxEy_[zz];
                int xx = static_cast<int>(zaxArr[0]);
                int yy = static_cast<int>(zaxArr[1]);
                int nZax = static_cast<int>(zaxArr[2]);
                vector<complex<double>> dystore(nZax, 0.0);

                zcopy_(nZax, &pmlArr_[kk].Dy_ -> point(xx-1,yy-1), stride, dystore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].Dy_ -> point(xx-1,yy-1), stride);
                zscal_(nZax, zaxArr[6],             &Ey_ -> point(xx,yy), stride_rel);

                zaxpy_(nZax, -1.0*zaxArr[5], &Hz_ -> point(xx  ,yy), stride_rel, &pmlArr_[kk].Dy_ -> point(xx-1,yy-1), stride);
                zaxpy_(nZax,      zaxArr[5], &Hz_ -> point(xx-1,yy), stride_rel, &pmlArr_[kk].Dy_ -> point(xx-1,yy-1), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].Dy_ -> point(xx-1,yy-1), stride, &Ey_ -> point(xx,yy), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], dystore.data()                  , 1     , &Ey_ -> point(xx,yy), stride_rel);
            }
            for(int zz = 0; zz < pmlArr_[kk].zaxEy_end_.size();zz++)
            {
                array<double,9> zaxArr = pmlArr_[kk].zaxEy_end_[zz];
                int xx = static_cast<int>(zaxArr[0]);
                int yy = static_cast<int>(zaxArr[1]);
                int xx_rel = ni + pow(-1, 1-d) * xx+2*(1-d); // If the pml is is the x direction xx_rel = nx - xx, else xrel = xx
                int yy_rel = nj + pow(-1, d)   * yy+2*d; // If the pml is is the y direction yy_rel = ny - yy, else yrel = yy
                int nZax = static_cast<int>(zaxArr[2]);

                vector<complex<double>> dystore(nZax, 0.0);
                zcopy_(nZax, &pmlArr_[kk].Dy_end_ -> point(xx-1,yy-1), stride, dystore.data(), 1);

                zscal_(nZax, zaxArr[4], &pmlArr_[kk].Dy_end_ -> point(xx-1,yy-1), stride);
                zscal_(nZax, zaxArr[6],                 &Ey_ -> point(xx_rel,yy_rel), stride_rel);

                zaxpy_(nZax, -1.0*zaxArr[5], &Hz_ -> point(xx_rel  ,yy_rel), stride_rel, &pmlArr_[kk].Dy_end_ -> point(xx-1,yy-1), stride);
                zaxpy_(nZax,      zaxArr[5], &Hz_ -> point(xx_rel-1,yy_rel), stride_rel, &pmlArr_[kk].Dy_end_ -> point(xx-1,yy-1), stride);

                zaxpy_(nZax,      zaxArr[7], &pmlArr_[kk].Dy_end_ -> point(xx-1,yy-1), stride, &Ey_ -> point(xx_rel,yy_rel), stride_rel);
                zaxpy_(nZax, -1.0*zaxArr[8], dystore.data()                      , 1     , &Ey_ -> point(xx_rel,yy_rel), stride_rel);
            }
        }
        if(pmlArr_.size() > 1)
        {
            // Setting everything such that the X-PML stores all the information for the corners
            int kx = 1;
            shared_ptr<vector<vector<array<double,5>>>> c_ex_0_n;
            shared_ptr<vector<vector<array<double,5>>>> c_ex_n_0;
            shared_ptr<vector<vector<array<double,5>>>> c_ey_0_n;
            shared_ptr<vector<vector<array<double,5>>>> c_ey_n_0;
            shared_ptr<vector<vector<array<double,5>>>> c_ex_0_0 = pmlArr_[kx].c_ex_0_0_;
            shared_ptr<vector<vector<array<double,5>>>> c_ex_n_n = pmlArr_[kx].c_ex_n_n_;
            shared_ptr<vector<vector<array<double,5>>>> c_ey_0_0 = pmlArr_[kx].c_ey_0_0_;
            shared_ptr<vector<vector<array<double,5>>>> c_ey_n_n = pmlArr_[kx].c_ey_n_n_;
            if(pmlArr_[1].d() == X)
            {
                c_ex_0_n = pmlArr_[1].c_ex_0_n_;
                c_ex_n_0 = pmlArr_[1].c_ex_n_0_;
                c_ey_0_n = pmlArr_[1].c_ey_0_n_;
                c_ey_n_0 = pmlArr_[1].c_ey_n_0_;
            }
            else
            {
                kx = 0;
                c_ex_0_n = pmlArr_[1].c_ex_n_0_;
                c_ex_n_0 = pmlArr_[1].c_ex_0_n_;
                c_ey_0_n = pmlArr_[1].c_ey_n_0_;
                c_ey_n_0 = pmlArr_[1].c_ey_0_n_;
            }

            complex<double> dxstore(0.0,0.0); complex<double> dystore(0.0,0.0);
            for(int ii = 0; ii < xPML_; ii++)
            {
                for(int jj = 0; jj < yPML_; jj++)
                {
                    //Bot Left
                    int xx = ii+1; int yy = jj+1;
                    dxstore = pmlArr_[kx].Dx_->point(ii,yy-1);
                    dystore = pmlArr_[kx].Dy_->point(ii,yy-1);
                    pmlArr_[kx].Dx_->point(ii,yy-1) = c_ex_0_0->at(ii).at(jj)[0] * pmlArr_[kx].Dx_->point(ii,yy-1) + c_ex_0_0->at(ii).at(jj)[1] * (Hz_->point(xx,yy)-Hz_->point(xx,yy-1));
                    pmlArr_[kx].Dy_->point(ii,yy-1) = c_ey_0_0->at(ii).at(jj)[0] * pmlArr_[kx].Dy_->point(ii,yy-1) - c_ey_0_0->at(ii).at(jj)[1] * (Hz_->point(xx,yy)-Hz_->point(xx-1,yy));
                    Ex_->point(xx,yy) = c_ex_0_0->at(ii).at(jj)[2] * Ex_->point(xx,yy) + c_ex_0_0->at(ii).at(jj)[3] * pmlArr_[kx].Dx_->point(ii,yy-1) - c_ex_0_0->at(ii).at(jj)[4] * dxstore;
                    Ey_->point(xx,yy) = c_ey_0_0->at(ii).at(jj)[2] * Ey_->point(xx,yy) + c_ey_0_0->at(ii).at(jj)[3] * pmlArr_[kx].Dy_->point(ii,yy-1) - c_ey_0_0->at(ii).at(jj)[4] * dystore;

                    //Top Left
                    xx = ii+1; yy = ny_-jj;
                    dxstore = pmlArr_[kx].Dx_->point(ii,yy-1);
                    dystore = pmlArr_[kx].Dy_->point(ii,yy-1);
                    pmlArr_[kx].Dx_->point(ii,yy-1) = c_ex_0_n->at(ii).at(jj)[0] * pmlArr_[kx].Dx_->point(ii,yy-1) + c_ex_0_n->at(ii).at(jj)[1] * (Hz_->point(xx,yy)-Hz_->point(xx,yy-1));
                    pmlArr_[kx].Dy_->point(ii,yy-1) = c_ey_0_n->at(ii).at(jj)[0] * pmlArr_[kx].Dy_->point(ii,yy-1) - c_ey_0_n->at(ii).at(jj)[1] * (Hz_->point(xx,yy)-Hz_->point(xx-1,yy));
                    Ex_->point(xx,yy) = c_ex_0_n->at(ii).at(jj)[2] * Ex_->point(xx,yy) + c_ex_0_n->at(ii).at(jj)[3] * pmlArr_[kx].Dx_->point(ii,yy-1) - c_ex_0_n->at(ii).at(jj)[4] * dxstore;
                    Ey_->point(xx,yy) = c_ey_0_n->at(ii).at(jj)[2] * Ey_->point(xx,yy) + c_ey_0_n->at(ii).at(jj)[3] * pmlArr_[kx].Dy_->point(ii,yy-1) - c_ey_0_n->at(ii).at(jj)[4] * dystore;

                    //Top Right
                    xx = nx_ -ii; yy = ny_-jj;
                    dxstore = pmlArr_[kx].Dx_end_->point(ii,yy-1);
                    dystore = pmlArr_[kx].Dy_end_->point(ii,yy-1);
                    pmlArr_[kx].Dx_end_->point(ii,yy-1) = c_ex_n_n->at(ii).at(jj)[0] * pmlArr_[kx].Dx_end_->point(ii,yy-1) + c_ex_n_n->at(ii).at(jj)[1] * (Hz_->point(xx,yy)-Hz_->point(xx,yy-1));
                    pmlArr_[kx].Dy_end_->point(ii,yy-1) = c_ey_n_n->at(ii).at(jj)[0] * pmlArr_[kx].Dy_end_->point(ii,yy-1) - c_ey_n_n->at(ii).at(jj)[1] * (Hz_->point(xx,yy)-Hz_->point(xx-1,yy));
                    Ex_->point(xx,yy) = c_ex_n_n->at(ii).at(jj)[2] * Ex_->point(xx,yy) + c_ex_n_n->at(ii).at(jj)[3] * pmlArr_[kx].Dx_end_->point(ii,yy-1) - c_ex_n_n->at(ii).at(jj)[4] * dxstore;
                    Ey_->point(xx,yy) = c_ey_n_n->at(ii).at(jj)[2] * Ey_->point(xx,yy) + c_ey_n_n->at(ii).at(jj)[3] * pmlArr_[kx].Dy_end_->point(ii,yy-1) - c_ey_n_n->at(ii).at(jj)[4] * dystore;

                    //Bot Right
                    xx = nx_-ii; yy = jj+1;
                    dxstore = pmlArr_[kx].Dx_end_->point(ii,yy-1);
                    dystore = pmlArr_[kx].Dy_end_->point(ii,yy-1);
                    pmlArr_[kx].Dx_end_->point(ii,yy-1) = c_ex_n_0->at(ii).at(jj)[0] * pmlArr_[kx].Dx_end_->point(ii,yy-1) + c_ex_n_0->at(ii).at(jj)[1] * (Hz_->point(xx,yy)-Hz_->point(xx,yy-1));
                    pmlArr_[kx].Dy_end_->point(ii,yy-1) = c_ey_n_0->at(ii).at(jj)[0] * pmlArr_[kx].Dy_end_->point(ii,yy-1) - c_ey_n_0->at(ii).at(jj)[1] * (Hz_->point(xx,yy)-Hz_->point(xx-1,yy));
                    Ex_->point(xx,yy) = c_ex_n_0->at(ii).at(jj)[2] * Ex_->point(xx,yy) + c_ex_n_0->at(ii).at(jj)[3] * pmlArr_[kx].Dx_end_->point(ii,yy-1) - c_ex_n_0->at(ii).at(jj)[4] * dxstore;
                    Ey_->point(xx,yy) = c_ey_n_0->at(ii).at(jj)[2] * Ey_->point(xx,yy) + c_ey_n_0->at(ii).at(jj)[3] * pmlArr_[kx].Dy_end_->point(ii,yy-1) - c_ey_n_0->at(ii).at(jj)[4] * dystore;
                }
            }
        }
    }
}

void FDTDField::updateDisp()
{
    if(Ez_)
    {
        for(int zz = 0; zz < zaxDz_.size(); zz++)
        {
            int xx = zaxDz_[zz][0]; int yy = zaxDz_[zz][1]; int nZax = zaxDz_[zz][2];
            double eps = objArr_[zaxDz_[zz][3]].dielectric();
            vector<double> alpha = objArr_[zaxDz_[zz][3]].alpha();
            vector<double> gamma = objArr_[zaxDz_[zz][3]].gamma();
            vector<double> zi    = objArr_[zaxDz_[zz][3]].zi();
            for(int pp = 0; pp < lorPz_.size(); pp ++)
            {
                vector<complex<double>> jzstore(nZax, 0.0);
                zcopy_(nZax, &lorPz_[pp]->point(xx,yy), 1, jzstore.data(), 1);

                zscal_(nZax, alpha[pp],     &lorPz_[pp]->point(xx,yy), 1);
                zaxpy_(nZax,    zi[pp], &prevLorPz_[pp]->point(xx,yy), 1, &lorPz_[pp] ->point(xx,yy), 1);
                zaxpy_(nZax, gamma[pp],            &Ez_->point(xx,yy), 1, &lorPz_[pp] ->point(xx,yy), 1);

                zcopy_(nZax, jzstore.data(),1, &prevLorPz_[pp]->point(xx,yy),1);
            }

            // zscal_(nZax,      c_eze, &Dz_->point(xx  ,yy  ), 1);
            zaxpy_(nZax,      dt_/dy_, &Hx_->point(xx  ,yy-1), 1, &Dz_ ->point(xx,yy), 1);
            zaxpy_(nZax, -1.0*dt_/dy_, &Hx_->point(xx  ,yy  ), 1, &Dz_ ->point(xx,yy), 1);
            zaxpy_(nZax, -1.0*dt_/dx_, &Hy_->point(xx-1,yy  ), 1, &Dz_ ->point(xx,yy), 1);
            zaxpy_(nZax,      dt_/dx_, &Hy_->point(xx  ,yy  ), 1, &Dz_ ->point(xx,yy), 1);

            zcopy_(nZax, &Dz_->point(xx,yy), 1,&Ez_->point(xx,yy), 1);
            zscal_(nZax, 1/eps, &Ez_->point(xx,yy), 1);
            for(int pp =0; pp < lorPz_.size(); pp ++)
                zaxpy_(nZax, -1.0/eps, &lorPz_[pp]->point(xx,yy), 1,&Ez_ ->point(xx,yy), 1);
        }
    }
    else
    {
        for(int zz = 0; zz < zaxDx_.size(); zz++)
        {
            int xx = zaxDx_[zz][0]; int yy = zaxDx_[zz][1]; int nZax = zaxDx_[zz][2];
            double eps = objArr_[zaxDx_[zz][3]].dielectric();
            vector<double> alpha = objArr_[zaxDx_[zz][3]].alpha();
            vector<double> gamma = objArr_[zaxDx_[zz][3]].gamma();
            vector<double> zi    = objArr_[zaxDx_[zz][3]].zi();
            for(int pp = 0; pp < lorPx_.size(); pp ++)
            {
                vector<complex<double>> jzstore(nZax, 0.0);
                zcopy_(nZax, &lorPx_[pp]->point(xx,yy), 1, jzstore.data(), 1);

                zscal_(nZax, alpha[pp],     &lorPx_[pp]->point(xx,yy), 1);
                zaxpy_(nZax,    zi[pp], &prevLorPx_[pp]->point(xx,yy), 1, &lorPx_[pp] ->point(xx,yy), 1);
                zaxpy_(nZax, gamma[pp],          &Ex_->point(xx,yy), 1, &lorPx_[pp] ->point(xx,yy), 1);

                zcopy_(nZax, jzstore.data(),1, &prevLorPx_[pp]->point(xx,yy),1);
            }

            // zscal_(nZax,      1.0    , &Dx_->point(xx  ,yy  ), 1);
            zaxpy_(nZax,      dt_/dy_, &Hz_->point(xx  ,yy  ), 1, &Dx_->point(xx,yy), 1);
            zaxpy_(nZax, -1.0*dt_/dy_, &Hz_->point(xx  ,yy-1), 1, &Dx_->point(xx,yy), 1);

            zcopy_(nZax, &Dx_->point(xx,yy), 1,&Ex_->point(xx,yy), 1);
            zscal_(nZax, 1/eps, &Ex_->point(xx,yy), 1);
            for(int pp =0; pp < lorPx_.size(); pp ++)
                zaxpy_(nZax, -1.0/eps, &lorPx_[pp]->point(xx,yy), 1,&Ex_ ->point(xx,yy), 1);
        }
        for(int zz = 0; zz < zaxDy_.size(); zz++)
        {
            int xx = zaxDy_[zz][0]; int yy = zaxDy_[zz][1]; int nZax = zaxDy_[zz][2];
            double eps = objArr_[zaxDy_[zz][3]].dielectric();
            vector<double> alpha = objArr_[zaxDy_[zz][3]].alpha();
            vector<double> gamma = objArr_[zaxDy_[zz][3]].gamma();
            vector<double> zi    = objArr_[zaxDy_[zz][3]].zi();
            for(int pp = 0; pp < lorPy_.size(); pp ++)
            {
                vector<complex<double>> jystore(nZax, 0.0);
                zcopy_(nZax, &lorPy_[pp]->point(xx,yy), 1, jystore.data(), 1);

                zscal_(nZax, alpha[pp],     &lorPy_[pp]->point(xx,yy), 1);
                zaxpy_(nZax,    zi[pp], &prevLorPy_[pp]->point(xx,yy), 1, &lorPy_[pp] ->point(xx,yy), 1);
                zaxpy_(nZax, gamma[pp],            &Ey_->point(xx,yy), 1, &lorPy_[pp] ->point(xx,yy), 1);

                zcopy_(nZax, jystore.data(),1, &prevLorPy_[pp]->point(xx,yy),1);
            }
            zaxpy_(nZax, -1.0*dt_/dx_, &Hz_->point(xx  ,yy  ), 1, &Dy_->point(xx,yy), 1);
            zaxpy_(nZax,      dt_/dx_, &Hz_->point(xx-1,yy  ), 1, &Dy_->point(xx,yy), 1);

            zcopy_(nZax, &Dy_->point(xx,yy), 1,&Ey_->point(xx,yy), 1);
            zscal_(nZax, 1/eps, &Ey_->point(xx,yy), 1);
            for(int pp =0; pp < lorPy_.size(); pp ++)
                zaxpy_(nZax, -1.0/eps, &lorPy_[pp]->point(xx,yy), 1,&Ey_ ->point(xx,yy), 1);
        }
    }
}
void FDTDField::applyPBC(shared_ptr<Grid2D<complex<double>>> fUp, int nx, int ny)
{
    vector<double> r = {0.0,0.0};
    for(int ii = 1; ii < nx; ii++)
    {
        r = {(ii-1)*dx_,(ny-1)*dy_};
        fUp->point(ii,0) = per_factor(r) * fUp->point(ii,ny-1);
        r = {(ii-1)*dx_,ny_*dy_};
        fUp->point(ii,ny) = per_factor(r) * fUp->point(ii,1);
    }
    for(int jj = 1; jj < ny; jj++)
    {
        r = {(nx-1)*dx_,(jj-1)*dy_};
        fUp->point(0,jj) = per_factor(r) * fUp->point(nx-1,jj);
        r = {0*dx_,(jj-1)*dy_};
        fUp->point(nx,jj) = per_factor(r) * fUp->point(1,jj);
    }
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
        int ii = srcArr_[kk].loc()[0];
        int jj = srcArr_[kk].loc()[1];
        switch ( srcArr_[kk].pol() )
        {
            case EZ: //if(srcArr[kk].pol() == EZ)
                if(abs(real(srcArr_[kk].prof().pulse(static_cast<double>(t_step_)))) > 1.0e-70)
                {
                    double eps = objArr_[phys_Ez_->point(ii,jj)].dielectric();
                    double c_ezj = dt_/(eps);
                    Ez_ -> point(ii,jj) += c_ezj * srcArr_[kk].prof().pulse(static_cast<double>(t_step_));
                }
                break;
            case HX: //else if(srcArr[kk].pol() == HX)
                if(abs(real(srcArr_[kk].prof().pulse(static_cast<double>(t_step_)))) > 1.0e-70)
                {
                    Hx_ -> point(ii,jj) += dt_ * srcArr_[kk].prof().pulse(static_cast<double>(t_step_));
                }
                break;
            case HY: //else if(srcArr[kk].pol() == HY)
                if(abs(real(srcArr_[kk].prof().pulse(static_cast<double>(t_step_)))) > 1.0e-70)
                {
                    Hy_ -> point(ii,jj) += dt_ * srcArr_[kk].prof().pulse(static_cast<double>(t_step_));
                }
                break;
            case HZ: //else if(srcArr[kk].pol() == HZ)
                if(abs(real(srcArr_[kk].prof().pulse(static_cast<double>(t_step_)))) > 1.0e-70)
                {
                    Hz_ -> point(ii,jj) += dt_ * srcArr_[kk].prof().pulse(static_cast<double>(t_step_));
                }
                break;
            case EX: //else if(srcArr[kk].pol() == EX)
                if(abs(real(srcArr_[kk].prof().pulse(static_cast<double>(t_step_)))) > 1.0e-70)
                {
                    double eps = objArr_[phys_Ex_->point(ii,jj)].dielectric();
                    double c_exj = dt_/(eps);
                    Ex_ -> point(ii,jj) += c_exj * srcArr_[kk].prof().pulse(static_cast<double>(t_step_));
                }
                break;
            case EY: //else if(srcArr[kk].pol() == EY)
                if(abs(real(srcArr_[kk].prof().pulse(static_cast<double>(t_step_)))) > 1.0e-70)
                {
                    double eps = objArr_[phys_Ey_->point(ii,jj)].dielectric();
                    double c_eyj = dt_/(eps);
                    Ey_ -> point(ii,jj) += c_eyj * srcArr_[kk].prof().pulse(static_cast<double>(t_step_));
                }
                break;
            default:
                throw logic_error("reached a default case in a switch state that should never happen!");
                break;
        }
    }
    updateH();
    if(periodic_)
    {
        if(Hz_)
            applyPBC(Hz_,nx_,ny_);
        else
        {
            applyPBC(Hx_,nx_+1,ny_);
            applyPBC(Hy_,nx_,ny_+1);
        }
    }
    updateE();
    updateDisp();
    if(periodic_)
    {
        if (Hz_)
        {
            applyPBC(Ex_,nx_,ny_+1);
            applyPBC(Ey_,nx_+1,ny_);
        }
        else
            applyPBC(Ez_,nx_+1,ny_+1);
    }
    for(int ii = 0; ii < dtcArr_.size(); ii ++)
        ouputField(dtcArr_[ii]);
    // if(abs(tcur_-floor(tcur_+0.5)) < 1e-7)
    if(false)
    {
        if(Hx_)
        {
            string fname("fout/Hx/HxField_t" + to_string(static_cast<int>(t_step_))+".dat");
            Hx_->gridOut(fname);
            fname = "fout/Hy/HyField_t" + to_string(static_cast<int>(t_step_))+".dat";
            Hy_->gridOut(fname);
            fname = "fout/Ez/EzField_t" + to_string(static_cast<int>(t_step_))+".dat";
            Ez_->gridOut(fname);
        }
        else
        {
            string fname("fout/Hx/HxField_t" + to_string(static_cast<int>(t_step_))+".dat");
            Ex_->gridOut(fname);
            fname = "fout/Hy/HyField_t" + to_string(static_cast<int>(t_step_))+".dat";
            Ey_->gridOut(fname);
            fname = "fout/Ez/EzField_t" + to_string(static_cast<int>(t_step_))+".dat";
            Hz_->gridOut(fname);
        }
    }
    tcur_ += dt_;
    t_step_ ++;
}