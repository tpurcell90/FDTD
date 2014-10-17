#ifndef FDTD_PML
#define FDTD_PML

#include "Source.hpp"
#include "Grid.hpp"
#include "Obj.hpp"
// #include <assert.h>
// #include <iomanip>
// #include <iostream>
#include <memory>
// #include <random>
// #include <stdexcept>
// #include <string>
#include <vector>
#include <map>
#include <tuple>

// #include <complex>
// typedef std::complex<double> cplx;
enum Direction {X,Y,Z};

template <typename T> class UPML
{
typedef  double (UPML<T>::*PMLMemFn)(double x, double y);
// typedef std::vector<std::vector<std::array<double,5>>>> corConsts;
protected:
    int thickness_;
    Direction d_;
    //PMLTopBot tb_;
    double m_;
    double R0_;
    double sigmaMax_;
    double kappaMax_;
    Polarization pol_;
    int nj_,ni_;
    double dx_,dy_,dt_;

public:
    std::shared_ptr<Grid2D<T>> Dx_,Dy_,Dz_,Bx_,By_,Bz_,Dx_end_,Dy_end_,Dz_end_,Bx_end_,By_end_,Bz_end_;
    std::shared_ptr<Grid2D<int>> phys_Hx_,phys_Hy_,phys_Hz_, phys_Hx_end_,phys_Hy_end_,phys_Hz_end_,phys_Ex_,phys_Ey_,phys_Ez_, phys_Ex_end_,phys_Ey_end_,phys_Ez_end_;

    std::shared_ptr<std::vector<std::vector<std::array<double,5>>>> c_hx_0_0_, c_hy_0_0_, c_ez_0_0_, c_hx_n_0_, c_hy_n_0_, c_ez_n_0_;
    std::shared_ptr<std::vector<std::vector<std::array<double,5>>>> c_ex_0_0_, c_ey_0_0_, c_hz_0_0_, c_ex_n_0_, c_ey_n_0_, c_hz_n_0_;
    std::shared_ptr<std::vector<std::vector<std::array<double,5>>>> c_hx_0_n_, c_hy_0_n_, c_ez_0_n_, c_hx_n_n_, c_hy_n_n_, c_ez_n_n_;
    std::shared_ptr<std::vector<std::vector<std::array<double,5>>>> c_ex_0_n_, c_ey_0_n_, c_hz_0_n_, c_ex_n_n_, c_ey_n_n_, c_hz_n_n_;

    std::vector<std::array<double,9>> zaxHx_, zaxHy_, zaxEz_, zaxHx_end_, zaxHy_end_, zaxEz_end_;
    std::vector<std::array<double,9>> zaxEx_, zaxEy_, zaxHz_, zaxEx_end_, zaxEy_end_, zaxHz_end_;

    int edgei_0_, edgei_n_, edgej_0_0_, edgej_0_n_, edgej_n_0_, edgej_n_n_;

    /**
     * @brief Constrcuts a PML for both ends of the cell
     * @details Uses the input to construct the functions and auxiliary fields for the PML calculations
     *
     * @param thickness Thickness of PML in number of unit cells
     * @param d Direction the PML is occupying
     * @param m Order of the polynomial grading for the kappa and sigma vectors
     * @param R0 Max allowed reflection for normally incident waves
     * @param nx Number of unit cells in the x direction
     * @param ny Number of unit cells in the y direction
     * @param dx unit cell spacing for the x direction
     * @param dy unit cell spacing for the y direction
     * @param xPML thickness of the x direction PML
     * @param yPML thickness of the y direction PML
     * @param pol A polarization so it can set up the right auxilliary fields
     */

    UPML(int thickness, Direction d, double m, double R0, int nx, int ny, double dx, double dy, double dt, int xPML, int yPML, Polarization pol) : thickness_(thickness), d_(d), m_(m), R0_(R0), dx_(dx), dy_(dy), dt_(dt), pol_(pol)
    {
        edgei_0_ = 0;
        edgei_n_ = 0;
        edgej_0_0_ = 0; edgej_0_n_ = 0;
        edgej_n_0_ = 0; edgej_n_n_ = 0;
        //sigmaMax_ = -(m_+1)*log(R0_)/(2*thickness_*dx); // eta should be included;
        sigmaMax_ = 0.8*(m_+1)/dx_;
        kappaMax_ = 1.0;
        zaxHx_ = {}; zaxHy_ = {}; zaxEz_ = {}; zaxHx_end_ = {}; zaxHy_end_ = {}; zaxEz_end_ = {};
        zaxEx_ = {}; zaxEy_ = {}; zaxHz_ = {}; zaxEx_end_ = {}; zaxEy_end_ = {}; zaxHz_end_ = {};
        int xmax; int ymax;
        if (d == X)
        {
            ni_ = nx;
            nj_ = ny;
            xmax = thickness_;
            ymax = ny;
        }
        else
        {
            xmax = nx;
            ymax = thickness_;
            nj_ = nx;
            ni_ = ny;
        }
        // Create fields and precostant containers
        if(pol == EX || pol == EY || pol == HZ)
        {
            Dx_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dy_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Bz_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dx_end_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dy_end_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Bz_end_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);

            Bx_ = nullptr;
            By_ = nullptr;
            Dz_ = nullptr;
            Bx_end_ = nullptr;
            By_end_ = nullptr;
            Dz_end_ = nullptr;

            phys_Hz_ = std::make_shared<Grid2D<int>>(xmax+2,ymax+2,dx,dy);
            phys_Hz_end_ = std::make_shared<Grid2D<int>>(xmax+2,ymax+2,dx,dy);
            phys_Ex_ = std::make_shared<Grid2D<int>>(xmax+2,ymax+2,dx,dy);
            phys_Ex_end_ = std::make_shared<Grid2D<int>>(xmax+2,ymax+2,dx,dy);
            phys_Ey_ = std::make_shared<Grid2D<int>>(xmax+2,ymax+2,dx,dy);
            phys_Ey_end_ = std::make_shared<Grid2D<int>>(xmax+2,ymax+2,dx,dy);

            phys_Hy_ = nullptr;
            phys_Hy_end_ = nullptr;
            phys_Hx_ = nullptr;
            phys_Hx_end_ = nullptr;
            phys_Ez_ = nullptr;
            phys_Ez_end_ = nullptr;
            if(xPML == 0)
            {
                c_hx_0_0_ = nullptr; c_hy_0_0_ = nullptr; c_ez_0_0_ = nullptr; c_hx_n_0_ = nullptr; c_hy_n_0_ = nullptr; c_ez_n_0_ = nullptr;
                c_hx_0_n_ = nullptr; c_hy_0_n_ = nullptr; c_ez_0_n_ = nullptr; c_hx_n_n_ = nullptr; c_hy_n_n_ = nullptr; c_ez_n_n_ = nullptr;
                c_ex_0_0_ = nullptr; c_ey_0_0_ = nullptr;  c_ex_n_0_ = nullptr; c_ey_n_0_ = nullptr;
                c_ex_0_n_ = nullptr; c_ey_0_n_ = nullptr;  c_ex_n_n_ = nullptr; c_ey_n_n_ = nullptr;
                c_hz_0_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(1, std::vector<std::array<double,5>>(yPML));
                c_hz_n_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(1, std::vector<std::array<double,5>>(yPML));
                c_hz_0_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(1, std::vector<std::array<double,5>>(yPML));
                c_hz_n_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(1, std::vector<std::array<double,5>>(yPML));
            }
            else if(yPML == 0)
            {
                c_hx_0_0_ = nullptr; c_hy_0_0_ = nullptr; c_ez_0_0_ = nullptr; c_hx_n_0_ = nullptr; c_hy_n_0_ = nullptr; c_ez_n_0_ = nullptr;
                c_hx_0_n_ = nullptr; c_hy_0_n_ = nullptr; c_ez_0_n_ = nullptr; c_hx_n_n_ = nullptr; c_hy_n_n_ = nullptr; c_ez_n_n_ = nullptr;
                c_ex_0_0_ = nullptr; c_ey_0_0_ = nullptr;  c_ex_n_0_ = nullptr; c_ey_n_0_ = nullptr;
                c_ex_0_n_ = nullptr; c_ey_0_n_ = nullptr;  c_ex_n_n_ = nullptr; c_ey_n_n_ = nullptr;
                c_hz_0_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(1));
                c_hz_n_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(1));
                c_hz_0_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(1));
                c_hz_n_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(1));
            }
            else
            {
                c_hx_0_0_ = nullptr; c_hy_0_0_ = nullptr; c_ez_0_0_ = nullptr; c_hx_n_0_ = nullptr; c_hy_n_0_ = nullptr; c_ez_n_0_ = nullptr;
                c_hx_0_n_ = nullptr; c_hy_0_n_ = nullptr; c_ez_0_n_ = nullptr; c_hx_n_n_ = nullptr; c_hy_n_n_ = nullptr; c_ez_n_n_ = nullptr;

                c_ex_0_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_ey_0_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_hz_0_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_ex_n_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_ey_n_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_hz_n_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_ex_0_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_ey_0_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_hz_0_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_ex_n_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_ey_n_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_hz_n_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
            }
        }
        else
        {
            Bx_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            By_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dz_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dx_ = nullptr;
            Dy_ = nullptr;
            Bz_ = nullptr;

            Bx_end_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            By_end_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dz_end_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dx_end_ = nullptr;
            Dy_end_ = nullptr;
            Bz_end_ = nullptr;

            phys_Hz_ = nullptr;
            phys_Hz_end_ = nullptr;
            phys_Hy_ = std::make_shared<Grid2D<int>>(xmax+2,ymax+2,dx,dy);
            phys_Hy_end_ = std::make_shared<Grid2D<int>>(xmax+2,ymax+2,dx,dy);
            phys_Hx_ = std::make_shared<Grid2D<int>>(xmax+2,ymax+2,dx,dy);
            phys_Hx_end_ = std::make_shared<Grid2D<int>>(xmax+2,ymax+2,dx,dy);

            phys_Ez_ = std::make_shared<Grid2D<int>>(xmax+2,ymax+2,dx,dy);
            phys_Ez_end_ = std::make_shared<Grid2D<int>>(xmax+2,ymax+2,dx,dy);
            phys_Ex_ = nullptr;
            phys_Ex_end_ = nullptr;
            phys_Ey_ = nullptr;
            phys_Ey_end_ = nullptr;
            // Corners v. no corners check. The preconstants for edges in the non corner case still need to be calculated for Ez
            if(yPML == 0)
            {
                c_ex_0_0_ = nullptr; c_ey_0_0_ = nullptr; c_hz_0_0_ = nullptr; c_ex_n_0_ = nullptr; c_ey_n_0_ = nullptr; c_hz_n_0_ = nullptr;
                c_ex_0_n_ = nullptr; c_ey_0_n_ = nullptr; c_hz_0_n_ = nullptr; c_ex_n_n_ = nullptr; c_ey_n_n_ = nullptr; c_hz_n_n_ = nullptr;
                c_hx_0_0_ = nullptr; c_hy_0_0_ = nullptr;  c_hx_n_0_ = nullptr; c_hy_n_0_ = nullptr;
                c_hx_0_n_ = nullptr; c_hy_0_n_ = nullptr;  c_hx_n_n_ = nullptr; c_hy_n_n_ = nullptr;
                c_ez_0_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(1));
                c_ez_n_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(1));
                c_ez_0_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(1));
                c_ez_n_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(1));
            }
            else if(xPML == 0)
            {
                c_ex_0_0_ = nullptr; c_ey_0_0_ = nullptr; c_hz_0_0_ = nullptr; c_ex_n_0_ = nullptr; c_ey_n_0_ = nullptr; c_hz_n_0_ = nullptr;
                c_ex_0_n_ = nullptr; c_ey_0_n_ = nullptr; c_hz_0_n_ = nullptr; c_ex_n_n_ = nullptr; c_ey_n_n_ = nullptr; c_hz_n_n_ = nullptr;
                c_hx_0_0_ = nullptr; c_hy_0_0_ = nullptr;  c_hx_n_0_ = nullptr; c_hy_n_0_ = nullptr;
                c_hx_0_n_ = nullptr; c_hy_0_n_ = nullptr;  c_hx_n_n_ = nullptr; c_hy_n_n_ = nullptr;
                c_ez_0_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(1, std::vector<std::array<double,5>>(yPML));
                c_ez_n_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(1, std::vector<std::array<double,5>>(yPML));
                c_ez_0_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(1, std::vector<std::array<double,5>>(yPML));
                c_ez_n_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(1, std::vector<std::array<double,5>>(yPML));
            }
            else
            {
                c_ex_0_0_ = nullptr; c_ey_0_0_ = nullptr; c_hz_0_0_ = nullptr; c_ex_n_0_ = nullptr; c_ey_n_0_ = nullptr; c_hz_n_0_ = nullptr;
                c_ex_0_n_ = nullptr; c_ey_0_n_ = nullptr; c_hz_0_n_ = nullptr; c_ex_n_n_ = nullptr; c_ey_n_n_ = nullptr; c_hz_n_n_ = nullptr;
                c_hx_0_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_hy_0_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_ez_0_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_hx_n_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_hy_n_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_ez_n_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_hx_0_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_hy_0_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_ez_0_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_hx_n_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_hy_n_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
                c_ez_n_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xPML, std::vector<std::array<double,5>>(yPML));
            }
        }
    }

    void initializeUPML(std::vector<Obj> objArr, int nx, int ny, double dx, double dy, double dt, int yPML, int xPML, UPML<T> *opp)
    {
        int jmax = 0;
        int pt_i = 0; int pt_j = 0;
        int ni   = 0; int nj   = 0;
        double di = 0.0; double dj = 0.0;
        int xmax; int ymax;
        int oppPML = 0;
        int delx =0; int dely = 0;
        int zaxJmax = 0; int zaxXJmax = 0; int zaxYJmax = 0;
        double eps=0.0;
        double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;

        double sigz = 0.0; double sigx = 0.0; double sigy = 0.0;
        double sigxx = 0.0; double sigxy = 0.0; double sigyx = 0.0; double sigyy = 0.0;

        int jj = 0; int ii = 0;
        int * xx; int * yy;

        PMLMemFn  sigmay;
        PMLMemFn  sigmax;

        UPML<T> *xpml; UPML<T> *ypml;
        if (d_ == X)
        {
            jmax = ny-1;
            pt_i = 0;  pt_j = 1;
            ni   = nx; nj   = ny;
            di   = dx; dj   = dy;
            xmax = thickness_; ymax = ny;
            sigmax = &UPML<T>::sigma;
            xpml = this;
            if (opp)
            {
                sigmay = opp ->sig_ptr();
                ypml = opp;
            }
            else
            {
                sigmay = &UPML<T>::sigmaopp;
                ypml = this;
            }
            xmax = thickness_; ymax = nj_;
            oppPML = yPML;
            delx = 0; dely = 1;
            zaxJmax = ny -yPML-1;
            zaxXJmax = zaxJmax;
            zaxYJmax = zaxJmax;
            xx = &ii;  yy = &jj;
        }
        else
        {
            jmax = nx-1;
            pt_i = 1; pt_j = 0;
            ni   = ny; nj   = nx;
            di   = dy; dj   = dx;
            xmax = nx; ymax = thickness_;
            sigmay = &UPML<T>::sigma;
            ypml = this;
            if (opp)
            {
                sigmax = opp ->sig_ptr();
                xpml = opp;
            }
            else
            {
                sigmax = &UPML<T>::sigmaopp;
                xpml = this;
            }
            ymax = thickness_; xmax = nj_;
            oppPML = xPML;
            delx = 1; dely=0;
            zaxJmax = nx -xPML-1;
            zaxXJmax = zaxJmax;
            zaxYJmax = zaxJmax;
            xx = &jj; yy = &ii;
        }
        int zmin = oppPML; int zmax = nj - oppPML - 1;
        if(oppPML == 0)
        {
            zaxXJmax -= dely;
            zaxYJmax -= delx;
        }
        for(int kk = 0; kk < objArr.size(); kk++)
        {
            std::vector<double> pt(2,0.0);
            if(pol_ == HZ || pol_ == EX || pol_ == EY)
            {
                for(ii =1; ii < thickness_+1; ii ++)
                {
                    for(jj = 1; jj < nj+1; jj++)
                    {
                        pt[0] = ((*xx)+0.5-(nx-1)/2.0)*dx;
                        pt[1] = ((*yy)-(ny-1)/2.0)*dy;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Ex_->point(*xx,*yy) = kk;

                        pt[1] += 0.5*dy;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Hz_->point(*xx,*yy) = kk;

                        pt[0] -= 0.5*dx;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Ey_->point(*xx,*yy) = kk;

                        pt[pt_i] = ((ni-ii)-(ni-1)/2.0)*di;
                        pt[pt_j] = (jj-(nj-1)/2.0)*dj;
                        pt[0] += 0.5*dx;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Ex_end_->point(*xx,*yy) = kk;
                        pt[1] += 0.5*dy;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Hz_end_->point(*xx,*yy) = kk;

                        pt[0] -= 0.5*dx;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Ey_end_->point(*xx,*yy) = kk;
                    }
                }
            }
            else
            {
                for(ii =1; ii < thickness_+1; ii ++)
                {
                    for(jj = 1; jj < nj+1; jj++)
                    {
                        pt[0] = ((*xx)+0.5-(nx-1)/2.0)*dx;
                        pt[1] = ((*yy)-(ny-1)/2.0)*dy;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Hy_->point(*xx,*yy) = kk;

                        pt[0] -= 0.5*dx;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Ez_->point(*xx,*yy) = kk;

                        pt[1]  += 0.5*dy;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Hx_->point(*xx,*yy) = kk;

                        pt[pt_i] = ((ni-ii)-(ni-1)/2.0)*di;
                        pt[pt_j] = (jj-(nj-1)/2.0)*dj;
                        pt[0] += 0.5*dx;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Hy_end_->point(*xx,*yy) = kk;

                        pt[0] -= 0.5*dx;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Ez_end_->point(*xx,*yy) = kk;

                        pt[1]  += 0.5*dy;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Hx_end_->point(*xx,*yy) = kk;
                    }
                }
            }
        }
        if(c_hx_0_0_)
        {
            for(ii = 0; ii < thickness_; ii++)
            {
                for(int kk = 0; kk < oppPML; kk++)
                {
                    jj = kk;
                    int hx = 0; int hy = 0;
                    //Update Hx factors nj_0 side
                    eps    = objArr[phys_Hx_->point(*xx+1,*yy+1)].dielectric();
                    sigxx  = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigyx  = (ypml->*sigmay)(static_cast<double>(*yy) + 0.5,eps);
                    c_hx_0_0_->at(*xx).at(*yy) = calcHPreConsts(sigxx, sigyx, sigz);
                    hx++;

                    eps    = objArr[phys_Hx_end_->point(*xx+1,*yy+1)].dielectric();
                    sigxx  = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigyx  = (ypml->*sigmay)(static_cast<double>(*yy) + pow(-1,floor(hx/(dely+1))) * 0.5,eps); //Switch to see if PML is on teh different side yet
                    c_hx_n_0_->at(*xx).at(*yy) = calcHPreConsts(sigxx, sigyx, sigz);
                    hx++;

                    //Update Hy factors nj_0 side
                    eps    = objArr[phys_Hy_->point(*xx+1,*yy+1)].dielectric();
                    sigxy = (xpml->*sigmax)(static_cast<double>(*xx) + 0.5,eps);
                    sigyy = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    c_hy_0_0_->at(*xx).at(*yy) = calcHPreConsts(sigyy, sigz, sigxy);
                    hy++;

                    eps    = objArr[phys_Hy_end_->point(*xx+1,*yy+1)].dielectric();
                    sigxy = (xpml->*sigmax)(static_cast<double>(*xx) + pow(-1,floor(hy/(delx+1))) * 0.5,eps); //Switch to see if PML is on teh different side yet
                    sigyy = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    c_hy_n_0_->at(*xx).at(*yy) = calcHPreConsts(sigyy, sigz, sigxy);
                    hy++;

                    //Update Ez factors nj_0 side
                    eps = objArr[phys_Ez_->point(*xx+1,*yy+1)].dielectric();
                    sigx = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigy = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    c_ez_0_0_->at(*xx).at(*yy) = calcEPreConsts(eps,sigz, sigx, sigy);

                    eps = objArr[phys_Ez_end_->point(*xx+1,*yy+1)].dielectric();
                    sigx = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigy = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    c_ez_n_0_->at(*xx).at(*yy) = calcEPreConsts(eps,sigz, sigx, sigy);

                    jj = nj -1 - kk;
                    //Update Hx factors nj_n side
                    eps    = objArr[phys_Hx_->point(*xx+1,*yy+1)].dielectric();
                    jj =kk;
                    sigxx  = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigyx  = (ypml->*sigmay)(static_cast<double>(*yy) + pow(-1,floor(hx/(dely+1))) * 0.5,eps); //Switch to see if PML is on teh different side yet
                    c_hx_0_n_->at(*xx).at(*yy) = calcHPreConsts(sigxx, sigyx, sigz);
                    hx++;

                    jj = nj -1 - kk;
                    eps    = objArr[phys_Hx_end_->point(*xx+1,*yy+1)].dielectric();
                    jj =kk;
                    sigxx  = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigyx  = (ypml->*sigmay)(static_cast<double>(*yy) - 0.5,eps);
                    c_hx_n_n_->at(*xx).at(*yy) = calcHPreConsts(sigxx, sigyx, sigz);
                    hx++;

                    //Update Hy factors nj_n side
                    jj = nj -1 - kk;
                    eps    = objArr[phys_Hy_->point(*xx+1,*yy+1)].dielectric();
                    jj =kk;
                    sigxy = (xpml->*sigmax)(static_cast<double>(*xx) + pow(-1,floor(hy/(delx+1))) * 0.5,eps); //Switch to see if PML is on teh different side yet
                    sigyy = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    c_hy_0_n_->at(*xx).at(*yy) = calcHPreConsts(sigyy, sigz, sigxy);
                    hy++;

                    jj = nj -1 - kk;
                    eps    = objArr[phys_Hy_end_->point(*xx+1,*yy+1)].dielectric();
                    jj =kk;
                    sigxy = (xpml->*sigmax)(static_cast<double>(*xx) - 0.5,eps);
                    sigyy = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    c_hy_n_n_->at(*xx).at(*yy) = calcHPreConsts(sigyy, sigz, sigxy);
                    hy++;

                    //Update Ez factors nj_n side
                    jj = nj -1 - kk;
                    eps = objArr[phys_Ez_->point(*xx+1,*yy+1)].dielectric();
                    jj =kk;
                    sigx = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigy = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    c_ez_0_n_->at(*xx).at(*yy) = calcEPreConsts(eps,sigz, sigx, sigy);

                    jj = nj -1 - kk;
                    eps = objArr[phys_Ez_end_->point(*xx+1,*yy+1)].dielectric();
                    jj =kk;
                    sigx = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigy = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    c_ez_n_n_->at(*xx).at(*yy) = calcEPreConsts(eps,sigz, sigx, sigy);
                }
            }
        }
        else if(c_ex_0_0_)
        {
            for(ii = 0; ii < thickness_; ii++)
            {
                for(int kk = 0; kk < oppPML; kk++)
                {
                    jj = kk;
                    int hx = 0; int hy = 0;
                    //Update Ex factors nj_0 side
                    eps    = objArr[phys_Ex_->point(*xx+1,*yy+1)].dielectric();
                    sigxx  = (xpml->*sigmax)(static_cast<double>(*xx) + 0.5,eps);
                    sigyx  = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    c_ex_0_0_->at(*xx).at(*yy) = calcEPreConsts(eps,sigxx, sigyx, sigz);

                    eps    = objArr[phys_Ex_end_->point(*xx+1,*yy+1)].dielectric();
                    sigxx  = (xpml->*sigmax)(static_cast<double>(*xx) + pow(-1,floor(1/(delx+1))) * 0.5,eps);
                    sigyx  = (ypml->*sigmay)(static_cast<double>(*yy),eps); //Switch to see if PML is on teh different side yet
                    c_ex_n_0_->at(*xx).at(*yy) = calcEPreConsts(eps,sigxx, sigyx, sigz);

                    //Update Ey factors nj_0 side
                    eps    = objArr[phys_Ey_->point(*xx+1,*yy+1)].dielectric();
                    sigxy = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigyy = (ypml->*sigmay)(static_cast<double>(*yy) + 0.5,eps);
                    c_ey_0_0_->at(*xx).at(*yy) = calcEPreConsts(eps,sigyy, sigz, sigxy);

                    eps    = objArr[phys_Ey_end_->point(*xx+1,*yy+1)].dielectric();
                    sigxy = (xpml->*sigmax)(static_cast<double>(*xx),eps); //Switch to see if PML is on teh different side yet
                    sigyy = (ypml->*sigmay)(static_cast<double>(*yy) + pow(-1,floor(1/(dely+1))) * 0.5,eps);
                    c_ey_n_0_->at(*xx).at(*yy) = calcEPreConsts(eps,sigyy, sigz, sigxy);

                    //Update Hz factors nj_0 side
                    eps = objArr[phys_Hz_->point(*xx+1,*yy+1)].dielectric();
                    sigx = (xpml->*sigmax)(static_cast<double>(*xx) + 0.5,eps);
                    sigy = (ypml->*sigmay)(static_cast<double>(*yy) + 0.5,eps);
                    c_hz_0_0_->at(*xx).at(*yy) = calcHPreConsts(sigz, sigx, sigy);

                    eps = objArr[phys_Hz_end_->point(*xx+1,*yy+1)].dielectric();
                    sigx = (xpml->*sigmax)(static_cast<double>(*xx)+pow(-1,floor(1/(delx+1))) * 0.5,eps);
                    sigy = (ypml->*sigmay)(static_cast<double>(*yy)+pow(-1,floor(1/(dely+1))) * 0.5,eps);
                    c_hz_n_0_->at(*xx).at(*yy) = calcHPreConsts(sigz, sigx, sigy);

                    jj = nj -1 - kk;
                    //Update Ex factors nj_n side
                    eps    = objArr[phys_Ex_->point(*xx+1,*yy+1)].dielectric();
                    jj =kk;
                    sigxx  = (xpml->*sigmax)(static_cast<double>(*xx) + pow(-1,floor(2/(delx+1))) * 0.5,eps);
                    sigyx  = (ypml->*sigmay)(static_cast<double>(*yy),eps); //Switch to see if PML is on teh different side yet
                    c_ex_0_n_->at(*xx).at(*yy) = calcEPreConsts(eps,sigxx, sigyx, sigz);

                    jj = nj -1 - kk;
                    eps    = objArr[phys_Ex_end_->point(*xx+1,*yy+1)].dielectric();
                    jj =kk;
                    sigxx  = (xpml->*sigmax)(static_cast<double>(*xx)-0.5,eps);
                    sigyx  = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    c_ex_n_n_->at(*xx).at(*yy) = calcEPreConsts(eps,sigxx, sigyx, sigz);

                    //Update Ey factors nj_n side
                    jj = nj -1 - kk;
                    eps    = objArr[phys_Ey_->point(*xx+1,*yy+1)].dielectric();
                    jj =kk;
                    sigxy = (xpml->*sigmax)(static_cast<double>(*xx),eps); //Switch to see if PML is on teh different side yet
                    sigyy = (ypml->*sigmay)(static_cast<double>(*yy) + pow(-1,floor(2/(dely+1))) * 0.5,eps);
                    c_ey_0_n_->at(*xx).at(*yy) = calcEPreConsts(eps,sigyy, sigz, sigxy);

                    jj = nj -1 - kk;
                    eps    = objArr[phys_Ey_end_->point(*xx+1,*yy+1)].dielectric();
                    jj =kk;
                    sigxy = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigyy = (ypml->*sigmay)(static_cast<double>(*yy)-0.5,eps);
                    c_ey_n_n_->at(*xx).at(*yy) = calcEPreConsts(eps,sigyy, sigz, sigxy);

                    //Update Hz factors nj_n side
                    jj = nj -1 - kk;
                    eps = objArr[phys_Hz_->point(*xx+1,*yy+1)].dielectric();
                    jj =kk;
                    sigx = (xpml->*sigmax)(static_cast<double>(*xx)+ pow(-1,floor(2/(delx+1))) * 0.5,eps);
                    sigy = (ypml->*sigmay)(static_cast<double>(*yy)+ pow(-1,floor(2/(dely+1))) * 0.5,eps);
                    c_hz_0_n_->at(*xx).at(*yy) = calcHPreConsts(sigz, sigx, sigy);

                    jj = nj -1 - kk;
                    eps = objArr[phys_Hz_end_->point(*xx+1,*yy+1)].dielectric();
                    jj =kk;
                    sigx = (xpml->*sigmax)(static_cast<double>(*xx)-0.5,eps);
                    sigy = (ypml->*sigmay)(static_cast<double>(*yy)-0.5,eps);
                    c_hz_n_n_->at(*xx).at(*yy) = calcHPreConsts(sigz, sigx, sigy);
                }
            }
        }
        if(pol_ == HZ || pol_ == EX || pol_ == EY)
        {
            for(ii= 0; ii < thickness_; ii++)
            {
                jj = zaxXJmax;
                while(jj > oppPML-1)
                {
                    int jjstore = jj;
                    while(jj > oppPML && phys_Ex_ -> point(*xx+1,*yy+1) == phys_Ex_ -> point(*xx+1-delx,*yy+1-dely)) //Fix
                        jj--;
                    std::array<double,9> tempArr = {static_cast<double>(*xx+1),static_cast<double>(*yy+1),static_cast<double>(jjstore - jj + 1),static_cast<double>(phys_Ex_->point(*xx+1,*yy+1))};
                    eps   = objArr[phys_Ex_->point(*xx+1,*yy+1)].dielectric();
                    sigxx = (xpml->*sigmax)(static_cast<double>(*xx) + 0.5,eps);
                    sigyx = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    std::array<double,5> preconsts = calcEPreConsts(eps,sigxx, sigyx, sigz);
                    std::copy_n(preconsts.begin(),5,tempArr.begin()+4);
                    zaxEx_.push_back(tempArr);
                    jj--;
                }
            }
            for(ii= delx; ii < thickness_; ii++)
            {
                jj = zaxXJmax;
                while(jj > oppPML-1)
                {
                    int jjstore = jj;
                    while(jj > oppPML && phys_Ex_end_ -> point(*xx+1,*yy+1) == phys_Ex_end_ -> point(*xx+1-delx,*yy+1-dely)) //Fix
                        jj--;
                    std::array<double,9> tempArr = {static_cast<double>(*xx+1),static_cast<double>(*yy+1),static_cast<double>(jjstore - jj + 1),static_cast<double>(phys_Ex_end_->point(*xx+1,*yy+1))};
                    eps   = objArr[phys_Ex_end_->point(*xx+1,*yy+1)].dielectric();
                    sigxx = (xpml->*sigmax)(static_cast<double>(*xx) + pow(-1,dely) * 0.5,eps); //Switch to see if PML is on teh different side yet
                    sigyx = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    std::array<double,5> preconsts = calcEPreConsts(eps,sigxx, sigyx, sigz);
                    std::copy_n(preconsts.begin(),5,tempArr.begin()+4);
                    zaxEx_end_.push_back(tempArr);
                    jj--;
                }
            }
            for(ii= 0; ii < thickness_; ii++)
            {
                jj = zaxYJmax;
                while(jj > oppPML-1)
                {
                    int jjstore = jj;
                    while(jj > oppPML && phys_Ey_ -> point(*xx+1,*yy+1) == phys_Ey_ -> point(*xx+1-delx,*yy+1-dely))
                        jj--;
                    std::array<double,9> tempArr = {static_cast<double>(*xx+1),static_cast<double>(*yy+1),static_cast<double>(jjstore - jj + 1),static_cast<double>(phys_Ey_->point(*xx+1,*yy+1))};
                    eps    = objArr[phys_Ey_->point(*xx+1,*yy+1)].dielectric();
                    sigxy = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigyy = (ypml->*sigmay)(static_cast<double>(*yy) + 0.5,eps);
                    std::array<double,5> preconsts = calcEPreConsts(eps,sigyy, sigz, sigxy);
                    std::copy_n(preconsts.begin(),5,tempArr.begin()+4);
                    zaxEy_.push_back(tempArr);
                    jj--;
                }
            }
            for(ii= dely; ii < thickness_; ii++)
            {
                jj = zaxYJmax;
                while(jj > oppPML-1)
                {
                    int jjstore = jj;
                    while(jj > oppPML && phys_Ey_end_ -> point(*xx+1,*yy+1) == phys_Ey_end_ -> point(*xx+1-delx,*yy+1-dely)) //Fix
                        jj--;
                    std::array<double,9> tempArr = {static_cast<double>(*xx+1),static_cast<double>(*yy+1),static_cast<double>(jjstore - jj + 1),static_cast<double>(phys_Ey_end_->point(*xx+1,*yy+1))};
                    eps    = objArr[phys_Ey_end_->point(*xx+1,*yy+1)].dielectric();
                    sigxy = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigyy = (ypml->*sigmay)(static_cast<double>(*yy) + pow(-1,delx) * 0.5,eps);//Switch to see if PML is on teh different side yet
                    std::array<double,5> preconsts = calcEPreConsts(eps,sigyy, sigz, sigxy);
                    std::copy_n(preconsts.begin(),5,tempArr.begin()+4);
                    zaxEy_end_.push_back(tempArr);
                    jj--;
                }
            }
            for(ii= 0; ii < thickness_; ii++)
            {
                jj = zmax;
                while(jj > zmin-1)
                {
                    int jjstore = jj;
                    while(jj > zmin && phys_Hz_ -> point(*xx+1,*yy+1) == phys_Hz_ -> point(*xx+1-delx,*yy+1-dely))
                        jj--;
                    std::array<double,9> tempArr = {static_cast<double>(*xx+1),static_cast<double>(*yy+1),static_cast<double>(jjstore - jj + 1),static_cast<double>(phys_Hz_->point(*xx+1,*yy+1))};
                    eps = objArr[phys_Hz_->point(*xx+1,*yy+1)].dielectric();
                    sigx = (xpml->*sigmax)(static_cast<double>(*xx)+0.5,eps);
                    sigy = (ypml->*sigmay)(static_cast<double>(*yy)+0.5,eps);
                    std::array<double,5> preconsts = calcHPreConsts(sigz, sigx, sigy);
                    std::copy_n(preconsts.begin(),5,tempArr.begin()+4);
                    zaxHz_.push_back(tempArr);
                    jj--;
                }
            }
            for(ii= 0; ii < thickness_; ii++)
            {
                jj = zmax;
                while(jj > zmin-1)
                {
                    int jjstore = jj;
                    while(jj > zmin && phys_Hz_end_ -> point(*xx+1,*yy+1) == phys_Hz_end_ -> point(*xx+1-delx,*yy+1-dely)) //Fix
                        jj--;
                    std::array<double,9> tempArr = {static_cast<double>(*xx+1),static_cast<double>(*yy+1),static_cast<double>(jjstore - jj + 1),static_cast<double>(phys_Hz_end_->point(*xx+1,*yy+1))};
                    eps = objArr[phys_Hz_end_->point(*xx+1,*yy+1)].dielectric();
                    sigx = (xpml->*sigmax)(static_cast<double>(*xx) + pow(-1,dely) * 0.5,eps);
                    sigy = (ypml->*sigmay)(static_cast<double>(*yy) + pow(-1,delx) * 0.5,eps);
                    std::array<double,5> preconsts = calcHPreConsts(sigz, sigx, sigy);
                    std::copy_n(preconsts.begin(),5,tempArr.begin()+4);
                    zaxHz_end_.push_back(tempArr);
                    jj--;
                }
            }
        }
        else
        {
            for(ii= 0; ii < thickness_; ii++)
            {
                jj = zaxXJmax;
                while(jj > oppPML-1)
                {
                    int jjstore = jj;
                    while(jj > oppPML && phys_Hx_ -> point(*xx+1,*yy+1) == phys_Hx_ -> point(*xx+1-delx,*yy+1-dely)) //Fix
                        jj--;
                    std::array<double,9> tempArr = {static_cast<double>(*xx+1),static_cast<double>(*yy+1),static_cast<double>(jjstore - jj + 1),static_cast<double>(phys_Hx_->point(*xx+1,*yy+1))};
                    eps   = objArr[phys_Hx_->point(*xx+1,*yy+1)].dielectric();
                    sigxx = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigyx = (ypml->*sigmay)(static_cast<double>(*yy) + 0.5,eps); //Switch to see if PML is on teh different side yet
                    std::array<double,5> preconsts = calcHPreConsts(sigxx, sigyx, sigz);
                    std::copy_n(preconsts.begin(),5,tempArr.begin()+4);
                    zaxHx_.push_back(tempArr);
                    jj--;
                }
            }
            for(ii= delx; ii < thickness_; ii++)
            {
                jj = zaxXJmax;
                while(jj > oppPML-1)
                {
                    int jjstore = jj;
                    while(jj > oppPML && phys_Hx_end_ -> point(*xx+1,*yy+1) == phys_Hx_end_ -> point(*xx+1-delx,*yy+1-dely)) //Fix
                        jj--;
                    std::array<double,9> tempArr = {static_cast<double>(*xx+1),static_cast<double>(*yy+1),static_cast<double>(jjstore - jj + 1),static_cast<double>(phys_Hx_end_->point(*xx+1,*yy+1))};
                    eps   = objArr[phys_Hx_end_->point(*xx+1,*yy+1)].dielectric();
                    sigxx = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigyx = (ypml->*sigmay)(static_cast<double>(*yy) + pow(-1,delx) * 0.5,eps); //Switch to see if PML is on teh different side yet
                    std::array<double,5> preconsts = calcHPreConsts(sigxx, sigyx, sigz);
                    std::copy_n(preconsts.begin(),5,tempArr.begin()+4);
                    zaxHx_end_.push_back(tempArr);
                    jj--;
                }
            }
            for(ii= 0; ii < thickness_; ii++)
            {
                jj = zaxYJmax;
                while(jj > oppPML-1)
                {
                    int jjstore = jj;
                    while(jj > oppPML && phys_Hy_ -> point(*xx+1,*yy+1) == phys_Hy_ -> point(*xx+1-delx,*yy+1-dely))
                        jj--;
                    std::array<double,9> tempArr = {static_cast<double>(*xx+1),static_cast<double>(*yy+1),static_cast<double>(jjstore - jj + 1),static_cast<double>(phys_Hy_->point(*xx+1,*yy+1))};
                    eps    = objArr[phys_Hy_->point(*xx+1,*yy+1)].dielectric();
                    sigxy = (xpml->*sigmax)(static_cast<double>(*xx) + 0.5,eps);
                    sigyy = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    std::array<double,5> preconsts = calcHPreConsts(sigyy, sigz, sigxy);
                    std::copy_n(preconsts.begin(),5,tempArr.begin()+4);
                    // std::cout << tempArr[0] << "\t" << tempArr[1] << "\t"  << std::setw(10)  <<tempArr[4] << "\t" << std::setw(10) << tempArr[5] << "\t" << std::setw(10) << tempArr[6] << "\t" << std::setw(10) << tempArr[7] << "\t" << std::setw(10) << tempArr[8] << std::endl;
                    zaxHy_.push_back(tempArr);
                    jj--;
                }
            }
            for(ii= dely; ii < thickness_; ii++)
            {
                jj = zaxYJmax;
                while(jj > oppPML-1)
                {
                    int jjstore = jj;
                    while(jj > oppPML && phys_Hy_end_ -> point(*xx+1,*yy+1) == phys_Hy_end_ -> point(*xx+1-delx,*yy+1-dely)) //Fix
                        jj--;
                    std::array<double,9> tempArr = {static_cast<double>(*xx+1),static_cast<double>(*yy+1),static_cast<double>(jjstore - jj + 1),static_cast<double>(phys_Hy_end_->point(*xx+1,*yy+1))};
                    eps    = objArr[phys_Hy_end_->point(*xx+1,*yy+1)].dielectric();
                    sigxy = (xpml->*sigmax)(static_cast<double>(*xx) + pow(-1,dely) * 0.5,eps); //Switch to see if PML is on teh different side yet
                    sigyy = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    std::array<double,5> preconsts = calcHPreConsts(sigyy, sigz, sigxy);
                    std::copy_n(preconsts.begin(),5,tempArr.begin()+4);
                    zaxHy_end_.push_back(tempArr);
                    jj--;
                }
            }
            for(ii= 0; ii < thickness_; ii++)
            {
                jj = zmax;
                while(jj > zmin-1)
                {
                    int jjstore = jj;
                    while(jj > zmin && phys_Ez_ -> point(*xx+1,*yy+1) == phys_Ez_ -> point(*xx+1-delx,*yy+1-dely))
                        jj--;
                    std::array<double,9> tempArr = {static_cast<double>(*xx+1),static_cast<double>(*yy+1),static_cast<double>(jjstore - jj + 1),static_cast<double>(phys_Ez_->point(*xx+1,*yy+1))};
                    eps = objArr[phys_Ez_->point(*xx+1,*yy+1)].dielectric();
                    sigx = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigy = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    std::array<double,5> preconsts = calcEPreConsts(eps,sigz, sigx, sigy);
                    std::copy_n(preconsts.begin(),5,tempArr.begin()+4);
                    zaxEz_.push_back(tempArr);
                    jj--;
                }
            }
            for(ii= 0; ii < thickness_; ii++)
            {
                jj = zmax;
                while(jj > zmin-1)
                {
                    int jjstore = jj;
                    while(jj > zmin && phys_Ez_end_ -> point(*xx+1,*yy+1) == phys_Ez_end_ -> point(*xx+1-delx,*yy+1-dely)) //Fix
                        jj--;
                    std::array<double,9> tempArr = {static_cast<double>(*xx+1),static_cast<double>(*yy+1),static_cast<double>(jjstore - jj + 1),static_cast<double>(phys_Ez_end_->point(*xx+1,*yy+1))};
                    eps = objArr[phys_Ez_end_->point(*xx+1,*yy+1)].dielectric();
                    sigx = (xpml->*sigmax)(static_cast<double>(*xx),eps);
                    sigy = (ypml->*sigmay)(static_cast<double>(*yy),eps);
                    std::array<double,5> preconsts = calcEPreConsts(eps,sigz, sigx, sigy);
                    std::copy_n(preconsts.begin(),5,tempArr.begin()+4);
                    zaxEz_end_.push_back(tempArr);
                    jj--;
                }
            }
        }
    }

    // Accessor Functions
    /**
     * @brief Reuturns the thickness of the PML
     */
    int thickness(){return thickness_;}
    /**
     * @brief Returns the direction of the PML
     */
    Direction d(){return d_;}
    /**
     * @brief Returns the alue of Kappa given a position
     * @details Will Return the value of Kappa at a given position once implimented, at the moment it's all zero
     *
     * @param x unit cell location of the point in question
     * @return [description]
     */
    double kappa(int x)
    {
        if(x <= thickness_-1 && x >=-0.25)
            return 1 + (kappaMax_ - 1) * pow((static_cast<double>(thickness_-1) - x) / static_cast<double>(thickness_-1) , m_);
        else if(ni_ -x <= thickness_-1  && x<ni_)
            return 1 + (kappaMax_ - 1) * pow((static_cast<double>(thickness_-1) - (ni_-1-x)) / static_cast<double>(thickness_-1) , m_) ;
        else
            return 0.0;
    }
    /**
     * @brief Returns the alue of sigma given a position
     * @details Will Return the value of sigma at a given position once implimented, at the moment it's all zero
     *
     * @param x unit cell location of the point in question
     * @return [description]
     */
    double sigma(double x,double eps)
    {
        if(x <= thickness_ && x >=-0.25)
            return sigmaMax_ * pow((static_cast<double>(thickness_-1) - x) / static_cast<double>(thickness_-1) , m_);// / sqrt(eps);
        else if(ni_ -x <= thickness_  && x<ni_)
            return sigmaMax_ * pow((static_cast<double>(thickness_-1) - (ni_-1-x)) / static_cast<double>(thickness_-1) , m_);// / sqrt(eps);
        else
            return 0.0;
    }
    /**
     * @brief lambda function that returns 0
     * @details FUnction with the same paramters as teh sigma fuction that returns 0 (used for oppisite pml clacs)
     *
     * @param x position
     * @param eps dielectric
     *
     * @return 0
     */
    double sigmaopp(double x,double eps){return 0.0;} // can't make lamda function
    /**
     * @brief Gets the member function pointer to sigma
     * @details  Gets the member function pointer to sigma
     * @return the member function pointer to sigma
     */
    PMLMemFn sig_ptr() {return &UPML<T>::sigma;}

    /**
     * @brief Finds the PML Preconstants
     * @details determines the H preconstants
     *
     * @param eps dielectric constant
     * @param sigi sigma value along the direction of the PML
     * @param sigj sigma value along next in the cylic version
     * @param sigk sigma value in the cyclic direction
     * @return H preconstants
     */
    std::array<double,5> calcHPreConsts(double sigi, double sigj, double sigk)
    {
        std::array<double,5> preFact = {0.0,0.0,0.0,0.0,0.0};
        double kapi = 1.0; double kapj = 1.0; double kapk = 1.0;
        preFact[0] = (2*kapj - sigj*dt_) / (2*kapj + sigj*dt_);
        preFact[1] = (2* dt_) / (dy_ * (2*kapj + sigj*dt_));
        preFact[2] = (2*kapk - sigk*dt_) / (2*kapk + sigk*dt_);
        preFact[3] = (2*kapi + sigi*dt_) / (2*kapk + sigk*dt_);
        preFact[4] = (2*kapi - sigi*dt_) / (2*kapk + sigk*dt_);
        return preFact;
    }
    /**
     * @brief Finds the PML Preconstants
     * @details determines the E preconstants
     *
     * @param eps dielectric constant
     * @param sigi sigma value along the direction of the PML
     * @param sigj sigma value along next in the cylic version
     * @param sigk sigma value in the cyclic direction
     * @return E preconstants
     */
    std::array<double,5> calcEPreConsts(double eps, double sigi, double sigj, double sigk)
    {
        std::array<double,5> preFact = {0.0,0.0,0.0,0.0,0.0};
        double kapi = 1.0; double kapj = 1.0; double kapk = 1.0;
        preFact[0] = (2*kapj - sigj*dt_) / (2*kapj + sigj*dt_);
        preFact[1] = (2* dt_) / (dy_ * (2*kapj + sigj*dt_));
        preFact[2] = (2*kapk - sigk*dt_) / (2*kapk + sigk*dt_);
        preFact[3] = (2*kapi + sigi*dt_) / (2*kapk + sigk*dt_) / eps;
        preFact[4] = (2*kapi - sigi*dt_) / (2*kapk + sigk*dt_) / eps;
        return preFact;
    }
};

#endif