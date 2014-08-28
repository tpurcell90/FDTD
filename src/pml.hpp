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
protected:
    int thickness_;
    Direction d_;
    //PMLTopBot tb_;
    double m_;
    double R0_;
    double sigmaMax_;
    double kappaMax_;
    bool precalc_;
    Polarization pol_;
    int nj_,ni_;
    double dx_,dy_,dt_;

public:
    std::shared_ptr<Grid2D<T>> Dx_,Dy,Dz_,Bx_,By_,Bz_,Dx_end_,Dyend_,Dz_end_,Bx_end_,By_end_,Bz_end_;
    std::shared_ptr<Grid2D<int>> phys_Hx_,phys_Hy_,phys_Hz_, phys_Hx_end_,phys_Hy_end_,phys_Hz_end_,phys_Ex_,phys_Ey_,phys_Ez_, phys_Ex_end_,phys_Ey_end_,phys_Ez_end_;

    std::shared_ptr<std::vector<std::vector<std::array<double,5>>>> c_hx_0_, c_hy_0_, c_ez_0_, c_hx_n_, c_hy_n_, c_ez_n_;
    std::shared_ptr<std::vector<std::vector<std::array<double,5>>>> c_ex_0_, c_ey_0_, c_hz_0_, c_ex_n_, c_ey_n_, c_hz_n_;

    /**
     * @brief Constrcuts a PML for both ends of the cell
     * @details Uses the input to construct the functions and auxiliary fields for the PML calculations
     *
     * @param thickness Thickness of PML in number of unit cells
     * @param d Direction the PML is occupying
     * @param m Order of the polynomial grading for the kappa and sigma vectors
     * @param R0 Max allowed reflection for normally incident waves
     * @param nx Number of unit cells in the normal direction of the PML
     * @param dx unit cell spacing for the x direction
     * @param dy unit cell spacing for the y direction
     * @param pol A polarization so it can set up the right auxilliary fields
     */
    UPML(int thickness, Direction d, double m, double R0, int nx, int ny, double dx, double dy, double dt, Polarization pol, bool precalc) : thickness_(thickness), d_(d), m_(m), R0_(R0), dx_(dx), dy_(dy), dt_(dt), precalc_(precalc), pol_(pol)
    {
        sigmaMax_ = -(m_+1)*log(R0_)/(2*thickness_*dx); // eta should be included;
        kappaMax_ = 1.0;
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
        if(pol == EX || pol == EY || pol == HZ)
        {
            Dx_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dy = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Bz_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dx_end_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dyend_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Bz_end_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);

            Bx_ = nullptr;
            By_ = nullptr;
            Dz_ = nullptr;
            Bx_end_ = nullptr;
            By_end_ = nullptr;
            Dz_end_ = nullptr;

            phys_Hz_ = std::make_shared<Grid2D<int>>(xmax,ymax,dx,dy);
            phys_Hz_end_ = std::make_shared<Grid2D<int>>(xmax,ymax,dx,dy);
            phys_Ex_ = std::make_shared<Grid2D<int>>(xmax,ymax,dx,dy);
            phys_Ex_end_ = std::make_shared<Grid2D<int>>(xmax,ymax,dx,dy);
            phys_Ey_ = std::make_shared<Grid2D<int>>(xmax,ymax,dx,dy);
            phys_Ey_end_ = std::make_shared<Grid2D<int>>(xmax,ymax,dx,dy);

            phys_Hy_ = nullptr;
            phys_Hy_end_ = nullptr;
            phys_Hx_ = nullptr;
            phys_Hx_end_ = nullptr;
            phys_Ez_ = nullptr;
            phys_Ez_end_ = nullptr;

            if(precalc_ == false)
            {
                c_hx_0_ = nullptr; c_hy_0_ = nullptr; c_ez_0_ = nullptr; c_hx_n_ = nullptr; c_hy_n_ = nullptr; c_ez_n_= nullptr;
                c_ex_0_ = nullptr; c_ey_0_ = nullptr; c_hz_0_ = nullptr; c_ex_n_ = nullptr; c_ey_n_ = nullptr; c_hz_n_ = nullptr;
            }
            else
            {
                c_hx_0_ = nullptr; c_hy_0_ = nullptr; c_ez_0_ = nullptr; c_hx_n_ = nullptr; c_hy_n_ = nullptr; c_ez_n_= nullptr;
                c_ex_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xmax, std::vector<std::array<double,5>>(ymax)); 
                c_ey_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xmax, std::vector<std::array<double,5>>(ymax));
                c_hz_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xmax, std::vector<std::array<double,5>>(ymax));
                c_ex_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xmax, std::vector<std::array<double,5>>(ymax));
                c_ey_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xmax, std::vector<std::array<double,5>>(ymax));
                c_hz_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xmax, std::vector<std::array<double,5>>(ymax));
            }
        }
        else
        {
            Bx_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            By_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dz_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dx_ = nullptr;
            Dy = nullptr;
            Bz_ = nullptr;

            Bx_end_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            By_end_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dz_end_ = std::make_shared<Grid2D<T>>(xmax,ymax,dx,dy);
            Dx_end_ = nullptr;
            Dyend_ = nullptr;
            Bz_end_ = nullptr;

            phys_Hz_ = nullptr;
            phys_Hz_end_ = nullptr;
            phys_Hy_ = std::make_shared<Grid2D<int>>(xmax,ymax,dx,dy);
            phys_Hy_end_ = std::make_shared<Grid2D<int>>(xmax,ymax,dx,dy);
            phys_Hx_ = std::make_shared<Grid2D<int>>(xmax,ymax,dx,dy);
            phys_Hx_end_ = std::make_shared<Grid2D<int>>(xmax,ymax,dx,dy);

            phys_Ez_ = std::make_shared<Grid2D<int>>(xmax,ymax,dx,dy);
            phys_Ez_end_ = std::make_shared<Grid2D<int>>(xmax,ymax,dx,dy);
            phys_Ex_ = nullptr;
            phys_Ex_end_ = nullptr;
            phys_Ey_ = nullptr;
            phys_Ey_end_ = nullptr;
           
            if(precalc_ == false)
            {
                c_hx_0_ = nullptr; c_hy_0_ = nullptr; c_ez_0_ = nullptr; c_hx_n_ = nullptr; c_hy_n_ = nullptr; c_ez_n_= nullptr;
                c_ex_0_ = nullptr; c_ey_0_ = nullptr; c_hz_0_ = nullptr; c_ex_n_ = nullptr; c_ey_n_ = nullptr; c_hz_n_ = nullptr;
            }
            else
            {
                c_ex_0_ = nullptr; c_ey_0_ = nullptr; c_hz_0_ = nullptr; c_ex_n_ = nullptr; c_ey_n_ = nullptr; c_hz_n_= nullptr;
                c_hx_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xmax, std::vector<std::array<double,5>>(ymax)); 
                c_hy_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xmax, std::vector<std::array<double,5>>(ymax));
                c_ez_0_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xmax, std::vector<std::array<double,5>>(ymax));
                c_hx_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xmax, std::vector<std::array<double,5>>(ymax));
                c_hy_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xmax, std::vector<std::array<double,5>>(ymax));
                c_ez_n_ = std::make_shared<std::vector<std::vector<std::array<double,5>>>>(xmax, std::vector<std::array<double,5>>(ymax));
            }
        }
    }

    void initializeUPML(std::vector<Obj> objArr, int nx, int ny, double dx, double dy, double dt, int yPML, int xPML, UPML<T> *opp)
    {
        //[](double x, double eps){return 0.0;}
        int jmax = 0;
        int pt_i = 0;
        int pt_j = 0;
        int ni   = 0;
        int nj   = 0;
        double di = 0.0;
        double dj = 0.0;
        int xmax; int ymax;
        if (d_ == X)
        {
            jmax = ny-1;
            pt_i = 0;
            pt_j = 1;
            ni   = nx;
            nj   = ny;
            di   = dx;
            dj   = dy;
            xmax = thickness_;
            ymax = ny;
        }
        else
        {
            jmax = nx-1;
            pt_i = 1;
            pt_j = 0;
            ni   = ny;
            nj   = nx;
            di   = dy;
            dj   = dx;
            xmax = nx;
            ymax = thickness_;
        }
        std::cout << 1 << std::endl;
        for(int kk = 0; kk < objArr.size(); kk++)
        {
            std::vector<double> pt(2,0.0);
            if(pol_ == HZ || pol_ == EX || pol_ == EY)
            {
                
            }
            else
            {
                if(objArr[kk].s() == sphere)
                {
                    
                }
                else if(objArr[kk].s() == block)
                {
                    
                    for(int ii =0; ii < xmax; ii ++)
                    {
                        for(int jj = 0; jj < ymax; jj++)
                        {
                            pt[pt_i] = (ii+0.5-(ni-1)/2.0)*di;
                            pt[pt_j]  = (jj-(nj-1)/2.0)*dj;
                            if(objArr[kk].isObj(pt)==true)
                                phys_Hy_->point(ii,jj) = kk;
                            pt[pt_i] -= 0.5*di;
                            if(objArr[kk].isObj(pt)==true)
                                phys_Ez_->point(ii,jj) = kk;
                            pt[pt_j]  += 0.5*dj;
                            if(objArr[kk].isObj(pt)==true && kk ==1)
                                phys_Hx_->point(ii,jj) = kk;
                            pt[pt_i] = ((ni-1-ii)+0.5-(ni-1)/2.0)*di;
                            pt[pt_j]  = (jj-(nj-1)/2.0)*dj;
                            if(objArr[kk].isObj(pt)==true)
                                phys_Hy_end_->point(ii,jj) = kk;
                            pt[pt_i] -= 0.5*di;
                            if(objArr[kk].isObj(pt)==true)
                                phys_Ez_end_->point(ii,jj) = kk;
                            pt[pt_j]  += 0.5*dj;
                            if(objArr[kk].isObj(pt)==true && kk ==1)
                                phys_Hx_end_->point(ii,jj) = kk;
                        }
                        pt[pt_i]=(ii+0.5-(ni-1)/2.0)*di;
                        pt[pt_j] =((nj-1)/2.0)*dj;

                        // if(objArr[kk].isObj(pt)==true)
                        //     phys_Hy_->point(ii,nj-1) = kk;
                        // pt[pt_i] -= 0.5*di;
                        // if(objArr[kk].isObj(pt)==true)
                        //     phys_Ez_->point(ii,nj-1) = kk;

                        // pt[pt_i]=((ni-1-ii)+0.5-(ni-1)/2.0)*di;
                        // pt[pt_j] =((nj-1)/2.0)*dj;
                        
                        // if(objArr[kk].isObj(pt)==true)
                        //     phys_Hy_end_->point(ii,nj-1) = kk;
                        // pt[pt_i] -= 0.5*di;
                        // if(objArr[kk].isObj(pt)==true)
                        //     phys_Ez_end_->point(ii,nj-1) = kk;
                    }
                }
            }
        }
        std::cout << 2 << std::endl;
        if(precalc_)
        {
            double eps=0.0;
            double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
           
            double sigz = 0.0; double sigx = 0.0; double sigy = 0.0;
            double sigxx = 0.0; double sigxy = 0.0; double sigyx = 0.0; double sigyy = 0.0;
            int xmax = 0; int ymax = 0;

            PMLMemFn  sigmay;
            PMLMemFn  sigmax;

            UPML<T> *xpml; UPML<T> *ypml; 
            if(d_==X)
            {
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
                xmax = thickness_;
                ymax = nj_;
            }
            else
            {
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
                ymax = thickness_;
                xmax = nj_;
            }
            if(pol_ == EZ || pol_ == HX || pol_ == HY)
            {
                for(int ii = 0; ii < xmax; ii++)
                {
                    for(int jj = 0; jj < ymax; jj++)
                    {
                        //Update Hx factors
                        eps    = objArr[phys_Hx_->point(ii,jj)].dielectric(1.0);
                        sigxx  = (xpml->*sigmax)(static_cast<double>(ii),eps);
                        sigyx  = (ypml->*sigmay)(static_cast<double>(jj) + 0.5,eps);
                        c_hx_0_->at(ii).at(jj) = calcPreConsts(eps,sigxx, sigyx, sigz);
                        
                        eps    = objArr[phys_Hx_end_->point(ii,jj)].dielectric(1.0);
                        sigxx  = (xpml->*sigmax)(static_cast<double>(ii),eps);
                        sigyx  = (ypml->*sigmay)(static_cast<double>(jj)-0.5,eps);
                        c_hx_n_->at(ii).at(jj) = calcPreConsts(eps,sigxx, sigyx, sigz);

                        //Update Hy factors
                        eps    = objArr[phys_Hy_->point(ii,jj)].dielectric(1.0);
                        sigxy = (xpml->*sigmax)(static_cast<double>(ii) + 0.5,eps);
                        sigyy = (ypml->*sigmay)(static_cast<double>(jj),eps);
                        c_hy_0_->at(ii).at(jj) = calcPreConsts(eps,sigyy, sigz, sigxy);

                        eps    = objArr[phys_Hy_end_->point(ii,jj)].dielectric(1.0);
                        sigxy = (xpml->*sigmax)(static_cast<double>(ii) - 0.5,eps);
                        sigyy = (ypml->*sigmay)(static_cast<double>(jj),eps);
                        c_hy_n_->at(ii).at(jj) = calcPreConsts(eps,sigyy, sigz, sigxy);

                        //Update Ez factors
                        eps = objArr[phys_Ez_->point(ii,jj)].dielectric(1.0);
                        sigx = (xpml->*sigmax)(static_cast<double>(ii),eps);
                        sigy = (ypml->*sigmay)(static_cast<double>(jj),eps);
                        c_ez_0_->at(ii).at(jj) = calcPreConsts(eps,sigz, sigx, sigy);
                        
                        eps = objArr[phys_Ez_end_->point(ii,jj)].dielectric(1.0);
                        sigx = (xpml->*sigmax)(static_cast<double>(ii),eps);
                        sigy = (ypml->*sigmay)(static_cast<double>(jj),eps);
                        c_ez_n_->at(ii).at(jj) = calcPreConsts(eps,sigz, sigx, sigy);
                    }
                }
            }
        }
    }

    // void updateH(FDTDField *FF)
    // {

    // }

    // void updateE(FDTDField *FF)
    // {

    // }
    // Accessor Functions
    /**
     * @brief Reuturns the thickness of the PML
     */
    int thickness(){return thickness_;}
    /**
     * @brief Returns the direction of the PML
     */
    Direction d(){return d_;}
    bool precalc() {return precalc_;}
    /**
     * @brief Returns the alue of Kappa given a position
     * @details Will Return the value of Kappa at a given position once implimented, at the moment it's all zero
     *
     * @param x unit cell location of the point in question
     * @return [description]
     */
    double kappa(int x){return kappaMax_;}
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
            return sigmaMax_ * sqrt(eps) * pow((static_cast<double>(thickness_-1) - x) / static_cast<double>(thickness_-1) , m_);
        else if(ni_ -x <= thickness_  && x<ni_)
            return sigmaMax_ * sqrt(eps) * pow((static_cast<double>(thickness_-1) - (ni_-1-x)) / static_cast<double>(thickness_-1) , m_);
        else
            return 0.0;
    }

    double sigmaopp(double x,double eps){return 0.0;} // can't make lamda function
    PMLMemFn sig_ptr() {return &UPML<T>::sigma;}

    std::array<double,5> calcPreConsts(double eps, double sigi, double sigj, double sigk)
    {
        std::array<double,5> preFact = {0.0,0.0,0.0,0.0,0.0};
        double kapi = 1.0; double kapj = 1.0; double kapk = 1.0;
        // c_bxb_0_ -> point(*ii,*jj) = (2*eps*kapj - sigj*dt_) / (2*eps*kapj + sigj*dt_);
        // c_bxe_0_ -> point(*ii,*jj) = (2 * eps * dt_) / (dy_ * (2*eps*kapj + sigj*dt_));
        // c_hxh_0_ -> point(*ii,*jj) = (2*eps*kapk - sigk*dt_) / (2*eps*kapk + sigk*dt_);
        // c_hxb0_0_-> point(*ii,*jj) = (2*eps*kapi - sigi*dt_) / (2*eps*kapk + sigk*dt_);
        // c_hxb1_0_-> point(*ii,*jj) = (2*eps*kapi + sigi*dt_) / (2*eps*kapk + sigk*dt_);
        preFact[0] = (2*eps*kapj - sigj*dt_) / (2*eps*kapj + sigj*dt_);
        preFact[1] = (2 * eps * dt_) / (dy_ * (2*eps*kapj + sigj*dt_));
        preFact[2] = (2*eps*kapk - sigk*dt_) / (2*eps*kapk + sigk*dt_);
        preFact[3] = (2*eps*kapi - sigi*dt_) / (2*eps*kapk + sigk*dt_);
        preFact[4] = (2*eps*kapi + sigi*dt_) / (2*eps*kapk + sigk*dt_);
        return preFact;
    }

};

#endif