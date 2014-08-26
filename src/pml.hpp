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
    int nj_;

public:
    std::shared_ptr<Grid2D<T>> Dx_,Dy,Dz_,Bx_,By_,Bz_,Dx_end_,Dyend_,Dz_end_,Bx_end_,By_end_,Bz_end_;
    std::shared_ptr<Grid2D<int>> phys_Hx_,phys_Hy_,phys_Hz_, phys_Hx_end_,phys_Hy_end_,phys_Hz_end_,phys_Ex_,phys_Ey_,phys_Ez_, phys_Ex_end_,phys_Ey_end_,phys_Ez_end_;

    std::shared_ptr<Grid2D<double>> c_bxb_0_, c_bxe_0_, c_hxh_0_, c_hxb0_0_, c_hxb1_0_,c_bxb_n_, c_bxe_n_, c_hxh_n_, c_hxb0_n_, c_hxb1_n_;
    std::shared_ptr<Grid2D<double>> c_byb_0_, c_bye_0_, c_hyh_0_, c_hyb0_0_, c_hyb1_0_,c_byb_n_, c_bye_n_, c_hyh_n_, c_hyb0_n_, c_hyb1_n_;
    std::shared_ptr<Grid2D<double>> c_bzb_0_, c_bze_0_, c_hzh_0_, c_hzb0_0_, c_hzb1_0_,c_bzb_n_, c_bze_n_, c_hzh_n_, c_hzb0_n_, c_hzb1_n_;
    std::shared_ptr<Grid2D<double>> c_dxd_0_, c_dxh_0_, c_exe_0_, c_exd0_0_, c_exd1_0_,c_dxd_n_, c_dxh_n_, c_exe_n_, c_exd0_n_, c_exd1_n_;
    std::shared_ptr<Grid2D<double>> c_dyd_0_, c_dyh_0_, c_eye_0_, c_eyd0_0_, c_eyd1_0_,c_dyd_n_, c_dyh_n_, c_eye_n_, c_eyd0_n_, c_eyd1_n_;
    std::shared_ptr<Grid2D<double>> c_dzd_0_, c_dzh_0_, c_eze_0_, c_ezd0_0_, c_ezd1_0_,c_dzd_n_, c_dzh_n_, c_eze_n_, c_ezd0_n_, c_ezd1_n_;


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
    UPML(int thickness, Direction d, double m, double R0, int nx, double dx, double dy, Polarization pol, bool precalc) : thickness_(thickness), d_(d), m_(m), R0_(R0), precalc_(precalc), pol_(pol),nj_(nx)
    {
        sigmaMax_ = -(m_+1)*log(R0_)/(2*thickness_*dx); // eta should be included;
        kappaMax_ = 1.0;
        if(pol == EX || pol == EY || pol == HZ)
        {
            Dx_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
            Dy = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
            Bz_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
            Dx_end_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
            Dyend_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
            Bz_end_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);

            Bx_ = nullptr;
            By_ = nullptr;
            Dz_ = nullptr;
            Bx_end_ = nullptr;
            By_end_ = nullptr;
            Dz_end_ = nullptr;

            phys_Hz_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
            phys_Hz_end_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
            phys_Ex_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
            phys_Ex_end_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
            phys_Ey_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
            phys_Ey_end_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);


            phys_Hy_ = nullptr;
            phys_Hy_end_ = nullptr;
            phys_Hx_ = nullptr;
            phys_Hx_end_ = nullptr;
            phys_Ez_ = nullptr;
            phys_Ez_end_ = nullptr;

            if(precalc_ == false)
            {
                c_bxb_0_ = nullptr; c_bxe_0_ = nullptr; c_hxh_0_ = nullptr; c_hxb0_0_ = nullptr; c_hxb1_0_ = nullptr;
                c_bxb_n_ = nullptr; c_bxe_n_ = nullptr; c_hxh_n_ = nullptr; c_hxb0_n_ = nullptr; c_hxb1_n_ = nullptr;
                c_byb_0_ = nullptr; c_bye_0_ = nullptr; c_hyh_0_ = nullptr; c_hyb0_0_ = nullptr; c_hyb1_0_ = nullptr;
                c_byb_n_ = nullptr; c_bye_n_ = nullptr; c_hyh_n_ = nullptr; c_hyb0_n_ = nullptr; c_hyb1_n_ = nullptr;
                c_bzb_0_ = nullptr; c_bze_0_ = nullptr; c_hzh_0_ = nullptr; c_hzb0_0_ = nullptr; c_hzb1_0_ = nullptr;
                c_bzb_n_ = nullptr; c_bze_n_ = nullptr; c_hzh_n_ = nullptr; c_hzb0_n_ = nullptr; c_hzb1_n_ = nullptr;
                c_dxd_0_ = nullptr; c_dxh_0_ = nullptr; c_exe_0_ = nullptr; c_exd0_0_ = nullptr; c_exd1_0_ = nullptr;
                c_dxd_n_ = nullptr; c_dxh_n_ = nullptr; c_exe_n_ = nullptr; c_exd0_n_ = nullptr; c_exd1_n_ = nullptr;
                c_dyd_0_ = nullptr; c_dyh_0_ = nullptr; c_eye_0_ = nullptr; c_eyd0_0_ = nullptr; c_eyd1_0_ = nullptr;
                c_dyd_n_ = nullptr; c_dyh_n_ = nullptr; c_eye_n_ = nullptr; c_eyd0_n_ = nullptr; c_eyd1_n_ = nullptr;
                c_dzd_0_ = nullptr; c_dzh_0_ = nullptr; c_eze_0_ = nullptr; c_ezd0_0_ = nullptr; c_ezd1_0_ = nullptr;
                c_dzd_n_ = nullptr; c_dzh_n_ = nullptr; c_eze_n_ = nullptr; c_ezd0_n_ = nullptr; c_ezd1_n_ = nullptr;
            }
            else
            {
                c_bxb_0_  = nullptr; c_bxe_0_ = nullptr; c_hxh_0_ = nullptr; c_hxb0_0_ = nullptr; c_hxb1_0_ = nullptr;
                c_bxb_n_  = nullptr; c_bxe_n_ = nullptr; c_hxh_n_ = nullptr; c_hxb0_n_ = nullptr; c_hxb1_n_ = nullptr;
                c_byb_0_  = nullptr; c_bye_0_ = nullptr; c_hyh_0_ = nullptr; c_hyb0_0_ = nullptr; c_hyb1_0_ = nullptr;
                c_byb_n_  = nullptr; c_bye_n_ = nullptr; c_hyh_n_ = nullptr; c_hyb0_n_ = nullptr; c_hyb1_n_ = nullptr;

                c_bzb_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_bze_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hzh_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hzb0_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hzb1_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_bzb_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_bze_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hzh_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hzb0_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hzb1_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                c_dxd_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_dxh_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_exe_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_exd0_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_exd1_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                c_dxd_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_dxh_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_exe_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_exd0_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_exd1_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                c_dyd_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_dyh_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_eye_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_eyd0_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_eyd1_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                c_dyd_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_dyh_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_eye_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_eyd0_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_eyd1_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                c_dzd_0_  = nullptr; c_dzh_0_ = nullptr; c_eze_0_ = nullptr; c_ezd0_0_ = nullptr; c_ezd1_0_ = nullptr;
                c_dzd_n_  = nullptr; c_dzh_n_ = nullptr; c_eze_n_ = nullptr; c_ezd0_n_ = nullptr; c_ezd1_n_ = nullptr;
            }
        }
        else
        {
            Bx_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
            By_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
            Dz_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
            Dx_ = nullptr;
            Dy = nullptr;
            Bz_ = nullptr;

            Bx_end_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
            By_end_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
            Dz_end_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
            Dx_end_ = nullptr;
            Dyend_ = nullptr;
            Bz_end_ = nullptr;

            phys_Hz_ = nullptr;
            phys_Hz_end_ = nullptr;
            phys_Hy_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
            phys_Hy_end_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
            phys_Hx_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
            phys_Hx_end_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);

            phys_Ez_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
            phys_Ez_end_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
            phys_Ex_ = nullptr;
            phys_Ex_end_ = nullptr;
            phys_Ey_ = nullptr;
            phys_Ey_end_ = nullptr;
           
            if(precalc_ == false)
            {
                c_bxb_0_ = nullptr; c_bxe_0_ = nullptr; c_hxh_0_ = nullptr; c_hxb0_0_ = nullptr; c_hxb1_0_ = nullptr;
                c_bxb_n_ = nullptr; c_bxe_n_ = nullptr; c_hxh_n_ = nullptr; c_hxb0_n_ = nullptr; c_hxb1_n_ = nullptr;
                c_byb_0_ = nullptr; c_bye_0_ = nullptr; c_hyh_0_ = nullptr; c_hyb0_0_ = nullptr; c_hyb1_0_ = nullptr;
                c_byb_n_ = nullptr; c_bye_n_ = nullptr; c_hyh_n_ = nullptr; c_hyb0_n_ = nullptr; c_hyb1_n_ = nullptr;
                c_bzb_0_ = nullptr; c_bze_0_ = nullptr; c_hzh_0_ = nullptr; c_hzb0_0_ = nullptr; c_hzb1_0_ = nullptr;
                c_bzb_n_ = nullptr; c_bze_n_ = nullptr; c_hzh_n_ = nullptr; c_hzb0_n_ = nullptr; c_hzb1_n_ = nullptr;
                c_dxd_0_ = nullptr; c_dxh_0_ = nullptr; c_exe_0_ = nullptr; c_exd0_0_ = nullptr; c_exd1_0_ = nullptr;
                c_dxd_n_ = nullptr; c_dxh_n_ = nullptr; c_exe_n_ = nullptr; c_exd0_n_ = nullptr; c_exd1_n_ = nullptr;
                c_dyd_0_ = nullptr; c_dyh_0_ = nullptr; c_eye_0_ = nullptr; c_eyd0_0_ = nullptr; c_eyd1_0_ = nullptr;
                c_dyd_n_ = nullptr; c_dyh_n_ = nullptr; c_eye_n_ = nullptr; c_eyd0_n_ = nullptr; c_eyd1_n_ = nullptr;
                c_dzd_0_ = nullptr; c_dzh_0_ = nullptr; c_eze_0_ = nullptr; c_ezd0_0_ = nullptr; c_ezd1_0_ = nullptr;
                c_dzd_n_ = nullptr; c_dzh_n_ = nullptr; c_eze_n_ = nullptr; c_ezd0_n_ = nullptr; c_ezd1_n_ = nullptr;
            }
            else
            {
                c_bxb_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_bxe_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hxh_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hxb0_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hxb1_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                c_bxb_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_bxe_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hxh_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hxb0_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hxb1_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                c_byb_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_bye_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hyh_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hyb0_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hyb1_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                c_byb_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_bye_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hyh_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hyb0_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_hyb1_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                c_bzb_0_  = nullptr; c_bze_0_ = nullptr; c_hzh_0_ = nullptr; c_hzb0_0_ = nullptr; c_hzb1_0_ = nullptr;
                c_bzb_n_  = nullptr; c_bze_n_ = nullptr; c_hzh_n_ = nullptr; c_hzb0_n_ = nullptr; c_hzb1_n_ = nullptr;
                c_dxd_0_  = nullptr; c_dxh_0_ = nullptr; c_exe_0_ = nullptr; c_exd0_0_ = nullptr; c_exd1_0_ = nullptr;
                c_dxd_n_  = nullptr; c_dxh_n_ = nullptr; c_exe_n_ = nullptr; c_exd0_n_ = nullptr; c_exd1_n_ = nullptr;
                c_dyd_0_  = nullptr; c_dyh_0_ = nullptr; c_eye_0_ = nullptr; c_eyd0_0_ = nullptr; c_eyd1_0_ = nullptr;
                c_dyd_n_  = nullptr; c_dyh_n_ = nullptr; c_eye_n_ = nullptr; c_eyd0_n_ = nullptr; c_eyd1_n_ = nullptr;

                c_dzd_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_dzh_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_eze_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_ezd0_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_ezd1_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                c_dzd_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_dzh_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_eze_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_ezd0_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                c_ezd1_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
            }
        }
        /*if(d == X)
        {
            if(pol == EX || pol == EY || pol == HZ)
            {
                Dx_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
                Dy = std::make_shared<Grid2D<T>>(thickness,nx-1,dx,dy);
                Bz_ = std::make_shared<Grid2D<T>>(thickness,nx-1,dx,dy);
                Dx_end_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
                Dyend_ = std::make_shared<Grid2D<T>>(thickness,nx-1,dx,dy);
                Bz_end_ = std::make_shared<Grid2D<T>>(thickness,nx-1,dx,dy);

                Bx_ = nullptr;
                By_ = nullptr;
                Dz_ = nullptr;
                Bx_end_ = nullptr;
                By_end_ = nullptr;
                Dz_end_ = nullptr;

                phys_Hz_ = std::make_shared<Grid2D<int>>(thickness,nx-1,dx,dy);
                phys_Hz_end_ = std::make_shared<Grid2D<int>>(thickness,nx-1,dx,dy);

                phys_Hy_ = nullptr;
                phys_Hy_end_ = nullptr;
                phys_Hx_ = nullptr;
                phys_Hx_end_ = nullptr;

                if(precalc_ == false)
                {
                    c_bxb_0_ = nullptr; c_bxe_0_ = nullptr; c_hxh_0_ = nullptr; c_hxb0_0_ = nullptr; c_hxb1_0_ = nullptr;
                    c_bxb_n_ = nullptr; c_bxe_n_ = nullptr; c_hxh_n_ = nullptr; c_hxb0_n_ = nullptr; c_hxb1_n_ = nullptr;
                    c_byb_0_ = nullptr; c_bye_0_ = nullptr; c_hyh_0_ = nullptr; c_hyb0_0_ = nullptr; c_hyb1_0_ = nullptr;
                    c_byb_n_ = nullptr; c_bye_n_ = nullptr; c_hyh_n_ = nullptr; c_hyb0_n_ = nullptr; c_hyb1_n_ = nullptr;
                    c_bzb_0_ = nullptr; c_bze_0_ = nullptr; c_hzh_0_ = nullptr; c_hzb0_0_ = nullptr; c_hzb1_0_ = nullptr;
                    c_bzb_n_ = nullptr; c_bze_n_ = nullptr; c_hzh_n_ = nullptr; c_hzb0_n_ = nullptr; c_hzb1_n_ = nullptr;
                    c_dxd_0_ = nullptr; c_dxh_0_ = nullptr; c_exe_0_ = nullptr; c_exd0_0_ = nullptr; c_exd1_0_ = nullptr;
                    c_dxd_n_ = nullptr; c_dxh_n_ = nullptr; c_exe_n_ = nullptr; c_exd0_n_ = nullptr; c_exd1_n_ = nullptr;
                    c_dyd_0_ = nullptr; c_dyh_0_ = nullptr; c_eye_0_ = nullptr; c_eyd0_0_ = nullptr; c_eyd1_0_ = nullptr;
                    c_dyd_n_ = nullptr; c_dyh_n_ = nullptr; c_eye_n_ = nullptr; c_eyd0_n_ = nullptr; c_eyd1_n_ = nullptr;
                    c_dzd_0_ = nullptr; c_dzh_0_ = nullptr; c_eze_0_ = nullptr; c_ezd0_0_ = nullptr; c_ezd1_0_ = nullptr;
                    c_dzd_n_ = nullptr; c_dzh_n_ = nullptr; c_eze_n_ = nullptr; c_ezd0_n_ = nullptr; c_ezd1_n_ = nullptr;
                }
                else
                {
                    c_bxb_0_  = nullptr; c_bxe_0_ = nullptr; c_hxh_0_ = nullptr; c_hxb0_0_ = nullptr; c_hxb1_0_ = nullptr;
                    c_bxb_n_  = nullptr; c_bxe_n_ = nullptr; c_hxh_n_ = nullptr; c_hxb0_n_ = nullptr; c_hxb1_n_ = nullptr;
                    c_byb_0_  = nullptr; c_bye_0_ = nullptr; c_hyh_0_ = nullptr; c_hyb0_0_ = nullptr; c_hyb1_0_ = nullptr;
                    c_byb_n_  = nullptr; c_bye_n_ = nullptr; c_hyh_n_ = nullptr; c_hyb0_n_ = nullptr; c_hyb1_n_ = nullptr;

                    c_bzb_0_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_bze_0_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_hzh_0_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_hzb0_0_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_hzb1_0_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_bzb_n_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_bze_n_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_hzh_n_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_hzb0_n_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_hzb1_n_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);

                    c_dxd_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_dxh_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_exe_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_exd0_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_exd1_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                    c_dxd_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_dxh_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_exe_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_exd0_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_exd1_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                    c_dyd_0_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_dyh_0_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_eye_0_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_eyd0_0_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_eyd1_0_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);

                    c_dyd_n_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_dyh_n_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_eye_n_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_eyd0_n_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_eyd1_n_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);

                    c_dzd_0_  = nullptr; c_dzh_0_ = nullptr; c_eze_0_ = nullptr; c_ezd0_0_ = nullptr; c_ezd1_0_ = nullptr;
                    c_dzd_n_  = nullptr; c_dzh_n_ = nullptr; c_eze_n_ = nullptr; c_ezd0_n_ = nullptr; c_ezd1_n_ = nullptr;
                }
            }
            else
            {
                Bx_ = std::make_shared<Grid2D<T>>(thickness,nx-1,dx,dy);
                By_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
                Dz_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
                Dx_ = nullptr;
                Dy = nullptr;
                Bz_ = nullptr;

                Bx_end_ = std::make_shared<Grid2D<T>>(thickness,nx-1,dx,dy);
                By_end_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
                Dz_end_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
                Dx_end_ = nullptr;
                Dyend_ = nullptr;
                Bz_end_ = nullptr;

                phys_Hz_ = nullptr;
                phys_Hz_end_ = nullptr;
                phys_Hy_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
                phys_Hy_end_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
                phys_Hx_ = std::make_shared<Grid2D<int>>(thickness,nx-1,dx,dy);
                phys_Hx_end_ = std::make_shared<Grid2D<int>>(thickness,nx-1,dx,dy);
                if(precalc_ == false)
                {
                    c_bxb_0_ = nullptr; c_bxe_0_ = nullptr; c_hxh_0_ = nullptr; c_hxb0_0_ = nullptr; c_hxb1_0_ = nullptr;
                    c_bxb_n_ = nullptr; c_bxe_n_ = nullptr; c_hxh_n_ = nullptr; c_hxb0_n_ = nullptr; c_hxb1_n_ = nullptr;
                    c_byb_0_ = nullptr; c_bye_0_ = nullptr; c_hyh_0_ = nullptr; c_hyb0_0_ = nullptr; c_hyb1_0_ = nullptr;
                    c_byb_n_ = nullptr; c_bye_n_ = nullptr; c_hyh_n_ = nullptr; c_hyb0_n_ = nullptr; c_hyb1_n_ = nullptr;
                    c_bzb_0_ = nullptr; c_bze_0_ = nullptr; c_hzh_0_ = nullptr; c_hzb0_0_ = nullptr; c_hzb1_0_ = nullptr;
                    c_bzb_n_ = nullptr; c_bze_n_ = nullptr; c_hzh_n_ = nullptr; c_hzb0_n_ = nullptr; c_hzb1_n_ = nullptr;
                    c_dxd_0_ = nullptr; c_dxh_0_ = nullptr; c_exe_0_ = nullptr; c_exd0_0_ = nullptr; c_exd1_0_ = nullptr;
                    c_dxd_n_ = nullptr; c_dxh_n_ = nullptr; c_exe_n_ = nullptr; c_exd0_n_ = nullptr; c_exd1_n_ = nullptr;
                    c_dyd_0_ = nullptr; c_dyh_0_ = nullptr; c_eye_0_ = nullptr; c_eyd0_0_ = nullptr; c_eyd1_0_ = nullptr;
                    c_dyd_n_ = nullptr; c_dyh_n_ = nullptr; c_eye_n_ = nullptr; c_eyd0_n_ = nullptr; c_eyd1_n_ = nullptr;
                    c_dzd_0_ = nullptr; c_dzh_0_ = nullptr; c_eze_0_ = nullptr; c_ezd0_0_ = nullptr; c_ezd1_0_ = nullptr;
                    c_dzd_n_ = nullptr; c_dzh_n_ = nullptr; c_eze_n_ = nullptr; c_ezd0_n_ = nullptr; c_ezd1_n_ = nullptr;
                }
                else
                {
                    c_bxb_0_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_bxe_0_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_hxh_0_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_hxb0_0_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_hxb1_0_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);

                    c_bxb_n_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_bxe_n_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_hxh_n_  = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_hxb0_n_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
                    c_hxb1_n_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);

                    c_byb_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_bye_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_hyh_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_hyb0_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_hyb1_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                    c_byb_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_bye_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_hyh_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_hyb0_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_hyb1_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                    c_bzb_0_  = nullptr; c_bze_0_ = nullptr; c_hzh_0_ = nullptr; c_hzb0_0_ = nullptr; c_hzb1_0_ = nullptr;
                    c_bzb_n_  = nullptr; c_bze_n_ = nullptr; c_hzh_n_ = nullptr; c_hzb0_n_ = nullptr; c_hzb1_n_ = nullptr;
                    c_dxd_0_  = nullptr; c_dxh_0_ = nullptr; c_exe_0_ = nullptr; c_exd0_0_ = nullptr; c_exd1_0_ = nullptr;
                    c_dxd_n_  = nullptr; c_dxh_n_ = nullptr; c_exe_n_ = nullptr; c_exd0_n_ = nullptr; c_exd1_n_ = nullptr;
                    c_dyd_0_  = nullptr; c_dyh_0_ = nullptr; c_eye_0_ = nullptr; c_eyd0_0_ = nullptr; c_eyd1_0_ = nullptr;
                    c_dyd_n_  = nullptr; c_dyh_n_ = nullptr; c_eye_n_ = nullptr; c_eyd0_n_ = nullptr; c_eyd1_n_ = nullptr;

                    c_dzd_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_dzh_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_eze_0_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_ezd0_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_ezd1_0_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);

                    c_dzd_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_dzh_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_eze_n_  = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_ezd0_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                    c_ezd1_n_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
                }
            }
        }
        else if (d == Y)
        {
            if(pol == EX || pol == EY || pol == HZ)
            {
                Dx_ = std::make_shared<Grid2D<T>>(nx-1,thickness,dx,dy);
                Dy = std::make_shared<Grid2D<T>>(nx,thickness,dx,dy);
                Bz_ = std::make_shared<Grid2D<T>>(nx-1,thickness,dx,dy);
                Bx_ = nullptr;
                By_ = nullptr;
                Dz_ = nullptr;
                Dx_end_ = std::make_shared<Grid2D<T>>(nx-1,thickness,dx,dy);
                Dyend_ = std::make_shared<Grid2D<T>>(nx,thickness,dx,dy);
                Bz_end_ = std::make_shared<Grid2D<T>>(nx-1,thickness,dx,dy);
                Bx_end_ = nullptr;
                By_end_ = nullptr;
                Dz_end_ = nullptr;
                phys_Hz_ = std::make_shared<Grid2D<int>>(nx-1,thickness,dx,dy);
                phys_Hz_end_ = std::make_shared<Grid2D<int>>(nx-1,thickness,dx,dy);
                phys_Hy_ = nullptr;
                phys_Hy_end_ = nullptr;
                phys_Hx_ = nullptr;
                phys_Hx_end_ = nullptr;
                if(precalc_ == true)
                {
                    c_bxb_0_ = nullptr; c_bxe_0_ = nullptr; c_hxh_0_ = nullptr; c_hxb0_0_ = nullptr; c_hxb1_0_ = nullptr;
                    c_bxb_n_ = nullptr; c_bxe_n_ = nullptr; c_hxh_n_ = nullptr; c_hxb0_n_ = nullptr; c_hxb1_n_ = nullptr;
                    c_byb_0_ = nullptr; c_bye_0_ = nullptr; c_hyh_0_ = nullptr; c_hyb0_0_ = nullptr; c_hyb1_0_ = nullptr;
                    c_byb_n_ = nullptr; c_bye_n_ = nullptr; c_hyh_n_ = nullptr; c_hyb0_n_ = nullptr; c_hyb1_n_ = nullptr;

                    c_bzb_0_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_bze_0_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_hzh_0_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_hzb0_0_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_hzb1_0_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);

                    c_bzb_n_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_bze_n_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_hzh_n_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_hzb0_n_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_hzb1_n_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);

                    c_dxd_0_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_dxh_0_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_exe_0_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_exd0_0_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_exd1_0_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);

                    c_dxd_n_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_dxh_n_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_exe_n_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_exd0_n_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_exd1_n_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);

                    c_dyd_0_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_dyh_0_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_eye_0_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_eyd0_0_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_eyd1_0_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);

                    c_dyd_n_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_dyh_n_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_eye_n_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_eyd0_n_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_eyd1_n_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);

                    c_dzd_0_ = nullptr; c_dzh_0_ = nullptr; c_eze_0_ = nullptr; c_ezd0_0_ = nullptr; c_ezd1_0_ = nullptr;
                    c_dzd_n_ = nullptr; c_dzh_n_ = nullptr; c_eze_n_ = nullptr; c_ezd0_n_ = nullptr; c_ezd1_n_ = nullptr;
                }
                else
                {
                    c_bxb_0_ = nullptr; c_bxe_0_ = nullptr; c_hxh_0_ = nullptr; c_hxb0_0_ = nullptr; c_hxb1_0_ = nullptr;
                    c_bxb_n_ = nullptr; c_bxe_n_ = nullptr; c_hxh_n_ = nullptr; c_hxb0_n_ = nullptr; c_hxb1_n_ = nullptr;
                    c_byb_0_ = nullptr; c_bye_0_ = nullptr; c_hyh_0_ = nullptr; c_hyb0_0_ = nullptr; c_hyb1_0_ = nullptr;
                    c_byb_n_ = nullptr; c_bye_n_ = nullptr; c_hyh_n_ = nullptr; c_hyb0_n_ = nullptr; c_hyb1_n_ = nullptr;
                    c_bzb_0_ = nullptr; c_bze_0_ = nullptr; c_hzh_0_ = nullptr; c_hzb0_0_ = nullptr; c_hzb1_0_ = nullptr;
                    c_bzb_n_ = nullptr; c_bze_n_ = nullptr; c_hzh_n_ = nullptr; c_hzb0_n_ = nullptr; c_hzb1_n_ = nullptr;
                    c_dxd_0_ = nullptr; c_dxh_0_ = nullptr; c_exe_0_ = nullptr; c_exd0_0_ = nullptr; c_exd1_0_ = nullptr;
                    c_dxd_n_ = nullptr; c_dxh_n_ = nullptr; c_exe_n_ = nullptr; c_exd0_n_ = nullptr; c_exd1_n_ = nullptr;
                    c_dyd_0_ = nullptr; c_dyh_0_ = nullptr; c_eye_0_ = nullptr; c_eyd0_0_ = nullptr; c_eyd1_0_ = nullptr;
                    c_dyd_n_ = nullptr; c_dyh_n_ = nullptr; c_eye_n_ = nullptr; c_eyd0_n_ = nullptr; c_eyd1_n_ = nullptr;
                    c_dzd_0_ = nullptr; c_dzh_0_ = nullptr; c_eze_0_ = nullptr; c_ezd0_0_ = nullptr; c_ezd1_0_ = nullptr;
                    c_dzd_n_ = nullptr; c_dzh_n_ = nullptr; c_eze_n_ = nullptr; c_ezd0_n_ = nullptr; c_ezd1_n_ = nullptr;
                }
            }
            else
            {
                Bx_ = std::make_shared<Grid2D<T>>(nx,thickness,dx,dy);
                By_ = std::make_shared<Grid2D<T>>(nx-1,thickness,dx,dy);
                Dz_ = std::make_shared<Grid2D<T>>(nx,thickness,dx,dy);
                Dx_ = nullptr;
                Dy = nullptr;
                Bz_ = nullptr;
                Bx_end_ = std::make_shared<Grid2D<T>>(nx,thickness,dx,dy);
                By_end_ = std::make_shared<Grid2D<T>>(nx-1,thickness,dx,dy);
                Dz_end_ = std::make_shared<Grid2D<T>>(nx,thickness,dx,dy);
                Dx_end_ = nullptr;
                Dyend_ = nullptr;
                Bz_end_ = nullptr;
                phys_Hz_ = nullptr;
                phys_Hz_end_ = nullptr;
                phys_Hy_ = std::make_shared<Grid2D<int>>(nx-1,thickness,dx,dy);
                phys_Hy_end_ = std::make_shared<Grid2D<int>>(nx-1,thickness,dx,dy);
                phys_Hx_ = std::make_shared<Grid2D<int>>(nx,thickness,dx,dy);
                phys_Hx_end_ = std::make_shared<Grid2D<int>>(nx,thickness,dx,dy);
                if(precalc_ == false)
                {
                    c_bxb_0_ = nullptr; c_bxe_0_ = nullptr; c_hxh_0_ = nullptr; c_hxb0_0_ = nullptr; c_hxb1_0_ = nullptr;
                    c_bxb_n_ = nullptr; c_bxe_n_ = nullptr; c_hxh_n_ = nullptr; c_hxb0_n_ = nullptr; c_hxb1_n_ = nullptr;
                    c_byb_0_ = nullptr; c_bye_0_ = nullptr; c_hyh_0_ = nullptr; c_hyb0_0_ = nullptr; c_hyb1_0_ = nullptr;
                    c_byb_n_ = nullptr; c_bye_n_ = nullptr; c_hyh_n_ = nullptr; c_hyb0_n_ = nullptr; c_hyb1_n_ = nullptr;
                    c_bzb_0_ = nullptr; c_bze_0_ = nullptr; c_hzh_0_ = nullptr; c_hzb0_0_ = nullptr; c_hzb1_0_ = nullptr;
                    c_bzb_n_ = nullptr; c_bze_n_ = nullptr; c_hzh_n_ = nullptr; c_hzb0_n_ = nullptr; c_hzb1_n_ = nullptr;
                    c_dxd_0_ = nullptr; c_dxh_0_ = nullptr; c_exe_0_ = nullptr; c_exd0_0_ = nullptr; c_exd1_0_ = nullptr;
                    c_dxd_n_ = nullptr; c_dxh_n_ = nullptr; c_exe_n_ = nullptr; c_exd0_n_ = nullptr; c_exd1_n_ = nullptr;
                    c_dyd_0_ = nullptr; c_dyh_0_ = nullptr; c_eye_0_ = nullptr; c_eyd0_0_ = nullptr; c_eyd1_0_ = nullptr;
                    c_dyd_n_ = nullptr; c_dyh_n_ = nullptr; c_eye_n_ = nullptr; c_eyd0_n_ = nullptr; c_eyd1_n_ = nullptr;
                    c_dzd_0_ = nullptr; c_dzh_0_ = nullptr; c_eze_0_ = nullptr; c_ezd0_0_ = nullptr; c_ezd1_0_ = nullptr;
                    c_dzd_n_ = nullptr; c_dzh_n_ = nullptr; c_eze_n_ = nullptr; c_ezd0_n_ = nullptr; c_ezd1_n_ = nullptr;
                }
                else
                {
                    c_bxb_0_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_bxe_0_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_hxh_0_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_hxb0_0_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_hxb1_0_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);

                    c_bxb_n_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_bxe_n_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_hxh_n_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_hxb0_n_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_hxb1_n_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);

                    c_byb_0_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_bye_0_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_hyh_0_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_hyb0_0_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_hyb1_0_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);

                    c_byb_n_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_bye_n_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_hyh_n_  = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_hyb0_n_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
                    c_hyb1_n_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);

                    c_bzb_0_ = nullptr; c_bze_0_ = nullptr; c_hzh_0_ = nullptr; c_hzb0_0_ = nullptr; c_hzb1_0_ = nullptr;
                    c_bzb_n_ = nullptr; c_bze_n_ = nullptr; c_hzh_n_ = nullptr; c_hzb0_n_ = nullptr; c_hzb1_n_ = nullptr;
                    c_dxd_0_ = nullptr; c_dxh_0_ = nullptr; c_exe_0_ = nullptr; c_exd0_0_ = nullptr; c_exd1_0_ = nullptr;
                    c_dxd_n_ = nullptr; c_dxh_n_ = nullptr; c_exe_n_ = nullptr; c_exd0_n_ = nullptr; c_exd1_n_ = nullptr;
                    c_dyd_0_ = nullptr; c_dyh_0_ = nullptr; c_eye_0_ = nullptr; c_eyd0_0_ = nullptr; c_eyd1_0_ = nullptr;
                    c_dyd_n_ = nullptr; c_dyh_n_ = nullptr; c_eye_n_ = nullptr; c_eyd0_n_ = nullptr; c_eyd1_n_ = nullptr;

                    c_dzd_0_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_dzh_0_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_eze_0_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_ezd0_0_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_ezd1_0_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);

                    c_dzd_n_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_dzh_n_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_eze_n_  = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_ezd0_n_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                    c_ezd1_n_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
                }
            }

        }
        else
             throw std::logic_error("While yes we could have a thrid dimension to run, I have yet to be implimented to do such a thing. So please accept this error as my sincerest appology.");
        */
    }

    void initializeUPML(std::vector<Obj> objArr, int nx, int ny, double dx, double dy, double dt, int yPML, int xPML, UPML* oppDir)
    {
        int jmax = 0;
        int pt_i = 0;
        int pt_j = 0;
        int ni   = 0;
        int nj   = 0;
        double di = 0.0;
        double dj = 0.0;
        if (d_ == X)
        {
            jmax = ny-1;
            pt_i = 0;
            pt_j = 1;
            di   = dx;
            dj   = dy;
            ni   = nx;
            nj   = ny;
        }
        else
        {
            jmax = nx-1;
            pt_i = 1;
            pt_j = 0;
            dj   = dx;
            di   = dy;
            nj   = nx;
            ni   = ny;
        }

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
                    
                    for(int ii =0; ii < thickness_; ii ++)
                    {
                        for(int jj = 0; jj < jmax; jj++)
                        {
                            pt[pt_i] = (ii+0.5-(ni-1)/2.0)*di;
                            pt[pt_j] = (jj-(nj-1)/2.0)*dj;
                            if(objArr[kk].isObj(pt)==true)
                                phys_Hy_->point(ii,jj) = kk;
                            pt[pt_i] -= 0.5*di;
                            if(objArr[kk].isObj(pt)==true)
                                phys_Ez_->point(ii,jj) = kk;
                            pt[pt_j] += 0.5*dj;
                            if(objArr[kk].isObj(pt)==true && kk ==1)
                                phys_Hx_->point(ii,jj) = kk;
                            pt[pt_i] = ((ni-1-ii)+0.5-(ni-1)/2.0)*di;
                            pt[pt_j] = (jj-(nj-1)/2.0)*dj;
                            if(objArr[kk].isObj(pt)==true)
                                phys_Hy_end_->point(ii,jj) = kk;
                            pt[pt_i] -= 0.5*di;
                            if(objArr[kk].isObj(pt)==true)
                                phys_Ez_->point(ni-1-ii,jj) = kk;
                            pt[pt_j] += 0.5*dj;
                            if(objArr[kk].isObj(pt)==true && kk ==1)
                                phys_Hx_end_->point(ii,jj) = kk;
                        }
                        pt[pt_i]=(ii+0.5-(ni-1)/2.0)*di;
                        pt[pt_j]=((nj-1)/2.0)*dj;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Hy_->point(ii,nj-1) = kk;
                        pt[pt_i] -= 0.5*di;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Ez_->point(ii,nj-1) = kk;

                        pt[pt_i]=((ni-1-ii)+0.5-(ni-1)/2.0)*di;
                        pt[pt_j]=((nj-1)/2.0)*dj;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Hy_end_->point(ii,nj-1) = kk;
                        pt[pt_i] -= 0.5*di;
                        if(objArr[kk].isObj(pt)==true)
                            phys_Ez_->point(ii,nj-1) = kk;
                    }
                }
            }
        }
        if(precalc_)
        {
            double eps=0.0;
            double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
           
            double sigz = 0.0; double sigx = 0.0; double sigy = 0.0;
            double sigxx = 0.0; double sigxy = 0.0; double sigyx = 0.0; double sigyy = 0.0;
            
            double *sigmax(double ,double);
            double *sigmay(double ,double);

            if(d_==X)
            {
                sigmax = sigma;
                if (oppDir)
                    sigmay = oppDir->sigma;
                else
                    sigmay = [] (double x,double eps) { return 0.0; };
            }
            else
            {
                sigmay = sigma;
                if (oppDir)
                    sigmax = oppDir->sigma;
                else
                    sigmax = [] (double x,double eps) { return 0.0; };
            }
            if(pol_ == EZ || pol_ == HX || pol_ == HY)
            {
                for(int ii = 0; ii < thickness_; ii++)
                {
                    for(int jj = 0; jj < nj; jj++)
                    {
                        //Update Hx factors
                        eps    = objArr[phys_Hx_->point(ii,jj)].dielectric(1.0);
                        sigxx  = sigmax(static_cast<double>(ii),eps);
                        sigyx  = sigmay(static_cast<double>(ii) + 0.5,eps);
                        c_bxb_0_ -> point(ii,jj) = (2*eps*kapy - sigyx*dt) / (2*eps*kapy + sigyx*dt);
                        c_bxe_0_ -> point(ii,jj) = (2 * eps * dt) / (dy * (2*eps*kapy + sigyx*dt));
                        c_hxh_0_ -> point(ii,jj) = (2*eps*kapz - sigz*dt) / (2*eps*kapz + sigz*dt);
                        c_hxb0_0_-> point(ii,jj) = (2*eps*kapx - sigxx*dt) / (2*eps*kapz + sigz*dt);
                        c_hxb1_0_-> point(ii,jj) = (2*eps*kapx + sigxx*dt) / (2*eps*kapz + sigz*dt);

                        eps    = objArr[phys_Hx_end_->point(ii,jj)].dielectric(1.0);
                        sigxx  = sigmax(static_cast<double>(ii),eps);
                        sigyx  = sigmay(static_cast<double>(ii)-0.5,eps);
                        c_bxb_n_ -> point(ii,jj) = (2*eps*kapy - sigyx*dt) / (2*eps*kapy + sigyx*dt);
                        c_bxe_n_ -> point(ii,jj) = (2 * eps * dt) / (dy * (2*eps*kapy + sigyx*dt));
                        c_hxh_n_ -> point(ii,jj) = (2*eps*kapz - sigz*dt) / (2*eps*kapz + sigz*dt);
                        c_hxb0_n_-> point(ii,jj) = (2*eps*kapx - sigxx*dt) / (2*eps*kapz + sigz*dt);
                        c_hxb1_n_-> point(ii,jj) = (2*eps*kapx + sigxx*dt) / (2*eps*kapz + sigz*dt);

                        //Update Hy factors
                        eps    = objArr[phys_Hy_->point(ii,jj)].dielectric(1.0);
                        sigxy = sigmax(static_cast<double>(ii) + 0.5,eps);
                        sigyy = sigmay(static_cast<double>(ii),eps);
                        c_byb_0_ -> point(ii,jj) = (2*eps*kapz - sigz*dt) / (2*eps*kapz + sigz*dt);
                        c_bye_0_ -> point(ii,jj) = (2 * eps * dt) / (dy * (2*eps*kapz + sigz*dt));
                        c_hyh_0_ -> point(ii,jj) = (2*eps*kapx - sigxy*dt) / (2*eps*kapx + sigxy*dt);
                        c_hyb0_0_-> point(ii,jj) = (2*eps*kapy - sigyy*dt) / (2*eps*kapx + sigxy*dt);
                        c_hyb1_0_-> point(ii,jj) = (2*eps*kapy + sigyy*dt) / (2*eps*kapx + sigxy*dt);

                        eps    = objArr[phys_Hy_end_->point(ii,jj)].dielectric(1.0);
                        sigxy = sigmax(static_cast<double>(ii) - 0.5,eps);
                        sigyy = sigmay(static_cast<double>(ii),eps);
                        c_byb_n_ -> point(ii,jj) = (2*eps*kapz - sigz*dt) / (2*eps*kapz + sigz*dt);
                        c_bye_n_ -> point(ii,jj) = (2 * eps * dt) / (dy * (2*eps*kapz + sigz*dt));
                        c_hyh_n_ -> point(ii,jj) = (2*eps*kapx - sigxy*dt) / (2*eps*kapx + sigxy*dt);
                        c_hyb0_n_-> point(ii,jj) = (2*eps*kapy - sigyy*dt) / (2*eps*kapx + sigxy*dt);
                        c_hyb1_n_-> point(ii,jj) = (2*eps*kapy + sigyy*dt) / (2*eps*kapx + sigxy*dt);

                        //Update Ez factors
                        eps = objArr[phys_Ez_->point(ii,jj)].dielectric(1.0);
                        sigx = sigmax(static_cast<double>(ii),eps);
                        sigy = sigmay(static_cast<double>(ii),eps);
                        c_dzd_0_ -> point(ii,jj) = (2*eps*kapx - sigx*dt) / (2*eps*kapx + sigx*dt);
                        c_dzh_0_ -> point(ii,jj) = (2 * eps * dt) / (dy * (2*eps*kapx + sigx*dt));
                        c_eze_0_ -> point(ii,jj) = (2*eps*kapy - sigy*dt) / (2*eps*kapy + sigy*dt);
                        c_ezd1_0_-> point(ii,jj) = (2*eps*kapz + sigz*dt) / (2*eps*kapy + sigy*dt) / eps;
                        c_ezd0_0_-> point(ii,jj) = (2*eps*kapz - sigz*dt) / (2*eps*kapy + sigy*dt) / eps;

                        eps = objArr[phys_Ez_->point(ni-1-ii,jj)].dielectric(1.0);
                        sigx = sigmax(static_cast<double>(ii),eps);
                        sigy = sigmay(static_cast<double>(ii),eps);
                        c_dzd_n_ -> point(ii,jj) = (2*eps*kapx - sigx*dt) / (2*eps*kapx + sigx*dt);
                        c_dzh_n_ -> point(ii,jj) = (2 * eps * dt) / (dy * (2*eps*kapx + sigx*dt));
                        c_eze_n_ -> point(ii,jj) = (2*eps*kapy - sigy*dt) / (2*eps*kapy + sigy*dt);
                        c_ezd1_n_-> point(ii,jj) = (2*eps*kapz + sigz*dt) / (2*eps*kapy + sigy*dt) / eps;
                        c_ezd0_n_-> point(ii,jj) = (2*eps*kapz - sigz*dt) / (2*eps*kapy + sigy*dt) / eps;
                    }
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
        if(x <= thickness_ -1 && x >=0)
            return sigmaMax_ * sqrt(eps) * pow((static_cast<double>(thickness_-1) - x) / static_cast<double>(thickness_-1) , m_);
        else
            return 0.0;
    }

};

#endif