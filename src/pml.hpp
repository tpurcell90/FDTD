#ifndef FDTD_PML
#define FDTD_PML

#include "Source.hpp"
#include "Grid.hpp"
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

public:
    std::shared_ptr<Grid2D<T>> Dx_,Dy_,Dz_,Bx_,By_,Bz_,Dx_end_,Dy_end_,Dz_end_,Bx_end_,By_end_,Bz_end_;
    std::shared_ptr<Grid2D<int>> phys_Hx_,phys_Hy_,phys_Hz_, phys_Hx_end_,phys_Hy_end_,phys_Hz_end_;

    std::shared_ptr<Grid2D<double>> c_bxb_0_, c_bxe_0_, c_hxh_0_, c_hxb0_0_, c_hxb1_0_,c_bxb_n_, c_bxe_n_, c_hxh_n_, c_hxb0_n_, c_hxb1_n_;
    std::shared_ptr<Grid2D<double>> c_byb_0_, c_bye_0_, c_hyh_0_, c_hyb0_0_, c_hyb1_0_,c_byb_n_, c_bye_n_, c_hyh_n_, c_hyb0_n_, c_hyb1_n_;
    std::shared_ptr<Grid2D<double>> c_bzb_0_, c_bze_0_, c_hzh_0_, c_hzb0_0_, c_hzb1_0_,c_bzb_n_, c_bze_n_, c_hzh_n_, c_hzb0_n_, c_hzb1_n_;
    std::shared_ptr<Grid2D<double>> c_dxd_0_, c_dxh_0_, c_exe_0_, c_exd0_0_, c_exd1_0_,c_dxd_n_, c_dxh_n_, c_exe_n_, c_exd0_n_, c_exd1_n_;
    std::shared_ptr<Grid2D<double>> c_dyd_0_, c_dyh_0_, c_eye_0_, c_eyd0_0_, c_eyd1_0_,c_dyd_n_, c_dyh_n_, c_eye_n_, c_eyd0_n_, c_eyd1_n_;
    std::shared_ptr<Grid2D<double>> c_dzd_0_, c_dzh_0_, c_eze_0_, c_ezd0_0_, c_ezd1_0_,c_dzd_n_, c_dzh_n_, c_eze_n_, c_ezd0_n_, c_ezd1_n_;
    std::vector<std::array<int,4>> zaxHxList_;
    std::vector<std::array<int,4>> zaxHyList_;
    std::vector<std::array<int,4>> zaxEzList_;
    std::vector<std::array<int,4>> zaxHxList_end_;
    std::vector<std::array<int,4>> zaxHyList_end_;
    std::vector<std::array<int,4>> zaxEzList_end_;

    int Hy0EdgeInd_, Hx0EdgeInd_, Ez0EdgeInd_, HynEdgeInd_, HxnEdgeInd_, EznEdgeInd_;
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
    UPML(int thickness, Direction d, double m, double R0, int nx, double dx, double dy, Polarization pol, bool precalc) : thickness_(thickness), d_(d), m_(m), R0_(R0)
    {
        sigmaMax_ = -(m_+1)*log(R0_)/(2*thickness_*dx); // eta should be included;
        kappaMax_ = 1.0;
        zaxEzList_  = {};
        zaxHxList_  = {};
        zaxHyList_  = {};
        Hy0EdgeInd_ = 0; Hx0EdgeInd_ = 0; Ez0EdgeInd_ = 0; HynEdgeInd_ = 0; HxnEdgeInd_ = 0; EznEdgeInd_ = 0;
        if(d == X)
        {
            if(pol == EX || pol == EY || pol == HZ)
            {
                Dx_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
                Dy_ = std::make_shared<Grid2D<T>>(thickness,nx-1,dx,dy);
                Bz_ = std::make_shared<Grid2D<T>>(thickness,nx-1,dx,dy);
                Dx_end_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
                Dy_end_ = std::make_shared<Grid2D<T>>(thickness,nx-1,dx,dy);
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

                if(precalc == false)
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
                Dy_ = nullptr;
                Bz_ = nullptr;

                Bx_end_ = std::make_shared<Grid2D<T>>(thickness,nx-1,dx,dy);
                By_end_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
                Dz_end_ = std::make_shared<Grid2D<T>>(thickness,nx,dx,dy);
                Dx_end_ = nullptr;
                Dy_end_ = nullptr;
                Bz_end_ = nullptr;

                phys_Hz_ = nullptr;
                phys_Hz_end_ = nullptr;
                phys_Hy_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
                phys_Hy_end_ = std::make_shared<Grid2D<int>>(thickness,nx,dx,dy);
                phys_Hx_ = std::make_shared<Grid2D<int>>(thickness,nx-1,dx,dy);
                phys_Hx_end_ = std::make_shared<Grid2D<int>>(thickness,nx-1,dx,dy);
                if(precalc == false)
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
                Dy_ = std::make_shared<Grid2D<T>>(nx,thickness,dx,dy);
                Bz_ = std::make_shared<Grid2D<T>>(nx-1,thickness,dx,dy);
                Bx_ = nullptr;
                By_ = nullptr;
                Dz_ = nullptr;
                Dx_end_ = std::make_shared<Grid2D<T>>(nx-1,thickness,dx,dy);
                Dy_end_ = std::make_shared<Grid2D<T>>(nx,thickness,dx,dy);
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
                if(precalc == true)
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
                Dy_ = nullptr;
                Bz_ = nullptr;
                Bx_end_ = std::make_shared<Grid2D<T>>(nx,thickness,dx,dy);
                By_end_ = std::make_shared<Grid2D<T>>(nx-1,thickness,dx,dy);
                Dz_end_ = std::make_shared<Grid2D<T>>(nx,thickness,dx,dy);
                Dx_end_ = nullptr;
                Dy_end_ = nullptr;
                Bz_end_ = nullptr;
                phys_Hz_ = nullptr;
                phys_Hz_end_ = nullptr;
                phys_Hy_ = std::make_shared<Grid2D<int>>(nx-1,thickness,dx,dy);
                phys_Hy_end_ = std::make_shared<Grid2D<int>>(nx-1,thickness,dx,dy);
                phys_Hx_ = std::make_shared<Grid2D<int>>(nx,thickness,dx,dy);
                phys_Hx_end_ = std::make_shared<Grid2D<int>>(nx,thickness,dx,dy);
                if(precalc == false)
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
        if(x <= thickness_ -1)
            return sigmaMax_ * sqrt(eps) * pow((static_cast<double>(thickness_-1) - x) / static_cast<double>(thickness_-1) , m_);
        else
            return 0.0;
    }

};

#endif