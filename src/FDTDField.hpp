#ifndef FDTD_FDTDFIELD
#define FDTD_FDTDFIELD

#include "enum.hpp"
#include "Detector.hpp"
#include "Grid.hpp"
#include "Inputs.hpp"
#include "pml.hpp"
// #include <assert.h>
// #include <iomanip>
// #include <iostream>
#include <memory>
// #include <random>
// #include <stdexcept>
// #include <string>
#include <vector>
#include <map>

// #include <complex>
// typedef std::complex<double> cplx;

class FDTDField
{
protected:
    size_t nx_,ny_;
    double dx_,dy_,dt_, tcur_;
    int res_, xPML_, yPML_, t_step_;
    std::shared_ptr<Grid2D<int>> phys_Ex_,phys_Ey_,phys_Ez_;//phys_Hx_,phys_Hy_,phys_Hz_;
    std::vector<Source<complex<double>>> srcArr_;
    std::vector<Obj> objArr_;
    std::vector<Detector<complex<double>>> dtcArr_;
    std::vector<UPML<complex<double>>> pmlArr_;
    bool periodic_;
    bool precalcPML_;
    std::vector<std::array<int,4>> zaxEzList_;
    int y0EdgeInd_, x0EdgeInd_, ynEdgeInd_, xnEdgeInd_;
    std::vector<double> k_point_;


public:
    std::shared_ptr<Grid2D<complex<double>>> Ex_,Ey_,Ez_,Hx_,Hy_,Hz_;
    FDTDField(programInputs &IP);
    void initializeGrid();

    // Access Functions
    size_t nx(){return nx_;}
    size_t ny(){return ny_;}
    double dx(){return dx_;}
    double dy(){return dy_;}
    double dt(){return dt_;}
    double getTime(){return tcur_;}
    Obj getObj(int n){return objArr_[n];}
    std::vector<Detector<complex<double>>> getDtcArr(){return dtcArr_;}
    int getRes(){return res_;}
    std::shared_ptr<Grid2D<int>> getPhysEz(){return phys_Ez_;}
    std::shared_ptr<Grid2D<int>> getPhysEx(){return phys_Ex_;}
    std::shared_ptr<Grid2D<int>> getPhysEy(){return phys_Ey_;}
    //std::shared_ptr<Grid2D<int>> getPhysHy(){return phys_Hy_;}
    //std::shared_ptr<Grid2D<int>> getPhysHx(){return phys_Hx_;}
    //std::shared_ptr<Grid2D<int>> getPhysHz(){return phys_Hz_;}



    void ouputField(Detector<complex<double>> d);
    void step();
    void updateH();
    void updateE();
    std::vector<double> pml(int npml, int m, int ma);
    std::complex<double> per_factor(std::vector<double> r);
};

#endif
