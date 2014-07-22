#ifndef FDTD_FDTDFIELD
#define FDTD_FDTDFIELD

#include "enum.hpp"
#include "Detector.hpp"
#include "Grid.hpp"
#include "Inputs.hpp"
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
    size_t nx,ny;
    double dx,dy,dt, t_cur;
    int res;
    std::shared_ptr<Grid2D<int>> phys_Ex,phys_Ey,phys_Ez,phys_Hx,phys_Hy,phys_Hz;
    std::vector<Source<double>> srcArr;
    std::vector<Obj> objArr;
    std::vector<Detector<double>> dtcArr;

public:
    std::shared_ptr<Grid2D<double>> Ex,Ey,Ez,Hx,Hy,Hz;
    FDTDField(programInputs &IP);
    void initializeGrid(programInputs &IP);

    // Access Functions
    size_t n_x();
    size_t n_y();
    double d_x();
    double d_y();
    double d_t();
    double getTime();
    std::vector<Detector<double>> getDtcArr();
    int getRes();
    std::shared_ptr<Grid2D<int>> getPhysEz();    

    void ouputField(Detector<double> d);
    void step();
    void inc_t();
    Obj makeSphere(vector<double> mater, double rad, vector<double> loc);
    Obj makeBlock(vector<double> mater, vector<double> dims, vector<double> loc);
    std::vector<double> pml(int npml, int m, int ma);
};

#endif
