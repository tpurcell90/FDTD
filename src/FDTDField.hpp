#ifndef FDTD_FDTDFIELD
#define FDTD_FDTDFIELD

#include "Grid.hpp"
#include "Source.hpp"
#include "Inputs.hpp"

// #include <assert.h>
// #include <iomanip>
// #include <iostream>
#include <memory>
// #include <random>
// #include <stdexcept>
// #include <string>
#include <vector>

// #include <complex>
// typedef std::complex<double> cplx;

class FDTDField
{
protected:
    size_t nx,ny,nt;
    double dx,dy,dt,res;
public:
    std::shared_ptr<Grid2D<double>> Ex,Ey,Ez,Hx,Hy,Hz,physGrid;
    std::vector<Source<double>> srcArr;

    FDTDField(programInputs *IP = NULL);
    void initializeGrid(programInputs *IP = NULL);
    void ouputField();
    void step();
};

#endif
