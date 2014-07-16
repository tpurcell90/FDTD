// How do I cite the Matrix class from Josh's Split Op code
#include <algorithm>
#include <assert.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <complex>
#include "Grid.hpp"
#include "Source.hpp"

typedef std::complex<double> cplx;

template <typename T> class FDTDField
{
protected:
    static unsigned int memSize;
    size_t nx,ny,nt;
    double dx,dy,dt,res;

public:
    Grid2D<double> Ex,Ey,Ez,Hx,Hy,Hz,physGrid;
    Source srcArr[];     
    FDTDField(programInputs &IP);
    void initializeGrid(programInputs &IP);
    void ouputField();
    void step()
    
}
