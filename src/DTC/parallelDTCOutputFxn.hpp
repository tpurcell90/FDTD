#ifndef FDTD_PARALLELDETECTORFXNS
#define FDTD_PARALLELDETECTORFXNS

#include <cmath>
#ifdef MKL
#include <UTIL/utilities_MKL.hpp>
#elif defined ACML
#include <UTIL/utilities_acml.hpp>
#else
#include <UTIL/utilities_MKL.hpp>
#endif
#include <UTIL/typedefs.hpp>

template <class T> void pwrOutputFunction(T* gridInBegin, T* gridInEnd, T* valsOut, double convFactor)
{
    std::transform(gridInBegin, gridInEnd, valsOut, valsOut, [&](T pt, T curVal){return curVal + convFactor*pow(std::abs( pt ), 2.0) ; } );
};

template <class T> void fieldOutputFunction(T* gridInBegin, T* gridInEnd, T* valsOut, double convFactor)
{
    std::transform(gridInBegin, gridInEnd, valsOut, [&](T pt){return convFactor*pt;} );
};

cplx fieldOutFreqFunction(int szFreq, cplx* gridInBegin, int szPwr, cplx* pwr, int nt, double convFactor);

cplx pwrOutFreqFunction(int szFreq, cplx* gridInBegin, int szPwr, cplx* pwr, int nt, double convFactor);

#endif