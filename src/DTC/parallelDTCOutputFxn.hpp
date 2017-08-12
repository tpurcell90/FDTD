#ifndef FDTD_PARALLELDETECTORFXNS
#define FDTD_PARALLELDETECTORFXNS
#include <UTIL/typedefs.hpp>

/**
 * @brief      Function used to output the power/intensity of a field as a function of time
 *
 * @param      gridInBegin  pointer the the first element of the field in array
 * @param      gridInEnd    pointer to the last element of the field in array
 * @param      valsOut      pointer to the output array
 * @param[in]  convFactor   The unit conversion factor
 *
 * @tparam     T            double or complex<double>
 */
template <class T> void pwrOutputFunction(T* gridInBegin, T* gridInEnd, T* valsOut, double convFactor)
{
    std::transform(gridInBegin, gridInEnd, valsOut, valsOut, [&](T pt, T curVal){return curVal + convFactor*pow(std::abs( pt ), 2.0) ; } );
};

/**
 * @brief      Function used to output a filed value as a function of time
 *
 * @param      gridInBegin  pointer the the first element of the field in array
 * @param      gridInEnd    pointer to the last element of the field in array
 * @param      valsOut      pointer to the output array
 * @param[in]  convFactor   The unit conversion factor
 *
 * @tparam     T            double or complex<double>
 */
template <class T> void fieldOutputFunction(T* gridInBegin, T* gridInEnd, T* valsOut, double convFactor)
{
    std::transform(gridInBegin, gridInEnd, valsOut, [&](T pt){return convFactor*pt;} );
};

/**
 * @brief      Function used to calculate the Fourier transformed field values
 *
 * @param[in]  szFreq       The size of the region the frequency detector is acting on
 * @param      gridInBegin  pointer the the first element of the field in array
 * @param[in]  szPwr        size of the vector storing the power once calculated (neglected for field output)
 * @param      pwr          The power temporary storage vector (neglected for field output)
 * @param[in]  nt           number of time steps
 * @param[in]  convFactor   The unit conversion factor
 *
 * @return     The Fourier transformed field value for a particular frequency
 */
cplx fieldOutFreqFunction(int szFreq, cplx* gridInBegin, int szPwr, cplx* pwr, int nt, double convFactor);

/**
 * @brief      Function used to calculate the Fourier transformed field power/intensity
 *
 * @param[in]  szFreq       The size of the region the frequency detector is acting on
 * @param      gridInBegin  pointer the the first element of the field in array
 * @param[in]  szPwr        size of the vector storing the power once calculated
 * @param      pwr          The power temporary storage vector
 * @param[in]  nt           number of time steps
 * @param[in]  convFactor   The unit conversion factor
 *
 * @return     The Fourier transformed field power/intensity for a particular frequency
 */
cplx pwrOutFreqFunction(int szFreq, cplx* gridInBegin, int szPwr, cplx* pwr, int nt, double convFactor);

#endif