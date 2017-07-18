#include <DTC/parallelDTCOutputFxn.hpp>
cplx fieldOutFreqFunction(int szFreq, cplx* gridInBegin, int szPwr, cplx* pwr, int nt, double convFactor)
{
    return convFactor * zdotc_(szFreq, std::vector<cplx>(szFreq, 1.0).data(), 1, gridInBegin, 1) / std::pow(static_cast<double>(nt), 1.0);
}

cplx pwrOutFreqFunction(int szFreq, cplx* gridInBegin, int szPwr, cplx* pwr, int nt, double convFactor)
{
    std::transform(gridInBegin, gridInBegin+szPwr, pwr, [](cplx field){ return field * std::conj(field) ; } );
    return convFactor * zdotc_(szFreq, std::vector<cplx>(szFreq, 1.0).data(), 1, pwr, 1) / std::pow(static_cast<double>(nt), 2.0);
}