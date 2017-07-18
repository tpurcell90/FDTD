#include <DTC/parallelDTC_FREQ.hpp>

parallelDetectorFREQReal::parallelDetectorFREQReal(std::string name, std::vector<pgrid_ptr> grids, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, DTCCLASSTYPE classType, int timeInt, int nfreq, double fcen, double fwidth, std::array<double,3> d, double dt, bool SI, double I0, double a) :
    parallelDetectorFREQ_Base(name, grids, loc, sz, type, classType, timeInt, nfreq, fcen, fwidth, d, dt, SI, I0, a)
{
    getIncdField_ = [](cplx a){return std::real(a); };
    for(auto& grid : grids)
    {
        gridsIn_.push_back(std::make_shared<parallelStorageFreqDTCReal>(masterProc_, grid, DIRECTION::NONE, loc_, sz_, freqList_) );
        if(gridComm_.rank() == masterProc_)
            freqFields_.push_back(std::make_shared<Grid<cplx>>(std::array<int,3>({{nfreq_, std::accumulate(sz_.begin(), sz_.end(), 1, std::multiplies<int>() ), 1 }}), std::array<double,3>({{dOmg_, d_[0], 1.0}}) ) );
    }
}

parallelDetectorFREQReal::parallelDetectorFREQReal(std::string name, std::vector<pgrid_ptr> grids, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, DTCCLASSTYPE classType, int timeInt, double lamL, double lamR, int nLam, std::array<double,3> d, double dt, bool SI, double I0, double a) :
    parallelDetectorFREQ_Base(name, grids, loc, sz, type, classType, timeInt, lamL, lamR, nLam, d, dt, SI, I0, a)
{
    getIncdField_ = [](cplx a){return std::real(a); };
    for(auto& grid : grids)
    {
        gridsIn_.push_back(std::make_shared<parallelStorageFreqDTCReal>(masterProc_, grid, DIRECTION::NONE, loc_, sz_, freqList_) );
        if(gridComm_.rank() == masterProc_)
            freqFields_.push_back(std::make_shared<Grid<cplx>>(std::array<int,3>({{nfreq_, std::accumulate(sz_.begin(), sz_.end(), 1, std::multiplies<int>() ), 1 }}), std::array<double,3>({{dLam_, d_[0], 1.0}}) ) );
    }
}

parallelDetectorFREQCplx::parallelDetectorFREQCplx(std::string name, std::vector<pgrid_ptr> grids, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, DTCCLASSTYPE classType, int timeInt, int nfreq, double fcen, double fwidth, std::array<double,3> d, double dt, bool SI, double I0, double a) :
    parallelDetectorFREQ_Base(name, grids, loc, sz, type, classType, timeInt, nfreq, fcen, fwidth, d, dt, SI, I0, a)
{
    getIncdField_ = [](cplx a){return a; };
    for(auto& grid : grids)
    {
        gridsIn_.push_back(std::make_shared<parallelStorageFreqDTCCplx>(masterProc_, grid, DIRECTION::NONE, loc_, sz_, freqList_) );
        if(gridComm_.rank() == masterProc_)
            freqFields_.push_back(std::make_shared<Grid<cplx>>(std::array<int,3>({{nfreq_, std::accumulate(sz_.begin(), sz_.end(), 1, std::multiplies<int>() ), 1 }}), std::array<double,3>({{dOmg_, d_[0], 1.0}}) ) );
    }
}

parallelDetectorFREQCplx::parallelDetectorFREQCplx(std::string name, std::vector<pgrid_ptr> grids, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, DTCCLASSTYPE classType, int timeInt, double lamL, double lamR, int nLam, std::array<double,3> d, double dt, bool SI, double I0, double a) :
    parallelDetectorFREQ_Base(name, grids, loc, sz, type, classType, timeInt, lamL, lamR, nLam, d, dt, SI, I0, a)
{
    getIncdField_ = [](cplx a){return a; };
    for(auto& grid : grids)
    {
        gridsIn_.push_back(std::make_shared<parallelStorageFreqDTCCplx>(masterProc_, grid, DIRECTION::NONE, loc_, sz_, freqList_) );
        if(gridComm_.rank() == masterProc_)
            freqFields_.push_back(std::make_shared<Grid<cplx>>(std::array<int,3>({{nfreq_, std::accumulate(sz_.begin(), sz_.end(), 1, std::multiplies<int>() ), 1 }}), std::array<double,3>({{dLam_, d_[0], 1.0}}) ) );
    }
}