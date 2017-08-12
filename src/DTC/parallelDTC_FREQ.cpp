#include <DTC/parallelDTC_FREQ.hpp>

parallelDetectorFREQReal::parallelDetectorFREQReal(std::string name, std::vector<pgrid_ptr> grids, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, int timeInt, std::vector<double> freqList, std::array<double,3> d, double dt, bool SI, double I0, double a) :
    parallelDetectorFREQ_Base(name, grids, loc, sz, type, timeInt, freqList, d, dt, SI, I0, a)
{
    // Se the getIncdField functions for both rea/complex values
    getIncdField_ = [](cplx a){return std::real(a); };
    for(auto& grid : grids)
    {
        // Construct the Storage fields
        gridsIn_.push_back(std::make_shared<parallelStorageFreqDTCReal>(masterProc_, grid, DIRECTION::NONE, loc_, sz_, freqList_) );
        // if Master make the frequency field masters
        if(gridComm_->rank() == masterProc_)
            freqFields_.push_back(std::make_shared<Grid<cplx>>(std::array<int,3>({{nfreq_, std::accumulate(sz_.begin(), sz_.end(), 1, std::multiplies<int>() ), 1 }}), std::array<double,3>({{dOmg_, d_[0], 1.0}}) ) );
    }
}

parallelDetectorFREQCplx::parallelDetectorFREQCplx(std::string name, std::vector<pgrid_ptr> grids, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, int timeInt, std::vector<double> freqList, std::array<double,3> d, double dt, bool SI, double I0, double a) :
    parallelDetectorFREQ_Base(name, grids, loc, sz, type, timeInt, freqList, d, dt, SI, I0, a)
{
    // Se the getIncdField functions for both rea/complex values
    getIncdField_ = [](cplx a){return a; };
    for(auto& grid : grids)
    {
        // Construct the Storage fields
        gridsIn_.push_back(std::make_shared<parallelStorageFreqDTCCplx>(masterProc_, grid, DIRECTION::NONE, loc_, sz_, freqList_) );
        // if Master make the frequency field masters
        if(gridComm_->rank() == masterProc_)
            freqFields_.push_back(std::make_shared<Grid<cplx>>(std::array<int,3>({{nfreq_, std::accumulate(sz_.begin(), sz_.end(), 1, std::multiplies<int>() ), 1 }}), std::array<double,3>({{dOmg_, d_[0], 1.0}}) ) );
    }
}