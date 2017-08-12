#include <DTC/parallelStorageFreqDTC.hpp>

parallelStorageFreqDTCReal::parallelStorageFreqDTCReal(int dtcNum, real_pgrid_ptr grid, DIRECTION propDir, std::array<int,3> loc, std::array<int,3> sz, std::vector<double> freqList) :
    parallelStorageFreqDTC(dtcNum, grid, propDir, loc, sz, freqList)
{}


void parallelStorageFreqDTCReal::fieldIn(cplx* fftFact)
{
    if(!fieldInFreq_)
        return;
    // Copy the field information into a vector
    for(int jj = 0; jj < fieldInFreq_->sz_[2]; ++jj)
    {
        for(int ii = 0; ii < fieldInFreq_->sz_[1]; ++ii)
        {
            dcopy_(fieldInFreq_->sz_[0], &grid_->point(fieldInFreq_->loc_[0]+ii*fieldInFreq_->addVec1_[0]+jj*fieldInFreq_->addVec2_[0], fieldInFreq_->loc_[1]+ii*fieldInFreq_->addVec1_[1]+jj*fieldInFreq_->addVec2_[1], fieldInFreq_->loc_[2]+ii*fieldInFreq_->addVec1_[2]+jj*fieldInFreq_->addVec2_[2]), fieldInFreq_->stride_, reinterpret_cast<double*>( &fIn_[ (ii*fieldInFreq_->sz_[2] + jj)*fieldInFreq_->sz_[0] ] ) , 2 );
        }
    }
    // Take an outer product of the prefactor vector and the field vectors to get the discrete Fourier Transform at all points
    zgerc_(nfreq_, fieldInFreq_->sz_[0]*fieldInFreq_->sz_[1]*fieldInFreq_->sz_[2], 1.0, fftFact, 1, fIn_.data(), 1, outGrid_->data(), nfreq_);
}

parallelStorageFreqDTCCplx::parallelStorageFreqDTCCplx(int dtcNum, cplx_pgrid_ptr grid, DIRECTION propDir, std::array<int,3> loc, std::array<int,3> sz, std::vector<double> freqList) :
    parallelStorageFreqDTC(dtcNum, grid, propDir, loc, sz, freqList)
{}


void parallelStorageFreqDTCCplx::fieldIn(cplx* fftFact)
{
    if(!fieldInFreq_)
        return;
    // Copy the field information into a vector
    for(int jj = 0; jj < fieldInFreq_->sz_[2]; ++jj)
    {
        for(int ii = 0; ii < fieldInFreq_->sz_[1]; ++ii)
        {
            zcopy_(fieldInFreq_->sz_[0], &grid_->point(fieldInFreq_->loc_[0]+ii*fieldInFreq_->addVec1_[0]+jj*fieldInFreq_->addVec2_[0], fieldInFreq_->loc_[1]+ii*fieldInFreq_->addVec1_[1]+jj*fieldInFreq_->addVec2_[1], fieldInFreq_->loc_[2]+ii*fieldInFreq_->addVec1_[2]+jj*fieldInFreq_->addVec2_[2]), fieldInFreq_->stride_, &fIn_[ (ii*fieldInFreq_->sz_[2] + jj)*fieldInFreq_->sz_[0]]  , 1 );
        }
    }
    // Take an outer product of the prefactor vector and the field vectors to get the discrete Fourier Transform at all points
    zgerc_(nfreq_, fieldInFreq_->sz_[0]*fieldInFreq_->sz_[1]*fieldInFreq_->sz_[2], ONE_, fftFact, 1, fIn_.data(), 1, outGrid_->data(), nfreq_);
}


