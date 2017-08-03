 #include <DTC/parallelStorageDTC.hpp>

parallelStorageDTCReal::parallelStorageDTCReal(int dtcNum, real_pgrid_ptr grid, std::array<int,3> loc, std::array<int,3> sz) :
    parallelStorageDTC(dtcNum, grid, loc, sz)
{}

void parallelStorageDTCReal::getField()
{
    grid_->gridComm().barrier();
    if(toOutGrid_)
    {
        for(int kk = 0; kk < toOutGrid_->opSz_[2]; ++kk )
        {
            for(int jj = 0; jj < toOutGrid_->opSz_[1]; ++jj)
            {
                dcopy_(toOutGrid_->opSz_[0], &grid_->point(toOutGrid_->loc_[0]+jj*toOutGrid_->addVec1_[0]+kk*toOutGrid_->addVec2_[0], toOutGrid_->loc_[1]+jj*toOutGrid_->addVec1_[1]+kk*toOutGrid_->addVec2_[1], toOutGrid_->loc_[2]+jj*toOutGrid_->addVec1_[2]+kk*toOutGrid_->addVec2_[2]), toOutGrid_->stride_, &outGrid_->point(toOutGrid_->locOutGrid_[0]+jj*toOutGrid_->addVec1_[0]+kk*toOutGrid_->addVec2_[0], toOutGrid_->locOutGrid_[1]+jj*toOutGrid_->addVec1_[1]+kk*toOutGrid_->addVec2_[1], toOutGrid_->locOutGrid_[2]+jj*toOutGrid_->addVec1_[2]+kk*toOutGrid_->addVec2_[2]), toOutGrid_->strideOutGrid_);
            }
        }
    }
    grid_->gridComm().barrier();
    if(slave_)
    {

        for(int kk = 0; kk < slave_->opSz_[2]; ++kk )
            for(int jj = 0; jj < slave_->opSz_[1]; ++jj)
                dcopy_(slave_->opSz_[0], &grid_->point(slave_->loc_[0]+jj*slave_->addVec1_[0]+kk*slave_->addVec2_[0], slave_->loc_[1]+jj*slave_->addVec1_[1]+kk*slave_->addVec2_[1], slave_->loc_[2]+jj*slave_->addVec1_[2]+kk*slave_->addVec2_[2]), slave_->stride_, &scratch_[ slave_->opSz_[0]*(jj + kk*slave_->opSz_[1]) ], 1);
        gridComm_.send(slave_->masterProc_, gridComm_.cantorTagGen(gridComm_.rank(), slave_->masterProc_, 1, 0), scratch_);
    }
    grid_->gridComm().barrier();
    if(masterBool_)
    {
        for(auto & slave : master_)
        {
            gridComm_.recv(slave->slaveProc_, gridComm_.cantorTagGen(slave->slaveProc_, gridComm_.rank(), 1, 0), scratch_);
            for(int kk = 0; kk < slave->sz_[2]; ++kk)
            {
                for(int jj = 0; jj < slave->sz_[1]; ++jj)
                {
                    dcopy_(slave->sz_[0], &scratch_[(jj + slave->sz_[1] * kk) * slave->sz_[0] ], 1, &outGrid_->point(slave->addVec1_[0]*jj+slave->addVec2_[0]*kk+slave->loc_[0], slave->addVec1_[1]*jj+slave->addVec2_[1]*kk+slave->loc_[1], slave->addVec1_[2]*jj+slave->addVec2_[2]*kk+slave->loc_[2]), slave->stride_);
                }
            }
        }
    }
    grid_->gridComm().barrier();
    return;
}

parallelStorageDTCCplx::parallelStorageDTCCplx(int dtcNum, cplx_pgrid_ptr grid, std::array<int,3> loc, std::array<int,3> sz) :
    parallelStorageDTC(dtcNum, grid, loc, sz)
{}

void parallelStorageDTCCplx::getField()
{
    grid_->gridComm().barrier();
    if(toOutGrid_)
    {
        for(int kk = 0; kk < toOutGrid_->opSz_[2]; ++kk )
        {
            for(int jj = 0; jj < toOutGrid_->opSz_[1]; ++jj)
            {
                zcopy_(toOutGrid_->opSz_[0], &grid_->point(toOutGrid_->loc_[0]+jj*toOutGrid_->addVec1_[0]+kk*toOutGrid_->addVec2_[0], toOutGrid_->loc_[1]+jj*toOutGrid_->addVec1_[1]+kk*toOutGrid_->addVec2_[1], toOutGrid_->loc_[2]+jj*toOutGrid_->addVec1_[2]+kk*toOutGrid_->addVec2_[2]), toOutGrid_->stride_, &outGrid_->point(toOutGrid_->locOutGrid_[0]+jj*toOutGrid_->addVec1_[0]+kk*toOutGrid_->addVec2_[0], toOutGrid_->locOutGrid_[1]+jj*toOutGrid_->addVec1_[1]+kk*toOutGrid_->addVec2_[1], toOutGrid_->locOutGrid_[2]+jj*toOutGrid_->addVec1_[2]+kk*toOutGrid_->addVec2_[2]), toOutGrid_->strideOutGrid_);
            }
        }
    }
    grid_->gridComm().barrier();
    if(slave_)
    {

        for(int kk = 0; kk < slave_->opSz_[2]; ++kk )
            for(int jj = 0; jj < slave_->opSz_[1]; ++jj)
                zcopy_(slave_->opSz_[0], &grid_->point(slave_->loc_[0]+jj*slave_->addVec1_[0]+kk*slave_->addVec2_[0], slave_->loc_[1]+jj*slave_->addVec1_[1]+kk*slave_->addVec2_[1], slave_->loc_[2]+jj*slave_->addVec1_[2]+kk*slave_->addVec2_[2]), slave_->stride_, &scratch_[ slave_->opSz_[0]*(jj + kk*slave_->opSz_[1]) ], 1);
        gridComm_.send(slave_->masterProc_, gridComm_.cantorTagGen(gridComm_.rank(), slave_->masterProc_, 1, 0), scratch_);
    }
    grid_->gridComm().barrier();
    if(masterBool_)
    {
        for(auto & slave : master_)
        {
            gridComm_.recv(slave->slaveProc_, gridComm_.cantorTagGen(slave->slaveProc_, gridComm_.rank(), 1, 0), scratch_);
            for(int kk = 0; kk < slave->sz_[2]; ++kk)
            {
                for(int jj = 0; jj < slave->sz_[1]; ++jj)
                {
                    zcopy_(slave->sz_[0], &scratch_[(jj + slave->sz_[1] * kk) * slave->sz_[0] ], 1, &outGrid_->point(slave->addVec1_[0]*jj+slave->addVec2_[0]*kk+slave->loc_[0], slave->addVec1_[1]*jj+slave->addVec2_[1]*kk+slave->loc_[1], slave->addVec1_[2]*jj+slave->addVec2_[2]*kk+slave->loc_[2]), slave->stride_);
                }
            }
        }
    }
    grid_->gridComm().barrier();
    return;
}