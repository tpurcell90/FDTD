#include <SOURCE/parallelSourceNormal.hpp>
parallelSourceNormalReal::parallelSourceNormalReal(mpiInterface comm, std::vector<std::shared_ptr<PulseBase>> pulse, real_pgrid_ptr grid, double dt, std::array<int,3> loc, std::array<int,3> sz) :
    parallelSourceNormalBase<double>(comm, pulse, grid, dt, loc, sz)
{}

void parallelSourceNormalReal::addPul(double t)
{
    if(slave_)
    {
        cplx pulVal = 0.0;
        for(auto& pul : pulse_)
            pulVal += pul->pulse(t);
        // std::cout << t << '\t' << pulVal << std::endl;
        std::fill_n(pulVec_.data(), slave_->sz_[0], std::real(pulVal));
        for(int kk = 0; kk < slave_->sz_[2]; ++kk)
        {
            for(int jj = 0; jj < slave_->sz_[1]; ++jj)
            {
                daxpy_(slave_->sz_[0], dt_, pulVec_.data(), 1, &grid_->point(slave_->loc_[0]+jj*slave_->addVec1_[0]+kk*slave_->addVec2_[0],   slave_->loc_[1]+jj*slave_->addVec1_[1]+kk*slave_->addVec2_[1],   slave_->loc_[2]+jj*slave_->addVec1_[2]+kk*slave_->addVec2_[2]), slave_->stride_);
            }
        }
    }
    grid_->transferDat();
}

parallelSourceNormalCplx::parallelSourceNormalCplx(mpiInterface comm, std::vector<std::shared_ptr<PulseBase>> pulse, cplx_pgrid_ptr grid, double dt, std::array<int,3> loc, std::array<int,3> sz) :
    parallelSourceNormalBase<cplx>(comm, pulse, grid, dt, loc, sz)
{}

void parallelSourceNormalCplx::addPul(double t)
{
    if(slave_)
    {
        cplx pulVal = 0.0;
        for(auto& pul : pulse_)
            pulVal += pul->pulse(t);

        std::fill_n(pulVec_.data(), slave_->sz_[0], pulVal);
        for(int kk = 0; kk < slave_->sz_[2]; ++kk)
        {
            for(int jj = 0; jj < slave_->sz_[1]; ++jj)
            {
                zaxpy_(slave_->sz_[0], dt_, pulVec_.data(), 1, &grid_->point(slave_->loc_[0]+jj*slave_->addVec1_[0]+kk*slave_->addVec2_[0],   slave_->loc_[1]+jj*slave_->addVec1_[1]+kk*slave_->addVec2_[1],   slave_->loc_[2]+jj*slave_->addVec1_[2]+kk*slave_->addVec2_[2]), slave_->stride_);
            }
        }
    }
    grid_->transferDat();
}