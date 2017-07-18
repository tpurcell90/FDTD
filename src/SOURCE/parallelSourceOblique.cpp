#include <SOURCE/parallelSourceOblique.hpp>


parallelSourceObliqueReal::parallelSourceObliqueReal(mpiInterface comm, std::vector<std::shared_ptr<PulseBase>> pulse, real_pgrid_ptr grid, POLARIZATION pol, double dt, std::array<int,3> loc, std::array<int,3> sz, double phi, double theta) :
    parallelSourceObliqueBase<double>(comm, pulse, grid, pol, dt, loc, sz, phi, theta)
{}

void parallelSourceObliqueReal::addPul(double t)
{
    for(auto& pul : pulse_)
        for(auto& param : updateSrcParams_)
            *param.loc_ += param.scalefact_ * std::real(pul->pulse(t - param.t_off_) );
    grid_->transferDat();

}
parallelSourceObliqueCplx::parallelSourceObliqueCplx(mpiInterface comm, std::vector<std::shared_ptr<PulseBase>> pulse, cplx_pgrid_ptr grid, POLARIZATION pol, double dt, std::array<int,3> loc, std::array<int,3> sz, double phi, double theta) :
    parallelSourceObliqueBase<cplx>(comm, pulse, grid, pol, dt, loc, sz, phi, theta)
{}

void parallelSourceObliqueCplx::addPul(double t)
{
    for(auto& pul : pulse_)
        for(auto& param : updateSrcParams_)
            *param.loc_ += param.scalefact_ * pul->pulse(t - param.t_off_);
    grid_->transferDat();
}
