#include <PML/parallelPML.hpp>

void pmlUpdateFxnReal::addPsi(std::vector<updateGridParams> &gridParamList, std::vector<updatePsiParams> &psiParamList, real_pgrid_ptr grid_i, real_pgrid_ptr psi, real_pgrid_ptr grid)
{
    updatePsiField(psiParamList, psi, grid);
    for(auto & param : gridParamList)
        daxpy_(param.nAx_, param.Db_, &psi->point(param.loc_[0], param.loc_[1], param.loc_[2]), param.stride_, &grid_i->point(param.loc_[0], param.loc_[1], param.loc_[2] ), param.stride_);
}

void pmlUpdateFxnReal::updatePsiField(std::vector<updatePsiParams> &paramList, real_pgrid_ptr psi, real_pgrid_ptr grid)
{
    for (auto & param : paramList)
    {
        dscal_(param.transSz_, param.b_   ,  &psi->point(param.loc_[0]   , param.loc_[1]   , param.loc_[2]   ), param.stride_);
        daxpy_(param.transSz_, param.c_   , &grid->point(param.loc_[0]   , param.loc_[1]   , param.loc_[2]   ), param.stride_, &psi->point(param.loc_[0], param.loc_[1], param.loc_[2]), param.stride_);
        daxpy_(param.transSz_, param.cOff_, &grid->point(param.locOff_[0], param.locOff_[1], param.locOff_[2]), param.stride_, &psi->point(param.loc_[0], param.loc_[1], param.loc_[2]), param.stride_);
    }
}

void pmlUpdateFxnCplx::addPsi(std::vector<updateGridParams> &gridParamList, std::vector<updatePsiParams> &psiParamList, cplx_pgrid_ptr grid_i, cplx_pgrid_ptr psi, cplx_pgrid_ptr grid)
{
    updatePsiField(psiParamList, psi, grid);
    for(auto & param : gridParamList)
        zaxpy_(param.nAx_, param.Db_, &psi->point(param.loc_[0], param.loc_[1], param.loc_[2]), param.stride_, &grid_i->point(param.loc_[0], param.loc_[1], param.loc_[2] ), param.stride_);
}

void pmlUpdateFxnCplx::updatePsiField(std::vector<updatePsiParams> &paramList, cplx_pgrid_ptr psi, cplx_pgrid_ptr grid)
{
    for (auto & param : paramList)
    {
        zscal_(param.transSz_, param.b_   ,  &psi->point(param.loc_[0]   , param.loc_[1]   , param.loc_[2]   ), param.stride_);
        zaxpy_(param.transSz_, param.c_   , &grid->point(param.loc_[0]   , param.loc_[1]   , param.loc_[2]   ), param.stride_, &psi->point(param.loc_[0], param.loc_[1], param.loc_[2]), param.stride_);
        zaxpy_(param.transSz_, param.cOff_, &grid->point(param.locOff_[0], param.locOff_[1], param.locOff_[2]), param.stride_, &psi->point(param.loc_[0], param.loc_[1], param.loc_[2]), param.stride_);
    }
}

parallelCPMLReal::parallelCPMLReal(std::shared_ptr<mpiInterface> gridComm, std::vector<real_grid_ptr> weights, real_pgrid_ptr grid_i, real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, POLARIZATION pol_i, std::array<int,3> n_vec, double m, double ma, double aMax, std::array<double,3> d, double dt, int_pgrid_ptr physGrid, std::vector<std::shared_ptr<Obj>> objArr) :
    parallelCPML<double>(gridComm, weights, grid_i, grid_j, grid_k, pol_i, n_vec, m, ma, aMax, d, dt, physGrid, objArr)
{
    if(psi_j_)
        upPsi_j_ = pmlUpdateFxnReal::addPsi;
    else
    {
        upPsi_j_ = [](std::vector<updateGridParams>&, std::vector<updatePsiParams>&, real_pgrid_ptr, real_pgrid_ptr, real_pgrid_ptr){return;};
    }
    if(psi_k_)
        upPsi_k_ = pmlUpdateFxnReal::addPsi;
    else
    {
        upPsi_k_ = [](std::vector<updateGridParams>&, std::vector<updatePsiParams>&, real_pgrid_ptr, real_pgrid_ptr, real_pgrid_ptr){return;};
    }
}
parallelCPMLCplx::parallelCPMLCplx(std::shared_ptr<mpiInterface> gridComm, std::vector<real_grid_ptr> weights, std::shared_ptr<parallelGrid<cplx > > grid_i, std::shared_ptr<parallelGrid<cplx > > grid_j, std::shared_ptr<parallelGrid<cplx > > grid_k, POLARIZATION pol_i, std::array<int,3> n_vec, double m, double ma, double aMax, std::array<double,3> d, double dt, int_pgrid_ptr physGrid, std::vector<std::shared_ptr<Obj>> objArr) :
    parallelCPML<cplx>(gridComm, weights, grid_i, grid_j, grid_k, pol_i, n_vec, m, ma, aMax, d, dt, physGrid, objArr)
{
    if(psi_j_)
        upPsi_j_ = pmlUpdateFxnCplx::addPsi;
    else
        upPsi_j_ = [](std::vector<updateGridParams>&, std::vector<updatePsiParams>&, cplx_pgrid_ptr, cplx_pgrid_ptr, cplx_pgrid_ptr){return;};
    if(psi_k_)
        upPsi_k_ = pmlUpdateFxnCplx::addPsi;
    else
        upPsi_k_ = [](std::vector<updateGridParams>&, std::vector<updatePsiParams>&, cplx_pgrid_ptr, cplx_pgrid_ptr, cplx_pgrid_ptr){return;};
}