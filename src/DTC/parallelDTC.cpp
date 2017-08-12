#include <DTC/parallelDTC.hpp>

parallelDetectorBaseReal::parallelDetectorBaseReal(std::vector<real_pgrid_ptr> grids, bool SI, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, double timeInterval, double a, double I0, double dt) :
    parallelDetectorBase( grids, SI, loc, sz, type, timeInterval, a, I0, dt)
{
    for(auto& grid : grids)
    {
        // Check if all the grids are the same size and then construct a storage object for it.
        if( grid->d().size() == grid->d().size() && grids[0]->dx() == grid->dx() && grids[0]->dy() == grid->dy() && grids[0]->dz() == grid->dz() )
            fields_.push_back(std::make_shared<parallelStorageDTCReal>(grid, loc, sz) );
        else
            throw std::logic_error("The step sizes of all the grids for a parallel dtc are not the same.");
    }
}
parallelDetectorBaseCplx::parallelDetectorBaseCplx(std::vector<cplx_pgrid_ptr> grids, bool SI, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, double timeInterval, double a, double I0, double dt) :
    parallelDetectorBase(grids, SI, loc, sz, type, timeInterval, a, I0, dt)
{
    for(auto& grid : grids)
    {
        // Check if all the grids are the same size and then construct a storage object for it.
        if( grid->d().size() == grid->d().size() && grids[0]->dx() == grid->dx() && grids[0]->dy() == grid->dy() && grids[0]->dz() == grid->dz() )
            fields_.push_back(std::make_shared<parallelStorageDTCCplx>(grid, loc, sz) );
        else
            throw std::logic_error("The step sizes of all the grids for a parallel dtc are not the same.");
    }
}