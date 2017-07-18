#include <DTC/parallelDTC.hpp>

parallelDetectorBaseReal::parallelDetectorBaseReal(int dtcNum, std::vector<real_pgrid_ptr> grids, bool SI, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt) :
    parallelDetectorBase(dtcNum, grids, SI, loc, sz, type, classType, timeInterval, firstComp, a, I0, dt)
{
    for(auto& grid : grids)
    {
        if( grid->d().size() == grid->d().size() && grids[0]->dx() == grid->dx() && grids[0]->dy() == grid->dy() && grids[0]->dz() == grid->dz() )
            fields_.push_back(std::make_shared<parallelStorageDTCReal>(dtcNum, grid, loc, sz) );
        else
            throw std::logic_error("The step sizes of all the grids for the parallel dtc number " + std::to_string(dtcNum) + " are not the same. This will definitely cause me some issues, I think I'll just stop here.");
    }
}
parallelDetectorBaseCplx::parallelDetectorBaseCplx(int dtcNum, std::vector<cplx_pgrid_ptr> grids, bool SI, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt) :
    parallelDetectorBase(dtcNum, grids, SI, loc, sz, type, classType, timeInterval, firstComp, a, I0, dt)
{
    for(auto& grid : grids)
    {
        if( std::equal(grid->d().begin(), grid->d().end(), grids[0]->d().begin() ) )
            fields_.push_back(std::make_shared<parallelStorageDTCCplx>(dtcNum, grid, loc, sz) );
        else
            throw std::logic_error("The step sizes of all the grids for the parallel dtc number " + std::to_string(dtcNum) + " are not the same. This will definitely cause me some issues, I think I'll just stop here.");
    }
}