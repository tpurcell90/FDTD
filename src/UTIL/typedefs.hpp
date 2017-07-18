#ifndef PARALLEL_FDTD_TYPEDEFS
#define PARALLEL_FDTD_TYPEDEFS

#include <array>
#include <GRID/Grid.hpp>
#include <GRID/parallelGrid.hpp>

typedef std::vector<std::pair<std::array<int,8>, std::array<double,2> > > upLists;
typedef std::complex<double> cplx;

typedef std::shared_ptr<Grid<cplx>> cplx_grid_ptr;
typedef std::shared_ptr<Grid<double>> real_grid_ptr;
typedef std::shared_ptr<Grid<int>> int_grid_ptr;

typedef std::shared_ptr<parallelGrid<cplx>> cplx_pgrid_ptr;
typedef std::shared_ptr<parallelGrid<double>> real_pgrid_ptr;
typedef std::shared_ptr<parallelGrid<int>> int_pgrid_ptr;
#endif