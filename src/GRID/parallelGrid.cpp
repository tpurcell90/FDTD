#include <src/GRID/parallelGrid.hpp>

parallelGridReal::parallelGridReal(mpiInterface & gridComm, bool PBC, std::array<int,3> n_vec, std::array<double,3> d, bool ylim) : parallelGrid<double>(gridComm, PBC, n_vec, d, ylim)
{
    copy_ = dcopy_;
}

parallelGridReal::parallelGridReal(mpiInterface & gridComm, bool PBC, std::vector<real_grid_ptr> weights, std::array<int,3> n_vec, std::array<double,3> d, bool ylim) : parallelGrid<double>(gridComm, PBC, weights, n_vec, d, ylim)
{
    copy_ = dcopy_;
}
// parallelGridReal::parallelGridReal(const parallelGridReal& o) :  parallelGrid<double>(o) {}

parallelGridCplx::parallelGridCplx(mpiInterface & gridComm, bool PBC, std::array<int,3> n_vec, std::array<double,3> d, bool ylim) : parallelGrid<std::complex<double>>(gridComm, PBC, n_vec, d, ylim)
{
    copy_ = zcopy_;
}

parallelGridCplx::parallelGridCplx(mpiInterface & gridComm, bool PBC, std::vector<real_grid_ptr> weights, std::array<int,3> n_vec, std::array<double,3> d, bool ylim) : parallelGrid<std::complex<double>>(gridComm, PBC, weights, n_vec, d, ylim)
{
    copy_ = zcopy_;
}
// parallelGridCplx::parallelGridCplx(const parallelGridCplx& o) :  parallelGrid<std::complex<double>>(o) {}

parallelGridInt::parallelGridInt(mpiInterface & gridComm, bool PBC, std::array<int,3> n_vec, std::array<double,3> d, bool ylim) : parallelGrid<int>(gridComm, PBC, n_vec, d, ylim)
{
    copy_ = scopy_;
}

parallelGridInt::parallelGridInt(mpiInterface & gridComm, bool PBC, std::vector<real_grid_ptr> weights, std::array<int,3> n_vec, std::array<double,3> d, bool ylim) : parallelGrid<int>(gridComm, PBC, weights, n_vec, d, ylim)
{
    copy_ = scopy_;
}