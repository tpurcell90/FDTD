#ifndef FDTD_SOURCE_NORMAL
#define FDTD_SOURCE_NORMAL
#include <SOURCE/parallelSource.hpp>

/**
 * @brief A parallel soft source for FDTD fields. Assuming phi is along a normal axis of the grid (0,90,270,360)
 *
 * @tparam T param for type of source, doulbe for real, complex<double> for complex field
 */
template <typename T> class parallelSourceNormalBase : public parallelSourceBase<T>
{
protected:
    using parallelSourceBase<T>::gridComm_; //!<  mpiInterface for the FDTD field
    using parallelSourceBase<T>::pulse_; //!< the pulse that the source is adding to the filed
    using parallelSourceBase<T>::grid_; //!< the grid that the source is adding the pulse to

    using parallelSourceBase<T>::dt_; //!< time step of the calculation

    using parallelSourceBase<T>::loc_; //!< location of the source's lower left corner
    using parallelSourceBase<T>::sz_; //!< location of the source's lower left corner

    using parallelSourceBase<T>::pulse; //!< accessor function for pulse_
    using parallelSourceBase<T>::grid; //!< accessor function for grid_
    using parallelSourceBase<T>::loc; //!< accessor function for loc_
    using parallelSourceBase<T>::sz; //!< accessor function for sz

    std::vector<T> pulVec_; //!< dummy vector for the pluse

    std::shared_ptr<SalveSource> slave_; //!< shared pointer to data structure with source adding parameters
public:
    /**
     * @brief Constructor for the parallel source
     *
     * @param[in]  comm       mpiInterface for the caclutlation
     * @param[in]  srcNum     index of the source in srcArr_
     * @param[in]  pulse      pulse of the calculation
     * @param[in]  grid       grid the source adds the pulse to
     * @param[in]  dt         time step of the calculation
     * @param[in]  loc        location of the lower left corner of source
     * @param[in]  sz         size of the soft source
     */
    parallelSourceNormalBase(mpiInterface comm, std::vector<std::shared_ptr<PulseBase>> pulse, std::shared_ptr<parallelGrid<T>> grid, double dt, std::array<int,3> loc, std::array<int,3> sz) :
        parallelSourceBase<T>(comm, pulse, grid, dt, loc, sz),
        pulVec_(std::max(sz_[0],sz_[1]), 0.0)
    {
        slave_ = nullptr;
        genDatStruct();
        if(slave_)
            pulse_ = pulse;
    }

    /**
     * @brief      generates the data structures necessary to add the pulse into the field
     *
     * @param[in]  srcNum  index of the source in the srcArr_
     */
    void genDatStruct()
    {
        SalveSource slaveSrc;
        // slaveSrc.sz_  = {0,0};
        slaveSrc.loc_ = {-1,-1,-1};

        // Checking  if any part of the source is in this location
        if(loc_[0] >= grid_->procLoc()[0] && loc_[0] < grid_->procLoc()[0] + grid_->local_x() -2) // Is the process in the process col that contains the left boundary of the source region?
        {
            slaveSrc.loc_[0] = loc_[0] - grid_->procLoc()[0] + 1;
            if(loc_[1] >= grid_->procLoc()[1] && loc_[1] < grid_->procLoc()[1] + grid_->local_y() -2) // Is the process in the process row that contains the lower boundary of the source region?
            {
                slaveSrc.loc_[1] = loc_[1] - grid_->procLoc()[1] + 1;
                if(grid_->local_z() == 1)
                    slaveSrc.loc_[2] = 0;
                else if(loc_[1] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                    slaveSrc.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                    slaveSrc.loc_[2] = 1;
            }
            else if(loc_[1] < grid_->procLoc()[1] && loc_[1] + sz_[1] > grid_->procLoc()[1]) // Is the process in a process row that the source region covers?
            {
                slaveSrc.loc_[1] = 1;
                if(grid_->local_z() == 1)
                    slaveSrc.loc_[2] = 0;
                else if(loc_[1] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                    slaveSrc.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                    slaveSrc.loc_[2] = 1;
            }
        }
        else if(loc_[0] < grid_->procLoc()[0] && loc_[0] + sz_[0] > grid_->procLoc()[0]) // Is the process in a process col that the source covers?
        {
            slaveSrc.loc_[0] = 1;

            if(loc_[1] >= grid_->procLoc()[1] && loc_[1] < grid_->procLoc()[1] + grid_->local_y() -2)
            {
                slaveSrc.loc_[1] = loc_[1] - grid_->procLoc()[1] + 1;
                if(grid_->local_z() == 1)
                    slaveSrc.loc_[2] = 0;
                else if(loc_[1] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                    slaveSrc.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                    slaveSrc.loc_[2] = 1;
            }
            else if(loc_[1] < grid_->procLoc()[1] && loc_[1] + sz_[1] > grid_->procLoc()[1])
            {
                slaveSrc.loc_[1] = 1;
                if(grid_->local_z() == 1)
                    slaveSrc.loc_[2] = 0;
                else if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                    slaveSrc.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                    slaveSrc.loc_[2] = 1;
            }
        }
        std::array<int, 3> sz = {{0,0,0}};
        if(slaveSrc.loc_[0] != -1 && slaveSrc.loc_[1] != -1 && slaveSrc.loc_[2] != -1 )
        {
            if(sz_[0] + loc_[0] > grid_->procLoc()[0] + grid_->local_x() - 2) // Does the detector go through the end of the process' grid?
                sz[0] = grid_->local_x() - slaveSrc.loc_[0] - 1;
            else
                sz[0] = loc_[0] + sz_[0] - (grid_->procLoc()[0] + slaveSrc.loc_[0] - 1);

            if(sz_[1] + loc_[1] > grid_->procLoc()[1] + grid_->local_y() - 2)
                sz[1] = grid_->local_y() - slaveSrc.loc_[1] - 1;
            else
                sz[1] = loc_[1] + sz_[1] - (grid_->procLoc()[1] + slaveSrc.loc_[1] - 1);

            if(grid_->local_z() == 1)
                sz[2] = 1;
            else if(sz_[2] + loc_[2] > grid_->procLoc()[2] + grid_->local_z() - 2)
                sz[2] = grid_->local_z() - slaveSrc.loc_[2] - 1;
            else
                sz[2] = loc_[2] + sz_[2] - (grid_->procLoc()[2] + slaveSrc.loc_[2] - 1);

            if( sz[0] >= sz[1] && sz[0] >= sz[2] )
            {
                slaveSrc.stride_  = 1;
                slaveSrc.sz_      = {{ sz[0], sz[1] , sz[2] }};
                slaveSrc.addVec1_ = {{ 0, 1, 0 }};
                slaveSrc.addVec2_ = {{ 0, 0, 1 }};
            }
            else if(sz[1] > sz[2] )
            {
                slaveSrc.stride_  = grid_->local_x() * grid_->local_z();
                slaveSrc.sz_      = {{ sz[1], sz[0] , sz[2] }};
                slaveSrc.addVec1_ = {{ 1, 0, 0 }};
                slaveSrc.addVec2_ = {{ 0, 0, 1 }};
            }
            else
            {
                slaveSrc.stride_  = grid_->local_x();
                slaveSrc.sz_      = {{ sz[2], sz[0] , sz[1] }};
                slaveSrc.addVec1_ = {{ 1, 0, 0 }};
                slaveSrc.addVec2_ = {{ 0, 1, 0 }};
            }
            slave_ = std::make_shared<SalveSource>(slaveSrc);
        }
    }
    /**
     * @brief      adds the pulse to the grid
     *
     * @param[in]  t     current time
     */
    virtual void addPul(double t) = 0;
};

class parallelSourceNormalReal : public parallelSourceNormalBase<double>
{
public:
    /**
     * @brief Constructor for the parallel source
     *
     * @param[in]  comm       mpiInterface for the caclutlation
     * @param[in]  srcNum     index of the source in srcArr_
     * @param[in]  pulse      pulse of the calculation
     * @param[in]  grid       grid the source adds the pulse to
     * @param[in]  dt         time step of the calculation
     * @param[in]  loc        location of the lower left corner of source
     * @param[in]  sz         size of the soft source
     */
    parallelSourceNormalReal(mpiInterface comm,  std::vector<std::shared_ptr<PulseBase>> pulse, real_pgrid_ptr grid, double dt, std::array<int,3> loc, std::array<int,3> sz);

    /**
     * @brief      adds the pulse to the grid
     *
     * @param[in]  t     current time
     */
    void addPul(double t);
};

class parallelSourceNormalCplx : public parallelSourceNormalBase<cplx>
{
public:
    /**
     * @brief Constructor for the parallel source
     *
     * @param[in]  comm       mpiInterface for the caclutlation
     * @param[in]  srcNum     index of the source in srcArr_
     * @param[in]  pulse      pulse of the calculation
     * @param[in]  grid       grid the source adds the pulse to
     * @param[in]  dt         time step of the calculation
     * @param[in]  loc        location of the lower left corner of source
     * @param[in]  sz         size of the soft source
     */
    parallelSourceNormalCplx(mpiInterface comm,  std::vector<std::shared_ptr<PulseBase>> pulse, cplx_pgrid_ptr grid, double dt, std::array<int,3> loc, std::array<int,3> sz);

    /**
     * @brief      adds the pulse to the grid
     *
     * @param[in]  t     current time
     */
    void addPul(double t);
};

#endif