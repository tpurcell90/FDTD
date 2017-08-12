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
     * @param[in]  gridComm   mpiInterface for the caclutlation
     * @param[in]  srcNum     index of the source in srcArr_
     * @param[in]  pulse      pulse of the calculation
     * @param[in]  grid       grid the source adds the pulse to
     * @param[in]  dt         time step of the calculation
     * @param[in]  loc        location of the lower left corner of source
     * @param[in]  sz         size of the soft source
     */
    parallelSourceNormalBase(std::shared_ptr<mpiInterface> gridComm, std::vector<std::shared_ptr<PulseBase>> pulse, std::shared_ptr<parallelGrid<T>> grid, double dt, std::array<int,3> loc, std::array<int,3> sz) :
        parallelSourceBase<T>(gridComm, pulse, grid, dt, loc, sz),
        pulVec_(std::max(sz_[0],sz_[1]), 0.0)
    {
        slave_ = nullptr;
        genDatStruct();
        if(slave_)
            pulse_ = pulse;
    }

    /**
     * @brief      Calculates the location of where to start the detector in this process
     *
     * @param[in]  ind  The index of the array to look at (0, 1, or 2)
     *
     * @return     The location of the detector's start in this process for the index ind
     */
    int getLocalLocEl(int ind)
    {
        if(loc_[ind] >= grid_->procLoc()[ind] && loc_[ind] < grid_->procLoc()[ind] + grid_->ln_vec()[ind] - 2) //!< Does this process hold the lower, left or back boundary of the detector
            return loc_[ind] - grid_->procLoc()[ind] + 1;
        else if(loc_[ind] < grid_->procLoc()[ind] && loc_[ind] + sz_[ind] > grid_->procLoc()[ind]) //!< Does this process start inside the detectors region?
            return 1;
        else
            return -1;
    }

    /**
     * @brief      Calculates the size of the detector on this process
     *
     * @param[in]  ind       The index of the array to look at (0, 1, or 2)
     * @param[in]  localLoc  The location of the detector's start in this process for the index ind
     *
     * @return     The local size el.
     */
    int getLocalSzEl(int ind, int localLoc)
    {
        if(sz_[ind] + loc_[ind] > grid_->procLoc()[ind] + grid_->ln_vec()[ind] - 2) // Does the detector go through the end of the process' grid?
            return grid_->ln_vec()[ind] - localLoc - 1;
        else
            return loc_[ind] + sz_[ind] - (grid_->procLoc()[ind] + localLoc - 1);
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

        // Get the location of the source within the detector
        for(int ii = 0; ii < 3; ++ii)
            slaveSrc.loc_[ii] = getLocalLocEl(ii);
        // if 2D calculation set loc[2] to 0
        if(grid_->local_z() == 1)
            slaveSrc.loc_[2] = 0;
        std::array<int, 3> sz = {{0,0,0}};
        if(slaveSrc.loc_[0] != -1 && slaveSrc.loc_[1] != -1 && slaveSrc.loc_[2] != -1 )
        {
            // get the size of the source within the process
            for(int ii = 0; ii < 3; ++ii)
                sz[ii] = getLocalSzEl(ii, slaveSrc.loc_[ii]);
            // set the size of the z to 1 if in 2D
            if(grid_->local_z() == 1)
                sz[2] = 1;

            // Set blass operations to be along the longest direction
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
     * @param[in]  gridComm       mpiInterface for the caclutlation
     * @param[in]  srcNum     index of the source in srcArr_
     * @param[in]  pulse      pulse of the calculation
     * @param[in]  grid       grid the source adds the pulse to
     * @param[in]  dt         time step of the calculation
     * @param[in]  loc        location of the lower left corner of source
     * @param[in]  sz         size of the soft source
     */
    parallelSourceNormalReal(std::shared_ptr<mpiInterface> gridComm,  std::vector<std::shared_ptr<PulseBase>> pulse, real_pgrid_ptr grid, double dt, std::array<int,3> loc, std::array<int,3> sz);

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
     * @param[in]  gridComm       mpiInterface for the caclutlation
     * @param[in]  srcNum     index of the source in srcArr_
     * @param[in]  pulse      pulse of the calculation
     * @param[in]  grid       grid the source adds the pulse to
     * @param[in]  dt         time step of the calculation
     * @param[in]  loc        location of the lower left corner of source
     * @param[in]  sz         size of the soft source
     */
    parallelSourceNormalCplx(std::shared_ptr<mpiInterface> gridComm,  std::vector<std::shared_ptr<PulseBase>> pulse, cplx_pgrid_ptr grid, double dt, std::array<int,3> loc, std::array<int,3> sz);

    /**
     * @brief      adds the pulse to the grid
     *
     * @param[in]  t     current time
     */
    void addPul(double t);
};

#endif