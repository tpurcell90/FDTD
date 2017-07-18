#ifndef FDTD_SOURCE
#define FDTD_SOURCE

// #include "enum.hpp"
#include <SOURCE/Pulse.hpp>
#include <src/GRID/parallelGrid.hpp>
/**
 * @brief data structure to know where to add the pulse
 *
 */
struct SalveSource
{
    int stride_; //!< stride for the mkl functions
    std::vector<int> loc_; //!< location of the lower, left, back corner of the source plane
    std::vector<int> sz_; //!< size in grid points of the source
    std::vector<int> addVec1_; //!< unit vector for looping through the source, loop over x->{1,0,0}; y->{0,1,0}; z->{0,0,1}
    std::vector<int> addVec2_; //!< unit vector for looping through the source a second time, loop over x->{1,0,0}; y->{0,1,0}; z->{0,0,1}

};

/**
 * @brief A parallel soft source for FDTD fields.
 *
 * @tparam T param for type of source, doulbe for real, complex<double> for complex field
 */
template <typename T> class parallelSourceBase
{
protected:
    mpiInterface gridComm_; //!<  mpiInterface for the FDTD field
    std::vector< std::shared_ptr<PulseBase> > pulse_; //!< the pulse that the source is adding to the filed
    std::shared_ptr<parallelGrid<T>> grid_; //!< the grid that the source is adding the pulse to

    double dt_; //!< time step of the calculation

    std::array<int,3> loc_; //!< location of the source's lower left corner
    std::array<int,3> sz_; //!< size of the source
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
     */
    parallelSourceBase(mpiInterface comm, std::vector<std::shared_ptr<PulseBase>> pulse, std::shared_ptr<parallelGrid<T>> grid, double dt, std::array<int,3> loc, std::array<int,3> sz) :
        gridComm_(comm),
        grid_(grid),
        dt_(dt),
        loc_(loc),
        sz_(sz)
    {}

    /**
     * @return     pulse_
     */
    inline std::vector<std::shared_ptr<PulseBase>> &pulse() {return pulse_;}
    /**
     * @return     grid_
     */
    inline std::shared_ptr<parallelGrid<T>> &grid() {return grid_;}
    /**
     * @return     loc_
     */
    inline std::array<int,3> &loc() {return loc_;}

    /**
     * @brief      generates the data structures necessary to add the pulse into the field
     *
     * @param[in]  srcNum  index of the source in the srcArr_
     */
    virtual void genDatStruct() = 0;

    /**
     * @brief      adds the pulse to the grid
     *
     * @param[in]  t     current time
     */
    virtual void addPul(double t) = 0;

    /**
     * @brief      Accessor function for source size
     *
     * @return     sz_
     */
    std::array<int,3> sz() {return sz_; }
};

#endif