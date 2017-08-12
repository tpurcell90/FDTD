#ifndef FDTD_PARALLELGRID
#define FDTD_PARALLELGRID

#include <MPI/mpiInterface.hpp>
#include <UTIL/enum.hpp>
#include <UTIL/mathUtils.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <iomanip>
#include <utility>
#include <numeric>

namespace mpi = boost::mpi;

typedef std::shared_ptr<Grid<double>> real_grid_ptr;
template <typename T>
class parallelGrid
{
protected:
    // The Grid's Communicator
    std::shared_ptr<mpiInterface> gridComm_; //!< The communicator for the processes that are storing the grid

    // global parameters
    std::array<int,3> n_vec_; //!< number of grid points of the full grid in all directions

    // Local Parameters
    std::array<int,3> ln_vec_; //!< number of grid points of the local (single process) gird in all directions

    int upProc_; //!< The process that stores the data directly above the current process, -1 for none
    int downProc_; //!< The process that stores the data directly below the current process, -1 for none

    int upSendTag_; //!< Tag used for sending data to the process above (positive y direction)
    int upRecvTag_; //!< Tag used for receiving data from the process above (positive y direction)
    int downSendTag_; //!< Tag used for sending data to the process below (negative y direction)
    int downRecvTag_; //!< Tag used for receiving data from the process below (negative y direction)

    int upSendIndex_; //!< Typically ln_vec_[1]-2, if periodic and filed is limited in y then ln_vec_[1]-3

    bool PBC_; //!< true if periodic boundary conditions

    std::array<double,3> d_; //!< step size in all directions
    // Scalapack specific
    std::array<int,3> procLoc_; //!< location of the lower left corner of the process in the total grid
    std::array<int,3> procR_C_; //!< location of the lower left corner of the process in the total grid including buffer space
    std::vector<int> xTrans_; //!< Vector storing the index of each process' starting location value
    std::vector<int> yTrans_; //!< Vector storing the index of each process' starting location value
    std::vector<int> zTrans_; //!< Vector storing the index of each process' starting location value

    std::array<PROC_DIR, 2> sendList_; //!< List of send command parameters for each process
    std::array<PROC_DIR, 2> recvList_; //!< List of recv command parameters for each process

    std::vector<mpi::request> reqs_; //!< vector to store all mpi requests for waiting at each transfer point

    // distributed parameters
    std::unique_ptr<T[]> local_; //!< Data array for the local grid

public:
    /**
     * @brief      constructs a parallelGrid without any information on the workload for each process
     *
     * @param      gridComm  The mpiInterface for communication
     * @param[in]  PBC       True if periodic
     * @param[in]  n_vec     Size of the grid in all directions
     * @param[in]  d         grid spacing in all directions
     * @param[in]  ylim      True if grid is for Ey, Hx or Hz fields.
     */
    parallelGrid(std::shared_ptr<mpiInterface> gridComm, bool PBC, std::array<int,3> n_vec, std::array<double, 3> d, bool ylim=false) :
        gridComm_(gridComm),
        n_vec_(n_vec),
        ln_vec_({{0,0,0}}),
        upProc_(-1),
        downProc_(-1),
        upSendTag_(-1),
        upRecvTag_(-1),
        downSendTag_(-1),
        downRecvTag_(-1),
        PBC_(PBC),
        d_(d),
        procLoc_(std::array<int,3>({0,0,0})),
        procR_C_(std::array<int,3>({0,0,0})),
        xTrans_(gridComm_->npX(),0),
        yTrans_(gridComm_->npY(),0),
        zTrans_(gridComm_->npZ(),0)
    {
        std::tie(ln_vec_[0],ln_vec_[1], ln_vec_[2]) = gridComm_->numroc(n_vec_[0], n_vec_[1], n_vec_[2]);
        if(n_vec_[0] == 1)
            ln_vec_[0] = 1;
        if(n_vec_[1] == 1)
            ln_vec_[1] = 1;
        if(n_vec_[2] == 1)
            ln_vec_[2] = 1;
        local_ = std::unique_ptr<T[]>(new T[size()]);
        zero();

        genProcSendRecv(PBC);
        reqs_ = std::vector<mpi::request>(sendList_.size()*2,mpi::request());
        upSendIndex_ = ln_vec_[1]-2;
        // Find the procLocation of lower left corner in the grid Assumes Cartesian Grid for the procs
        determineProcLoc();
        if(PBC && (gridComm_->rank() == gridComm_->size() - 1) && ylim)
            upSendIndex_ -= 1;
    }
    /**
     * @brief      constructs a parallelGrid with information on the workload for each process
     *
     * @param      gridComm  The mpiInterface for communication
     * @param[in]  PBC       True if periodic
     * @param[in]  weights   Vector containing information of the work load for each point
     * @param[in]  n_vec     Size of the grid in all directions
     * @param[in]  d         grid spacing in all directions
     * @param[in]  ylim      True if grid is for Ey, Hx or Hz fields.
     */
    parallelGrid(std::shared_ptr<mpiInterface> gridComm, bool PBC, std::vector<real_grid_ptr> weights, std::array<int,3> n_vec, std::array<double,3> d, bool ylim=false) :
        gridComm_(gridComm),
        n_vec_(n_vec),
        ln_vec_({{0,0,0}}),
        upProc_(-1),
        downProc_(-1),
        upSendTag_(-1),
        upRecvTag_(-1),
        downSendTag_(-1),
        downRecvTag_(-1),
        PBC_(PBC),
        d_(d),
        procLoc_(std::array<int,3>({0,0,0})),
        procR_C_(std::array<int,3>({0,0,0})),
        xTrans_(gridComm_->npX(),0),
        yTrans_(gridComm_->npY(),0),
        zTrans_(gridComm_->npZ(),0)
    {
        std::tie(ln_vec_[0],ln_vec_[1], ln_vec_[2]) = gridComm_->getLocxLocyLocz(weights);
        if(n_vec_[0] == 1)
            ln_vec_[0] = 1;
        if(n_vec_[1] == 1)
            ln_vec_[1] = 1;
        if(n_vec_[2] == 1)
            ln_vec_[2] = 1;
        local_ = std::unique_ptr<T[]>(new T[size()]);
        zero();

        genProcSendRecv(PBC);

        reqs_ = std::vector<mpi::request>(sendList_.size()*2,mpi::request());

        // Find the procLocation of lower left corner in the grid Assumes Cartesian Grid for the procs
        determineProcLoc();
        upSendIndex_ = ln_vec_[1]-2;
        if(PBC && (gridComm_->rank() == gridComm_->size() - 1) && ylim)
            upSendIndex_ -= 1;
    }


    /**
     * @brief     etermines which processors borders the current process
     *
     * @param[in]  PBC   True if using periodic boundary conditions
     */
    void genProcSendRecv(bool PBC)
    {
        int ii = gridComm_->rank();
        // if the process is at the lower boundary down proc only exists with PBC
        if(ii % gridComm_->npY() != 0)
        {
            downProc_ =  ii - 1;
            downSendTag_ = gridComm_->cantorTagGen(gridComm_->rank(), downProc_, 4, 0);
            downRecvTag_ = gridComm_->cantorTagGen(downProc_, gridComm_->rank(), 4, 1);
        }
        else if(PBC)
        {
            downProc_ =  ii + gridComm_->npY()-1;
            downSendTag_ = gridComm_->cantorTagGen(gridComm_->rank(), downProc_, 4, 0);
            downRecvTag_ = gridComm_->cantorTagGen(downProc_, gridComm_->rank(), 4, 1);
        }

        // if the process is at the upper boundary up proc only exists with PBC
        if(ii % gridComm_->npY() != gridComm_->npY() - 1)
        {
            upProc_ =  ii + 1;
            upSendTag_ = gridComm_->cantorTagGen(gridComm_->rank(), upProc_, 4, 1);
            upRecvTag_ = gridComm_->cantorTagGen(upProc_, gridComm_->rank(), 4, 0);
        }
        else if(PBC)
        {
            upProc_ =  ii - gridComm_->npY()+1;
            upSendTag_ = gridComm_->cantorTagGen(gridComm_->rank(), upProc_, 4, 1);
            upRecvTag_ = gridComm_->cantorTagGen(upProc_, gridComm_->rank(), 4, 0);
        }
        // if upProc_ is still -1 PROC_DIR == NONE otherwise its UP
        if(upProc_ != -1)
        {
            sendList_[0] = PROC_DIR::UP;
            recvList_[1] = PROC_DIR::UP;
        }
        else
        {
            sendList_[0] = PROC_DIR::NONE;
            recvList_[1] = PROC_DIR::NONE;
        }
        // if downProc_ is still -1 PROC_DIR == NONE otherwise its DOWN
        if(downProc_ != -1)
        {
            sendList_[1] = PROC_DIR::DOWN;
            recvList_[0] = PROC_DIR::DOWN;
        }
        else
        {
            sendList_[1] = PROC_DIR::NONE;
            recvList_[0] = PROC_DIR::NONE;
        }
    }

    /**
     * @brief      Determines what grid point represents the lower left corner of the grid
     */
    void determineProcLoc()
    {
        // Grids are split only in the y direction so these are always 0
        procLoc_[0] = 0;
        procLoc_[2] = 0;

        procR_C_[0] = 0;
        procR_C_[2] = 0;

        int sumY = 0, tempY;
        for(int ii = 0; ii < gridComm_->npY(); ii++)
        {
            // yTrans stores all processes procR_C[1] value
            yTrans_[ii] = sumY + 2*ii;
            if(gridComm_->rank() % gridComm_->npY() == ii)
            {
                // procLoc does not include boundary PBC storage regions, procR_C does
                procLoc_[1] = sumY;
                procR_C_[1] = sumY + ii*2;
            }
            // tempY does not include boundary
            if(gridComm_->rank() == ii)
                tempY = ln_vec_[1]-2;
            // send tempY for the process to all processes
            broadcast(*gridComm_, tempY, ii);
            // add tempY to sumY
            sumY += tempY;
        }

    }

    /**
     * @brief      Acessor function for  xTrans_
     *
     * @return      xTrans_
     */
    inline std::vector<int> x_trans() {return xTrans_;};
    /**
     * @brief      Acessor function for  yTrans_
     *
     * @return      yTrans_
     */
    inline std::vector<int> y_trans() {return yTrans_;};
    /**
     * @brief      Acessor function for  zTrans_
     *
     * @return      zTrans_
     */
    inline std::vector<int> z_trans() {return zTrans_;};

    /**
     * @brief      returns the value at a point x_val,y_val,0
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     *
     * @return     reference to the data point at (x,y,0)
     */
    inline T& point(const int x, const int y)
    {
        assert(0 <= x && x < ln_vec_[0] && 0 <= y && y < ln_vec_[1] );
        return local_[y*ln_vec_[0]*ln_vec_[2] + x];
    }

    /**
     * @brief      returns the value at a point x_val,y_val,0
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     *
     * @return     reference to the data point at (x,y,0)
     */
    inline const T& point(const int x, const int y) const  {
        assert(0 <= x && x < ln_vec_[0] && 0 <= y && y < ln_vec_[1] );
        return local_[y*ln_vec_[0]*ln_vec_[2] + x];
    }

    /**
     * @brief      returns the value at a point x_val,y_val,0
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     *
     * @return     reference to the data point at (x,y,0)
     */
    inline T& operator()(const int x, const int y) { return point(x,y); }

    /**
     * @brief      returns the value at a point x_val,y_val,0
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     *
     * @return     reference to the data point at (x,y,0)
     */
    inline const T& operator()(const int x, const int y) const { return point(x,y); }

    /**
     * @brief      returns the value at a point x_val,y_val,z_val
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     * @param[in]  z_val  grid coordinate in the z direction
     *
     * @return     reference to the data point at (x,y,z)
     */
    inline T& point(const int x, const int y, const int z)
    {
        assert(0 <= x && x < ln_vec_[0] && 0 <= y && y < ln_vec_[1] && 0 <= z && z < ln_vec_[2] );
        return local_[y*ln_vec_[0]*ln_vec_[2] + z*ln_vec_[0] + x];
    }

    /**
     * @brief      returns the value at a point x_val,y_val,z_val
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     * @param[in]  z_val  grid coordinate in the z direction
     *
     * @return     reference to the data point at (x,y,z)
     */
    inline const T& point(const int x, const int y, const int z) const  {
        assert(0 <= x && x < ln_vec_[0] && 0 <= y && y < ln_vec_[1] && 0 <= z && z < ln_vec_[2] );
        return local_[y*ln_vec_[0]*ln_vec_[2] + z*ln_vec_[0] + x];
    }

    /**
     * @brief      returns the value at a point x_val,y_val,z_val
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     * @param[in]  z_val  grid coordinate in the z direction
     *
     * @return     reference to the data point at (x,y,z)
     */
    inline T& operator()(const int x, const int y, const int z) { return point(x,y,z); }

    /**
     * @brief      returns the value at a point x_val,y_val,z_val
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     * @param[in]  z_val  grid coordinate in the z direction
     *
     * @return     reference to the data point at (x,y,z)
     */
    inline const T& operator()(const int x, const int y, const int z) const { return point(x,y,z); }

    /**
     * @brief Return the Data storage of the local grid
     * @return the Data storage of the local grid
     */
    /**
     *
     * @return     a pointer to the data vector
     */
    inline const std::unique_ptr<T[]>& local() const { return local_; }

    /**
     * @return     total size of the storage vector
     */
    inline int size() const { return ln_vec_[0]*ln_vec_[1]*ln_vec_[2]; }

    /**
     * @return     the global number of grid points in the x direction
     */
    inline int x() const { return n_vec_[0]; }

    /**
     * @return     the global number of grid points in the y direction
     */
    inline int y() const { return n_vec_[1]; }

    /**
     * @return     the global number of grid points in the z direction
     */
    inline int z() const { return n_vec_[2]; }

    /**
     * @return     an array containing the global number of grid points in all directions
     */
    inline std::array<int,3> n_vec() const {return n_vec_;}

    /**
     * @return     an array containing the global number of grid points in all directions
     */
    inline std::array<int,3> ln_vec() const {return ln_vec_;}

    /**
     * @return     the local number of grid points in the x direction
     */
    inline int local_x() const { return ln_vec_[0]; }

    /**
     * @return     the local number of grid points in the y direction
     */
    inline int local_y() const { return ln_vec_[1]; }

    /**
     * @return     the local number of grid points in the z direction
     */
    inline int local_z() const { return ln_vec_[2]; }

    /**
     * @return     the grid spacing in the x direction
     */
    inline double dx() const { return d_[0]; }

    /**
     * @return     the grid spacing in the y direction
     */
    inline double dy() const { return d_[1]; }

    /**
     * @return     the grid spacing in the z direction
     */
    inline double dz() const { return d_[2]; }

    /**
     * @return     an array storing grid spacing in all directions
     */
    inline std::array<double,3> d() const {return d_;}

    /**
     * @return     an array storing the first point this process holds in the global storage vector neglecting transfer boundaries
     */
    inline std::array<int,3> procLoc() {return procLoc_;}

    /**
     * @return     an array storing the first point this process holds in the global storage vector including transfer boundaries
     */
    inline std::array<int,3> procR_C() {return procR_C_;}


    /**
     * @brief      Fills the local grid with a specified value
     *
     * @param[in]  a     value to fill vector with
     */
    inline void fill(const T a) { std::fill_n(local_.get(), ln_vec_[0]*ln_vec_[1]*ln_vec_[2], a); }

    /**
     * @brief Fills the local grid with zeros
     *
     */
    inline void zero() { fill(static_cast<T>(0.0)); }

    /**
     * @return     true if periodic boundary conditions are used
     */
    inline bool PBC() {return PBC_;}

    /**
     * @return     shared_ptr to the mpiInterface
     */
    inline std::shared_ptr<mpiInterface> gridComm(){return gridComm_;}

    /**
     * @brief      Transfer data from one process to another
     */
    void transferDat()
    {
        for(int ii = 0; ii < sendList_.size(); ii++)
        {
            switch(sendList_[ii])
            {
                case PROC_DIR::UP:
                    reqs_[ii*2] = gridComm_->isend(  upProc_,   upSendTag_, &point(0,upSendIndex_,0), ln_vec_[0]*ln_vec_[2]);
                    break;
                case PROC_DIR::DOWN:
                    reqs_[ii*2] = gridComm_->isend(downProc_, downSendTag_, &point(0,1,0), ln_vec_[0]*ln_vec_[2]);
                    break;
                case PROC_DIR::NONE:
                    break;
                default:
                    throw std::logic_error("The PROC_DIR default has been hit, see transferDat() function to add the direction here.");
            }
            switch(recvList_[ii])
            {
                case PROC_DIR::UP:
                    reqs_[ii*2+1] = gridComm_->irecv(  upProc_,   upRecvTag_, &point(0, upSendIndex_+1,0), ln_vec_[0]*ln_vec_[2]);
                    break;
                case PROC_DIR::DOWN:
                    reqs_[ii*2+1] = gridComm_->irecv(downProc_, downRecvTag_, &point(0,0,0), ln_vec_[0]*ln_vec_[2]);
                    break;
                case PROC_DIR::NONE:
                    break;
                default:
                    throw std::logic_error("The PROC_DIR default has been hit, see transferDat() function to add the direction here.");
            }
            mpi::wait_all(reqs_.data(), reqs_.data() + reqs_.size() );
        }
    }


    /**
     * @brief      Finds the x coordinate in the process grid that has the value of i in the x coordinate of the real grid
     *
     * @param[in]  i     value of the grid's x coordinate you are looking to find
     *
     * @return     The x coordinate in the process grid, and the offset of where i is relative to the process Grid's x value
     */
    std::pair<int, int> locate_x(const int i)
    {
        int ii = 0;
        while (ii < xTrans_.size()-1 && i >= xTrans_[ii+1])
            ++ii;
        return {ii, i-xTrans_[ii]};
    }

    /**
     * @brief      Finds the y coordinate in the process grid that has the value of i in the y coordinate of the real grid
     *
     * @param[in]  j     value of the grid's y coordinate you are looking to find
     *
     * @return     The y coordinate in the process grid, and the offset of where i is relative to the process Grid's y value
     */
    std::pair<int, int> locate_y(const int j)
    {
        int jj = 0;
        while (jj < yTrans_.size()-1 && j >= yTrans_[jj+1])
            jj++;
        return {jj, j-yTrans_[jj]};
    }

    /**
     * @brief      Finds the z coordinate in the process grid that has the value of i in the z coordinate of the real grid
     *
     * @param[in]  k     value of the grid's z coordinate you are looking to find
     *
     * @return     The z coordinate in the process grid, and the offset of where i is relative to the process Grid's z value
     */
    std::pair<int, int> locate_z(const int k)
    {
        int kk = 0;
        while (kk < zTrans_.size()-1 && k >= zTrans_[kk+1])
            kk++;
        return {kk, k-zTrans_[kk]};
    }

    /**
     * @brief      Finds the x coordinate in the process grid that has the value of i in the x coordinate of the real grid (neglecting boundaries)
     *
     * @param[in]  i     value of the grid's x coordinate you are looking to find
     *
     * @return     The x coordinate in the process grid, and the offset of where i is relative to the process Grid's x value (neglecting boundaries)
     */
    std::pair<int, int> locate_x_no_boundaries(const int i)
    {
        int ii = 0;
        while (ii < xTrans_.size()-1 && i >= xTrans_[ii+1]-2*(ii+1) )
            ii++;
        return {ii, i-xTrans_[ii]};
    }

    /**
     * @brief      Finds the y coordinate in the process grid that has the value of i in the y coordinate of the real grid (neglecting boundaries)
     *
     * @param[in]  j     value of the grid's y coordinate you are looking to find
     *
     * @return     The y coordinate in the process grid, and the offset of where i is relative to the process Grid's y value (neglecting boundaries)
     */
    std::pair<int, int> locate_y_no_boundaries(const int j)
    {
        int jj = 0;
        while (jj < yTrans_.size()-1 && j >= yTrans_[jj+1]-2*(jj+1) )
            jj++;
        return {jj, j-yTrans_[jj]};
    }

    /**
     * @brief      Finds the z coordinate in the process grid that has the value of i in the z coordinate of the real grid (neglecting boundaries)
     *
     * @param[in]  k     value of the grid's z coordinate you are looking to find
     *
     * @return     The z coordinate in the process grid, and the offset of where i is relative to the process Grid's z value (neglecting boundaries)
     */
    std::pair<int, int> locate_z_no_boundaries(const int k)
    {
        int kk = 0;
        while (kk < yTrans_.size()-1 && k >= yTrans_[kk+1]-2*(kk+1) )
            kk++;
        return {kk, k-yTrans_[kk]};
    }

    /**
     * @brief      Returns the process rank storing a particular point
     *
     * @param[in]  xx    the point coordinate in the x direction
     * @param[in]  yy    the point coordinate in the y direction
     * @param[in]  zz    the point coordinate in the z direction
     *
     * @return     The process that stores that point (xx, yy, zz) including boundaries
     */
    int getLocsProc(int xx, int yy, int zz)
    {
        int py, off;
        std::tie(py, off) = locate_y(yy);

        return py;
    }

    /**
     * @brief      Returns the process rank storing a particular point (neglecting boundaries)
     *
     * @param[in]  xx    the point coordinate in the x direction
     * @param[in]  yy    the point coordinate in the y direction
     * @param[in]  zz    the point coordinate in the z direction
     *
     * @return     The process that stores the point (xx, yy, zz).
     */
    int getLocsProc_no_boundaries(int xx, int yy, int zz)
    {
        int py, off;

        std::tie(py, off) = locate_y_no_boundaries(yy);

        return py;
    }

    /**
     * @brief      Returns the process rank storing a particular point
     *
     * @param[in]  xx    the point coordinate in the x direction
     * @param[in]  yy    the point coordinate in the y direction
     *
     * @return     The process that stores that point (xx, yy, 0) including boundaries
     */
    int getLocsProc(int xx, int yy)
    {
        int py, off;
        std::tie(py, off) = locate_y(yy);

        return py;
    }

    /**
     * @brief      Returns the process rank storing a particular point (neglecting boundaries)
     *
     * @param[in]  xx    the point coordinate in the x direction
     * @param[in]  yy    the point coordinate in the y direction
     *
     * @return     The process that stores the point (xx,yy).
     */
    int getLocsProc_no_boundaries(int xx, int yy)
    {
        int py, off;

        std::tie(py, off) = locate_y_no_boundaries(yy);

        return py;
    }

    /**
     * @brief      Gets the XZ plane at global point y=j.
     *
     * @param[in]  j     global value of y to take the plane at
     *
     * @return     The XZ plane at global point y=j.
     */
    std::vector<T> getPlaneXZ(const int j)
    {
        int py, off;
        std::tie(py, off) = locate_y(j);
        std::vector<T> planeXZ(n_vec_[0]*n_vec_[2], 0.0);
        if(gridComm_->rank() == py)
            std::copy_n(&point(0,off,0), ln_vec_[0]*ln_vec_[2], planeXZ.data());

        mpi::broadcast(*gridComm_, planeXZ, py);
        return planeXZ;
    }

    /**
     * @brief      Gets the XY plane at global point z=k.
     *
     * @param[in]  k     global value of z to take the plane at
     *
     * @return     The XZ plane at global point z=k.
     */
    std::vector<T> getPlaneXY(const int k)
    {
        int pz, off;
        std::tie(pz, off) = locate_z(k);
        std::vector<T> planeXY(n_vec_[0]*n_vec_[1], 0.0);

        if(gridComm_->rank() != pz)
        {
            std::vector<T> toSend(ln_vec_[1]*ln_vec_[0], 0.0);
            for(int yy = 0; yy < ln_vec_[1]; ++yy)
                std::copy_n( &point(0, yy, off), ln_vec_[0], &toSend[ln_vec_[0]*yy] );
            gridComm_->send(pz, gridComm_->cantorTagGen(gridComm_->rank(), pz, 2, 0), toSend);
        }
        else
        {
            for(int yy = 0; yy < ln_vec_[1]; ++yy)
                std::copy_n( &point(0, yy, off), ln_vec_[0], &planeXY[ln_vec_[0]*yy] );
            int spot = ln_vec_[0]*ln_vec_[1];
            for(int ii = 0; ii < gridComm_->size(); ii ++)
            {
                if(ii != pz )
                {
                    std::vector<T> tempVec;
                    gridComm_->recv(ii, gridComm_->cantorTagGen(ii, gridComm_->rank(), 2, 0), tempVec);
                    std::copy_n( tempVec.data(), tempVec.size(), &planeXY[spot] );
                    spot += tempVec.size();
                }
            }
        }
        mpi::broadcast(*gridComm_, planeXY, pz);
        return planeXY;
    }

    /**
     * @brief      Gets the YZ plane at global point x=i.
     *
     * @param[in]  i     global value of x to take the plane at
     *
     * @return     The YZ plane at global point x=i.
     */
    std::vector<T> getPlaneYZ(const int i)
    {
        gridComm_->barrier();
        int px, off;
        std::tie(px, off) = locate_x(i);
        std::vector<T> planeYZ(n_vec_[1]*n_vec_[2], 0.0);

        gridComm_->barrier();
        if(gridComm_->rank() != px)
        {
            // gridComm_->barrier();
            std::vector<T> toSend(ln_vec_[1]*ln_vec_[2], 0.0);
            // gridComm_->barrier();
            for(int yy = 0; yy < ln_vec_[1]; ++yy)
                for(int zz = 0; zz < ln_vec_[2]; ++zz)
                    toSend[ln_vec_[2]*yy + zz] = point(off, yy, zz);
            // gridComm_->barrier();
            gridComm_->send(px, gridComm_->cantorTagGen(gridComm_->rank(), px, 2, 0), toSend);
            // gridComm_->barrier();
        }
        else
        {
            // gridComm_->barrier();
            for(int yy = 0; yy < ln_vec_[1]; ++yy)
                for(int zz = 0; zz < ln_vec_[2]; ++zz)
                    planeYZ[ln_vec_[2]*yy+zz] = point(off, yy, zz);
            // gridComm_->barrier();
            int spot = ln_vec_[2]*ln_vec_[1];
            for(int ii = 0; ii < gridComm_->size(); ii ++)
            {
                if(ii != px )
                {
                    std::vector<T> tempVec;
                    gridComm_->recv(ii, gridComm_->cantorTagGen(ii, gridComm_->rank(), 2, 0), tempVec);
                    std::copy_n(tempVec.begin(), tempVec.size(), &planeYZ[spot]);
                    spot += tempVec.size();
                }
            }
            // gridComm_->barrier();
        }
        gridComm_->barrier();
        mpi::broadcast(*gridComm_, planeYZ, px);
        gridComm_->barrier();
        return planeYZ;
    }

    template <typename U> friend std::ostream &operator<<(std::ostream &out, const parallelGrid <U> &o);
};

#endif