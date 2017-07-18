#ifndef FDTD_PARALLELGRID
#define FDTD_PARALLELGRID

#include <MPI/mpiInterface.hpp>
#include <GRID/Grid.hpp>
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
    mpiInterface gridComm_; //!< The communicator for the processes that are storing the grid

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
    void (*copy_)(const int a, const T* b, const int c, T* d, const int e); //!< copy command d or zcopy_ from BLAS
    /**
     * @brief      constructs a parallelGrid without any information on the workload for each process
     *
     * @param      gridComm  The mpiInterface for communication
     * @param[in]  PBC       True if periodic
     * @param[in]  n_vec     Size of the grid in all directions
     * @param[in]  d         grid spacing in all directions
     * @param[in]  ylim      True if grid is for Ey, Hx or Hz fields.
     */
    parallelGrid(mpiInterface & gridComm, bool PBC, std::array<int,3> n_vec, std::array<double, 3> d, bool ylim=false) :
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
        xTrans_(gridComm_.npX(),0),
        yTrans_(gridComm_.npY(),0),
        zTrans_(gridComm_.npZ(),0)
    {
        std::tie(ln_vec_[0],ln_vec_[1], ln_vec_[2]) = gridComm_.numroc(n_vec_[0], n_vec_[1], n_vec_[2]);
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
        if(PBC && (gridComm_.rank() == gridComm_.size() - 1) && ylim)
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
    parallelGrid(mpiInterface & gridComm, bool PBC, std::vector<real_grid_ptr> weights, std::array<int,3> n_vec, std::array<double,3> d, bool ylim=false) :
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
        xTrans_(gridComm_.npX(),0),
        yTrans_(gridComm_.npY(),0),
        zTrans_(gridComm_.npZ(),0)
    {
        std::tie(ln_vec_[0],ln_vec_[1], ln_vec_[2]) = gridComm_.getLocxLocyLocz(weights);
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
        if(PBC && (gridComm_.rank() == gridComm_.size() - 1) && ylim)
            upSendIndex_ -= 1;
    }

    /**
     * @brief Determines which processors borders the current process
     * @details Figures out where to send the boundary information in each process
     *
     * @param PBC boolean that sets if PBC are present
     */
    void genProcSendRecv(bool PBC)
    {
        int ii = gridComm_.rank();
        if(ii % gridComm_.npY() != 0)
        {
            downProc_ =  ii - 1;
            downSendTag_ = gridComm_.cantorTagGen(gridComm_.rank(), downProc_, 4, 0);
            downRecvTag_ = gridComm_.cantorTagGen(downProc_, gridComm_.rank(), 4, 1);
        }
        else if(PBC)
        {
            downProc_ =  ii + gridComm_.npY()-1;
            downSendTag_ = gridComm_.cantorTagGen(gridComm_.rank(), downProc_, 4, 0);
            downRecvTag_ = gridComm_.cantorTagGen(downProc_, gridComm_.rank(), 4, 1);
        }

        if(ii % gridComm_.npY() != gridComm_.npY() - 1)
        {
            upProc_ =  ii + 1;
            upSendTag_ = gridComm_.cantorTagGen(gridComm_.rank(), upProc_, 4, 1);
            upRecvTag_ = gridComm_.cantorTagGen(upProc_, gridComm_.rank(), 4, 0);
        }
        else if(PBC)
        {
            upProc_ =  ii - gridComm_.npY()+1;
            upSendTag_ = gridComm_.cantorTagGen(gridComm_.rank(), upProc_, 4, 1);
            upRecvTag_ = gridComm_.cantorTagGen(upProc_, gridComm_.rank(), 4, 0);
        }

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
     * @brief Determines what grid point represents the lower left corner of the grid
     * @details Stores the information for where the processor is in the real grid and the boundary included grid
     */
    void determineProcLoc()
    {
        procLoc_[0] = 0;
        procLoc_[2] = 0;

        procR_C_[0] = 0;
        procR_C_[2] = 0;

        int sumY = 0, tempY;
        for(int ii = 0; ii < gridComm_.npY(); ii++)
        {
            yTrans_[ii] = sumY + 2*ii;
            if(gridComm_.rank() % gridComm_.npY() == ii)
            {
                procLoc_[1] = sumY;
                procR_C_[1] = sumY + ii*2;
            }
            if(gridComm_.rank() == ii)
                tempY = ln_vec_[1]-2;
            broadcast(gridComm_, tempY, ii);
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
     * @brief Get grid element at specified point
     *
     * @param x x grid point
     * @param y y grid point
     *
     * @return Reference to data point
     */
    inline T& point(const int x, const int y)
    {
        assert(0 <= x && x < ln_vec_[0] && 0 <= y && y < ln_vec_[1] );
        return local_[y*ln_vec_[0]*ln_vec_[2] + x];
    }

    /**
     * @brief Get grid element at specified point
     *
     * @param x x grid point
     * @param y y grid point
     *
     * @return Reference to data point
     */
    inline const T& point(const int x, const int y) const  {
        assert(0 <= x && x < ln_vec_[0] && 0 <= y && y < ln_vec_[1] );
        return local_[y*ln_vec_[0]*ln_vec_[2] + x];
    }

    /**
     * @brief Get grid element at specified point
     *
     * @param x x grid point
     * @param y y grid point
     *
     * @return Reference to data point
     */
    inline T& operator()(const int x, const int y) { return point(x,y); }

    /**
     * @brief Get grid element at specified point
     *
     * @param x x grid point
     * @param y y grid point
     *
     * @return Reference to data point
     */
    inline const T& operator()(const int x, const int y) const { return point(x,y); }

    /**
     * @brief Get grid element at specified point
     *
     * @param x x grid point
     * @param y y grid point
     * @param z z grid point
     *
     * @return Reference to data point
     */
    inline T& point(const int x, const int y, const int z)
    {
        assert(0 <= x && x < ln_vec_[0] && 0 <= y && y < ln_vec_[1] && 0 <= z && z < ln_vec_[2] );
        return local_[y*ln_vec_[0]*ln_vec_[2] + z*ln_vec_[0] + x];
    }

    /**
     * @brief Get grid element at specified point
     *
     * @param x x grid point
     * @param y y grid point
     * @param z z grid point
     *
     * @return Reference to data point
     */
    inline const T& point(const int x, const int y, const int z) const  {
        assert(0 <= x && x < ln_vec_[0] && 0 <= y && y < ln_vec_[1] && 0 <= z && z < ln_vec_[2] );
        return local_[y*ln_vec_[0]*ln_vec_[2] + z*ln_vec_[0] + x];
    }

    /**
     * @brief Get grid element at specified point
     *
     * @param x x grid point
     * @param y y grid point
     * @param z z grid point
     *
     * @return Reference to data point
     */
    inline T& operator()(const int x, const int y, const int z) { return point(x,y,z); }

    /**
     * @brief Get grid element at specified point
     *
     * @param x x grid point
     * @param y y grid point
     * @param z z grid point
     *
     * @return Reference to data point
     */
    inline const T& operator()(const int x, const int y, const int z) const { return point(x,y,z); }

    /**
     * @brief Return the Data storage of the local grid
     * @return the Data storage of the local grid
     */
    inline const std::unique_ptr<T[]>& local() const { return local_; }

    /**
     * @brief Get the size of the data storage
     * @return number of points in the grid
     */
    inline int size() const { return ln_vec_[0]*ln_vec_[1]*ln_vec_[2]; }

    /**
     * @brief Get number of global grid points in the x direction
     * @return number of global grid points in the x direction
     */
    inline int x() const { return n_vec_[0]; }
    /**
     * @brief Get number of global grid points in the y direction
     * @return number of global grid points in the y direction
     */
    inline int y() const { return n_vec_[1]; }
    /**
     * @brief Get number of global grid points in the z direction
     * @return number of global grid points in the z direction
     */
    inline int z() const { return n_vec_[2]; }

    /**
     * @brief      returns the total field's size array including buffer space
     *
     * @return     n_vec_
     */
    inline std::array<int,3> n_vec() const {return n_vec_;}
    /**
     * @brief      returns the local process' field's size array including buffer space
     *
     * @return     ;n_vec_
     */
    inline std::array<int,3> ln_vec() const {return ln_vec_;}
    /**
     * @brief Get number of local grid points in the x direction
     * @return number of local grid points in the x direction
     */
    inline int local_x() const { return ln_vec_[0]; }

    /**
     * @brief Get number of local grid points in the y direction
     * @return number of local grid points in the y direction
     */
    inline int local_y() const { return ln_vec_[1]; }

    /**
     * @brief Get number of local grid points in the z direction
     * @return number of local grid points in the z direction
     */
    inline int local_z() const { return ln_vec_[2]; }

    /**
     * @brief Assessor to step size in x direction
     * @return d_[0]
     */
    inline double dx() const { return d_[0]; }
    /**
     * @brief Assessor to step size in y direction
     * @return d_[1]
     */
    inline double dy() const { return d_[1]; }
    /**
     * @brief Assessor to step size in z direction
     * @return d_[2]
     */
    inline double dz() const { return d_[2]; }

    inline std::array<double,3> d() const {return d_;}

    /**
     * @brief Accessor to procLoc
     * @return procLoc
     */
    inline std::array<int,3> procLoc() {return procLoc_;}

    /**
     * @brief Accessor to procR_C_
     * @return procR_C_
     */
    inline std::array<int,3> procR_C() {return procR_C_;}


    /**
     * @brief Fills the local grid with a specified value
     *
     * @param a value to fill the grid with
     */
    inline void fill(const T a) { std::fill_n(local_.get(), ln_vec_[0]*ln_vec_[1]*ln_vec_[2], a); }
    /**
     * @brief Fills the local grid with the process rank
     *
     */
    inline void fill_rank() { fill(static_cast<T>(gridComm_.rank())); }
    /**
     * @brief Fills the local grid with zeros
     *
     */
    inline void zero() { fill(static_cast<T>(0.0)); }

    /**
     * @brief      returns true if periodic boundary conditions are used
     *
     * @return     PBC_
     */
    inline bool PBC() {return PBC_;}
    /**
     * @brief      Accessor function to the grid's MPIInterface
     *
     * @return     gridComm_
     */
    inline mpiInterface gridComm(){return gridComm_;}
    /**
     * @brief Send boundary data between local grids
     * @details Sends data from the current process to update the boundary of a neighboring process
     *
     * @param sendDir Direction of Communication
     */
    void sendDat(PROC_DIR sendDir)
    {
        switch(sendDir)
        {
            case PROC_DIR::UP:
                gridComm_.send(upProc_, upSendTag_, &point(0,upSendIndex_, 0), (ln_vec_[0]) * (ln_vec_[2]));
                break;
            case PROC_DIR::DOWN:
                gridComm_.send(downProc_, downSendTag_, &point(0,1,0), (ln_vec_[0]) * (ln_vec_[2]));
                break;
            case PROC_DIR::NONE:
                break;
            default:
                throw std::logic_error("You made it to the default, you added directions but did not add the cases in sendDat in ParlellGrid.hpp");
        }
    }
    /**
     * @brief Recv boundary data between local grids
     * @details Receives data from a neighboring process to update the boundary of its process
     *
     * @param sendDir Direction of Communication
     */
    void recvDat(PROC_DIR recvDir)
    {
        switch(recvDir)
        {
            case PROC_DIR::UP:
                gridComm_.recv(upProc_, upRecvTag_, &point(0, upSendIndex_+1,0), ln_vec_[0]*ln_vec_[2]);
                break;
            case PROC_DIR::DOWN:
                gridComm_.recv(downProc_, downRecvTag_, &point(0,0,0), ln_vec_[0]*ln_vec_[2]);
                break;
            case PROC_DIR::NONE:
                break;
            default:
                throw std::logic_error("You made it to the default, you added directions but did not add the cases in sendDat in ParlellGrid.hpp");
        }
    }
    /**
     * @brief Update all the local boundaries
     * @details Update each boundary direction
     */
    void transferDat()
    {
        for(int ii = 0; ii < sendList_.size(); ii++)
        {
            switch(sendList_[ii])
            {
                case PROC_DIR::UP:
                    reqs_[ii*2] = gridComm_.isend(  upProc_,   upSendTag_, &point(0,upSendIndex_,0), ln_vec_[0]*ln_vec_[2]);
                    break;
                case PROC_DIR::DOWN:
                    reqs_[ii*2] = gridComm_.isend(downProc_, downSendTag_, &point(0,1,0), ln_vec_[0]*ln_vec_[2]);
                    break;
                case PROC_DIR::NONE:
                    break;
                default:
                    throw std::logic_error("You made it to the default, you added directions but did not add the cases in sendDat in ParlellGrid.hpp");
            }
            switch(recvList_[ii])
            {
                case PROC_DIR::UP:
                    reqs_[ii*2+1] = gridComm_.irecv(  upProc_,   upRecvTag_, &point(0, upSendIndex_+1,0), ln_vec_[0]*ln_vec_[2]);
                    break;
                case PROC_DIR::DOWN:
                    reqs_[ii*2+1] = gridComm_.irecv(downProc_, downRecvTag_, &point(0,0,0), ln_vec_[0]*ln_vec_[2]);
                    break;
                case PROC_DIR::NONE:
                    break;
                default:
                    throw std::logic_error("You made it to the default, you added directions but did not add the cases in sendDat in ParlellGrid.hpp");
            }
            mpi::wait_all(reqs_.data(), reqs_.data() + reqs_.size() );
        }
    }

    /**
     * @brief Returns prow and local x offset for ith x coordinate
     *
     * @param i x coordinate value
     * @return mypX and local x coordinate offset for ith x coordinate
     */
    std::pair<int, int> locate_x(const int i)
    {
        int ii = 0;
        while (ii < xTrans_.size()-1 && i >= xTrans_[ii+1])
            ++ii;
        return {ii, i-xTrans_[ii]};
    }

    /**
     * @brief Returns prow and local y offset for jth y coordinate
     *
     * @param j y coordinate value
     * @return mypY and local y coordinate offset for jth y coordinate
     */
    std::pair<int, int> locate_y(const int j)
    {
        int jj = 0;
        while (jj < yTrans_.size()-1 && j >= yTrans_[jj+1])
            jj++;
        return {jj, j-yTrans_[jj]};
    }

    /**
     * @brief Returns prow and local z offset for kth z coordinate
     *
     * @param k z coordinate value
     * @return mypZ and local z coordinate offset for kth z coordinate
     */
    std::pair<int, int> locate_z(const int k)
    {
        int kk = 0;
        while (kk < zTrans_.size()-1 && k >= zTrans_[kk+1])
            kk++;
        return {kk, k-zTrans_[kk]};
    }

    /**
     * @brief Returns prow and local x offset for ith x coordinate
     *
     * @param i x coordinate value
     * @return mypX and local x coordinate offset for ith x coordinate
     */
    std::pair<int, int> locate_x_no_boundaries(const int i)
    {
        int ii = 0;
        while (ii < xTrans_.size()-1 && i >= xTrans_[ii+1]-2*(ii+1) )
            ii++;
        return {ii, i-xTrans_[ii]};
    }

    /**
     * @brief Returns prow and local y offset for jth y coordinate
     *
     * @param j y coordinate value
     * @return mypY and local y coordinate offset for jth y coordinate
     */
    std::pair<int, int> locate_y_no_boundaries(const int j)
    {
        int jj = 0;
        while (jj < yTrans_.size()-1 && j >= yTrans_[jj+1]-2*(jj+1) )
            jj++;
        return {jj, j-yTrans_[jj]};
    }

    /**
     * @brief Returns prow and local z offset for kth z coordinate
     *
     * @param k z coordinate value
     * @return mypY and local z coordinate offset for kth z coordinate
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
     * @return     The locs proc.
     */
    int getLocsProc(int xx, int yy, int zz)
    {
        int py, off;
        std::tie(py, off) = locate_y(yy);

        return py;
    }

    /**
     * @brief      Returns the process rank storing a particular point excluding all boundary buffers
     *
     * @param[in]  xx    the point coordinate in the x direction
     * @param[in]  yy    the point coordinate in the y direction
     * @param[in]  zz    the point coordinate in the z direction
     *
     * @return     The locs proc.
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
     * @return     The locs proc.
     */
    int getLocsProc(int xx, int yy)
    {
        int py, off;
        std::tie(py, off) = locate_y(yy);

        return py;
    }

    /**
     * @brief      Returns the process rank storing a particular point excluding all boundary buffers
     *
     * @param[in]  xx    the point coordinate in the x direction
     * @param[in]  yy    the point coordinate in the y direction
     *
     * @return     The locs proc.
     */
    int getLocsProc_no_boundaries(int xx, int yy)
    {
        int py, off;

        std::tie(py, off) = locate_y_no_boundaries(yy);

        return py;
    }

    /**
     * @brief returns values of the grid along a particular row
     *
     * @param j The y coordinate for where to take the plane
     * @return A vector of all the values in a particular row
     */
    std::vector<T> getPlaneXZ(const int j)
    {
        int py, off;
        std::tie(py, off) = locate_y(j);
        std::vector<T> planeXZ(n_vec_[0]*n_vec_[2], 0.0);
        if(gridComm_.rank() == py)
            std::copy_n(&point(0,off,0), ln_vec_[0]*ln_vec_[2], planeXZ.data());

        mpi::broadcast(gridComm_, planeXZ, py);
        return planeXZ;
    }

    /**
     * @brief returns values of the grid along a particular column
     *
     * @param k the z coordinate where to take the plane
     * @return A vector of all the values in a particular column
     */
    std::vector<T> getPlaneXY(const int k)
    {
        int pz, off;
        std::tie(pz, off) = locate_z(k);
        std::vector<T> planeXY(n_vec_[0]*n_vec_[1], 0.0);

        if(gridComm_.rank() != pz)
        {
            std::vector<T> toSend(ln_vec_[1]*ln_vec_[0], 0.0);
            for(int yy = 0; yy < ln_vec_[1]; ++yy)
                std::copy_n( &point(0, yy, off), ln_vec_[0], &toSend[ln_vec_[0]*yy] );
                // copy_(ln_vec_[0], &point(0, yy, off), 1, &toSend[ln_vec_[0]*yy], 1);
            gridComm_.send(pz, gridComm_.cantorTagGen(gridComm_.rank(), pz, 2, 0), toSend);
        }
        else
        {
            for(int yy = 0; yy < ln_vec_[1]; ++yy)
                std::copy_n( &point(0, yy, off), ln_vec_[0], &planeXY[ln_vec_[0]*yy] );
                // copy_(ln_vec_[0], &point(0, yy, off), 1, &planeXY[ln_vec_[0]*yy], 1);
            int spot = ln_vec_[0]*ln_vec_[1];
            for(int ii = 0; ii < gridComm_.size(); ii ++)
            {
                if(ii != pz )
                {
                    std::vector<T> tempVec;
                    gridComm_.recv(ii, gridComm_.cantorTagGen(ii, gridComm_.rank(), 2, 0), tempVec);
                    std::copy_n( tempVec.data(), tempVec.size(), &planeXY[spot] );
                    // copy_(tempVec.size(), tempVec.data(), 1, &planeXY[spot], 1);
                    spot += tempVec.size();
                }
            }
        }
        mpi::broadcast(gridComm_, planeXY, pz);
        return planeXY;
    }

    /**
     * @brief returns values of the grid along a particular column
     *
     * @param i the x coordinate of where to get the plane
     * @return A vector of all the values in a particular column
     */
    std::vector<T> getPlaneYZ(const int i)
    {
        gridComm_.barrier();
        int px, off;
        std::tie(px, off) = locate_x(i);
        std::vector<T> planeYZ(n_vec_[1]*n_vec_[2], 0.0);

        gridComm_.barrier();
        if(gridComm_.rank() != px)
        {
            // gridComm_.barrier();
            std::vector<T> toSend(ln_vec_[1]*ln_vec_[2], 0.0);
            // gridComm_.barrier();
            for(int yy = 0; yy < ln_vec_[1]; ++yy)
                copy_(ln_vec_[2], &point(off, yy, 0), ln_vec_[0], &toSend[ln_vec_[2]*yy], 1);
            // gridComm_.barrier();
            gridComm_.send(px, gridComm_.cantorTagGen(gridComm_.rank(), px, 2, 0), toSend);
            // gridComm_.barrier();
        }
        else
        {
            // gridComm_.barrier();
            for(int yy = 0; yy < ln_vec_[1]; ++yy)
                copy_(ln_vec_[2], &point(off, yy, 0), ln_vec_[0], &planeYZ[ln_vec_[2]*yy], 1);
            // gridComm_.barrier();
            int spot = ln_vec_[2]*ln_vec_[1];
            for(int ii = 0; ii < gridComm_.size(); ii ++)
            {
                if(ii != px )
                {
                    std::vector<T> tempVec;
                    gridComm_.recv(ii, gridComm_.cantorTagGen(ii, gridComm_.rank(), 2, 0), tempVec);
                    std::copy_n(tempVec.begin(), tempVec.size(), &planeYZ[spot]);
                    // copy_(tempVec.size(), tempVec.data(), 1, &planeYZ[spot], 1);
                    spot += tempVec.size();
                }
            }
            // gridComm_.barrier();
        }
        gridComm_.barrier();
        mpi::broadcast(gridComm_, planeYZ, px);
        gridComm_.barrier();
        return planeYZ;
    }

    template <typename U> friend std::ostream &operator<<(std::ostream &out, const parallelGrid <U> &o);
};

class parallelGridReal : public parallelGrid<double>
{
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
    parallelGridReal(mpiInterface & gridComm, bool PBC, std::array<int,3> n_vec, std::array<double,3> d, bool ylim=false);
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
    parallelGridReal(mpiInterface & gridComm, bool PBC, std::vector<real_grid_ptr> weights, std::array<int,3> n_vec, std::array<double,3> d, bool ylim=false);
};

class parallelGridCplx : public parallelGrid<std::complex<double>>
{
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
    parallelGridCplx(mpiInterface & gridComm, bool PBC, std::array<int,3> n_vec, std::array<double,3> d, bool ylim=false);
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
    parallelGridCplx(mpiInterface & gridComm, bool PBC, std::vector<real_grid_ptr> weights, std::array<int,3> n_vec, std::array<double,3> d, bool ylim=false);
};

class parallelGridInt : public parallelGrid<int>
{
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
    parallelGridInt(mpiInterface & gridComm, bool PBC, std::array<int,3> n_vec, std::array<double,3> d, bool ylim=false);
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
    parallelGridInt(mpiInterface & gridComm, bool PBC, std::vector<real_grid_ptr> weights, std::array<int,3> n_vec, std::array<double,3> d, bool ylim=false);
};



/**
 * @details stream printer
 *
 * @param out output stream
 * @param o grid to be outputted
 */
template <typename T>
std::ostream &operator<<(std::ostream &out, const parallelGrid <T> &o)
{
    for (int yy = 0; yy < std::min(40,int(o.ln_vec_[1])); yy++)
    {
        for (int xx = 0; xx < std::min(40,int(o.ln_vec_[0])); xx++)
        {
            out << std::setprecision(3) << o(xx,yy) << "\t";
        }
        out << "\n";
    }
    out << std::endl;
    return out;
};

#endif