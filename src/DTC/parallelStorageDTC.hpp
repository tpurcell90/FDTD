#ifndef FDTD_PARALLELDETECTORSTORAGE
#define FDTD_PARALLELDETECTORSTORAGE

#include <DTC/parallelStorageDTCSructs.hpp>
#include <src/GRID/parallelGrid.hpp>
#include <src/UTIL/ml_consts.hpp>
#include <cstdio>
#include <UTIL/typedefs.hpp>

namespace mpi = boost::mpi;

/**
* @brief The data storage class for the detectors
* @details Creates the output grid on a single processor that then can be used to output the full data set into a file in the correct order
*
*/
template <typename T> class parallelStorageDTC
{
protected:
    bool masterBool_;
    std::array<int,3> loc_; //!< Grid point location of the lower left corner of the detector (full grid space)
    std::array<int,3> sz_; //!< number of grid points in each direction the detector is storing
    std::vector<std::shared_ptr<slaveProcInfo>> master_; //!< A shared pointer that is used to generate the output data grod. If the process is not master this is set to a nullptr.
    std::shared_ptr<slaveProcDtc> slave_; //!< A shared pointer for the slaveProcDtc struct, set to a nullptr if the detector area does not include the process
    std::shared_ptr<copyProcDtc> toOutGrid_; //!< used to transfer from actual grid to the outGrid
    std::shared_ptr<parallelGrid<T>> grid_; //!< A shared pointer to a grid that the detector will need to generate the correct output
    mpiInterface gridComm_; //!< The communicator for the grid and detector (must be the same communicator)

    std::vector<T> scratch_; //!< vector for scratch space

    std::shared_ptr<Grid<T>> outGrid_; //!< A shared pointer to the serial output grid (what is used to make the file)
public:

    /**
     * @brief Constructor of parallelDTC
     * @param dtcNum the detector number
     * @param grid pointer to output grid
     * @param loc location of lower left back corner of the dtc
     * @param sz size in grid points for the dtc
     */
     parallelStorageDTC(int dtcNum, std::shared_ptr<parallelGrid<T>> grid, std::array<int,3> loc, std::array<int,3> sz) :
        masterBool_(false),
        loc_(loc),
        sz_(sz),
        grid_(grid),
        gridComm_(grid_->gridComm()),
        scratch_(std::accumulate(sz_.begin(), sz_.end(), 1, std::multiplies<int>()),0.0)
    {
        genDatStruct(dtcNum);

        /// Output grid only made in the master process
        if(masterBool_ || toOutGrid_ )
            outGrid_ = std::make_shared<Grid<T>>( sz_, grid_->d()  );
        else
            outGrid_ = nullptr;
    }
    /**
     * @return  reference to loc_
     */
    inline std::vector<int> &loc() {return loc_;}

    /**
     * @return reference to sz_
     */
    inline std::vector<int> &sz() {return sz_;}

    /**
     * @return reference to outGrid_
     */
    inline std::shared_ptr<Grid<T>> &outGrid() { return outGrid_; }

    inline std::shared_ptr<parallelGrid<T>> grid() {return grid_; }

    /**
     * @brief Generates the slave and master process data structures so the output grid is constructed correctly.
     * @details use size and location of the dtc to determine what information to send where
     *
     * @param the index of the detector in the FDTDField dtcArr_
     */
    void genDatStruct(int dtcNum)
    {
        /// initialize slave and master to null pointers
        slave_ = nullptr;
        master_ = {};

        // Generate the slaveProcDtc and initialize all variables in the slave process
        slaveProcDtc slaveProc;
        slaveProc.masterProc_ = grid_->getLocsProc(loc_[0], loc_[1], loc_[2]);
        // if(gridComm_.rank() == slaveProc.masterProc_ )
        slaveProc.stride_ = 0;
        slaveProc.sz_ = {0,0,0};
        slaveProc.loc_ = {-1, -1, -1};
        if(gridComm_.rank() != slaveProc.masterProc_)
        {
            if(loc_[0] >= grid_->procLoc()[0] && loc_[0] < grid_->procLoc()[0] + grid_->local_x() -2) //!< Is the process in the process col that contains the left boundary of the detector region?
            {
                slaveProc.loc_[0] = loc_[0] - grid_->procLoc()[0] + 1;

                if(loc_[1] >= grid_->procLoc()[1] && loc_[1] < grid_->procLoc()[1] + grid_->local_y() -2) //!< Is the process in the process row that contains the lower boundary of the detector region?
                {
                    slaveProc.loc_[1] = loc_[1] - grid_->procLoc()[1] + 1;
                    if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                    {
                        slaveProc.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                    }
                    else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                    {
                        slaveProc.loc_[2] = 1;
                    }
                }
                else if(loc_[1] < grid_->procLoc()[1] && loc_[1] + sz_[1] > grid_->procLoc()[1]) //!< Is the process in a process row that the detector region covers?
                {
                    slaveProc.loc_[1] = 1;
                    if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                    {
                        slaveProc.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                    }
                    else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                    {
                        slaveProc.loc_[2] = 1;
                    }
                }
            }
            else if(loc_[0] < grid_->procLoc()[0] && loc_[0] + sz_[0] > grid_->procLoc()[0]) //!< Is the process in a process col that the detector covers?
            {
                slaveProc.loc_[0] = 1;

                if(loc_[1] >= grid_->procLoc()[1] && loc_[1] < grid_->procLoc()[1] + grid_->local_y() -2)
                {
                    slaveProc.loc_[1] = loc_[1] - grid_->procLoc()[1] + 1;
                    if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                    {
                        slaveProc.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                    }
                    else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                    {
                        slaveProc.loc_[2] = 1;
                    }
                }
                else if(loc_[1] < grid_->procLoc()[1] && loc_[1] + sz_[1] > grid_->procLoc()[1])
                {
                    slaveProc.loc_[1] = 1;
                    if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                    {
                        slaveProc.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                    }
                    else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                    {
                        slaveProc.loc_[2] = 1;
                    }
                }
            }
            if(grid_->local_z() == 1)
            {
                slaveProc.loc_[2] = 0;
            }
            // If either slaveProc.loc_ values is 0 then the detector is not in the process
            if(slaveProc.loc_[0] != -1 && slaveProc.loc_[1] != -1 && slaveProc.loc_[2] != -1  )
            {
                if(sz_[0] + loc_[0] > grid_->procLoc()[0] + grid_->local_x() - 2) // Does the detector go through the end of the process' grid?
                {
                    slaveProc.sz_[0] = grid_->local_x() - slaveProc.loc_[0] - 1;
                }
                else
                {
                    slaveProc.sz_[0] = loc_[0] + sz_[0] - (grid_->procLoc()[0] + slaveProc.loc_[0] - 1);
                }

                if(sz_[1] + loc_[1] > grid_->procLoc()[1] + grid_->local_y() - 2)
                {
                    slaveProc.sz_[1] = grid_->local_y() - slaveProc.loc_[1] - 1;
                }
                else
                {
                    slaveProc.sz_[1] = loc_[1] + sz_[1] - (grid_->procLoc()[1] + slaveProc.loc_[1] - 1);
                }

                if(grid_->local_z() != 1)
                {
                    if(sz_[2] + loc_[2] > grid_->procLoc()[2] + grid_->local_z() - 2)
                        slaveProc.sz_[2] = grid_->local_z() - slaveProc.loc_[2] - 1;
                    else
                        slaveProc.sz_[2] = loc_[2] + sz_[2] - (grid_->procLoc()[2] + slaveProc.loc_[2] - 1);
                }
                else
                {
                    slaveProc.sz_[2] = 1;
                }

                if(isamax_(slaveProc.sz_.size(), slaveProc.sz_.data(), 1) - 1 == 0 )
                {
                    slaveProc.stride_ =  1;
                    slaveProc.opSz_ = {{ slaveProc.sz_[0] , slaveProc.sz_[1], slaveProc.sz_ [2] }};
                    if(slaveProc.loc_[2] == -1)
                    {
                        slaveProc.opSz_[2] = 1;
                        slaveProc.loc_[2] = 0;
                    }
                    slaveProc.addVec1_ = {{ 0, 1, 0 }};
                    slaveProc.addVec2_ = {{ 0, 0, 1 }};
                }
                else if(isamax_(slaveProc.sz_.size(), slaveProc.sz_.data(), 1) - 1 == 1 )
                {
                    slaveProc.stride_ =  grid_->local_x() * grid_->local_z();
                    slaveProc.opSz_ = {{ slaveProc.sz_[1] , slaveProc.sz_[0], slaveProc.sz_ [2] }};
                    if(slaveProc.loc_[2] == -1)
                    {
                        slaveProc.opSz_[2] = 1;
                        slaveProc.loc_[2] = 0;
                    }
                    slaveProc.addVec1_ = {{ 1, 0, 0 }};
                    slaveProc.addVec2_ = {{ 0, 0, 1 }};
                }
                else
                {
                    slaveProc.stride_ =  grid_->local_x();
                    slaveProc.opSz_ = {{ slaveProc.sz_[2] , slaveProc.sz_[0], slaveProc.sz_ [1] }};
                    if(slaveProc.loc_[2] == -1)
                        throw std::logic_error("The max size is in a direction that is not defined, for 2D calcs set sz[2] to 0");
                    slaveProc.addVec1_ = {{ 1, 0, 0 }};
                    slaveProc.addVec2_ = {{ 0, 1, 0 }};
                }
                slave_ = std::make_shared<slaveProcDtc>(slaveProc); //!< only make slave active if necessary
            }
        }
        else
        {
            copyProcDtc toOutGrid;
            toOutGrid.strideOutGrid_ = 0;
            toOutGrid.stride_ = 0;

            toOutGrid.locOutGrid_ = {0,0};
            toOutGrid.sz_ = {0,0,0};
            toOutGrid.loc_ = {-1,-1,-1};
            if(loc_[0] >= grid_->procLoc()[0] && loc_[0] < grid_->procLoc()[0] + grid_->local_x() -2) //!< Is the process in the process col that contains the left boundary of the detector region?
            {
                toOutGrid.loc_[0] = loc_[0] - grid_->procLoc()[0] + 1;

                if(loc_[1] >= grid_->procLoc()[1] && loc_[1] < grid_->procLoc()[1] + grid_->local_y() -2) //!< Is the process in the process row that contains the lower boundary of the detector region?
                {
                    toOutGrid.loc_[1] = loc_[1] - grid_->procLoc()[1] + 1;
                    if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                        toOutGrid.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                    else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                        toOutGrid.loc_[2] = 1;
                }
                else if(loc_[1] < grid_->procLoc()[1] && loc_[1] + sz_[1] > grid_->procLoc()[1]) //!< Is the process in a process row that the detector region covers?
                {
                    toOutGrid.loc_[1] = 1;
                    if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                        toOutGrid.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                    else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                        toOutGrid.loc_[2] = 1;
                }
            }
            else if(loc_[0] < grid_->procLoc()[0] && loc_[0] + sz_[0] > grid_->procLoc()[0]) //!< Is the process in a process col that the detector covers?
            {
                toOutGrid.loc_[0] = 1;

                if(loc_[1] >= grid_->procLoc()[1] && loc_[1] < grid_->procLoc()[1] + grid_->local_y() -2)
                {
                    toOutGrid.loc_[1] = loc_[1] - grid_->procLoc()[1] + 1;
                    if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                        toOutGrid.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                    else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                        toOutGrid.loc_[2] = 1;
                }
                else if(loc_[1] < grid_->procLoc()[1] && loc_[1] + sz_[1] > grid_->procLoc()[1])
                {
                    toOutGrid.loc_[1] = 1;
                    if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                        toOutGrid.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                    else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                        toOutGrid.loc_[2] = 1;
                }
            }
            if(grid_->local_z() == 1)
            {
                toOutGrid.loc_[2] = 0;
            }
            // If either toOutGrid.loc_ values is 0 then the detector is not in the process
            if(toOutGrid.loc_[0] != -1 && toOutGrid.loc_[1] != -1 && toOutGrid.loc_[2] != -1)
            {
                if(sz_[0] + loc_[0] > grid_->procLoc()[0] + grid_->local_x() - 2) // Does the detector go through the end of the process' grid?
                    toOutGrid.sz_[0] = grid_->local_x() - toOutGrid.loc_[0] - 1;
                else
                    toOutGrid.sz_[0] = loc_[0] + sz_[0] - (grid_->procLoc()[0] + toOutGrid.loc_[0] - 1);

                if(sz_[1] + loc_[1] > grid_->procLoc()[1] + grid_->local_y() - 2)
                    toOutGrid.sz_[1] = grid_->local_y() - toOutGrid.loc_[1] - 1;
                else
                    toOutGrid.sz_[1] = loc_[1] + sz_[1] - (grid_->procLoc()[1] + toOutGrid.loc_[1] - 1);

                if(grid_->local_z() != 1)
                {
                    if(sz_[2] + loc_[2] > grid_->procLoc()[2] + grid_->local_z() - 2)
                        toOutGrid.sz_[2] = grid_->local_z() - toOutGrid.loc_[2] - 1;
                    else
                        toOutGrid.sz_[2] = loc_[2] + sz_[2] - (grid_->procLoc()[2] + toOutGrid.loc_[2] - 1);
                }
                else
                    toOutGrid.sz_[2] = 1;


                toOutGrid.locOutGrid_ = {grid_->procLoc()[0] + toOutGrid.loc_[0] - 1 - loc_[0], grid_->procLoc()[1] + toOutGrid.loc_[1] - 1 - loc_[1], grid_->local_z() == 1 ? 0 : grid_->procLoc()[2] + toOutGrid.loc_[2] - 1 - loc_[2]}; //!< make loc relative to process grid
                if(isamax_(toOutGrid.sz_.size(), toOutGrid.sz_.data(), 1) - 1 == 0 || (toOutGrid.sz_[0] == toOutGrid.sz_[1] && toOutGrid.sz_[0] <= toOutGrid.sz_[2]) || (toOutGrid.sz_[0] == toOutGrid.sz_[2] && toOutGrid.sz_[0] <= toOutGrid.sz_[1]) )
                {
                    toOutGrid.stride_ =  1;
                    toOutGrid.opSz_ = {{ toOutGrid.sz_[0] , toOutGrid.sz_[1], toOutGrid.sz_ [2] }};
                    if(toOutGrid.loc_[2] == -1)
                    {
                        toOutGrid.loc_[2] = 0;
                        toOutGrid.opSz_[2] = 1;
                    }
                    toOutGrid.addVec1_ = {{ 0, 1, 0 }};
                    toOutGrid.addVec2_ = {{ 0, 0, 1 }};

                    toOutGrid.szOutGrid_ = {toOutGrid.sz_[0], toOutGrid.sz_[1], toOutGrid.sz_[2]};
                    toOutGrid.strideOutGrid_ = 1;
                }
                else if(isamax_(toOutGrid.sz_.size(), toOutGrid.sz_.data(), 1) - 1 == 1 || (toOutGrid.sz_[1] == toOutGrid.sz_[2]) )
                {
                    toOutGrid.stride_ =  grid_->local_x() * grid_->local_z();
                    toOutGrid.opSz_ = {{ toOutGrid.sz_[1] , toOutGrid.sz_[0], toOutGrid.sz_ [2] }};
                    if(toOutGrid.loc_[2] == -1)
                    {
                        toOutGrid.loc_[2] = 0;
                        toOutGrid.opSz_[2] = 1;
                    }
                    toOutGrid.addVec1_ = {{ 1, 0, 0 }};
                    toOutGrid.addVec2_ = {{ 0, 0, 1 }};

                    toOutGrid.szOutGrid_ = {toOutGrid.sz_[1], toOutGrid.sz_[0], toOutGrid.sz_[2]};
                    toOutGrid.strideOutGrid_ = sz_[0]*sz_[2];
                }
                else
                {
                    toOutGrid.stride_ =  grid_->local_x();
                    toOutGrid.opSz_ = {{ toOutGrid.sz_[2] , toOutGrid.sz_[0], toOutGrid.sz_ [1] }};
                    if(toOutGrid.loc_[2] == -1)
                        throw std::logic_error("The max size is in a direction that is not defined, for 2D calcs set sz[2] to 0");
                    toOutGrid.addVec1_ = {{ 1, 0, 0 }};
                    toOutGrid.addVec2_ = {{ 0, 1, 0 }};

                    toOutGrid.szOutGrid_ = {toOutGrid.sz_[2], toOutGrid.sz_[0], toOutGrid.sz_[1]};
                    toOutGrid.strideOutGrid_ = sz_[0];
                }
                toOutGrid_ = std::make_shared<copyProcDtc>(toOutGrid);
            }
        }
        slaveProcInfo toMaster;
        if(slave_)
        {
            toMaster.slaveProc_ = gridComm_.rank();
            toMaster.loc_ = {grid_->procLoc()[0] + slave_->loc_[0] - 1 - loc_[0], grid_->procLoc()[1] + slave_->loc_[1] - 1 - loc_[1], grid_->procLoc()[2] + slave_->loc_[2] - 1 - loc_[2]}; //!< make loc relative to process grid
            if(grid_->local_z() == 1)
                toMaster.loc_[2]  = -1;

            if(isamax_(slave_->sz_.size(), slave_->sz_.data(), 1) == 0 )
            {
                toMaster.stride_ =  1;
                toMaster.sz_ = {{ slave_->sz_[0] , slave_->sz_[1], slave_->sz_[2] }};
                if(slave_->loc_[2] == -1)
                {
                    toMaster.loc_[2] = 0;
                    toMaster.sz_[2] = 1;
                }
                toMaster.addVec1_ = {{ 0, 1, 0 }};
                toMaster.addVec2_ = {{ 0, 0, 1 }};
            }
            else if(isamax_(slave_->sz_.size(), slave_->sz_.data(), 1) == 1 )
            {
                toMaster.stride_ =  grid_->local_x() * grid_->local_z();
                toMaster.sz_ = {{ slave_->sz_[1] , slave_->sz_[0], slave_->sz_ [2] }};
                if(toMaster.loc_[2] == -1)
                {
                    toMaster.loc_[2] = 0;
                    toMaster.sz_[2] = 1;
                }
                toMaster.addVec1_ = {{ 0, 0, 1 }};
                toMaster.addVec2_ = {{ 1, 0, 0 }};
            }
            else
            {
                toMaster.stride_ =  grid_->local_x();
                toMaster.sz_ = {{ slave_->sz_[2] , slave_->sz_[0], slave_->sz_ [1] }};
                if(toMaster.loc_[2] == -1)
                    throw std::logic_error("The max size is in a direction that is not defined, for 2D calcs set sz[2] to 0");
                toMaster.addVec1_ = {{ 1, 0, 0 }};
                toMaster.addVec2_ = {{ 0, 1, 0 }};
            }
        }
        else
            toMaster.slaveProc_ = -1;
        gridComm_.barrier();
        if(gridComm_.rank() == slaveProc.masterProc_)
        {
            masterBool_ = true;
            std::vector<slaveProcInfo> allProcs;
            mpi::gather(gridComm_, toMaster, allProcs, slaveProc.masterProc_);
            for(auto & proc : allProcs)
            {
                if(proc.slaveProc_ != -1)
                    master_.push_back(std::make_shared<slaveProcInfo>(proc) );
            }

            // master_ = std::make_shared<std::vector<slaveProcInfo>>(masterProc);
        }
        else
            mpi::gather(gridComm_, toMaster,  slaveProc.masterProc_);

        return;
    }
    /**
     * @brief Communicate the fields to the detector output
     * @details The slave processes send their data to the master and the master moves it to the output grid
     *
     */
    virtual void getField() = 0;

    inline bool master() { return masterBool_; }
};

class parallelStorageDTCReal : public parallelStorageDTC<double>
{
public:
    /**
     * @brief Constructor of parallelStorageDTC
     * @param dtcNum the detector number
     * @param grid pointer to output grid
     * @param loc location of lower left back corner of the dtc
     * @param sz size in grid points for the dtc
     */
    parallelStorageDTCReal(int dtcNum, real_pgrid_ptr grid, std::array<int,3> loc, std::array<int,3> sz);
    void getField();
};

/**
 * Complex version of FluxDTC see base class for more descriptions
 */
class parallelStorageDTCCplx : public parallelStorageDTC<cplx>
{
public:
    /**
     * @brief Constructor of parallelStorageDTC
     * @param dtcNum the detector number
     * @param grid pointer to output grid
     * @param loc location of lower left back corner of the dtc
     * @param sz size in grid points for the dtc
     */
    parallelStorageDTCCplx(int dtcNum, cplx_pgrid_ptr grid, std::array<int,3> loc, std::array<int,3> sz);

    /**
     * @brief Communicate the fields to the detector output
     * @details The slave processes send their data to the master and the master moves it to the output grid
     *
     */
    void getField();
};

#endif