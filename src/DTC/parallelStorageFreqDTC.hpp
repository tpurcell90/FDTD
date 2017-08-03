#ifndef FDTD_PARALLELDETECTORSTORAGEFREQ
#define FDTD_PARALLELDETECTORSTORAGEFREQ

#include <DTC/parallelStorageDTCSructs.hpp>
#include <src/GRID/parallelGrid.hpp>
#include <src/UTIL/ml_consts.hpp>
#include <UTIL/typedefs.hpp>
#include <cstdio>

namespace mpi = boost::mpi;

template <typename T> class parallelStorageFreqDTC
{
protected:
    char noTranspose_; //!< char for mkl functions
    char transpose_; //!< char for mkl functions
    mpiInterface gridComm_; //!< mpi interface for all mpi calls
    int nfreq_; //!<  number of frequencies to detect
    int zgemmK_; //!< k value for mkl functions
    int shift_j_; //!< 0 if j = outGrid dir 1; 1 if j is outGrid dir 2
    int shift_k_; //!< 0 if k = outGrid dir 1; 1 if k is outGrid dir 2
    cplx ONE_; //!< 1 for mkl functions
    std::array<int,3> loc_; //!< location  of lower left back corner of detection region
    std::array<int,3> sz_; //!< size of the detection region in grid points
    std::vector<double> freqList_; //!< list of all frequencies
    std::vector<cplx> fftFact_; //!< vector storing the exp(i $\omg$ t) values each time step
    std::vector<cplx> fIn_; //!< vector storing the field input values
    std::shared_ptr<std::vector<slaveProcInfo>> master_; //!< parameters for the master process to take in all the slave processes info and combine it
    std::shared_ptr<slaveProcDtc> slave_; //!< parameters for slave processes to get and send info to master
    std::shared_ptr<copyProcDtc> toOutGrid_; //!< a copy param set if master also needs to get info
    std::shared_ptr<parallelGrid<T>> grid_; //!< grid pointer to the field that is being stored

    std::vector<T> scratch_; //!< scratch space

    cplx_grid_ptr outGrid_; //!< the output gird storage

    std::shared_ptr<fInParam> fieldInFreq_; //!< parameters to take in field
public:
    /**
     * @brief      construct the frequency storage dtc
     *
     * @param[in]  dtcNum    the detector number
     * @param[in]  grid      pointer to output grid
     * @param[in]  propDir   The direction of propagation
     * @param[in]  loc       location of lower left corner of the dtc
     * @param[in]  sz        size in grid points for the dtc
     * @param[in]  freqList  The frequency list
     */
    parallelStorageFreqDTC(int dtcNum, std::shared_ptr<parallelGrid<T>> grid, DIRECTION propDir, std::array<int,3> loc, std::array<int,3> sz, std::vector<double> freqList) :
        gridComm_(grid->gridComm()),
        grid_(grid),
        noTranspose_('N'),
        transpose_('T'),
        freqList_(freqList),
        nfreq_(freqList.size()),
        zgemmK_(1),
        shift_j_(-1),
        shift_k_(-1),
        ONE_(1.0,0.0),
        loc_(loc),
        sz_(sz),
        fftFact_(nfreq_, 0.0)
    {
        genDatStruct(dtcNum, propDir);
        if(fieldInFreq_)
        {
            int szProd = std::accumulate(fieldInFreq_->sz_.begin(), fieldInFreq_->sz_.end(), 1, std::multiplies<int>() );
            fIn_ = std::vector<cplx>(szProd, 0.0);
            scratch_ = std::vector<T>(szProd,0.0);
            if(propDir == DIRECTION::X )
            {
                outGrid_ = std::make_shared<Grid<cplx>>( std::array<int,3>( {{nfreq_, fieldInFreq_->sz_[1], fieldInFreq_->sz_[0] }}) ,  std::array<double,3>( {{freqList_[1] - freqList_[0], grid_->dz(), grid_->dy() }} ) );
            }
            else if(propDir == DIRECTION::Y )
            {
                outGrid_ = std::make_shared<Grid<cplx>>( std::array<int,3>( {{nfreq_, fieldInFreq_->sz_[1], fieldInFreq_->sz_[0] }}) ,  std::array<double,3>( {{freqList_[1] - freqList_[0], grid_->dz(), grid_->dx() }} ) );
            }
            else if(propDir == DIRECTION::Z )
            {
                outGrid_ = std::make_shared<Grid<cplx>>( std::array<int,3>( {{nfreq_, fieldInFreq_->sz_[1], fieldInFreq_->sz_[0] }}) ,  std::array<double,3>( {{freqList_[1] - freqList_[0], grid_->dy(), grid_->dx() }} ) );
            }
            else if(propDir == DIRECTION::NONE )
            {
                outGrid_ = std::make_shared<Grid<cplx>>( std::array<int,3>( {{nfreq_, szProd, 1 }}) ,  std::array<double,3>( {{freqList_[1] - freqList_[0], grid_->dx(), 1.0 }} ) );
            }

        }
        else
            outGrid_ = nullptr;
    }
    /**
     * @return  reference to loc_
     */
    inline std::array<int,3> &loc() {return loc_;}

    /**
     * @return reference to sz_
     */
    inline std::array<int,3> &sz() {return sz_;}

    /**
     * @return reference to outGrid_
     */
    inline cplx_grid_ptr &outGrid() { return outGrid_; }

    /**
     * @brief Generates the slave and master process data structures so the output grid is constructed correctly.
     * @details use size and location of the dtc to determine what information to send where
     *
     * @param the index of the detector in the FDTDField dtcArr_
     */
    void genDatStruct(int dtcNum, DIRECTION propDir)
    {
        fInParam fieldInFreq;
        fieldInFreq.loc_ = {-2, -2, -2};
        fieldInFreq.sz_ = {0,0,0};
        fieldInFreq.stride_ = 0;

        if(loc_[0] >= grid_->procLoc()[0] && loc_[0] < grid_->procLoc()[0] + grid_->local_x() -2) //!< Is the process in the process col that contains the left boundary of the detector region?
        {
            fieldInFreq.loc_[0] = loc_[0] - grid_->procLoc()[0] + 1;

            if(loc_[1] >= grid_->procLoc()[1] && loc_[1] < grid_->procLoc()[1] + grid_->local_y() -2) //!< Is the process in the process row that contains the lower boundary of the detector region?
            {
                fieldInFreq.loc_[1] = loc_[1] - grid_->procLoc()[1] + 1;
                if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                {
                    fieldInFreq.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                }
                else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                {
                    fieldInFreq.loc_[2] = 1;
                }
            }
            else if(loc_[1] < grid_->procLoc()[1] && loc_[1] + sz_[1] > grid_->procLoc()[1]) //!< Is the process in a process row that the detector region covers?
            {
                fieldInFreq.loc_[1] = 1;
                if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                {
                    fieldInFreq.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                }
                else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                {
                    fieldInFreq.loc_[2] = 1;
                }
            }
        }
        else if(loc_[0] < grid_->procLoc()[0] && loc_[0] + sz_[0] > grid_->procLoc()[0]) //!< Is the process in a process col that the detector covers?
        {
            fieldInFreq.loc_[0] = 1;

            if(loc_[1] >= grid_->procLoc()[1] && loc_[1] < grid_->procLoc()[1] + grid_->local_y() -2)
            {
                fieldInFreq.loc_[1] = loc_[1] - grid_->procLoc()[1] + 1;
                if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                {
                    fieldInFreq.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                }
                else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                {
                    fieldInFreq.loc_[2] = 1;
                }
            }
            else if(loc_[1] < grid_->procLoc()[1] && loc_[1] + sz_[1] > grid_->procLoc()[1])
            {
                fieldInFreq.loc_[1] = 1;
                if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
                {
                    fieldInFreq.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
                }
                else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
                {
                    fieldInFreq.loc_[2] = 1;
                }
            }
        }

        // if(loc_[0] >= grid_->procLoc()[0] && loc_[0] < grid_->procLoc()[0] + grid_->local_x() -2) //!< Is the process in the process col that contains the left boundary of the detector region?
        // {
        //     fieldInFreq.loc_[0] = loc_[0] - grid_->procLoc()[0] + 1;
        //     if(loc_[1] >= grid_->procLoc()[1] && loc_[1] < grid_->procLoc()[1] + grid_->local_y() -2) //!< Is the process in the process row that contains the lower boundary of the detector region?
        //     {
        //         fieldInFreq.loc_[1] = loc_[1] - grid_->procLoc()[1] + 1;
        //         if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
        //         {
        //             fieldInFreq.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
        //         }
        //         else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
        //         {
        //             fieldInFreq.loc_[2] = 1;
        //         }
        //     }
        //     else if(loc_[1] < grid_->procLoc()[1] && loc_[1] + sz_[1] > grid_->procLoc()[1]) //!< Is the process in a process row that the detector region covers?
        //     {
        //         fieldInFreq.loc_[1] = 1;
        //         if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
        //         {
        //             fieldInFreq.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
        //         }
        //         else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
        //         {
        //             fieldInFreq.loc_[2] = 1;
        //         }
        //     }
        // }
        // else if(loc_[0] < grid_->procLoc()[0] && loc_[0] + sz_[0] > grid_->procLoc()[0]) //!< Is the process in a process col that the detector covers?
        // {
        //     fieldInFreq.loc_[0] = 1;
        //     if(loc_[1] >= grid_->procLoc()[1] && loc_[1] < grid_->procLoc()[1] + grid_->local_y() -2)
        //     {
        //         fieldInFreq.loc_[1] = loc_[1] - grid_->procLoc()[1] + 1;
        //         if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
        //         {
        //             fieldInFreq.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
        //         }
        //         else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
        //         {
        //             fieldInFreq.loc_[2] = 1;
        //         }
        //     }
        //     else if(loc_[1] < grid_->procLoc()[1] && loc_[1] + sz_[1] > grid_->procLoc()[1])
        //     {
        //         fieldInFreq.loc_[1] = 1;
        //         if(loc_[2] >= grid_->procLoc()[2] && loc_[2] < grid_->procLoc()[2] + grid_->local_z() -2)
        //         {
        //             fieldInFreq.loc_[2] = loc_[2] - grid_->procLoc()[2] + 1;
        //         }
        //         else if(loc_[2] < grid_->procLoc()[2] && loc_[2] + sz_[2] > grid_->procLoc()[2])
        //         {
        //             fieldInFreq.loc_[2] = 1;
        //         }
        //     }
        // }
        if(grid_->local_z() == 1)
        {
            fieldInFreq.loc_[2] = 0;
        }
        // If either fieldInFreq.loc_ values is 0 then the detector is not in the process
        if(fieldInFreq.loc_[0] != -2 && fieldInFreq.loc_[1] != -2 && fieldInFreq.loc_[2] != -2)
        {
            if( propDir == DIRECTION::X || (propDir == DIRECTION::NONE  && isamax_(sz_.size(), sz_.data(), 1) - 1 == 1) )
            {
                fieldInFreq.stride_ = grid_->local_x();
                // 2D use y to copy else use z
                if( grid_->local_z() == 1 )
                {
                    shift_j_ = 0;
                    shift_k_ = 1;
                    if( (sz_[1] + loc_[1] > grid_->procLoc()[1] + grid_->local_y() - 2)) // Does the detector go through the end of the process' grid?
                        fieldInFreq.sz_[0] = grid_->local_y() - fieldInFreq.loc_[1] - 1;
                    else
                        fieldInFreq.sz_[0] = loc_[1] + sz_[1] - (grid_->procLoc()[1] + fieldInFreq.loc_[1] - 1);

                    fieldInFreq.addVec1_ = {0, 0, 1};
                    fieldInFreq.addVec2_ = {1, 0, 0};

                    if(sz_[0] + loc_[0] > grid_->procLoc()[0] + grid_->local_x() - 2)
                        fieldInFreq.sz_[2] = grid_->local_x() - fieldInFreq.loc_[0] - 1;
                    else
                        fieldInFreq.sz_[2] = loc_[0] + sz_[0] - (grid_->procLoc()[0] + fieldInFreq.loc_[0] - 1);

                    if(grid_->local_z() == 1)
                         fieldInFreq.sz_[1] = 1;
                    else if(sz_[2] + loc_[2] > grid_->procLoc()[2] + grid_->local_z() - 2)
                    {
                        fieldInFreq.sz_[1] = grid_->local_z() - fieldInFreq.loc_[2] - 1;
                    }
                    else
                    {
                        fieldInFreq.sz_[1] = loc_[2] + sz_[2] - (grid_->procLoc()[2] + fieldInFreq.loc_[2] - 1);
                    }
                }
                else
                {
                    shift_j_ = 1;
                    shift_k_ = 0;
                    if( (sz_[1] + loc_[1] > grid_->procLoc()[1] + grid_->local_y() - 2)) // Does the detector go through the end of the process' grid?
                        fieldInFreq.sz_[1] = grid_->local_y() - fieldInFreq.loc_[1] - 1;
                    else
                        fieldInFreq.sz_[1] = loc_[1] + sz_[1] - (grid_->procLoc()[1] + fieldInFreq.loc_[1] - 1);

                    fieldInFreq.addVec1_ = {0, 1, 0};
                    fieldInFreq.addVec2_ = {1, 0, 0};

                    if(sz_[0] + loc_[0] > grid_->procLoc()[0] + grid_->local_x() - 2)
                        fieldInFreq.sz_[2] = grid_->local_x() - fieldInFreq.loc_[0] - 1;
                    else
                        fieldInFreq.sz_[2] = loc_[0] + sz_[0] - (grid_->procLoc()[0] + fieldInFreq.loc_[0] - 1);

                    if(sz_[2] + loc_[2] > grid_->procLoc()[2] + grid_->local_z() - 2)
                    {
                        fieldInFreq.sz_[0] = grid_->local_z() - fieldInFreq.loc_[2] - 1;
                    }
                    else
                    {
                        fieldInFreq.sz_[0] = loc_[2] + sz_[2] - (grid_->procLoc()[2] + fieldInFreq.loc_[2] - 1);
                    }
                }
            }
            else if( propDir == DIRECTION::Y || (propDir == DIRECTION::NONE  && isamax_(sz_.size(), sz_.data(), 1) - 1 == 0) )
            {

                shift_j_ = 1;
                shift_k_ = 0;
                fieldInFreq.stride_ = 1;
                if(sz_[0] + loc_[0] > grid_->procLoc()[0] + grid_->local_x() - 2)
                    fieldInFreq.sz_[0] = grid_->local_x() - fieldInFreq.loc_[0] - 1;
                else
                    fieldInFreq.sz_[0] = loc_[0] + sz_[0] - (grid_->procLoc()[0] + fieldInFreq.loc_[0] - 1);

                fieldInFreq.addVec1_ = {0, 0, 1};
                fieldInFreq.addVec2_ = {0, 1, 0};

                if( (sz_[1] + loc_[1] > grid_->procLoc()[1] + grid_->local_y() - 2)) // Does the detector go through the end of the process' grid?
                    fieldInFreq.sz_[2] = grid_->local_y() - fieldInFreq.loc_[1] - 1;
                else
                    fieldInFreq.sz_[2] = loc_[1] + sz_[1] - (grid_->procLoc()[1] + fieldInFreq.loc_[1] - 1);

                if(grid_->local_z() == 1)
                {
                    fieldInFreq.sz_[1] = 1;
                }
                else if(sz_[2] + loc_[2] > grid_->procLoc()[2] + grid_->local_z() - 2)
                {
                    fieldInFreq.sz_[1] = grid_->local_z() - fieldInFreq.loc_[2] - 1;
                }
                else
                {
                    fieldInFreq.sz_[1] = loc_[2] + sz_[2] - (grid_->procLoc()[2] + fieldInFreq.loc_[2] - 1);
                }
            }
            else
            {
                shift_j_ = 0;
                shift_k_ = 1;
                fieldInFreq.stride_ = 1;
                if(sz_[0] + loc_[0] > grid_->procLoc()[0] + grid_->local_x() - 2)
                    fieldInFreq.sz_[0] = grid_->local_x() - fieldInFreq.loc_[0] - 1;
                else
                    fieldInFreq.sz_[0] = loc_[0] + sz_[0] - (grid_->procLoc()[0] + fieldInFreq.loc_[0] - 1);

                fieldInFreq.addVec1_ = {0, 1, 0};
                fieldInFreq.addVec2_ = {0, 0, 1};

                if( (sz_[1] + loc_[1] > grid_->procLoc()[1] + grid_->local_y() - 2)) // Does the detector go through the end of the process' grid?
                    fieldInFreq.sz_[1] = grid_->local_y() - fieldInFreq.loc_[1] - 1;
                else
                    fieldInFreq.sz_[1] = loc_[1] + sz_[1] - (grid_->procLoc()[1] + fieldInFreq.loc_[1] - 1);
                if(sz_[2] + loc_[2] > grid_->procLoc()[2] + grid_->local_z() - 2)
                    fieldInFreq.sz_[2] = grid_->local_z() - fieldInFreq.loc_[2] - 1;
                else
                    fieldInFreq.sz_[2] = loc_[2] + sz_[2] - (grid_->procLoc()[2] + fieldInFreq.loc_[2] - 1);
                // gridComm_.barrier();
                // if(gridComm_.rank() == 0)
                //     std::cout << gridComm_.rank() << '\t' << fieldInFreq.loc_[0] << '\t' << fieldInFreq.loc_[1] << '\t' << fieldInFreq.loc_[2] << '\t' << fieldInFreq.sz_[0] << '\t' << fieldInFreq.sz_[1] << '\t' << fieldInFreq.sz_[2] << std::endl;
                // gridComm_.barrier();
                // if(gridComm_.rank() == 1)
                //     std::cout << gridComm_.rank() << '\t' << fieldInFreq.loc_[0] << '\t' << fieldInFreq.loc_[1] << '\t' << fieldInFreq.loc_[2] << '\t' << fieldInFreq.sz_[0] << '\t' << fieldInFreq.sz_[1] << '\t' << fieldInFreq.sz_[2] << std::endl;
                // gridComm_.barrier();
            }
            fieldInFreq_ = std::make_shared<fInParam>(fieldInFreq); //!< only make slave active if necessary
        }
        return;
    }

    /**
     * @brief      take in fields
     *
     * @param      fftFact  The fourier transform factors
     */
    virtual void fieldIn(cplx* fftFact) = 0;

    inline std::shared_ptr<std::vector<slaveProcInfo>> master() { return master_; }
};


class parallelStorageFreqDTCReal : public parallelStorageFreqDTC<double>
{
public:
    /**
     * @brief      construct the frequency storage dtc
     *
     * @param[in]  dtcNum    the detector number
     * @param[in]  grid      pointer to output grid
     * @param[in]  propDir   The direction of propagation
     * @param[in]  loc       location of lower left corner of the dtc
     * @param[in]  sz        size in grid points for the dtc
     * @param[in]  freqList  The frequency list
     */
    parallelStorageFreqDTCReal(int dtcNum, real_pgrid_ptr grid, DIRECTION propDir, std::array<int,3> loc, std::array<int,3> sz, std::vector<double> freqList);
    /**
     * @brief      take in fields
     *
     * @param      fftFact  The fourier transform factors
     */
    void fieldIn(cplx* fftFact);
};

/**
 * Complex version of FluxDTC see base class for more descriptions
 */
class parallelStorageFreqDTCCplx : public parallelStorageFreqDTC<cplx>
{
public:
    /**
     * @brief      construct the frequency storage dtc
     *
     * @param[in]  dtcNum    the detector number
     * @param[in]  grid      pointer to output grid
     * @param[in]  propDir   The direction of propagation
     * @param[in]  loc       location of lower left corner of the dtc
     * @param[in]  sz        size in grid points for the dtc
     * @param[in]  freqList  The frequency list
     */
    parallelStorageFreqDTCCplx(int dtcNum, cplx_pgrid_ptr grid, DIRECTION propDir, std::array<int,3> loc, std::array<int,3> sz, std::vector<double> freqList);
    /**
     * @brief      take in fields
     *
     * @param      fftFact  The fourier transform factors
     */
    void fieldIn(cplx* fftFact);
};

#endif