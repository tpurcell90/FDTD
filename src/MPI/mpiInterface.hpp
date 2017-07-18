#ifndef SRC_MPIINTERFACE_HPP
#define SRC_MPIINTERFACE_HPP

#include <boost/mpi.hpp>
#include <GRID/Grid.hpp>
#include "parallelUtilities.hpp"
#ifdef MKL
#include <UTIL/utilities_MKL.hpp>
#elif defined ACML
#include <UTIL/utilities_acml.hpp>
#else
#include <UTIL/utilities_MKL.hpp>
#endif

// const static int blocksize__ = 1 ;
typedef std::complex<double> cplx;
typedef std::shared_ptr<Grid<double>> real_grid_ptr;

extern boost::mpi::environment env;
namespace mpi = boost::mpi;

static std::tuple<int,int,int> numgrid(int);

// Augment the boost mpi communicator class with some other useful data
/**
 * @brief Augment the boost MPI communicator class with some other useful data
 * @details MPI communicator used to transfer data throughout the cell.
 *
 */
class mpiInterface : public boost::mpi::communicator
{
protected:
    int npX_; //!< number of process in the X direction
    int npY_; //!< number of process in the Y direction
    int npZ_; //!< number of process in the Z direction
    int mypX_; //!< the local process' process X
    int mypY_; //!< the local Process' process Y
    int mypZ_; //!< the local Process' process Z

public:
    /**
     * @brief      Constructor for the mpiInterface
     */
    mpiInterface();

    /**
     * @brief Accessors to npX_
     * @return npX_
     */
    inline int npX() {return npX_;}

    /**
     * @brief Accessors to npY_
     * @return npY_
     */
    inline int npY() {return npY_;}

    /**
     * @brief Accessors to npZ_
     * @return npZ_
     */
    inline int npZ() {return npZ_;}

    /**
     * @brief Accessors to mypX_
     * @return mypX_
     */
    inline int mypX() {return mypX_;}

    /**
     * @brief Accessors to mypY_
     * @return mypY_
     */
    inline int mypY() {return mypY_;}
    /**
     * @brief Accessors to mypZ_
     * @return mypZ_
     */
    inline int mypZ() {return mypZ_;}
    /**
     * @brief numroc caller
     * @details runs numroc with some of the local parameters
     *
     * @param nx number of grid points in the x direction
     * @param ny number of grid points in the y direction
     * @param nz number of grid points in the z direction
     *
     * @return numroc
     */
    std::tuple<int,int, int> numroc(const int nx, const int ny, const int nz) const;


    /**
     * @brief      Finds the lower, left, back corner of the processor (proc loc for Grids)
     *
     * @param[in]  weight maps of the grids
     *
     * @return     The procLoc for the Grids
     */
    std::tuple<int,int, int> getLocxLocyLocz(std::vector<real_grid_ptr> weights) const;

    /**
     * @brief Unique int tag generator
     * @details uses a modified cantor function to get a unique int tag for each communication, based on process number and what type it is
     *
     * @param procSend sending process
     * @param procRecv receiving process
     * @param maxOffset number of different communication processes possible between two processes within the same operation
     * @param offset the assigned offset corresponding to a single communication within the same operation
     * @return [description]
     */
    int cantorTagGen(unsigned int procSend, unsigned int procRecv, unsigned int maxOffest, unsigned int offest) { return (int((procSend + procRecv) * (procSend + procSend +1) / 2) + procRecv) * maxOffest + offest; }


};
/**
 * @brief      numgrid fuction from scalapack
 *
 * @param[in]  numproc  The ramnk of the process
 *
 * @return     the procLoc
 */
static std::tuple<int, int, int> numgrid(int numproc);

extern mpiInterface world;

#endif