#ifndef FDTD_pARALLELDETECTOR_BIN
#define FDTD_pARALLELDETECTOR_BIN

#include <src/DTC/parallelDTC.hpp>

class parallelDetectorBINReal : public parallelDetectorBaseReal
{
protected:
    using parallelDetectorBaseReal::timeInterval_; //!< The stride of the time (How often should the detector print?)
    using parallelDetectorBaseReal::tConv_; //!< conversion factor for t to get it in the correct units
    using parallelDetectorBaseReal::convFactor_; //!< Conversion factor for the type of output (to SI units from FDTD)
    using parallelDetectorBaseReal::sz_; //!< the location of the detector's lower left corner
    using parallelDetectorBaseReal::loc_; //!< the size in grid points of the detector
    using parallelDetectorBaseReal::realSpaceLoc_; //!< Location of lower, left, back corner in real spaceZ
    using parallelDetectorBaseReal::fields_; //!< A vector of shared pointers to each of the grids associated with the detector
    using parallelDetectorBaseReal::outputFunction_; //!< function to take grids and output to file in the correct manner

    std::string outFile_; //!< output file name
    std::vector<double> gridVals_; //!< temporary storage of grid values before transfer into a file
public:
    /**
     * @brief      Constructs a detector that outputs to a binary file
     *
     * @param[in]  grid          a vector of pointers to output grids
     * @param[in]  SI            True if using SI units
     * @param[in]  loc           The location in grid points
     * @param[in]  sz            The size in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  type          The type: output type of dtc: fields or power
     * @param[in]  timeInterval  The time interval
     * @param[in]  a             unit length
     * @param[in]  I0            unit current
     * @param[in]  dt            time step
     */
    parallelDetectorBINReal(std::vector<real_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, double timeInterval, double a, double I0, double dt) ;
    /**
     * @brief Outputs to a binary file
     *
     * @param[in] t time of the simulation
     */
    void output(double t);
    /**
     * @brief returns the output file name
     */
    inline std::string outfile() {return outFile_;}
};

class parallelDetectorBINCplx : public parallelDetectorBaseCplx
{
protected:
    using parallelDetectorBaseCplx::timeInterval_; //!< The stride of the time (How often should the detector print?)
    using parallelDetectorBaseCplx::tConv_; //!< conversion factor for t to get it in the correct units
    using parallelDetectorBaseCplx::convFactor_; //!< Conversion factor for the type of output (to SI units from FDTD)
    using parallelDetectorBaseCplx::sz_; //!< the location of the detector's lower left corner
    using parallelDetectorBaseCplx::loc_; //!< the size in grid points of the detector
    using parallelDetectorBaseCplx::realSpaceLoc_; //!< Location of lower, left, back corner in real spaceZ
    using parallelDetectorBaseCplx::fields_; //!< A vector of shared pointers to each of the grids associated with the detector
    using parallelDetectorBaseCplx::outputFunction_; //!< function to take grids and output to file in the correct manner

    std::string outFile_; //!< output file name
    std::vector<cplx> gridVals_; //!< temporary storage of grid values before transfer into a file
public:
    /**
     * @brief      Constructs a detector that outputs to a binary file
     *
     * @param[in]  grid          a vector of pointers to output grids
     * @param[in]  SI            True if using SI units
     * @param[in]  loc           The location in grid points
     * @param[in]  sz            The size in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  type          The type: output type of dtc: fields or power
     * @param[in]  timeInterval  The time interval
     * @param[in]  a             unit length
     * @param[in]  I0            unit current
     * @param[in]  dt            time step
     */
    parallelDetectorBINCplx(std::vector<cplx_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, double timeInterval, double a, double I0, double dt) ;

   /**
     * @brief Outputs to a binary file
     *
     * @param[in] t time of the simulation
     */
    void output(double t);

    /**
     * @brief returns the output file name
     */
    inline std::string outfile() {return outFile_;}
};

#endif
