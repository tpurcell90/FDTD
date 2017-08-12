#ifndef FDTD_pARALLELDETECTOR_TXT
#define FDTD_pARALLELDETECTOR_TXT

#include "parallelDTC.hpp"
/**
 * @brief parallelDetector class for printing out field information to a text file
 * @details Copies the field magnitude to a text file
 * @param fieldConv unit conversion for field magnitude
 * @param outFile_
 * @param outFileStream_
 *
 */
class parallelDetectorTXTReal : public parallelDetectorBaseReal
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

    int t_percision_; //!< precison used to store time (easier reading of files)
    std::string outFile_; //!< output file name
    std::shared_ptr<std::ofstream> outFileStream_; //!< the output file stream

public:
    /**
     * @brief      Constructs a detector that outputs to a text file
     *
     * @param[in]  grid          a vector of pointers to output grids
     * @param[in]  SI            bool to determine if SI units are used
     * @param[in]  loc           location of lower left corner of the dtc
     * @param[in]  sz            size in grid points for the dtc
     * @param[in]  out_name      output file name
     * @param[in]  type          The type: output type of dtc: fields or power
     * @param[in]  timeInterval  Time interval for how often to record data
     * @param[in]  a             unit length
     * @param[in]  I0            unit current
     * @param[in]  dt            time step
     */
    parallelDetectorTXTReal(std::vector<real_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, double timeInterval, double a, double I0, double dt);
    /**
     * @brief Output the fields to a text file
     *
     * @param[in] t Time of the output.
     */
    void output (double t);
    /**
     * @brief returns the output file name
     */
    inline std::string outfile() {return outFile_;}
};

class parallelDetectorTXTCplx : public parallelDetectorBaseCplx
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

    int t_percision_; //!< precison used to store time (easier reading of files)
    std::string outFile_; //!< output file name
    std::shared_ptr<std::ofstream> outFileStream_; //!< the output file stream

public:
    /**
     * @brief      Constructs a detector that outputs to a text file
     *
     * @param[in]  grid          a vector of pointers to output grids
     * @param[in]  SI            bool to determine if SI units are used
     * @param[in]  loc           location of lower left corner of the dtc
     * @param[in]  sz            size in grid points for the dtc
     * @param[in]  out_name      output file name
     * @param[in]  type          The type: output type of dtc: fields or power
     * @param[in]  timeInterval  Time interval for how often to record data
     * @param[in]  a             unit length
     * @param[in]  I0            unit current
     * @param[in]  dt            time step
     */
    parallelDetectorTXTCplx(std::vector<cplx_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, double timeInterval, double a, double I0, double dt );
    /**
     * @brief Output the fields to a text file
     *
     * @param[in] t Time of the output.
     */
    void output (double t);
    /**
     * @brief returns the output file name
     */
    inline std::string outfile() {return outFile_;}
};

#endif