#ifndef FDTD_pARALLELDETECTOR_BIN
#define FDTD_pARALLELDETECTOR_BIN

#include <src/DTC/parallelDTC.hpp>

class parallelDetectorBINReal : public parallelDetectorBaseReal
{
protected:
    using parallelDetectorBaseReal::fields_;
    using parallelDetectorBaseReal::sz_;
    using parallelDetectorBaseReal::loc_;
    using parallelDetectorBaseReal::firstComp_;
    using parallelDetectorBaseReal::timeInterval_;
    using parallelDetectorBaseReal::dirConv_;
    using parallelDetectorBaseReal::tConv_;
    using parallelDetectorBaseReal::convFactor_;
    using parallelDetectorBaseReal::outputFunction_;

    std::string outFile_; //!< output file name
    std::vector<double> gridVals_; //!< temporary storage of grid values before transfer into a file
public:
    /**
     * @brief      Constructs a detector that outputs to a binary file
     *
     * @param[in]  dtcNum        The dtc number
     * @param[in]  grid          The grid pointer
     * @param[in]  SI            True if using SI units
     * @param[in]  loc           The location in grid points
     * @param[in]  sz            The size in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  type          The type: output type of dtc: fields or power
     * @param[in]  classType     The class type: Is it output a base field, polarization or power
     * @param[in]  timeInterval  The time interval
     * @param[in]  firstComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length
     * @param[in]  I0            unit current
     * @param[in]  dt            time step
     */
    parallelDetectorBINReal(int dtcNum,  std::vector<real_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt) ;
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
    using parallelDetectorBaseCplx::fields_;
    using parallelDetectorBaseCplx::sz_;
    using parallelDetectorBaseCplx::loc_;
    using parallelDetectorBaseCplx::firstComp_;
    using parallelDetectorBaseCplx::timeInterval_;
    using parallelDetectorBaseCplx::dirConv_;
    using parallelDetectorBaseCplx::tConv_;
    using parallelDetectorBaseCplx::convFactor_;
    using parallelDetectorBaseCplx::outputFunction_;

    std::string outFile_; //!< output file name
    std::vector<cplx> gridVals_; //!< temporary storage of grid values before transfer into a file
public:
    /**
     * @brief      Constructs a detector that outputs to a binary file
     *
     * @param[in]  dtcNum        The dtc number
     * @param[in]  grid          The grid pointer
     * @param[in]  SI            True if using SI units
     * @param[in]  loc           The location in grid points
     * @param[in]  sz            The size in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  type          The type: output type of dtc: fields or power
     * @param[in]  classType     The class type: Is it output a base field, polarization or power
     * @param[in]  timeInterval  The time interval
     * @param[in]  firstComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length
     * @param[in]  I0            unit current
     * @param[in]  dt            time step
     */
    parallelDetectorBINCplx(int dtcNum,  std::vector<cplx_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt) ;

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
