#ifndef FDTD_pARALLELDETECTOR_COUT
#define FDTD_pARALLELDETECTOR_COUT

#include "parallelDTC.hpp"

class parallelDetectorCOUTReal: public parallelDetectorBaseReal
{
protected:
    using parallelDetectorBaseReal::fields_;
    using parallelDetectorBaseReal::loc_;
    using parallelDetectorBaseReal::sz_;
    using parallelDetectorBaseReal::firstComp_;
    using parallelDetectorBaseReal::timeInterval_;
    using parallelDetectorBaseReal::dirConv_;
    using parallelDetectorBaseReal::tConv_;

    using parallelDetectorBaseReal::convFactor_;
    using parallelDetectorBaseReal::outputFunction_;

    std::vector<int> dtc2cart_; //!< converts the ijk of the detector to the x, y , z cart coordinates
public:
    /**
     * @brief      Constructs a detector that outputs to console
     *
     * @param[in]  dtcNum        The dtc number
     * @param[in]  grid          pointer to output grid
     * @param[in]  SI            bool to determine if SI units are used
     * @param[in]  loc           location of lower left corner of the dtc
     * @param[in]  sz            size in grid points for the dtc
     * @param[in]  out_name      output file name
     * @param[in]  type          The type: output type of dtc: fields or power
     * @param[in]  classType     The class type: Is it output a base field, polarization or power
     * @param[in]  timeInterval  Time interval for how often to record data
     * @param[in]  firstComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length
     * @param[in]  I0            unit current
     * @param[in]  dt            time step
     */
    parallelDetectorCOUTReal(int dtcNum,  std::vector<real_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt);
    /**
     * @brief Outputs the information to console
     *
     * @param[in] t Time of the output.
     */
    void output(double t);
};

class parallelDetectorCOUTCplx: public parallelDetectorBaseCplx
{
protected:
    using parallelDetectorBaseCplx::fields_;
    using parallelDetectorBaseCplx::loc_;
    using parallelDetectorBaseCplx::sz_;
    using parallelDetectorBaseCplx::firstComp_;
    using parallelDetectorBaseCplx::timeInterval_;
    using parallelDetectorBaseCplx::dirConv_;
    using parallelDetectorBaseCplx::tConv_;
    using parallelDetectorBaseCplx::convFactor_;
    using parallelDetectorBaseCplx::outputFunction_;

    std::vector<int> dtc2cart_; //!< converts the ijk of the detector to the x, y , z cart coordinates
public:
    /**
     * @brief      Constructs a detector that outputs to console
     *
     * @param[in]  dtcNum        The dtc number
     * @param[in]  grid          pointer to output grid
     * @param[in]  SI            bool to determine if SI units are used
     * @param[in]  loc           location of lower left corner of the dtc
     * @param[in]  sz            size in grid points for the dtc
     * @param[in]  out_name      output file name
     * @param[in]  type          The type: output type of dtc: fields or power
     * @param[in]  classType     The class type: Is it output a base field, polarization or power
     * @param[in]  timeInterval  Time interval for how often to record data
     * @param[in]  firstComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length
     * @param[in]  I0            unit current
     * @param[in]  dt            time step
     */
    parallelDetectorCOUTCplx(int dtcNum,  std::vector<cplx_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt);
    /**
     * @brief Outputs the information to console
     *
     * @param[in] t Time of the output.
     */
    void output(double t);
};
#endif