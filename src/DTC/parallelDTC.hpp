#ifndef FDTD_PARALLELDETECTOR
#define FDTD_PARALLELDETECTOR

#include <DTC/parallelDTCOutputFxn.hpp>
#include <DTC/parallelStorageDTC.hpp>

template <typename T> class parallelDetectorBase
{
protected:
    DTCTYPE type_; //!< The type of the detector: EX,EY,EZ,HX,HY,HZ,EPWR,HPWR
    DIRECTION firstComp_; //!< loc_[0] corresponds to this direction
    int timeInterval_; //!< The stride of the time (How often should the detector print?)
    double tConv_; //!< conversion factor for t to get it in the correct units
    double convFactor_; //!< Conversion factor for the type of output (to SI units from FDTD)
    std::array<int,3> loc_; //!< the location of the detector's lower left corner
    std::array<int,3> sz_; //!< the size in grid points of the detector
    std::array<double,3> dirConv_; //!< Conversion factor for directions
    std::vector<std::shared_ptr<parallelStorageDTC<T>>> fields_; //!< A vector of shared pointers to each of the grids associated with the detector
    std::function<void(T* gridInBegin, T* gridInEnd, T* valsOut, double convFactor)> outputFunction_; //!< function to take grids and output to file in the correct manner
public:
    /**
     * @brief      Construct a basic detector
     *
     * @param[in]  dtcNum        The dtc number
     * @param[in]  grids         The fields that need to be output
     * @param[in]  SI            True if using SI units
     * @param[in]  loc           The location in grid point position
     * @param[in]  sz            The size in number of grid points
     * @param[in]  type          The type: output type of dtc: fields or power
     * @param[in]  classType     The class type: Is it output a base field, polarization or power
     * @param[in]  timeInterval  The time interval
     * @param[in]  firstComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length
     * @param[in]  I0            unit current
     * @param[in]  dt            time step
     */
    parallelDetectorBase(int dtcNum, std::vector<std::shared_ptr<parallelGrid<T>>> grids, bool SI, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt)  :
        type_(type),
        firstComp_(firstComp),
        timeInterval_(static_cast<int>(std::floor(timeInterval/dt) ) ),
        dirConv_(grids[0]->d()),
        tConv_(1.0),
        convFactor_(1.0),
        loc_(loc),
        sz_(sz)
    {
        if(timeInterval_ == 0)
            throw std::logic_error("Well it looks like you are trying to record the fields at a time less than the actual time step, if I were you I would fix that. Anyway I'll just give you this error as a reward for your valiant attempts to break time discretization.");
        if(static_cast<double>(timeInterval_) * dt != timeInterval)
        if(SI)
        {
            tConv_ *= a / SPEED_OF_LIGHT;
            dscal_(dirConv_.size(), a, dirConv_.data(), 1);
        }

        if(classType == DTCCLASSTYPE::FIELD )
        {
            if(SI)
            {
                convFactor_ = (I0/a);
                if(type == DTCTYPE::EX || type == DTCTYPE::EY || type == DTCTYPE::EZ)
                    convFactor_ /= EPS0*SPEED_OF_LIGHT;
            }
            outputFunction_ = fieldOutputFunction<T>;
        }
        else if(classType == DTCCLASSTYPE::POW )
        {
            if(SI)
            {
                convFactor_ = pow(I0/a, 2.0);
                if(type == DTCTYPE::EX || type == DTCTYPE::EY || type == DTCTYPE::EZ)
                    convFactor_ /= pow(EPS0*SPEED_OF_LIGHT, 2.0);
            }
            outputFunction_ = pwrOutputFunction<T>;
        }
        else if(classType == DTCCLASSTYPE::POL )
        {
            if(SI)
            {
                convFactor_ = ( I0/(a*SPEED_OF_LIGHT) );
            }
            outputFunction_ = fieldOutputFunction<T>;
        }
        else
            throw std::logic_error("DTCCLASSTYPE is not defined.");
    }
    /**
     * @return timeInterval_
     */
    inline int &timeInt() {return timeInterval_;}

    /**
     * @return location of dtc lower left corner
     */
    inline std::array<int,3> &loc() {return fields_[0]->loc();}

    /**
     * @return size of dtc lower left corner
     */
    inline std::array<int,3> &sz() {return fields_[0]->sz();}
    /**
     * @brief Output the fields
     * @details At time t output the fields
     *
     * @param t Time of the output.
     */
    virtual void output(double t) = 0;
};

class parallelDetectorBaseReal : public parallelDetectorBase<double>
{
public:
    /**
     * @brief      Construct a basic detector for real fields
     *
     * @param[in]  dtcNum        The dtc number
     * @param[in]  grids         The fields that need to be output
     * @param[in]  SI            True if using SI units
     * @param[in]  loc           The location in grid point position
     * @param[in]  sz            The size in number of grid points
     * @param[in]  type          The type: output type of dtc: fields or power
     * @param[in]  classType     The class type: Is it output a base field, polarization or power
     * @param[in]  timeInterval  The time interval
     * @param[in]  firstComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length
     * @param[in]  I0            unit current
     * @param[in]  dt            time step
     */
    parallelDetectorBaseReal(int dtcNum, std::vector<real_pgrid_ptr> grids, bool SI, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt);
    /**
     * @brief      outputs the data to the field
     *
     * @param[in]  t     current simulation time
     */
    virtual void output(double t) = 0;
};

/**
 * Complex version of FluxDTC see base class for more descriptions
 */
class parallelDetectorBaseCplx : public parallelDetectorBase<cplx>
{
public:
    /**
     * @brief      Construct a basic detector for complex fields
     *
     * @param[in]  dtcNum        The dtc number
     * @param[in]  grids         The fields that need to be output
     * @param[in]  SI            True if using SI units
     * @param[in]  loc           The location in grid point position
     * @param[in]  sz            The size in number of grid points
     * @param[in]  type          The type: output type of dtc: fields or power
     * @param[in]  classType     The class type: Is it output a base field, polarization or power
     * @param[in]  timeInterval  The time interval
     * @param[in]  firstComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length
     * @param[in]  I0            unit current
     * @param[in]  dt            time step
     */
    parallelDetectorBaseCplx(int dtcNum, std::vector<cplx_pgrid_ptr> grids, bool SI, std::array<int,3> loc, std::array<int,3> sz, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt);
    /**
     * @brief      outputs the data to the field
     *
     * @param[in]  t     current simulation time
     */
    virtual void output(double t) = 0;
};

#endif
