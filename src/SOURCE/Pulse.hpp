#ifndef FDTD_PULSE
#define FDTD_PULSE

#include <UTIL/typedefs.hpp>

class PulseBase
{
protected:
    std::vector<double> param_; //!< pulse function parameters
    double E0_; //!< peak intensity of the pulse
    cplx phaseFact_; //!< complex value to offset the phase of the pulse
    double dt_; //!< time step

public:
    /**
     * @brief Constructor
     *
     * @param param functional parameters of the pulse
     * @param E0_ peak intensity of the pulse
     * @param dt of the time step
     * @param alpha phase of the pulse
     */
    PulseBase(std::vector<double> param, double E0, double dt, double alpha = 0.0) :
        param_(param),
        E0_(E0),
        phaseFact_( exp( cplx(0.0, alpha) ) ),
        dt_(dt)
    {};
    //Acessor Functions
    /**
     * @return the parameters of the pulse
     */
    std::vector<double> param() {return param_;}

    /**
     * @brief Creates the pulse
     * @details Returns the value of the pulse at a given time
     *
     * @param t current time
     * @return the pulse value
     */
    virtual const cplx pulse(double t) = 0;

    /**
     * @brief      modifies the E0 of a pulse by a factor fact
     *
     * @param[in]  fact  The factor to modify the E0 of the pulse
     */
    virtual void modE0(double fact) = 0 ;

};
class PulseGauss : public PulseBase
{
protected:
    using PulseBase::param_;
    using PulseBase::dt_;
    using PulseBase::E0_;
    cplx omg_; //!< center frequency of the pulse
    double t0_; //!< center time of the pulse
    double cutoff_; //!< width * cutoff = when to stop calculating the pulse
public:
    /**
     * @brief Constructs a Gaussian Pulse
     *
     * @param param Mathematical Parameters of the PulseGauss
     * @param E0 Peak intensity of the pulse
     * @param dt time step of the FDTD simulation
     * @param alpha phase of the pulse
     */
    PulseGauss(std::vector<double> param, double E0, double dt, double alpha = 0.0);
    /**
     * @brief Calculates the value of the pulse at time t
     *
     * @param t current time of the simulation
     * @return pulse value
     */
    const cplx pulse(double t);
        /**
     * @brief      modifies the E0 of a pulse by a factor fact
     *
     * @param[in]  fact  The factor to modify the E0 of the pulse
     */
    inline void modE0(double fact){E0_ *= fact;}

};
class PulseCont : public PulseBase
{
protected:
    using PulseBase::param_;
    using PulseBase::dt_;
    using PulseBase::E0_;
    cplx omg_; //!< center frequency of the pulse

public:
    /**
     * @brief Constructs a Continuous Pulse
     *
     * @param param Mathematical Parameters of the PulseCont
     * @param E0 Peak intensity of the pulse
     * @param dt time step of the FDTD simulation
     * @param alpha phase of the pulse
     */
    PulseCont(std::vector<double> param, double E0, double dt, double alpha = 0.0) ;
    /**
     * @brief Calculates the value of the pulse at time t
     *
     * @param t current time of the simulation
     * @return pulse value
     */
    const cplx pulse(double t);
    /**
     * @brief      modifies the E0 of a pulse by a factor fact
     *
     * @param[in]  fact  The factor to modify the E0 of the pulse
     */
    inline void modE0(double fact){E0_ *= fact;}
};

class PulseRicker : public PulseBase
{
protected:
    using PulseBase::param_;
    using PulseBase::dt_;
    using PulseBase::E0_;
public:
    /**
     * @brief Constructs a Ricker Pulse
     * @details Constructs a PulseRicker with the base pulse constructor
     *
     * @param param Mathematical Parameters of the PulseRicker
     * @param E0 Peak intensity of the pulse
     * @param dt time step of the FDTD simulation
     * @param alpha phase of the pulse
     */
    PulseRicker(std::vector<double> param, double E0, double dt, double alpha = 0.0);
    /**
     * @brief Calculates the value of the pulse at time t
     *
     * @param t current time of the simulation
     * @return pulse value
     */
    const cplx pulse(double t);
        /**
     * @brief      modifies the E0 of a pulse by a factor fact
     *
     * @param[in]  fact  The factor to modify the E0 of the pulse
     */
    inline void modE0(double fact){E0_ *= fact;}
};

class PulseRampCont : public PulseBase
{
protected:
    using PulseBase::param_;
    using PulseBase::dt_;
    using PulseBase::E0_;
    cplx omg_; //!< center frequency of the pulse
public:
    /**
     * @brief Constructs a Ramped Continuous pulse
     *
     * @param param Mathematical Parameters of the PulseRampCont
     * @param E0 Peak intensity of the pulse
     * @param dt time step of the FDTD simulation
     */
    PulseRampCont(std::vector<double> param, double E0, double dt, double alpha = 0.0);
    /**
     * @brief Calculates the value of the pulse at time t
     *
     * @param t current time of the simulation
     * @return pulse value
     * @param alpha phase of the pulse
     */
    const cplx pulse(double t);
        /**
     * @brief      modifies the E0 of a pulse by a factor fact
     *
     * @param[in]  fact  The factor to modify the E0 of the pulse
     */
    inline void modE0(double fact){E0_ *= fact;  if(fact < 0.0) param_[1] *= -1.0;}
};

class PulseBH : public PulseBase
{
protected:
    using PulseBase::param_;
    using PulseBase::dt_;
    using PulseBase::E0_;
    cplx omg_; //!< central frequency of the pulse
    double t0_; //!< time of the center of the pulse
    double startPulse_; //!< point to start calculating the pulse t0-width/2.0
    double endPulse_; //!< point to stop calculating the pulse t0+width/2.0
public:
    /**
     * @brief Constructs a Blackman-Harris Pulse
     *
     * @param param Mathematical Parameters of the
     * @param E0 Peak intensity of the pulse
     * @param dt time step of the FDTD simulation
     * @param alpha phase of the pulse
     */
    PulseBH(std::vector<double> param, double E0, double dt, double alpha = 0.0);
    /**
     * @brief Calculates the value of the pulse at time t
     *
     * @param t current time of the simulation
     * @return pulse value
     */
    const cplx pulse(double t);
        /**
     * @brief      modifies the E0 of a pulse by a factor fact
     *
     * @param[in]  fact  The factor to modify the E0 of the pulse
     */
    inline void modE0(double fact){E0_ *= fact;}
};

class PulseRect : public PulseBase
{
protected:
    using PulseBase::param_;
    using PulseBase::dt_;
    using PulseBase::E0_;
    cplx omg_; //!< central frequency of the pulse
    double tau_; //!< width of the pulse
    double t0_; //!< time of the peak of the pulse
    double n_; //!< value controlling shape of the rectangle approx

public:
    /**
     * @brief Constructs a Rectangular pulse
     *
     * @param param Mathematical Parameters of the PulseRect
     * @param E0 Peak intensity of the pulse
     * @param dt time step of the FDTD simulation
     * @param alpha phase of the pulse
     */
    PulseRect(std::vector<double> param, double E0, double dt, double alpha = 0.0);

    /**
     * @brief Calculates the value of the pulse at time t
     *
     * @param t current time of the simulation
     * @return pulse value
     */
    const cplx pulse(double t);
    /**
     * @brief      modifies the E0 of a pulse by a factor fact
     *
     * @param[in]  fact  The factor to modify the E0 of the pulse
     */
    inline void modE0(double fact){E0_ *= fact;}
};

#endif