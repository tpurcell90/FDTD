#include "Pulse.hpp"
#include <iostream>
PulseCont::PulseCont(std::vector<double> param, double E0, double dt, double alpha) :
    PulseBase(param, E0, dt, alpha),
    omg_(0.0, -1.0*param_[0] * 2.0 * M_PI)
{}

const cplx PulseCont::pulse(double t)
{
    return E0_ * exp(omg_ * t) * phaseFact_;
}

PulseRicker::PulseRicker(std::vector<double> param, double E0, double dt, double alpha) :
    PulseBase(param, E0, dt, alpha)
{}

const cplx PulseRicker::pulse(double t)
{
    if (t < param_[2] * param_[1]/param_[0])
        return E0_ * (1.0-2.0*pow(M_PI*(t*param_[0] - param_[1]),2.0))*exp(-pow(M_PI*(t*param_[0] - param_[1]),2.0)) * phaseFact_;
    else
        return cplx(0.0,0.0);
}

PulseRampCont::PulseRampCont(std::vector<double> param, double E0, double dt, double alpha) :
    PulseBase(param, E0, dt, alpha),
    omg_(0.0, -1.0*param_[0] * 2.0 * M_PI)
{
    if(E0 < 0.0)
        param_[1] *= -1.0;
    // if(std::abs(std::imag(phaseFact_)) <1e-15)
    //     phaseFact_ = cplx(std::real(phaseFact_), 0.0);
    // if(std::real(phaseFact_) <1e-15)
    //     phaseFact_ = cplx(0.0, std::imag(phaseFact_) );
}
const cplx PulseRampCont::pulse(double t)
{
    if(std::abs(param_[1]*t) <= std::abs(E0_) )
    {
        return param_[1] * t * exp(omg_ * t) * phaseFact_;
    }
    else
    {
        return E0_ * exp(omg_ * t) * phaseFact_;
    }
}

PulseGauss::PulseGauss(std::vector<double> param, double E0, double dt, double alpha) :
    PulseBase(param, E0, dt, alpha),
    omg_(0.0, -1.0 * param_[0] * 2 * M_PI),
    t0_(param_[3]),
    cutoff_(param_[1] * param_[2] + t0_)
{}

const cplx PulseGauss::pulse(double t)
{
    if (t <= cutoff_ + dt_)
    {
        return E0_ * exp(-1.0*pow((t-t0_)/param_[1],2.0)/2.0) * exp(omg_*t) * phaseFact_;
    }
    else
        return 0.0;
}

PulseBH::PulseBH(std::vector<double> param, double E0, double dt, double alpha) :
    PulseBase(param, E0, dt, alpha),
    omg_(0.0,-1.0*param_[0] * 2 * M_PI),
    t0_(param_[2]),
    startPulse_(t0_ ),
    endPulse_(t0_ +param_[1])
{}


const cplx PulseBH::pulse(double t)
{
    if (startPulse_ <= t && t <=  endPulse_)
        return E0_ * exp( omg_*(t-t0_) ) * ( param_[3] + param_[4] * cos( 2.0*M_PI*(t-t0_) / param_[1] ) + param_[5] * cos(4.0*M_PI*(t-t0_) / param_[1]) + param_[6] * cos(6.0*M_PI*(t-t0_) / param_[1]) ) * phaseFact_;
    else
        return 0.0;
}

PulseRect::PulseRect(std::vector<double> param, double E0, double dt, double alpha) :
    PulseBase(param, E0, dt, alpha),
    omg_(0.0,-1.0*param_[3] * 2 * M_PI),
    tau_(param_[0]/2.0  ),
    t0_(param_[1] + param_[0]/2.0),
    n_(param_[2])
{}

const cplx PulseRect::pulse(double t)
{
    return E0_/(pow((t - t0_)/tau_, n_) + 1) * ( exp(omg_*(t-t0_) ) ) * phaseFact_;
}
