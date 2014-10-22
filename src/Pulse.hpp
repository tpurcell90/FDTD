#ifndef FDTD_PULSE
#define FDTD_PULSE

#include <string>
#include <vector>
#include <complex>
#include "enum.hpp"
#include <iostream>

enum plsShape {gaussian, continuous, ricker};

template <typename T> class Pulse
{
protected:
    std::vector<double> param_;
    plsShape plsType_;
    double dt_;

public:
    /**
     * @brief Constructor
     * @details Constructs a Pulse shape given the correct parameters
     *
     * @param param functional parameters of the pulse
     * @param type Shape of the pulse
     */
    Pulse(std::vector<double> param, plsShape type, double dt) : param_(param), plsType_(type), dt_(dt) {}
    /**
     * @brief Copy Constructor
     * @details Constructs a Pulse from another Pulse
     *
     * @param o The pulseshpe to be copied
     */
    Pulse(const Pulse& o) : param_(o.param_), plsType_(o.plsType_), dt_(o.dt_) {}
    //Acessor Functions
    /**
     * @brief Returns the parameters of the pulse
     */
    std::vector<T> param() {return param_;}
    /**
     * @brief Retuns the pulse shape
     * @return The Pulse shape
     */
    plsShape type() {return plsType_;}
    /**
     * @brief Creates the pulse
     * @details Returns the value of the pulse at a given time
     *
     * @param t current time
     * @return the pulse value
     */
    const T pulse(double t)
    {
        T pul;
        switch ( plsType_ )
        {
            case gaussian: //if(srcArr[kk].pol() == EZ)
                pul = gauss_pulse(t);
                break;
            case continuous: //else if(srcArr[kk].pol() == HX)
                pul = const_pulse(t);
                break;
            case ricker:
                pul = ricker_pulse(t);
                break;
            default:
                pul =  T(0.0);
                break;
        }
        return pul;
        //if(plsType_ == gaussian)
        //    return gauss_pulse(t);
        //else if
        //return T(0.0);
    }
    //Pulse functions
    /**
     * @brief Creates a gaussian pulse
     * @details Returns the value of the gaussian pulse at a given time
     *
     * @param t current time
     * @return the gaussian pulse value
     */
    const T gauss_pulse(double n)
    {
        double t = n*dt_;
        std::complex<double> i(0.0,1.0);
        double omg = param_[0] * 2 * M_PI;
        double t0 = param_[1] * param_[2];
        //if (t < param_[1]*param_[3])
   //return real(-1.0 / (imag*param_[0]) * (-1*param_[0]*imag + (param_[2]-t) / pow(2*param_[1],2)) * exp(-1*param_[0]*imag - pow(((param_[2]-t)/pow(2*param_[1],2.0)),2.0)));
        if (t <= 2*param_[1] * param_[2] + dt_)
            return -100.0/(omg*i*dt_) * (exp(-i*omg*t - pow((t-t0)/param_[1],2.0)/2.0) - exp(-i*omg*(t-dt_) - pow(((t-dt_)-t0)/param_[1],2.0)/2.0));
            // return exp(-1.0 * pow((n - param_[1]*param_[2]/2.0)/(sqrt(2.0)*param_[1]),2.0))*exp(i * 2.0*M_PI * param_[0] * (n- param_[1]*param_[2]/2.0));
        else
            return T(0.0);
        // look for the best way to calculate gaussian pulse
    }
    /**
     * @brief Creates a continuous pulse
     * @details Returns the value of the continuous pulse at a given time
     *
     * @param t current time
     * @return the continuous pulse value
     */
    const T const_pulse(double t)
    {
        std::complex<double> i(0.0,1.0);
        if (t *param_[2] < param_[1] / (param_[0]/ param_[2]))
            return sin(2.0*M_PI * param_[0] * t);
        else
            return T(0.0);
    }
    const T ricker_pulse(double t)
    {
        if (t < param_[2] * param_[1]/param_[0])
            return (1.0-2.0*pow(M_PI*(t*param_[0] - param_[1]),2.0))*exp(-pow(M_PI*(t*param_[0] - param_[1]),2.0));
        else
            return T(0.0);
    }
};
#endif