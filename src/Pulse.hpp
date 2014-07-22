#ifndef FDTD_PULSE
#define FDTD_PULSE

#include <string>
#include <vector>
#include <complex>
#include "enum.hpp"

enum ProfType {gaussian, continuous};

template <typename T> class Pulse
{
protected:
    std::vector<T> Param;
    ProfType Type;

public:
    // Constructor
    Pulse(std::vector<T> param, ProfType type) : Type(type)
    {
        for (int ii = 0 ; ii < param.size(); ii++)
            Param.push_back(param[ii]);
    }
    // Copy Constructor
    Pulse(const Pulse& o) : Param(o.Param), Type(o.Type) {}
    //Acessor Functions
    std::vector<T> param() {return Param;}
    ProfType type() {return Type;}

    const T pulse(double t)
    {
        if(Type == gaussian)
            return gauss_pulse(t);
        return T(0.0);
    }
    //Pulse functions
    const T gauss_pulse(double t)
    {
        std::complex<double> imag(0.0,1.0);
        //if (t < Param[1]*Param[3])
        //return real(-1.0 / (imag*Param[0]) * (-1*Param[0]*imag + (Param[2]-t) / pow(2*Param[1],2)) * exp(-1*Param[0]*imag - pow(((Param[2]-t)/pow(2*Param[1],2.0)),2.0)));
        if (t < 3.00)
            return t*exp(-t*t);
        return T(0.0);
        // look for the best way to calculate gaussian pulse
    }
};
#endif
