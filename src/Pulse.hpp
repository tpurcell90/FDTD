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
    Pulse(const Pulse& o) :  Type(o.Type)
    {
        for (int ii = 0 ; ii < o.Param.size(); ii++)
            Param.push_back(o.Param[ii]);
    }
    //Acessor Functions
    std::vector<T> param() {return Param;}
    ProfType type() {return Type;}

    const T pulse(double t)
    {
        T pul;
        switch ( Type )
        {
            case gaussian: //if(srcArr[kk].pol() == EZ)
                pul = gauss_pulse(t);
                break;
            case continuous: //else if(srcArr[kk].pol() == HX)
                pul = const_pulse(t);
                break;
            default:
                pul =  T(0.0);
                break;
        }
        return pul;
        //if(Type == gaussian)
        //    return gauss_pulse(t);
        //else if 
        //return T(0.0);
    }
    //Pulse functions
    const T gauss_pulse(double t)
    {
        std::complex<double> imag(0.0,1.0);
        //if (t < Param[1]*Param[3])
   //return real(-1.0 / (imag*Param[0]) * (-1*Param[0]*imag + (Param[2]-t) / pow(2*Param[1],2)) * exp(-1*Param[0]*imag - pow(((Param[2]-t)/pow(2*Param[1],2.0)),2.0)));
        return t*t*exp(-1 * pow(t-2,2));
        // look for the best way to calculate gaussian pulse
    }
    
    const T const_pulse(double t)
    {
        return T(0.0);
    }
};
#endif
