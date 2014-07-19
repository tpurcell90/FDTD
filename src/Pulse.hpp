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
    Pulse(std::vector<T> param, ProfType type) : Param(param), Type(type) {}
    // Copy Constructor
    Pulse(const Pulse& o) : Param(o.Param), Type(o.Type) {}
    //Acessor Functions
    std::vector<T> param() {return Param;}
    ProfType type() {return Type;}

    T& pulse(double t)
    {
        if(Type == gaussian)
            return gauss_pulse(t);
        return T(0.0);
    }
    //Pulse functions
    T& gauss_pulse(double t)
    {
        // look for the best way to calculate gaussian pulse
    }    
};
#endif
