#ifndef FDTD_PULSE
#define FDTD_PULSE

#include <string>
#include <vector>
#include <complex>

template <typename T> class Pulse
{
protected:
    std::vector<T> Param;
    std::string Type;

public:
    // Constructor
    Pulse(std::vector<T> param, std::string type) : Param(param), Type(type) {}
    // Copy Constructor
    Pulse(const Pulse& o) : Param(o.Param), Type(o.Type) {}
    //Acessor Functions
    std::vector<T> param() {return Param;}
    std::string type() {return Type;}

    //Pulse functions
    T& gauss_pulse(double t)
    {
        // look for the best way to calculate gaussian pulse
    }    
};
#endif
