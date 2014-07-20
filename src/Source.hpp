#ifndef FDTD_SOURCE
#define FDTD_SOURCE

#include "enum.hpp"
#include "Pulse.hpp"
#include <string>
#include <complex>
#include <vector>
#include <functional>


enum Polarization {EX,EY,EZ,HX,HY,HZ};

// Needs modification for complex I think
template <typename T> class Source
{
protected:
    Pulse<T> profile;
    Polarization polarization;
    std::vector<int> location;

public:
    // Constructor
    Source(Pulse<double> prof, Polarization pol,std::vector<int> loc) : profile(prof), polarization(pol), location(loc) {}

    // Copy Constructor
    Source(const Source& o) : profile(o.profile), polarization(o.polarization), location(o.location) {}
   
    //Access Functions
    Polarization pol() const {return polarization;}
    std::vector<double> loc() const {return location;}
    Pulse<T> prof() const {return profile;}
};

#endif
