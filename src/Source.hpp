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
    Polarization polarization_;
    std::vector<int> location_;

public:
    // Constructor
    Source(Pulse<double> prof, Polarization pol,std::vector<int> loc) : profile(prof), polarization_(pol), location_(loc) {}
    // Copy Constructor
    Source(const Source& o) : profile(o.profile), polarization_(o.polarization_), location_(o.location_) {}
    //Access Functions
    Polarization pol() {return polarization_;}
    std::vector<int> loc()  {return location_;}
    Pulse<T> prof() {return profile;}
};

#endif
