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
    Source(Pulse<double> prof, Polarization pol,std::vector<int> loc) : profile(prof), polarization(pol)
    {
        for (int ii = 0 ; ii < loc.size(); ii++)
            location.push_back(loc[ii]);

    }
    // Copy Constructor
    Source(const Source& o) : profile(o.profile), polarization(o.polarization)
    {
        for (int ii = 0 ; ii < o.location.size(); ii++)
            location.push_back(o.location[ii]);
    }
    //Access Functions
    Polarization pol() {return polarization;}
    std::vector<int> loc()  {return location;}
    Pulse<T> prof() {return profile;}
};

#endif
