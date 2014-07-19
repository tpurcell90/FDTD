#ifndef FDTD_SOURCE
#define FDTD_SOURCE

#include "ProfileFunctions.hpp"
#include <string>
#include <complex>
#include <vector>
#include <functional>

enum Polarization {EX,EY,EZ,HX,HY,HZ};

// Needs modification for complex I think
template <typename T> class Source
{
protected:
    ProfileFunctions<double> profile;
    Polarization polarization;
    std::vector<double> location;
    std::vector<double> size;
    std::vector<double> direction;
    std::vector<double> fxnParam;

public:
    // Constructor
    Source(ProfileFunctions<double> prof, Polarization pol,std::vector<double> loc, std::vector<double> sz,std::vector<double> dir,std::vector<double> fxn) : profile(prof), polarization(pol), location(loc), size(sz), direction(dir), fxnParam(fxn) {}

    // Copy Constructor
    Source(const Source& o) : profile(o.profile), polarization(o.polarization), location(o.location), size(o.size), direction(o.direction), fxnParam(o.fxnParam) {}

    //Access Functions
    Polarization pol() const {return polarization;}
    T field(double t) const {return profile(t);}
};

#endif
