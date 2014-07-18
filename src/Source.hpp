#ifndef FDTD_SOURCE
#define FDTD_SOURCE

#include <string>
#include <complex>
#include <vector>

enum Polarization {EX,EY,EZ,HX,HY,HZ};

// Needs modification for complex I think
template <typename T> class Source
{
protected:
    std::function<T(double)> profile;
    Polarization polarization;
    std::vector<double> location;
    std::vector<double> size;
    std::vector<double> direction;

public:
    // Constructor
    Source(const std::function<T(T)> prof, Polarization pol,std::vector<double> loc, std::vector<double> sz,std::vector<double> dir) : profile(prof), polarization(pol), location(loc), size(sz), direction(dir) {}

    // Copy Constructor
    Source(const Source& o) : profile(o.profile), polarization(o.polarization), location(o.location), size(o.size), direction(o.direction) {}

    //Access Functions
    Polarization pol() const {return polarization;}
    T field(double t) const {return profile(t);}
};

#endif
