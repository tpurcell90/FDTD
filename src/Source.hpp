#ifndef FDTD_SOURCE
#define FDTD_SOURCE

#include <string>
#include <complex>

enum Polarization {EX,EY,EZ,HX,HY,HZ};

// Needs modification for complex I think
template <typename T> class Source
{
protected:
    std::function<T(double)> profile;
    Polarization polarization;

public:
    // Constructor
    Source(const std::function<T(T)> prof, Polarization pol) : profile(prof), polarization(pol) {}

    // Copy Constructor
    Source(const Source& o) : profile(o.profile), polarization(o.polarization) {}

    //Access Functions
    Polarization pol() const {return polarization;}
    T field(double t) const {return profile(t);}
};

#endif
