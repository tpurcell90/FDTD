#ifndef FDTD_DETECTOR
#define FDTD_DETECTOR

#include <string>
#include <vector>
#include <complex>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include "Grid.hpp"

enum OupuptsData  {field, flux};

template <typename T> class Detector
{
protected:
    std::vector<int> location;
    OupuptsData Type;
public:
    // Constructor
    Detector(std::vector<int> loc, OupuptsData type) : location(loc), Type(type){}
    // Copy Constructor
    Detector(const Detector& o) : location(o.location), Type(o.Type){}

    // Accessor Function
    std::vector<int> loc() {return location;}
    std::string type() {return Type;}

    // Output functions
    T output (std::shared_ptr<Grid2D<T>> field_in)
    {
        T out;
        (Type == field) ? out = field_in->point(location[0],location[1]) : out = T(0);
        return out;
    }
};
#endif
