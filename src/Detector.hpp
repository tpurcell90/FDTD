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

template <typename T> class Detector
{
protected:
    std::vector<double> location;
    std::vector<double> size;
    std::string Type;
    T out;
public:
    // Constructor
    Detector(std::vector<double> loc, std::vector<double> sz, std::string type) : location(loc), size(sz), Type(type), out(T(0.0)){}
    // Copy Constructor
    Detector(const Detector& o) : location(o.location), size(o.size), Type(o.Type), out(o.out){}

    // Accessor Function
    std::vector<double> loc() {return location;}
    std::vector<double> sz() {return size;}
    std::string type() {return Type;}

    // Output functions
    T& outField (Grid2D<T> field)
    {

    }
};
#endif
