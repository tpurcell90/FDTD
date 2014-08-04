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
#include "Source.hpp"

enum dtcOutType  {field, flux};

template <typename T> class Detector
{
protected:
    std::vector<int> location_;
    dtcOutType dtcType_;
    std::string outFile_;
    Polarization pol_;
public:
    // Constructor
    Detector(std::vector<int> loc, dtcOutType type, std::string out_name, Polarization pol) : location_(loc),dtcType_(type), outFile_(out_name), pol_(pol) {}
    // Copy Constructor
    Detector(const Detector& o) : location_(o.location_),dtcType_(o.dtcType_), outFile_(o.outFile_), pol_(o.pol_) {}

    // Accessor Function
    std::vector<int> loc() {return location_;}
    dtcOutType type() {return dtcType_;}
    std::string outfile() {return outFile_;}
    Polarization pol() {return pol_;}
    // Output functions
    T output (std::shared_ptr<Grid2D<T>> field_in, double eps)
    {
        T out;
        switch (dtcType_)
        {
            case field:
                out = field_in->point(location_[0],location_[1]);
                break;
            case flux:
                out = field_in->flux(location_,eps);
                break;
            default:
                out =  T(0.0);
                break;
        }
        //(Type == field) ? out = field_in->point(location[0],location[1]) : out = T(0);
        return out;
    }
};
#endif
