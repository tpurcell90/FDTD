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
    /**
     * @brief Constructs a detector
     * @details Constructs a detector for the FDTD field
     *
     * @param loc Grid coordinates for the detector
     * @param type Whether it will record flux of field data
     * @param out_name output file name
     * @param pol polarization of the filed you are detecting
     */
    Detector(std::vector<int> loc, dtcOutType type, std::string out_name, Polarization pol) : location_(loc),dtcType_(type), outFile_(out_name), pol_(pol) {}
    /**
     * @brief Copy Constructor
     * @details Inputs a detector and makes a copy of it
     *
     * @param o The detector to be copied
     */
    Detector(const Detector& o) : location_(o.location_),dtcType_(o.dtcType_), outFile_(o.outFile_), pol_(o.pol_) {}

    // Accessor Function
    /**
     * @brief Returns the location of a dectector
     * @return the location
     */
    std::vector<int> loc() {return location_;}
    /**
     * @brief Returns the type of detector
     */
    dtcOutType type() {return dtcType_;}
    /**
     * @brief returns the output file name
     */
    std::string outfile() {return outFile_;}
    /**
     * @brief returns the Polarzation
     */
    Polarization pol() {return pol_;}
    // Output functions
    /**
     * @brief Ouputs the field info for the detector
     * @details Field accepts a field input and dielectric constant value and returns that value at it's locaton
     *
     * @param field_in The field that is being looked at
     * @param eps The dielectric constatn at that frequency
     *
     * @return [description]
     */
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
