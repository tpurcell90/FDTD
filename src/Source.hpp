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
    /**
     * @brief Constructs a Source
     * @details Constructs a source given a pulse and Polarization and location
     *
     * @param loc location of the source
     * @param pol polarization of the pulse
     * @param prof the Pulse for the source
     */
    Source(Pulse<double> prof, Polarization pol,std::vector<int> loc) : profile(prof), polarization_(pol), location_(loc) {}
    // Copy Constructor
    /**
     * @brief Copy Constructor
     * @details Creates a new source from the input source
     *
     * @param o The source to be copied
     * @param _ [description]
     * @param _ [description]
     */
    Source(const Source& o) : profile(o.profile), polarization_(o.polarization_), location_(o.location_) {}
    //Access Functions
    /**
     * @brief Get the source polarization
     * @return Source polarization
     */
    Polarization pol() {return polarization_;}
    /**
     * @brief Get the Source location
     * @return The location of the source
     */
    std::vector<int> loc()  {return location_;}
    /**
     * @brief Get the Source's pulse
     * @return The Pulse of teh source
     */
    Pulse<T> prof() {return profile;}
};

#endif
