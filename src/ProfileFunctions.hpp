#ifndef FDTD_PROFILEFUNCTION
#define FDTD_PROFILEFUNCTION

#include <string>
#include <vector>

template <typename T> class ProfileFunctions
{
protected:
    std::vector<T> Param;
    std::string Type;

public:
    // Constructor
    ProfileFunctions(std::vector<T> param, std::string type) : Param(param), Type(type) {}
    //  Copy Constructor
    ProfileFunctions(const ProfileFunctions& o) : Param(o.Param), Type(o.Type) {}
    
    //Acessor Functions
    std::vector<T> param() {return Param;}
    std::string type() {return Type;}

    //Pulse functions
    T gauss_pulse(double t)
    {
        // look for the best way to calculate gaussian pulse
    }

};
#endif
