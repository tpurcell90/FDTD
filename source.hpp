#include <string>
#include <complex>

using namespace std;

// Needs modification for complex I think
template <typename T> class Source
{
protected:
    std::function<double(double)> profile
    string polarization
    
public:
    // Constructor
    Source(const function<double(double)> prof, const string pol) : profile(prof), polarization(pol) {}
    
    // Copy Constructor
    Source(const Source o) : profile(o.profile), polarization(o.polarization) {}
    
    // Destructor
    ~Source(){}
    
    //Access Functions
    string pol() const {return polarization;}
    double field(double t) const {return profile(t);}

}

