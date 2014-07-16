#include <string>
#include <complex>
#include <function>

using namespace std;

template <typename T> class Source
{
protected:
    function profile
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

