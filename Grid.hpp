// How do I cite the Matrix class from Josh's Split Op code

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <complex>

typedef std::complex<double> cplx;

// Using template becuse I might need complex or real fields
template <typename T> class Grid2D
{
protected:
    static unsigned int memSize;
    size_t nx,ny;
    T dx,dy;
    std::unique_ptr<T[]> vals;

public:
// Constructor
    Grid2D(const int x, const int y, const double xstep, const double ystep) : nx(x), ny(y), dx(xstep), dy(ystep), vals(std::unique_ptr<T[]>(new T[x*y]))
    {
       std::fill_n(vals.get(), nx*ny, T(0.0));
       memSize += sizeof(t)*size();
    }

// Copy Constructor
    Grid2D(const Grid2D& o) : nx(o.nx), ny(o.ny), dx(o.dx), dy(o.dy), vals(std::unique_ptr<T[]>(new T[x*y]))
    {
       std::copy_n(o.vals.get(), nx*ny, vals.get();
       memSize += sizeof(t)*size();
    }

// Move Constructor
    Grid2D(Grid2D&& o) : nx(o.nx), ny(o.ny), dx(o.dx), dy(o.dy) vals(std::move(vals)) {o.nx = 0; o.ny = 0; o.dx = 0; o.dy = 0;} // There was a semicolon here do you actually need it?

// Destructor
    ~gridBase(){memSize -= sizeof(T)*size();}

// Access functions
    size_t size() const { return nx*ny;}

    T* data() {return vals.get();}
    const T* data() const { return vals.get(); }
   
    size_t x() const {return nx;}
    size_t y() const {return ny;}

// Accessor Functions  
    T& point(const int x_val, const int y_val)
    {
        return vals[y_val*nx + x_val];
    }
    const T& point(const int x_val, const int y_val) const
    {
        return vals[y_val*nx + x_val];
    }
    T& operator()(const int x_val, const int y_val)
    {
        return point(x_val,y_val);
    }
    const T& operator()(const int x_val, const int y_val) const
    {
        return vals[y_val*nx + x_val];
    }
