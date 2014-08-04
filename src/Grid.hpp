#ifndef FDTD_GRID
#define FDTD_GRID

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <complex>
#include <algorithm>
#include <assert.h>

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
       memSize += sizeof(T)*size();
    }

// Copy Constructor
    Grid2D(const Grid2D& o) : nx(o.nx), ny(o.ny), dx(o.dx), dy(o.dy), vals(std::unique_ptr<T[]>(new T[nx*ny]))
    {
       std::copy_n(o.vals.get(), nx*ny, vals.get());
       memSize += sizeof(T)*size();
    }

// Move Constructor
    Grid2D(Grid2D&& o) : nx(o.nx), ny(o.ny), dx(o.dx), dy(o.dy), vals(std::move(vals)) {o.nx = 0; o.ny = 0;}

// Destructor
    ~Grid2D(){memSize -= sizeof(T)*size();}

// Access functions
    size_t size() const { return nx*ny;}

    T* data() {return vals.get();}
    const T* data() const { return vals.get(); }

    size_t x() const {return nx;}
    size_t y() const {return ny;}

// Accessor Functions
    /*T& point(const int x_val, const int y_val) { return vals[(y_val%ny)*nx + x_val%nx];}
    const T& point(const int x_val, const int y_val) const{ return vals[(y_val%ny)*nx + x_val%nx];}
    T& operator()(const int x_val, const int y_val) { return point(x_val,y_val); }
    const T& operator()(const int x_val, const int y_val) const { return vals[(y_val%ny)*nx + x_val%nx]; }*/

    T& point(const int x_val, const int y_val)
    {
        assert(0 <= x_val && x_val < nx && 0 <= y_val && y_val < ny );
        return vals[y_val*nx + x_val];
    }
    const T& point(const int x_val, const int y_val) const{ assert(0 <= x_val && x_val < nx && 0 <= y_val && y_val < ny );return vals[y_val*nx + x_val];}
    T& operator()(const int x_val, const int y_val) {assert(0 <= x_val && x_val < nx && 0 <= y_val && y_val < ny ); return point(x_val,y_val); }
    const T& operator()(const int x_val, const int y_val) const { assert(0 <= x_val && x_val < nx && 0 <= y_val && y_val < ny );return vals[(y_val)*nx + x_val]; }

// Flux calculation
    double flux(std::vector<int> loc, double eps)
    {
        T& pt = point(loc[0],loc[1]);
        // NEed to do something to determine if complex or not
        //return pow(point(loc[0],loc[1]),2.00) * eps / 2.00;
        return std::real(pt*std::conj(pt) * eps / 2.0);
    }
};

template <typename T> unsigned int Grid2D<T>::memSize = 0;

#endif
