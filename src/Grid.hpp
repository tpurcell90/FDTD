#ifndef FDTD_GRID
#define FDTD_GRID

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <complex>
#include <algorithm>
#include <assert.h>
#include <fstream>

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
    /**
     * @brief Constructs a grid
     * @details Constructs a grid given its size and spacing
     *
     * @param x Real size of the grid in the x direction
     * @param y Real size of the grid in the y direction
     * @param xstep grid step size in the x direction
     * @param ystep grid step size in the y direction
     */
    Grid2D(const int x, const int y, const double xstep, const double ystep) : nx(x), ny(y), dx(xstep), dy(ystep), vals(std::unique_ptr<T[]>(new T[x*y]))
    {
       std::fill_n(vals.get(), nx*ny, T(0.0));
       memSize += sizeof(T)*size();
    }

// Copy Constructor
    /**
     * @brief Copies a grid into a new grid object
     * @details Consturcts a new grid from the parameters of the copied grid
     *
     * @param o [description]
     */
    Grid2D(const Grid2D& o) : nx(o.nx), ny(o.ny), dx(o.dx), dy(o.dy), vals(std::unique_ptr<T[]>(new T[nx*ny]))
    {
       std::copy_n(o.vals.get(), nx*ny, vals.get());
       memSize += sizeof(T)*size();
    }

// Move Constructor
    /**
     * @brief Move constructor
     * @details Copies the grid into a new object and destructs the old grid
     */
    Grid2D(Grid2D&& o) : nx(o.nx), ny(o.ny), dx(o.dx), dy(o.dy), vals(std::move(vals)) {o.nx = 0; o.ny = 0;}

// Destructor
    /**
     * @brief Grid destructor
     * @details Destructs the Grid
     */
    ~Grid2D(){memSize -= sizeof(T)*size();}

// Access functions
    /**
     * @brief returns teh size of teh grid
     */
    size_t size() const { return nx*ny;}
    /**
     * @brief returns the data stored in teh grid
     */
    T* data() {return vals.get();}
    /**
     * @brief Returns a const form of the data
     */
    const T* data() const { return vals.get(); }
    /**
     * @brief returns the number of grid points in teh x direction
     */
    size_t x() const {return nx;}
    /**
     * @brief Reutns the number of grid points in the y direction
     */
    size_t y() const {return ny;}

// Accessor Functions
    /*T& point(const int x_val, const int y_val) { return vals[(y_val%ny)*nx + x_val%nx];}
    const T& point(const int x_val, const int y_val) const{ return vals[(y_val%ny)*nx + x_val%nx];}
    T& operator()(const int x_val, const int y_val) { return point(x_val,y_val); }
    const T& operator()(const int x_val, const int y_val) const { return vals[(y_val%ny)*nx + x_val%nx]; }*/
    /**
     * @brief returns the value at a point x_val,y_val
     *
     * @param x_val [gird coordinate in the x direction]
     * @param y_val [grid coordinate in the y direction]
     *
     */
    T& point(const int x_val, const int y_val)
    {
        assert(0 <= x_val && x_val < nx && 0 <= y_val && y_val < ny );
        return vals[y_val*nx + x_val];
    }
    /**
     * @brief returns the value at a point x_val,y_val
     *
     * @param x_val [gird coordinate in the x direction]
     * @param y_val [grid coordinate in the y direction]
     *
     */
    const T& point(const int x_val, const int y_val) const{ assert(0 <= x_val && x_val < nx && 0 <= y_val && y_val < ny );return vals[y_val*nx + x_val];}
    /**
     * @brief returns the value at a point x_val,y_val
     *
     * @param x_val [gird coordinate in the x direction]
     * @param y_val [grid coordinate in the y direction]
     *
     */
    T& operator()(const int x_val, const int y_val) {assert(0 <= x_val && x_val < nx && 0 <= y_val && y_val < ny ); return point(x_val,y_val); }
    /**
     * @brief returns the value at a point x_val,y_val
     *
     * @param x_val [gird coordinate in the x direction]
     * @param y_val [grid coordinate in the y direction]
     *
     */
    const T& operator()(const int x_val, const int y_val) const { assert(0 <= x_val && x_val < nx && 0 <= y_val && y_val < ny );return vals[(y_val)*nx + x_val]; }
    /**
     * @brief Prints out a grid snapshot
     * @details Prints the values for the entire grid into the file filename
     *
     * @param filename name of the file for the grid to be outputted to
     */
    void gridOut(std::string filename)
    {
        std::ofstream outFile;
        outFile.open(filename,std::ios_base::app);
        for(int ii = 0; ii < nx; ii++)
            for(int jj =0; jj < ny; jj++)
                outFile<< ii << "\t" << jj << "\t" << point(ii,jj) <<std::endl;
    }
// Flux calculation
    /**
     * @brief Calculates the flux at a point
     * @details cacluates the flux of the field at a point using the power as the basis
     *
     * @param loc location to calculate the flux
     * @param eps dielectric function
     *
     * @return [description]
     */
    double flux(std::vector<int> loc, double eps)
    {
        T& pt = point(loc[0],loc[1]);
        // NEed to do something to determine if complex or not
        //return pow(point(loc[0],loc[1]),2.00) * eps / 2.00;
        return std::real(pt*std::conj(pt) * sqrt(eps) / 2.0); // sqrt(eps) from eps*c
    }
};

template <typename T> unsigned int Grid2D<T>::memSize = 0;

#endif
