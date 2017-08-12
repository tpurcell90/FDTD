#ifndef FDTD_GRID
#define FDTD_GRID

#include <memory>
#include <string>
#include <vector>
#include <complex>
#include <algorithm>
#include <assert.h>
#include <fstream>
#include <functional>
#include <iostream>

template <typename T> class Grid
{
protected:
    static unsigned int memSize; //!< total memory of the grid
    std::array<int,3> n_vec_; //!< the number of grid points in all directions
    std::array<double,3>  d_; //!< the grid point spacing in all directions
    std::unique_ptr<T[]> vals_; //!< the field values in the grid

public:
    /**
     * @brief      Constructs a grid
     *
     * @param[in]  n_vec  grid size in terms of number of grid points in all directions
     * @param[in]  d      the gird spacing in all directions
     */
    Grid(std::array<int,3> n_vec, std::array<double,3> d) :
        n_vec_(n_vec),
        d_(d),
        vals_(std::unique_ptr<T[]>(new T[std::accumulate (n_vec.begin(), n_vec.end(), 1, std::multiplies<int>())]))
    {
        std::fill_n(vals_.get(), std::accumulate(n_vec_.begin(), n_vec_.end(), 1, std::multiplies<int>()), T(0.0));
        memSize += sizeof(T)*size();
    }

    /**
     * @brief      Copies a grid into a new grid object
     *
     * @param[in]  o     Grid to be copied
     */
    Grid(const Grid& o) :
        n_vec_(o.n_vec_),
        d_(o.d_),
        vals_(std::unique_ptr<T[]>(new T[std::accumulate (o.n_vec_.begin(), o.n_vec_.end(), 1, std::multiplies<int>())]))
    {
       std::copy_n(o.vals_.get(), std::accumulate (n_vec_.begin(), n_vec_.end(), 1, std::multiplies<int>()), vals_.get());
       memSize += sizeof(T)*size();
    }

    /**
     * @brief Move constructor
     */
    /**
     * @brief      Moves a grid into a new grid object
     *
     * @param[in]  o  Grid to be moved
     */
    Grid(Grid&& o) :
        n_vec_(o.n_vec_),
        d_(o.d_),
        vals_(std::move(o.vals_))
    {
        o.n_vec_ = std::array<int,3>(3,0);
    }

    /**
     * @brief Grid destructor
     */
    ~Grid(){memSize -= sizeof(T)*size();}

    /**
     * @return the size of the grid
     */
    inline int size() const { return std::accumulate(n_vec_.begin(), n_vec_.end(), 1, std::multiplies<int>());}

    /**
     * @return the data stored in the grid
     */
    T* data() {return vals_.get();}

    /**
     * @return a const form of the data
     */
    const T* data() const { return vals_.get(); }

    /**
     * @return the number of grid points in the x direction
     */
    inline int x() const {return n_vec_[0];}

    /**
     * @return the number of grid points in the y direction
     */
    inline int y() const {return n_vec_[1];}

    /**
     * @return the number of grid points in the z direction
     */
    inline int z() const {return n_vec_[2];}

    /**
     * @return the number of grid points in all directions
     */
    inline std::array<int,3> n_vec(){return n_vec_;}

    /**
     * @return the grid spacing in the x direction
     */
    inline double dx() const {return d_[0];}

    /**
     * @return the grid spacing in the y direction
     */
    inline double dy() const {return d_[1];}

    /**
     * @return the grid spacing in the z direction
     */
    inline double dz() const {return d_[2];}

    /**
     * @return the grid spacing in all directions
     */
    inline std::array<double,3> d(){return d_;}


    /**
     * @brief      returns the value at a point x_val,y_val,0
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     *
     * @return     reference to the data point at (x,y,0)
     */
    T& point(const int x_val, const int y_val)
    {
        assert(0 <= x_val && x_val < n_vec_[0] && 0 <= y_val && y_val < n_vec_[1] );
        return vals_[y_val*n_vec_[0]*n_vec_[2] + x_val];
    }

    /**
     * @brief      returns the value at a point x_val,y_val,0
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     *
     * @return     reference to the data point at (x,y,0)
     */
    const T& point(const int x_val, const int y_val) const{ assert(0 <= x_val && x_val < n_vec_[0] && 0 <= y_val && y_val < n_vec_[1] );return vals_[y_val*n_vec_[0]*n_vec_[2] + x_val];}

    /**
     * @brief      returns the value at a point x_val,y_val,0
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     *
     * @return     reference to the data point at (x,y,0)
     */
    T& operator()(const int x_val, const int y_val) {assert(0 <= x_val && x_val < n_vec_[0] && 0 <= y_val && y_val < n_vec_[1] ); return vals_[y_val*n_vec_[0]*n_vec_[2] + x_val]; }

    /**
     * @brief      returns the value at a point x_val,y_val,0
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     *
     * @return     reference to the data point at (x,y,0)
     */
    const T& operator()(const int x_val, const int y_val) const { assert(0 <= x_val && x_val < n_vec_[0] && 0 <= y_val && y_val < n_vec_[1] );return vals_[(y_val)*n_vec_[0]*n_vec_[2] + x_val]; }

    /**
     * @brief      returns the value at a point x_val,y_val,z_val
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     * @param[in]  z_val  grid coordinate in the z direction
     *
     * @return     reference to the data point at (x,y,z)
     */
    T& point(const int x_val, const int y_val, const int z_val)
    {
        assert(0 <= x_val && x_val < n_vec_[0] && 0 <= y_val && y_val < n_vec_[1] && 0 <= z_val && z_val < n_vec_[2] );
        return vals_[z_val*n_vec_[0] + y_val*n_vec_[2]*n_vec_[0] + x_val];
    }

    /**
     * @brief      returns the value at a point x_val,y_val,z_val
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     * @param[in]  z_val  grid coordinate in the z direction
     *
     * @return     reference to the data point at (x,y,z)
     */
    const T& point(const int x_val, const int y_val, const int z_val) const
    {
        assert(0 <= x_val && x_val < n_vec_[0] && 0 <= y_val && y_val < n_vec_[1] && 0 <= z_val && z_val < n_vec_[2] );
        return vals_[z_val*n_vec_[0] + y_val*n_vec_[2]*n_vec_[0] + x_val];
    }

    /**
     * @brief      returns the value at a point x_val,y_val,z_val
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     * @param[in]  z_val  grid coordinate in the z direction
     *
     * @return     reference to the data point at (x,y,z)
     */
    T& operator()(const int x_val, const int y_val, const int z_val)
    {
        assert(0 <= x_val && x_val < n_vec_[0] && 0 <= y_val && y_val < n_vec_[1] && 0 <= z_val && z_val < n_vec_[2] );
        return vals_[z_val*n_vec_[0] + y_val*n_vec_[2]*n_vec_[0] + x_val];
    }

    /**
     * @brief      returns the value at a point x_val,y_val,z_val
     *
     * @param[in]  x_val  gird coordinate in the x direction
     * @param[in]  y_val  grid coordinate in the y direction
     * @param[in]  z_val  grid coordinate in the z direction
     *
     * @return     reference to the data point at (x,y,z)
     */
    const T& operator()(const int x_val, const int y_val, const int z_val) const
    {
        assert(0 <= x_val && x_val < n_vec_[0] && 0 <= y_val && y_val < n_vec_[1] && 0 <= z_val && z_val < n_vec_[2] );
    return vals_[z_val*n_vec_[0] + (y_val)*n_vec_[2]*n_vec_[0] + x_val];
}

    /**
     * @brief      Prints out a field snapshot in the area specified by lower, left, back corner loc and size of sz, of the field values in a grid, in a box format separate XY planes separated by a line
     *
     * @param[in]  filename  output file name
     * @param[in]  loc       location of the lower left corner of the field area desired
     * @param[in]  sz        size in grid points of the region you want to print out
     * @param[in]  fxn       function defining what to print out from the field, real, imaginary, magnitude intensity
     */
    void gridOutBox(std::string filename, std::array<int,3> loc, std::array<int,3> sz, std::function<double(T)> fxn)
    {
        std::ofstream outFile;
        outFile.open(filename,std::ios_base::out);
        for(int kk = loc[2]; kk < loc[2]+sz[2]; ++kk)
        {
            for(int jj = loc[1]; jj < loc[1]+sz[1]; ++jj)
            {
                for(int ii = loc[0]; ii < loc[0]+sz[0]; ++ii)
                {
                    outFile << fxn(point(ii,jj,kk)) << "\t";
                }
                outFile << "\n";
            }
            outFile << "\n";
        }
    }

    /**
     * @brief      Prints out a field snapshot of the field values in a grid, in a box format separate XY planes separated by a line
     *
     * @param[in]  filename  output file name
     * @param[in]  fxn       function defining what to print out from the field, real, imaginary, magnitude intensity
     */
    void gridOutBox(std::string filename, std::function<double(T)> fxn)
    {
        std::ofstream outFile;
        outFile.open(filename,std::ios_base::out);
        for(int kk = 0; kk < n_vec_[2]; ++kk)
        {
            for(int jj = 0; jj < n_vec_[1]; ++jj)
            {
                for(int ii = 0; ii < n_vec_[0]; ++ii)
                {
                    outFile<< fxn(point(ii,jj,kk)) << "\t";
                }
                outFile << "\n";
            }
            outFile << "\n";
        }
    }

    /**
     * @brief      Prints out a field snapshot in the area specified by lower, left, back corner loc and size of sz, of the field values in a grid, in a coordinate list format
     *
     * @param[in]  filename  output file name
     * @param[in]  loc       location of the lower left corner of the field area desired
     * @param[in]  sz        size in grid points of the region you want to print out
     * @param[in]  fxn       function defining what to print out from the field, real, imaginary, magnitude intensity
     */
    void gridOutList(std::string filename, std::array<int,3> loc, std::array<int,3> sz, std::function<double(T)> fxn)
    {
        std::ofstream outFile;
        outFile.open(filename,std::ios_base::out);
        for(int ii = loc[0]; ii < loc[0]+sz[0]; ++ii)
        {
            for(int jj = loc[1]; jj < loc[1]+sz[1]; ++jj)
            {
                for(int kk = loc[2]; kk < loc[2]+sz[2]; ++kk)
                {
                    outFile << ii << "\t" << jj << '\t' << kk << '\t' << fxn( point(ii,jj,kk) ) <<  "\n";
                }
            }
        }
    }

    /**
     * @brief      Prints out a field snapshot in the area specified by lower, left, back corner loc and size of sz, of the field values in a grid, in a coordinate list format
     *
     * @param[in]  filename  output file name
     * @param[in]  fxn       function defining what to print out from the field, real, imaginary, magnitude intensity
     */
    void gridOutList(std::string filename, std::function<double(T)> fxn)
    {
        std::ofstream outFile;
        outFile.open(filename,std::ios_base::out);
        for(int ii = 0; ii < n_vec_[0]; ++ii)
        {
            for(int jj = 0; jj < n_vec_[1]; ++jj)
            {
                for(int kk = 0; kk < n_vec_[2]; ++kk)
                {
                    outFile << ii << "\t" << jj << '\t' << kk << '\t' << fxn( point(ii,jj,kk) ) << "\n";
                }
            }
        }
    }

    /**
     * @brief      Prints out a field snapshot in the area specified by lower, left, back corner loc and size of sz, of the field values in a grid, in a box format separate XY planes separated by a line
     *
     * @param[in]  filename     output file name
     * @param[in]  otherFields  A vector storing other field grid pointers to need to be added to the output
     * @param[in]  loc          location of the lower left corner of the field area desired
     * @param[in]  sz           size in grid points of the region you want to print out
     * @param[in]  fxn          function defining what to print out from the field, real, imaginary, magnitude intensity
     */
    void gridOutBox(std::string filename, std::vector<std::shared_ptr<Grid<T>>> otherFields, std::array<int,3> loc, std::array<int,3> sz, std::function<double(T)> fxn)
    {
        std::ofstream outFile;
        outFile.open(filename,std::ios_base::out);
        T pt = 0.0;
        for(int kk = loc[2]; kk < loc[2]+sz[2]; ++kk)
        {
            for(int jj = loc[1]; jj < loc[1]+sz[1]; ++jj)
            {
                for(int ii = loc[0]; ii < loc[0]+sz[0]; ++ii)
                {
                    pt = fxn(point(ii,jj,kk));
                    for(auto& grid : otherFields)
                        pt += fxn(grid->point(ii,jj,kk));
                    outFile << pt << "\t";
                }
                outFile<< "\n";
            }
           outFile<< "\n";
        }
    }

    /**
     * @brief      Prints out a field snapshot in the area specified by lower, left, back corner loc and size of sz, of the field values in a grid, in a coordinate list format
     *
     * @param[in]  filename     output file name
     * @param[in]  otherFields  A vector storing other field grid pointers to need to be added to the output
     * @param[in]  loc          location of the lower left corner of the field area desired
     * @param[in]  sz           size in grid points of the region you want to print out
     * @param[in]  fxn          function defining what to print out from the field, real, imaginary, magnitude intensity
     */
    void gridOutList(std::string filename, std::array<int,3> loc, std::vector<std::shared_ptr<Grid<T>>> otherFields, std::array<int,3> sz, std::function<double(T)> fxn)
    {
        std::ofstream outFile;
        outFile.open(filename,std::ios_base::out);
        T pt = 0.0;
        for(int ii = loc[0]; ii < loc[0]+sz[0]; ++ii)
        {
            for(int jj = loc[1]; jj < loc[1]+sz[1]; ++jj)
            {
                for(int kk = loc[2]; kk < loc[2]+sz[2]; ++kk)
                {
                    pt = fxn(point(ii,jj,kk));
                    for(auto& grid : otherFields)
                        pt += fxn(grid->point(ii,jj,kk));
                    outFile << ii << "\t" << jj < '\t' << kk << '\t' << pt << "\n";
                }
                outFile << "\n";
            }
        }
    }
};

template <typename T> unsigned int Grid<T>::memSize = 0;

#endif