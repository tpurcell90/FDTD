#include "mpiInterface.hpp"
#include <iostream>
#include <iomanip>

mpiInterface::mpiInterface() : boost::mpi::communicator()
{
    std::tie(npX_, npY_, npZ_) = numgrid(size());
    if (rank() == 0)
        std::cout << " A processor grid of dimension (" << npX_ << ", " << npY_ << ") will be used" << std::endl;
    mypY_  =  rank();
    mypX_  = 0;
    mypZ_ = 0;
    // sl_init_(context_, npX_, npY_);
    // blacs_gridinfo_(context_, npX_, npY_, mypRow_, mypCol_);
}

// std::vector<int> mpiInterface::descinit(const int ndim, const int ncol) const
// {
//     // std::cout << "contex     " << context_ << std::endl;
//     std::vector<int> desc(9);
//     const int localrow = numroc_(ndim, 1, mypRow_, 0, npX_);
//     int info;
//     descinit_(desc.data(), ndim, ncol, 1, 1, 0, 0, context_, std::max(1,localrow), info);
//     return desc;
// }

std::tuple<int,int,int> mpiInterface::getLocxLocyLocz(std::vector<real_grid_ptr> weights) const
{
    // std::cout << "in" << std::endl;
    std::vector<double> yE_weights( weights[0]->y(), 0 );
    // std::vector<double> yH_weights( weights[0]->y(), 0 );
    for(auto& grid : weights)
    {
        for(int rr = 0; rr < weights[0]->y(); rr++)
        {
            yE_weights[rr] += std::accumulate(&grid->point(0,rr,0), &grid->point(0,rr,0)+grid->x()*grid->z(), 0);
            // std::cout << yE_weights[rr] << std::endl;
        }
    }

    // if(weights[2]->point(0,0) == 10.2)
    // {
    //     // std::cout << "HI!!!!!" << std::endl;
    //     // TE Mode
    //     for(int rr = 0; rr < weights[0]->y(); rr++)
    //     {
    //         yE_weights[rr]-> = std::accumulate(&weights[0]->point(0,rr,0), &weights[0]->point(0,rr,0)+weights[0]->x()*weights[0]->z(), 0);
    //         yE_weights[rr]->+= std::accumulate(&weights[1]->point(0,rr,0), &weights[1]->point(0,rr,0)+weights[1]->x()*weights[1]->z(), 0);
    //         yH_weights[rr]-> = std::accumulate(&weights[2]->point(0,rr,0), &weights[2]->point(0,rr,0)+weights[2]->x()*weights[2]->z(), 0);
    //     }

    //     // if(rank() == 0)
    //     //     for(auto& val : colE_weights )
    //     //         std::cout << val << std::endl;
    // }
    // else
    // {
    //     //TM Mode
    //     for(int rr = 0; rr < weights[0]->y(); rr++)
    //     {
    //         yH_weights[rr]-> = std::accumulate(&weights[0]->point(0,rr,0), &weights[0]->point(0,rr,0)+weights[0]->x()*weights[0]->z(), 0);
    //         yH_weights[rr]->+= std::accumulate(&weights[1]->point(0,rr,0), &weights[1]->point(0,rr,0)+weights[1]->x()*weights[1]->z(), 0);
    //         yE_weights[rr]-> = std::accumulate(&weights[2]->point(0,rr,0), &weights[2]->point(0,rr,0)+weights[2]->x()*weights[2]->z(), 0);
    //     }
    // }

    double yWeightSum_E = std::accumulate(yE_weights.begin(), yE_weights.end(), 0.0);
    // double yWeightSum_H = std::accumulate(yH_weights.begin(), yH_weights.end(), 0.0);

    double yWeightAvg_E = static_cast<double>(yWeightSum_E) / static_cast<double>(npY_);
    // double yWeightAvg_H = static_cast<double>(yWeightSum_H) / static_cast<double>(npY_);
    double valWeight = 0.0;

    std::vector<int> yStartBound(npY_, 0);
    std::vector<int> yEndBound(npY_, weights[0]->y()-1);
    int curY = 0;
    for(int cc = 0; cc < yEndBound.size()-1; cc++)
    {
        valWeight = 0;
        curY = yStartBound[cc];
        while(valWeight + yE_weights[curY] / 2.0 < yWeightAvg_E )
        {
            valWeight += yE_weights[curY];
            curY++;
        }
        yEndBound[cc]     = curY-1;
        yStartBound[cc+1] = curY;
    }
    barrier();
    return std::make_tuple(weights[0]->x()+2, yEndBound[mypY_] - yStartBound[mypY_] + 3, weights[0]->z()+2 );

}

// std::vector<int> getLocSz(Grid2d<int> weights)
// {
//     double weightAvg = 0;
//     int loc = 0;
//     for(int dd = 0; dd < dimMax; dd++)
//         rowWeightAvg += sasum(weights.x(), weights.point(0,rr), 1);
// }
std::tuple<int,int,int> mpiInterface::numroc(const int ndim, const int ncol, const int nface) const
{
    return std::make_tuple(ndim + 2, numroc_(ncol, 1, mypY_, 0, npY_), nface + 2 );
}

static std::tuple<int, int, int> numgrid(int numproc)
{
    // int sq = static_cast<int>(std::sqrt(static_cast<double>(numproc)))+1;
    // std::pair<int,int> returnVals = {-1,-1};
    // for (int i = sq; i != 0; i--)
    //     if (numproc % i == 0)
    //     {
    //         returnVals = {i, numproc/i};
    //         break;
    //     };
    // if (returnVals.first == -1)
    //     throw std::logic_error("Error determining MPI Cartesian topology.");
    // std::pair<int,int> returnVals = {numproc,1};
    std::tuple<int,int,int> returnVals = std::make_tuple( 1,numproc,1 );
    return returnVals;
}


