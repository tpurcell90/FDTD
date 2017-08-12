#ifndef FDTD_PARALLELDETECTORSTORAGESTRUCTS
#define FDTD_PARALLELDETECTORSTORAGESTRUCTS

#include <vector>
#include <boost/serialization/vector.hpp>
#include <array>

struct slaveProcDtc
{
    int masterProc_; //!< masterProc_ where to send the data
    int stride_; //!< stride of the copy
    std::array<int,3> loc_; //!< loc_ lower left corner of the region inside the slave process to send the data to the master processor
    std::array<int,3> sz_; //!< sz_ size of the region in grid points of the detection region inside the slave process
    std::array<int,3> opSz_; //!< size of the field listed in order needed for the copy
    std::array<int,3> addVec1_; //!< vector to determine what component to iterate over for the first copy loop
    std::array<int,3> addVec2_; //!< vector to determine what component to iterate over for the second copy loop
};

struct slaveProcInfo
{
    int slaveProc_; //!< masterProc_ where to send the data
    int stride_; //!< stride of the copy
    std::array<int,3> loc_; //!< loc_ lower left corner of the region inside the slave process to send the data to the master processor
    std::array<int,3> sz_; //!< sz_ size of the region in grid points of the detection region inside the slave process
    std::array<int,3> addVec1_; //!< vector to determine what component to iterate over for the first copy loop
    std::array<int,3> addVec2_; //!< vector to determine what component to iterate over for the second copy loop
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & slaveProc_;
        ar & stride_;
        ar & loc_;
        ar & sz_;
        ar & addVec1_;
        ar & addVec2_;
    }
};

struct copyProcDtc
{
    int stride_; //!< stride of the copy
    int strideOutGrid_; //!< stride inside the output grid
    std::array<int,3> loc_; //!< loc_ lower left corner of the region inside the slave process to send the data to the master processor
    std::array<int,3> sz_; //!< sz_ size of the region in grid points of the detection region inside the slave process
    std::array<int,3> opSz_; //!< size of the field listed in order needed for the copy
    std::array<int,3> addVec1_; //!< vector to determine what component to iterate over for the first copy loop
    std::array<int,3> addVec2_; //!< vector to determine what component to iterate over for the second copy loop
    std::array<int,3> locOutGrid_; //!< loction of lower, left, back corner in the output grid
};

struct fInParam
{
    int stride_; //!< stride of the copy
    std::array<int,3> addVec1_; //!< vector to determine what component to iterate over for the first copy loop
    std::array<int,3> addVec2_; //!< vector to determine what component to iterate over for the second copy loop
    std::array<int,3> sz_; //!< sz_ size of the region in grid points of the detection region inside the slave process
    std::array<int,3> loc_; //!< loc_ lower left corner of the region inside the slave process to send the data to the master processor
};

#endif