#include <src/DTC/parallelDTC_BIN.hpp>

parallelDetectorBINReal::parallelDetectorBINReal(int dtcNum,  std::vector<real_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt) :
    parallelDetectorBaseReal(dtcNum, grid, SI, loc, sz, type, classType, timeInterval, firstComp, a, I0, dt),
    outFile_(out_name),
    gridVals_(sz_[0], 0.0)
{
    std::ofstream dat(outFile_.c_str(), std::ios::out | std::ios::binary);
    if(firstComp_ == DIRECTION::X)
        dat.write("x", sizeof(char) );
    else if(firstComp_ == DIRECTION::Y)
        dat.write("y", sizeof(char) );
    else
        dat.write("z", sizeof(char) );
    dat.write(reinterpret_cast<const char*>(&sz_[0]), sz_.size()*sizeof(int) );
    dat.write(reinterpret_cast<const char*>(&loc_[0]), loc_.size()*sizeof(int) );
    dat.close();
}
void parallelDetectorBINReal::output(double t)
{
    for(auto & field : fields_)
        field->getField();
    if(fields_[0]->master())
    {
        double tt = t*tConv_;
        std::ofstream dat(outFile_, std::ios::app | std::ios::binary);

        dat.write(reinterpret_cast<char *>(&tt), sizeof(tt));
        for(int kk = 0; kk < sz_[2]; ++kk)
        {
            for(int jj = 0; jj < sz_[1]; ++jj)
            {
                std::fill_n(gridVals_.begin(), gridVals_.size(), 0.0);
                for(auto & field : fields_)
                    outputFunction_(&field->outGrid()->point(0,jj,kk), &field->outGrid()->point(0,jj,kk)+sz_[0], gridVals_.data(), convFactor_);
                dat.write(reinterpret_cast<char *>(&gridVals_[0]), gridVals_.size()*sizeof(double));
            }
        }
        dat.close();
    }
}

parallelDetectorBINCplx::parallelDetectorBINCplx(int dtcNum,  std::vector<cplx_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt) :
    parallelDetectorBaseCplx(dtcNum, grid, SI, loc, sz, type, classType, timeInterval, firstComp, a, I0, dt),
    outFile_(out_name),
    gridVals_(sz_[0], 0.0)
{
    std::ofstream dat(outFile_.c_str(), std::ios::out | std::ios::binary);
    if(firstComp_ == DIRECTION::X)
        dat.write("x", sizeof(char) );
    else if(firstComp_ == DIRECTION::Y)
        dat.write("y", sizeof(char) );
    else
        dat.write("z", sizeof(char) );
    dat.write(reinterpret_cast<const char*>(&sz_[0]), sz_.size()*sizeof(int) );
    dat.write(reinterpret_cast<const char*>(&loc_[0]), loc_.size()*sizeof(int) );
    dat.close();
}
void parallelDetectorBINCplx::output(double t)
{
    for(auto & field : fields_)
        field->getField();
    if(fields_[0]->master())
    {
        double tt = t*tConv_;
        std::ofstream dat(outFile_, std::ios::app | std::ios::binary);

        dat.write(reinterpret_cast<char *>(&tt), sizeof(tt));
        for(int kk = 0; kk < sz_[2]; ++kk)
        {
            for(int jj = 0; jj < sz_[1]; ++jj)
            {
                std::fill_n(gridVals_.begin(), gridVals_.size(), 0.0);
                for(auto & field : fields_)
                    outputFunction_(&field->outGrid()->point(0,jj,kk), &field->outGrid()->point(0,jj,kk)+sz_[0], gridVals_.data(), convFactor_);
                dat.write(reinterpret_cast<char *>(&gridVals_[0]), gridVals_.size()*sizeof(double));
            }
        }
        dat.close();
    }
}
