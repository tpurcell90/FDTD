#include <src/DTC/parallelDTC_COUT.hpp>

parallelDetectorCOUTReal::parallelDetectorCOUTReal(int dtcNum,  std::vector<real_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt) :
    parallelDetectorBaseReal(dtcNum, grid, SI, loc, sz, type, classType, timeInterval, firstComp, a, I0, dt),
    dtc2cart_(3, 0)
{
    if(firstComp == DIRECTION::X)
        dtc2cart_ = {{ 0, 1, 2 }};
    else if(firstComp == DIRECTION::Y )
        dtc2cart_ = {{ 1, 2, 0 }};
    else
        dtc2cart_ = {{ 2, 0, 1 }};
}

void parallelDetectorCOUTReal::output(double t)
{
    fields_[0]->getField();
    if(fields_[0]->master() )
    {
        std::cout << t*tConv_ << "\t" << loc_[dtc2cart_[0]]*dirConv_[0] << "\t" << loc_[dtc2cart_[1]] * dirConv_[1] << '\t' << loc_[dtc2cart_[2]]*dirConv_[2] << '\t' << std::endl;
        double point = 0.0;
        for(int kk = fields_[0]->outGrid()->z()-1; kk >= 0; --kk )
        {
            for(int jj = fields_[0]->outGrid()->y()-1; jj >= 0; --jj)
            {
                for(int ii = 0; ii < fields_[0]->outGrid()->x(); ++ii)
                {
                    point = 0.0;
                    for(auto & field :fields_)
                        outputFunction_(&field->outGrid()->point(ii,jj,kk), &field->outGrid()->point(ii,jj,kk)+1, &point, convFactor_);
                    std::cout << "\t" << point;
                }
                std::cout << std::endl;
            }
        }
    }
}

parallelDetectorCOUTCplx::parallelDetectorCOUTCplx(int dtcNum,  std::vector<cplx_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, DTCCLASSTYPE classType, double timeInterval, DIRECTION firstComp, double a, double I0, double dt) :
    parallelDetectorBaseCplx(dtcNum, grid, SI, loc, sz, type, classType, timeInterval, firstComp, a, I0, dt),
    dtc2cart_(3, 0)
{
    if(firstComp == DIRECTION::X)
        dtc2cart_ = {{ 0, 1, 2 }};
    else if(firstComp == DIRECTION::Y )
        dtc2cart_ = {{ 1, 2, 0 }};
    else
        dtc2cart_ = {{ 2, 0, 1 }};
}

void parallelDetectorCOUTCplx::output(double t)
{
    fields_[0]->getField();
    if(fields_[0]->master() )
    {
        std::cout << t*tConv_ << "\t" << loc_[dtc2cart_[0]]*dirConv_[0] << "\t" << loc_[dtc2cart_[1]] * dirConv_[1] << '\t' << loc_[dtc2cart_[2]]*dirConv_[2] << std::endl;
        cplx point = 0.0;
        for(int kk = fields_[0]->outGrid()->z()-1; kk >= 0; --kk )
        {
            for(int jj = fields_[0]->outGrid()->y()-1; jj >= 0; --jj)
            {
                for(int ii = 0; ii < fields_[0]->outGrid()->x(); ++ii)
                {
                    point = 0.0;
                    for(auto & field :fields_)
                        outputFunction_(&field->outGrid()->point(ii,jj,kk), &field->outGrid()->point(ii,jj,kk)+1, &point, convFactor_);
                    std::cout << "\t" << point;
                }
                std::cout << std::endl;
            }
        }
    }
}