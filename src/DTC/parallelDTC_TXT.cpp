#include <src/DTC/parallelDTC_TXT.hpp>

parallelDetectorTXTReal::parallelDetectorTXTReal(std::vector<real_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, double timeInterval, double a, double I0, double dt) :
    parallelDetectorBaseReal(grid, SI, loc, sz, type, timeInterval, a, I0, dt),
    outFile_(out_name)
{
    // Construct output file stream
    outFileStream_ = std::make_shared<std::ofstream>();
    if(fields_[0]->master())
        outFileStream_->open(outFile_);
}
void parallelDetectorTXTReal::output(double t)
{
    // Import fields from outside processes
    for(auto & field :fields_)
        field->getField();
    if(fields_[0]->master())
    {
        double point = 0.0;
        // output time/location
        *outFileStream_ << t*tConv_ << "\t" << realSpaceLoc_[0] << "\t" << realSpaceLoc_[1] << "\t" << realSpaceLoc_[2];// << "\n";
        for(int kk = fields_[0]->outGrid()->z()- 1; kk >= 0; --kk)
        {
            for(int jj = fields_[0]->outGrid()->y()-1; jj >= 0; --jj)
            {
                for(int ii = 0; ii < fields_[0]->outGrid()->x(); ++ii)
                {
                    // Calculate and output point for all fields
                    point = 0.0;
                    for(auto & field :fields_)
                    {
                        outputFunction_(&field->outGrid()->point(ii,jj,kk), &field->outGrid()->point(ii,jj,kk)+1, &point, convFactor_);
                    }
                    *outFileStream_ << "\t" <<  point;
                    // Format files properly
                }
            }
        }
        *outFileStream_ << '\n';
    }
}
parallelDetectorTXTCplx::parallelDetectorTXTCplx(std::vector<cplx_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, double timeInterval, double a, double I0, double dt) :
    parallelDetectorBaseCplx(grid, SI, loc, sz, type, timeInterval, a, I0, dt),
    outFile_(out_name)
{
    // Construct output file stream
    outFileStream_ = std::make_shared<std::ofstream>();
    if(fields_[0]->master())
        outFileStream_->open(outFile_);
}
void parallelDetectorTXTCplx::output(double t)
{
    // Import fields from outside processes
    for(auto & field :fields_)
        field->getField();
    if(fields_[0]->master())
    {
        cplx point;
        // output time/location
        *outFileStream_ << t*tConv_ << "\t" << realSpaceLoc_[0] << "\t" << realSpaceLoc_[1] << "\t" << realSpaceLoc_[2];// << "\n";
        for(int kk = fields_[0]->outGrid()->z()- 1; kk >= 0; --kk)
        {
            for(int jj = fields_[0]->outGrid()->y()-1; jj >= 0; --jj)
            {
                for(int ii = 0; ii < fields_[0]->outGrid()->x(); ++ii)
                {
                    // Calculate and output point for all fields
                    point = 0.0;
                    for(auto & field :fields_)
                        outputFunction_(&field->outGrid()->point(ii,jj,kk), &field->outGrid()->point(ii,jj,kk)+1, &point, convFactor_);
                    *outFileStream_ << "\t" <<  std::real(point);
                }
            }
            // Format files properly
        }
        *outFileStream_ << "\n";
    }
}
