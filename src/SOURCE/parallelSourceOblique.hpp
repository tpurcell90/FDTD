#ifndef FDTD_SOURCE_OBLIQUE
#define FDTD_SOURCE_OBLIQUE

#include <SOURCE/parallelSource.hpp>
/**
 * @brief A parallel soft source for FDTD fields. Assuming phi is not along a normal axis of the grid (0,90,270,360)
 *
 * @tparam T param for type of source, doulbe for real, complex<double> for complex field
 */
template <typename T> class parallelSourceObliqueBase : public parallelSourceBase<T>
{
public:
    public:
    struct PulseAddParams
    {
        T* loc_; //!< reference to gird point
        double t_off_; //!< factor to offset the time of the pulse
        double scalefact_; //!< factor to scale the pulse based on angle of incidence
    };
protected:
    using parallelSourceBase<T>::gridComm_; //!<  mpiInterface for the FDTD field
    using parallelSourceBase<T>::pulse_; //!< the pulse that the source is adding to the filed
    using parallelSourceBase<T>::grid_; //!< the grid that the source is adding the pulse to
    using parallelSourceBase<T>::dt_; //!< time step of the calculation
    using parallelSourceBase<T>::loc_; //!< location of the source's lower left corner
    using parallelSourceBase<T>::sz_; //!< location of the source's lower left corner

    using parallelSourceBase<T>::pulse; //!< accessor function for pulse_
    using parallelSourceBase<T>::grid; //!< accessor function for grid_
    using parallelSourceBase<T>::loc; //!< accessor function for loc_
    using parallelSourceBase<T>::sz; //!< accessor function for sz

    std::vector<PulseAddParams> updateSrcParams_; //!< pulse update parameters
    double phi_; //<! azimuthal angle of light propagation
    double theta_; //!< polar angle of light propagation
    POLARIZATION pol_; //!< Polarization of the EM field
public:
    /**
     * @brief Constructor for the parallel source
     *
     * @param[in]  gridComm    mpiInterface for the caclutlation
     * @param[in]  srcNum     index of the source in srcArr_
     * @param[in]  pulse      pulse of the calculation
     * @param[in]  grid       grid the source adds the pulse to
     * @param[in]  dt         time step of the calculation
     * @param[in]  loc        location of the lower left corner of source
     * @param[in]  sz         size of the soft source
     */
    parallelSourceObliqueBase(std::shared_ptr<mpiInterface> gridComm, std::vector<std::shared_ptr<PulseBase>> pulse, std::shared_ptr<parallelGrid<T>> grid, POLARIZATION pol, double dt, std::array<int,3> loc, std::array<int,3> sz, double phi, double theta) :
        parallelSourceBase<T>(gridComm, pulse, grid, dt, loc, sz),
        phi_( fmod(phi,360.0) * M_PI/180.0 ),
        theta_( fmod(theta,360.0) * M_PI/180.0 - M_PI / 2.0 ),
        pol_(pol)
    {

        if( std::abs(theta_) > M_PI )
            throw std::logic_error("The polar angle has to be between 0 and 180 degrees");
        if(grid_->local_z() == 1)
            theta_ = M_PI/2.0;
        genDatStruct();
        if(updateSrcParams_.size() > 0)
            pulse_ = pulse;
    }

    /**
     * @brief      Calculates the cotanget of x
     *
     * @param[in]  x     value to get cot of
     *
     * @return     cot(x)
     */
    inline double cot(double x) { return tan(M_PI_2 - x); }

    /**
     * @brief      generates the data structures necessary to add the pulse into the field
     *
     * @param[in]  srcNum  index of the source in the srcArr_
     */
    void genDatStruct()
    {
        // std::vector<std::array<int,3>> allLocs_;
        double phi = fmod( phi_ + M_PI/2.0, 2.0*M_PI );
        double theta;
        if(theta_ >= M_PI/2)
            theta = theta_ - M_PI / 2.0;
        else
            theta = theta_ + M_PI / 2.0;
        std::vector<double> normVec = {{ sin(theta_)*cos(phi_), sin(theta_)*sin(phi_), cos(theta_) }}; //!< defined using light angles since the plane should be perpendicular to dir of propagation
        dscal_(normVec.size(), 1.0/sqrt(pow(normVec[0], 2.0) + pow(normVec[1], 2.0) + pow(normVec[2], 2.0)), normVec.data(), 1 ); //!< normalize the normal vector
        double dx = grid_->dx();
        double dy = grid_->dy();
        double dz = 1.0;
        if(grid_->local_z() != 1)
            dz = grid_->dz();

        int xMin = loc_[0] - fabs( int( round(sz_[0] * cos(phi) / 2.0 ) ) ); //!< will check the effect of the polar angle inside the point check (szgths along sz_[0])
        xMin = std::min(xMin, int( fabs( round( xMin + sz_[1]*cos(phi)*sin(theta) ) ) ) ); //!< Does the projection of the polar tilt on the x/y plane affect the results?
        int yMin = loc_[1] - fabs( int( round(sz_[0] * sin(phi) / 2.0 ) ) ); //!< will check the effect of the polar angle inside the point check (szgths along sz_[0])
        yMin = std::min(yMin, int( fabs( round( yMin + sz_[1]*sin(phi)*sin(theta) ) ) ) ); //!< Does the projection of the polar tilt on the x/y plane affect the results?
        int zMin = 0;
        if(grid_->local_z() != 1)
             zMin = loc_[2] - fabs( int( round(sz_[1] * cos(theta) / 2.0 ) ) );

        int xMax = loc_[0] + fabs( int( round(sz_[0] * cos(phi) / 2.0 ) ) ); //!< will check the effect of the polar angle inside the point check (szgths along sz_[0])
        xMax = std::min(xMax, int( fabs( round( xMax + sz_[1]*cos(phi)*sin(theta) ) ) ) ); //!< Does the projection of the polar tilt on the x/y plane affect the results?
        int yMax = loc_[1] + fabs( int( round(sz_[0] * sin(phi) / 2.0 ) ) ); //!< will check the effect of the polar angle inside the point check (szgths along sz_[0])
        yMax = std::min(yMax, int( fabs( round( yMax + sz_[1]*sin(phi)*sin(theta) ) ) ) ); //!< Does the projection of the polar tilt on the x/y plane affect the results?
        int zMax = 0;
        if(grid_->local_z() != 1)
             zMax = loc_[2] + fabs( int( round(sz_[1] * cos(theta) / 2.0 ) ) );


        std::vector<double>off(3,0.0);
        if(pol_ == POLARIZATION::EX || pol_ == POLARIZATION::HY || pol_ == POLARIZATION::HZ)
            off[0] = 0.5;
        if(pol_ == POLARIZATION::HX || pol_ == POLARIZATION::EY || pol_ == POLARIZATION::HZ)
            off[1] = 0.5;
        if(grid_->local_z() != 1 && (pol_ == POLARIZATION::HX || pol_ == POLARIZATION::HY || pol_ == POLARIZATION::EZ) )
            off[2] = 0.5;

        for( int xx = xMin; xx <= xMax; ++xx)
        {
            for( int yy = yMin; yy <= yMax; ++yy)
            {
                for( int zz = zMin; zz <= zMax; ++zz)
                {
                    std::vector<double> pt = {{(xx-loc_[0] + off[0])*dx, (yy-loc_[1] + off[1])*dy, (zz-loc_[2] + off[2])*dz}};
                    if( fabs( ddot_(3, pt.data(), 1, normVec.data(), 1 ) ) <= dx/sqrt(2.0) )
                    {
                        std::vector<int> gridPt =  {{xx, yy, zz}};
                        if(grid_->getLocsProc_no_boundaries(gridPt[0], gridPt[1], gridPt[2]) == gridComm_->rank())
                        {
                            PulseAddParams srcEl;
                            if(grid_->local_z() == 1)
                                srcEl.loc_ = &grid_->point(gridPt[0] - grid_->procLoc()[0] + 1, gridPt[1] - grid_->procLoc()[1] + 1, 0 );
                            else
                                srcEl.loc_ = &grid_->point(gridPt[0] - grid_->procLoc()[0] + 1, gridPt[1] - grid_->procLoc()[1] + 1, gridPt[2] - grid_->procLoc()[2] + 1 );

                            if(pol_ == POLARIZATION::EX || pol_ == POLARIZATION::HX)
                                srcEl.scalefact_ = cos(phi) * sin(theta_);
                            else if(pol_ == POLARIZATION::EY || pol_ == POLARIZATION::HY)
                                srcEl.scalefact_ = sin(phi) * sin(theta_);
                            else
                                srcEl.scalefact_ = cos(theta_);

                            srcEl.t_off_ = ddot_(3, normVec.data(), 1, pt.data(), 1);

                            updateSrcParams_.push_back(srcEl);
                        }
                    }
                }
            }
        }
    }
    /**
     * @brief      adds the pulse to the grid
     *
     * @param[in]  t     current time
     */
    virtual void addPul(double t) = 0;
};

class parallelSourceObliqueReal : public parallelSourceObliqueBase<double>
{
public:
    /**
     * @brief Constructor for the parallel source
     *
     * @param[in]  gridComm       mpiInterface for the caclutlation
     * @param[in]  srcNum     index of the source in srcArr_
     * @param[in]  pulse      pulse of the calculation
     * @param[in]  grid       grid the source adds the pulse to
     * @param[in]  dt         time step of the calculation
     * @param[in]  loc        location of the lower left corner of source
     * @param[in]  sz         size of the soft source
     */
    parallelSourceObliqueReal(std::shared_ptr<mpiInterface> gridComm,  std::vector<std::shared_ptr<PulseBase>> pulse, real_pgrid_ptr grid, POLARIZATION pol, double dt, std::array<int,3> loc, std::array<int,3> sz, double phi, double theta);
    /**
     * @brief      adds the pulse to the grid
     *
     * @param[in]  t     current time
     */
    void addPul(double t);
};

class parallelSourceObliqueCplx : public parallelSourceObliqueBase<cplx>
{
public:
    /**
     * @brief Constructor for the parallel source
     *
     * @param[in]  gridComm       mpiInterface for the caclutlation
     * @param[in]  srcNum     index of the source in srcArr_
     * @param[in]  pulse      pulse of the calculation
     * @param[in]  grid       grid the source adds the pulse to
     * @param[in]  dt         time step of the calculation
     * @param[in]  loc        location of the lower left corner of source
     * @param[in]  sz         size of the soft source
     */
    parallelSourceObliqueCplx(std::shared_ptr<mpiInterface> gridComm,  std::vector<std::shared_ptr<PulseBase>> pulse, cplx_pgrid_ptr grid, POLARIZATION pol, double dt, std::array<int,3> loc, std::array<int,3> sz, double phi, double theta);
    /**
     * @brief      adds the pulse to the grid
     *
     * @param[in]  t     current time
     */
    void addPul(double t);
};

#endif