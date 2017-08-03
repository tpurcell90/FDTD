#ifndef FDTD_PARALLEL_FLUX
#define FDTD_PARALLEL_FLUX

#include <DTC/parallelStorageFreqDTC.hpp>


template <typename T> class parallelFluxDTC
{
public:

    struct masterImportDat
    {
        std::vector<int> addIndex_; //!< index of the first point on the surface
        std::vector< std::vector< std::array<int,9> > > szProcOffsetEj_; //!< A vector of the same size as number of offsets needed for spatial averaging storing an array storing information on the size, slave process rank, and location offset for each slave process for teh Ej detectors
        std::vector< std::vector< std::array<int,9> > > szProcOffsetEk_; //!< A vector of the same size as number of offsets needed for spatial averaging storing an array storing information on the size, slave process rank, and location offset for each slave process for teh Ek detectors
        std::vector< std::vector< std::array<int,9> > > szProcOffsetHj_; //!< A vector of the same size as number of offsets needed for spatial averaging storing an array storing information on the size, slave process rank, and location offset for each slave process for teh Hj detectors
        std::vector< std::vector< std::array<int,9> > > szProcOffsetHk_; //!< A vector of the same size as number of offsets needed for spatial averaging storing an array storing information on the size, slave process rank, and location offset for each slave process for teh HK detectors
        // Needed to send the struct from the slave process to the master
        template <typename Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar & addIndex_;
            ar & szProcOffsetEj_;
            ar & szProcOffsetEk_;
            ar & szProcOffsetHj_;
            ar & szProcOffsetHk_;
        }
    };
    struct FieldInputParamsFlux
    {
        double weight_; //!< weight factor for each surface
        int sz_;       //!< length of surface
        std::vector<int> addIndex_; //!< the lower, left, back corner of the flux region on any given process (if region goes across multiple processses)
        std::vector<std::shared_ptr<parallelStorageFreqDTC<T>>> Ej_dtc_; //!< shared_ptr to the Ej field where i is the direction of flux detection; vector to handle offset for averaging
        std::vector<std::shared_ptr<parallelStorageFreqDTC<T>>> Ek_dtc_; //!< shared_ptr to the Ek field where i is the direction of flux detection; vector to handle offset for averaging
        std::vector<std::shared_ptr<parallelStorageFreqDTC<T>>> Hj_dtc_; //!< shared_ptr to the Hj field where i is the direction of flux detection; vector to handle offset for averaging
        std::vector<std::shared_ptr<parallelStorageFreqDTC<T>>> Hk_dtc_; //!< shared_ptr to the Hk field where i is the direction of flux detection; vector to handle offset for averaging
        std::vector<std::shared_ptr<masterImportDat>> combineEjFields_; //!< Structs needed to combine the field info from each process' Ej_dtc_ (size of this vector should be the number of processes that has info for the Ej field)
        std::vector<std::shared_ptr<masterImportDat>> combineEkFields_; //!< Structs needed to combine the field info from each process' Ek_dtc_ (size of this vector should be the number of processes that has info for the Ek field)
        std::vector<std::shared_ptr<masterImportDat>> combineHjFields_; //!< Structs needed to combine the field info from each process' Hj_dtc_ (size of this vector should be the number of processes that has info for the Hj field)
        std::vector<std::shared_ptr<masterImportDat>> combineHkFields_; //!< Structs needed to combine the field info from each process' Hk_dtc_ (size of this vector should be the number of processes that has info for the Hk field)
    };
protected:
    typedef std::shared_ptr<parallelGrid<T>> pgrid_ptr;
    mpiInterface gridComm_; //!< mpiInterface for all communication of the flux detector
    bool save_; //!< if save is true it will save the final field data in a binary file

    int outProc_; //!< process rank of the collector/outputting process (process holding lower, left, back corner of the flux region)
    int t_step_; //!< the number of times the fields have been inputted
    int timeInt_; //!< the number of time steps in the main calculation for each field input step
    int nfreq_; //!< number of frequencies that the fulx will be calculated

    double dt_; //!< the time step size for the flux detector (dt$_{\text{main_calc} * timeInt)
    double fluxConv_; //!< Convesion factor for the final flux calculation (unit conversions from FDTD standard)
    double fluxIncdConv_; //!< Conversion factor for the incident flux calculation (unit converstions from FDTD standard)
    double freqConv_; //!< Conversion factor for frequency values (unit conversion from FDTD standard)
    double dOmg_; //!< step size of the frequency list (-1 if region defined by wavelength step)
    double dLam_; //!< step size of the frequency list in terms of wavelength (-1 if region defined by frequency step)

    cplx incd_prefactor_Ej_; //!< prefactor for the Ej field for incd flux calculation
    cplx incd_prefactor_Ek_; //!< prefactor for the Ek field for incd flux calculation
    cplx incd_prefactor_Hj_; //!< prefactor for the Hj field for incd flux calculation
    cplx incd_prefactor_Hk_; //!< prefactor for the Hk field for incd flux calculation

    std::vector<FieldInputParamsFlux> fInParam_; //!< struct for inputting fields

    std::array<int,3> loc_; //!< location of the lower, left, back corner of the flux region in grid point
    std::array<int,3> sz_; //!< size of the flux region in grid points

    std::array<double,3> d_; //!< step size of all grids inputted

    std::vector<cplx> fftFact_; //!< frequency factor for fourier transform for each frequency (find $\exp(i\omega t)$)

    std::string fname_; //!< file name for the flux detector
    std::string incd_fields_file_; //!< file name of the incident field files for the region if necessary

    std::vector<cplx_grid_ptr> Ej_freq_; //!< grid storing the Ej fields at points defined to be the same points as Hi where i is the direction of flux
    std::vector<cplx_grid_ptr> Ek_freq_; //!< grid storing the Ek fields at points defined to be the same points as Hi where i is the direction of flux
    std::vector<cplx_grid_ptr> Hj_freq_; //!< grid storing the Hj fields at points defined to be the same points as Hi where i is the direction of flux
    std::vector<cplx_grid_ptr> Hk_freq_; //!< grid storing the Hk fields at points defined to be the same points as Hi where i is the direction of flux
    std::vector<double> freqList_; //!< list of all frequencies to be monitered

    std::function<cplx(cplx)> getIncdField_; //!< if real fields return (real, 0.0) else return cplx field
public:

    /**
     * @brief Constructs a parallel flux detector with constant frequency step
     *
     * @param[in]  gridComm      mpiCommunicator used for the current calculation
     * @param[in]  name          file name for the outputted flux
     * @param[in]  weight        additional factor to adjust the final value of the detector (typically +/- 1 and set to the unit normal of the region)
     * @param[in]  Ex            grid pointer to the Ex grid in the main calculation
     * @param[in]  Ey            grid pointer to the Ey grid in the main calculation
     * @param[in]  Ez            grid pointer to the Ez grid in the main calculation
     * @param[in]  Hx            grid pointer to the Hx grid in the main calculation
     * @param[in]  Hy            grid pointer to the Hy grid in the main calculation
     * @param[in]  Hz            grid pointer to the Hz grid in the main calculation
     * @param[in]  location      of the lower, left and back corner of the flux region in number grid points
     * @param[in]  sz            size of the flux region in number of grid points
     * @param[in]  cross_sec     true if output should give optical cross-sections
     * @param[in]  save          true if the fourier transformed fields need to be outputed to files
     * @param[in]  load          true if previously saved fourier transformed fields need to be inputted to the calculation
     * @param[in]  timeInt       number of time steps of the main calcualtion per filed input step
     * @param[in]  nfreq         number of frequencies to be monitered
     * @param[in]  fwidth        width of the frequency region
     * @param[in]  fcen          center of teh frequency region
     * @param[in]  propDir       direction of propagation of the incident light (from TFSF surface)
     * @param[in]  d             vector describing the step size of the grids
     * @param[in]  dt            time step of the main calculation
     * @param[in]  incd_file     file name storing the inicident fields for the flux region
     * @param[in]  SI            true if outputting in SI units
     * @param[in]  I0            unit current of the calculation
     * @param[in]  a             unit distance of the calculation
     */
    parallelFluxDTC(mpiInterface gridComm, std::string name, double weight, pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, std::array<int,3> loc, std::array<int,3> sz, bool cross_sec, bool save, bool load, int timeInt, int nfreq, double fwidth, double fcen, DIRECTION propDir, std::array<double,3> d, double dt, double theta, double phi, double psi, double alpha, std::string incd_file, bool SI, double I0, double a) :
        gridComm_( gridComm ),
        save_(save),
        t_step_(0),
        timeInt_(timeInt),
        nfreq_(nfreq),
        dt_(dt*static_cast<double>(timeInt)),
        d_(d),
        fluxConv_(weight),
        fluxIncdConv_(1.0),
        freqConv_(1.0),
        dOmg_(fwidth / static_cast<double>(nfreq-1) * 2 * M_PI),
        dLam_(-1.0),
        incd_prefactor_Ej_(0.0),
        incd_prefactor_Ek_(0.0),
        incd_prefactor_Hj_(0.0),
        incd_prefactor_Hk_(0.0),
        loc_(loc),
        sz_(sz),
        fftFact_(nfreq_,0.0),
        fname_(name),
        incd_fields_file_(incd_file),
        freqList_(nfreq_, 0.0)
    {

        // set outProc to be in the process with the lower left back corner of the region
        if(Ex)
            outProc_ = Ex->getLocsProc( loc_[0], loc_[1], loc_[2] );
        else
            outProc_ = Hx->getLocsProc( loc_[0], loc_[1], loc_[2] );
        // if outputting in SI do the unit vconversions
        if(SI)
        {
            freqConv_ = SPEED_OF_LIGHT / a;
            fluxConv_ = I0/a * I0/(EPS0*a*SPEED_OF_LIGHT);
            fluxIncdConv_ = I0/a * I0/(EPS0*a*SPEED_OF_LIGHT);
        }
        // if not calculating cross-sections the iciident flux needs to be multiplied by the area of the detector
        if(!cross_sec)
        {
            if(propDir == DIRECTION::Y)
                fluxIncdConv_ *= std::real( simps2DNorm( sz_[0], sz_[2], d_[0], d_[2] ) );
            else if(propDir == DIRECTION::X)
                fluxIncdConv_ *= std::real( simps2DNorm( sz_[1], sz_[2], d_[1], d_[2] ) );
            else if(propDir == DIRECTION::Z)
                fluxIncdConv_ *= std::real( simps2DNorm( sz_[0], sz_[1], d_[0], d_[1] ) );
            else
                throw std::logic_error("Need a direction of propgation to correctly normalize flux" );
        }

        if(propDir == DIRECTION::X)
        {
            // Coordinates corresponding to how it is in terms of adding from main girds to detector grids
            if(Hz && Ez)
            {
                incd_prefactor_Ej_ = weight * -1.0 * ( cos(psi)*cos(phi) + sin(psi)*cos(theta)*sin(phi) );
                incd_prefactor_Ek_ = weight * sin(psi)*sin(theta) * std::exp( cplx(0, alpha) );
                incd_prefactor_Hj_ = weight * -1.0 * ( sin(psi)*cos(phi) - cos(psi)*cos(theta)*sin(phi) ) * std::exp( cplx(0, alpha) );
                incd_prefactor_Hk_ = weight * -1.0 * cos(psi)*sin(theta);
            }
            else if(Hz)
            {
                incd_prefactor_Ej_ = weight * std::cos(phi);
                incd_prefactor_Hk_ = weight * 1.0;
            }
            else
            {
                incd_prefactor_Hj_ = weight * -1.0*std::cos(phi);
                incd_prefactor_Ek_ = weight * 1.0;
            }
        }
        else if(propDir == DIRECTION::Y)
        {
            if(Hz && Ez)
            {
                incd_prefactor_Ej_ = weight * sin(psi)*sin(theta);
                incd_prefactor_Ek_ = weight * ( cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi) ) * std::exp( cplx(0, alpha) );
                incd_prefactor_Hj_ = weight * -1.0 * cos(psi)*sin(theta) * std::exp( cplx(0, alpha) );
                incd_prefactor_Hk_ = weight * sin(psi)*sin(phi) + cos(psi)*cos(theta)*cos(phi);
            }
            else if(Hz)
            {
                incd_prefactor_Ek_ = weight * -1.0*std::sin(phi);
                incd_prefactor_Hj_ = weight * 1.0;
            }
            else
            {
                incd_prefactor_Hk_ = weight * std::sin(phi);
                incd_prefactor_Ej_ = weight * 1.0;
            }
        }
        else
        {
            incd_prefactor_Ej_ = weight * cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi);
            incd_prefactor_Ek_ = weight * -1.0 * ( cos(psi)*cos(phi) + sin(psi)*cos(theta)*sin(phi) ) * std::exp( cplx(0, alpha) );
            incd_prefactor_Hj_ = weight * ( sin(psi)*sin(phi) + cos(psi)*cos(theta)*cos(phi) ) * std::exp( cplx(0, alpha) );
            incd_prefactor_Hk_ = weight * -1.0 * ( sin(psi)*cos(phi) - cos(psi)*cos(theta)*sin(phi) );
        }

        freqConv_ /= (M_PI*2.0);
        for(int ii = 0; ii < nfreq; ii ++)
            freqList_[ii] = (fcen - fwidth/2.0) * 2.0 * M_PI + ii*dOmg_;
        if(load)
            loadFields(-1.0);
    }

    /**
     * @brief Constructs a parallel flux detector with constant frequency step
     *
     * @param[in]  gridComm      mpiCommunicator used for the current calculation
     * @param[in]  name          file name for the outputted flux
     * @param[in]  weight        additional factor to adjust the final value of the detector (typically +/- 1 and set to the unit normal of the region)
     * @param[in]  Ex            grid pointer to the Ex grid in the main calculation
     * @param[in]  Ey            grid pointer to the Ey grid in the main calculation
     * @param[in]  Ez            grid pointer to the Ez grid in the main calculation
     * @param[in]  Hx            grid pointer to the Hx grid in the main calculation
     * @param[in]  Hy            grid pointer to the Hy grid in the main calculation
     * @param[in]  Hz            grid pointer to the Hz grid in the main calculation
     * @param[in]  location      of the lower, left and back corner of the flux region in number grid points
     * @param[in]  sz            size of the flux region in number of grid points
     * @param[in]  cross_sec     true if output should give optical cross-sections
     * @param[in]  save          true if the fourier transformed fields need to be outputed to files
     * @param[in]  load          true if previously saved fourier transformed fields need to be inputted to the calculation
     * @param[in]  timeInt       number of time steps of the main calcualtion per filed input step
     * @param[in]  lamL          Left end point of the wavelength region to be monitered
     * @param[in]  lamR          Right end point of the wavelength region to be monitered by the region
     * @param[in]  nLam          number of wavelengths to be monitered
     * @param[in]  propDir       direction of propagation of the incident light (from TFSF surface)
     * @param[in]  d             vector describing the step size of the grids
     * @param[in]  dt            time step of the main calculation
     * @param[in]  incd_file     file name storing the inicident fields for the flux region
     * @param[in]  SI            true if outputting in SI units
     * @param[in]  I0            unit current of the calculation
     * @param[in]  a             unit distance of the calculation
     */
    parallelFluxDTC(mpiInterface gridComm, std::string name, double weight, pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, std::array<int,3> loc, std::array<int,3> sz, bool cross_sec, bool save, bool load, int timeInt, double lamL, double lamR, int nLam, DIRECTION propDir, std::array<double,3> d, double dt, double theta, double phi, double psi, double alpha, std::string incd_file, bool SI, double I0, double a) :
        gridComm_( gridComm ),
        save_(save),
        t_step_(0),
        timeInt_(timeInt),
        nfreq_(nLam),
        dt_(dt*static_cast<double>(timeInt)),
        d_(d),
        fluxConv_(weight),
        fluxIncdConv_(1.0),
        freqConv_(1.0),
        dOmg_(-1.0),
        dLam_((lamR - lamL) / static_cast<double>(nLam-1)),
        loc_(loc),
        sz_(sz),
        fftFact_(nfreq_,0.0),
        fname_(name),
        incd_fields_file_(incd_file),
        freqList_(nfreq_, 0.0)
    {
        if(Ex)
            outProc_ = Ex->getLocsProc( loc_[0], loc_[1], loc_[2] );
        else
            outProc_ = Hx->getLocsProc( loc_[0], loc_[1], loc_[2] );
        // if SI do teh unit conversions
        if(SI)
        {
            freqConv_ = SPEED_OF_LIGHT / a;
            fluxConv_ = I0/a * I0/(EPS0*a*SPEED_OF_LIGHT);
            fluxIncdConv_ = I0/a * I0/(EPS0*a*SPEED_OF_LIGHT);
        }
        if(!cross_sec)
        {
            if(propDir == DIRECTION::Y)
                fluxIncdConv_ *= std::real( simps2DNorm( sz_[0], sz_[2], d_[0], d_[2] ) );
            else if(propDir == DIRECTION::X)
                fluxIncdConv_ *= std::real( simps2DNorm( sz_[1], sz_[2], d_[1], d_[2] ) );
            else if(propDir == DIRECTION::Z)
                fluxIncdConv_ *= std::real( simps2DNorm( sz_[0], sz_[1], d_[0], d_[1] ) );
            else
                throw std::logic_error("Need a direction of propgation to correctly normalize flux" );
        }
        freqConv_ /= (M_PI*2.0);

        if(propDir == DIRECTION::X)
        {
            // Coordinates corresponding to how it is in terms of adding from main girds to detector grids
            if(Hz && Ez)
            {
                incd_prefactor_Ej_ = -1.0 * ( cos(psi)*cos(phi) + sin(psi)*cos(theta)*sin(phi) );
                incd_prefactor_Ek_ = sin(psi)*sin(theta) * std::exp( cplx(0, alpha) );

                incd_prefactor_Hj_ = -1.0 * ( sin(psi)*cos(phi) - cos(psi)*cos(theta)*sin(phi) ) * std::exp( cplx(0, alpha) );
                incd_prefactor_Hk_ = -1.0 * cos(psi)*sin(theta);
            }
            else if(Hz)
            {
                incd_prefactor_Ej_ = std::cos(phi);
                incd_prefactor_Hj_ = std::cos(phi);
            }
            else
            {
                incd_prefactor_Ej_ = -1.0*std::cos(phi);
                incd_prefactor_Hj_ = -1.0*std::cos(phi);
            }
        }
        else if(propDir == DIRECTION::Y)
        {
            if(Hz && Ez)
            {
                incd_prefactor_Ej_ = sin(psi)*sin(theta);
                incd_prefactor_Ek_ = ( cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi) ) * std::exp( cplx(0, alpha) );

                incd_prefactor_Hj_ = -1.0 * cos(psi)*sin(theta) * std::exp( cplx(0, alpha) );
                incd_prefactor_Hk_ = sin(psi)*sin(phi) + cos(psi)*cos(theta)*cos(phi);
            }
            else if(Hz)
            {
                incd_prefactor_Ej_ = -1.0*std::sin(phi);
                incd_prefactor_Hj_ = -1.0*std::sin(phi);
            }
            else
            {
                incd_prefactor_Ej_ = std::sin(phi);
                incd_prefactor_Hj_ = std::sin(phi);
            }
        }
        else
        {
            incd_prefactor_Ej_ = cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi);
            incd_prefactor_Ek_ = -1.0 * ( cos(psi)*cos(phi) + sin(psi)*cos(theta)*sin(phi) ) * std::exp( cplx(0, alpha) );

            incd_prefactor_Hj_ = ( sin(psi)*sin(phi) + cos(psi)*cos(theta)*cos(phi) ) * std::exp( cplx(0, alpha) );
            incd_prefactor_Hk_ = -1.0 * ( sin(psi)*cos(phi) - cos(psi)*cos(theta)*sin(phi) );
        }
        for(int ii = 0; ii < nLam; ii ++)
            freqList_[ii] = 1.0 / (lamL + ii*dLam_) * 2.0 * M_PI;

        if(load)
            loadFields(-1.0);
    }

    std::vector< std::vector< std::array<int,9> > > constructSzProcOffsetLists( std::vector<std::shared_ptr<parallelStorageFreqDTC<T>>> dtcArr, std::vector<int> addIndex, pgrid_ptr act_field, pgrid_ptr trans_field, int corJK, int cor, int transCor1 )
    {
        std::vector< std::vector< std::array<int,9> > > szProcOffset;
        for(auto& dtc : dtcArr)
        {
            if( act_field && trans_field )
            {
                szProcOffset.push_back( std::vector<std::array<int, 9> >(2) );
                // If it does not need to send anything send a blank
                if(dtc->outGrid() )
                {
                    szProcOffset.back()[0] = {{ gridComm_.rank(), dtc->outGrid()->y(), dtc->outGrid()->z(), 0, 0, 0, 0, 0, 0 }} ;
                    szProcOffset.back()[1] = {{ gridComm_.rank(), dtc->outGrid()->y(), dtc->outGrid()->z(), 0, 0, 0, 0, 0, 0 }} ;
                }
                else
                {
                    szProcOffset.back()[0] = {{ -1, 0, 0, 0, 0, 0, 0, 0, 0 }} ;
                    szProcOffset.back()[1] = {{ -1, 0, 0, 0, 0, 0, 0, 0, 0 }} ;
                }
                if(cor == corJK)
                {
                    if(addIndex[1] == 0)
                    {
                        szProcOffset.back()[1][6] = 1;
                        szProcOffset.back()[1][4] = 1;
                    }
                    else
                    {
                        szProcOffset.back()[0][8] = 1;
                    }
                    // if( (dtc->loc()[corJK] >= act_field->procLoc()[corJK]) && (dtc->loc()[corJK] + dtc->sz()[corJK]-1 < act_field->procLoc()[corJK] + act_field->ln_vec()[corJK]) )
                    if( (act_field->procLoc()[corJK] <= dtc->loc()[corJK] + dtc->sz()[corJK]-1 ) && (dtc->loc()[corJK] + dtc->sz()[corJK]-1 < act_field->procLoc()[corJK] + act_field->ln_vec()[corJK] - 2) )
                        szProcOffset.back()[0][4] = 1;

                }
                else if(transCor1 == corJK)
                {
                    if(addIndex[0] == 0)
                    {
                        szProcOffset.back()[1][5] = 1;
                        szProcOffset.back()[1][3] = 1;
                    }
                    else
                    {
                        szProcOffset.back()[0][7] = 1;
                    }
                    // if( (dtc->loc()[corJK] >= act_field->procLoc()[corJK]) && (dtc->loc()[corJK] + dtc->sz()[corJK]-1 < act_field->procLoc()[corJK] + act_field->ln_vec()[corJK]) )
                    if( (act_field->procLoc()[corJK] <= dtc->loc()[corJK] + dtc->sz()[corJK]-1 ) && (dtc->loc()[corJK] + dtc->sz()[corJK]-1 < act_field->procLoc()[corJK] + act_field->ln_vec()[corJK] - 2) )
                        szProcOffset.back()[0][3] = 1;
                }
            }
            else
            {
                szProcOffset.push_back( std::vector<std::array<int, 9> >(1) );
                // If it does not need to send anything send a blank
                if(dtc->outGrid() )
                {
                    szProcOffset.back()[0] = {{ gridComm_.rank(), dtc->outGrid()->y(), dtc->outGrid()->z(), 0, 0, 0, 0, 0, 0 }} ;
                    szProcOffset.back()[1] = {{ gridComm_.rank(), dtc->outGrid()->y(), dtc->outGrid()->z(), 0, 0, 0, 0, 0, 0 }} ;
                }
                else
                {
                    szProcOffset.back()[0] = {{ -1, 0, 0, 0, 0, 0, 0, 0, 0 }} ;
                    szProcOffset.back()[1] = {{ -1, 0, 0, 0, 0, 0, 0, 0, 0 }} ;
                }
            }

        }
        return szProcOffset;
    }

    /**
     * @brief constrcut a parameter strcut for the slave and master processes
     * @details [long description]
     *
     * @param Ex grid pointer to the Ex grid in the main calculation
     * @param Ey grid pointer to the Ey grid in the main calculation
     * @param Ez grid pointer to the Ez grid in the main calculation
     * @param Hx grid pointer to the Hx grid in the main calculation
     * @param Hy grid pointer to the Hy grid in the main calculation
     * @param Hz grid pointer to the Hz grid in the main calculation
     * @param dir direction of the flux calculation
     * @param pl true if surface is the up, front or right surface, or if it is a stand alone surface
     * @return struct for the slave process for the surface described by the flux detector parameters
     */
    virtual FieldInputParamsFlux makeParamIn(pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, DIRECTION dir, bool pl) = 0;

    /**
     * @return fname_
     */
    inline std::string &fname(){return fname_;}
    /**
     * @return refl_
     */
    inline bool save() {return save_;}
    /**
     * @return location of the flux dtector
     */
    inline std::array<int,3> loc() {return loc_;}
    /**
     * @return dize of flux detector
     */
    inline std::array<int,3> sz() {return sz_;}

    /**
     * @return the number of time steps of the main calculation for each field input step
     */
    inline int & timeInt(){return timeInt_;}

    /**
     * @brief take in the field information
     * @details takes in teh electromagnetic filed infomratino at each time step
     *
     * @param tt current time in the simulation
     */
    void fieldIn(double& tt)
    {
        std::transform(freqList_.begin(), freqList_.end(), fftFact_.begin(), [&tt](double freq){return std::exp(cplx(0.0,-1.0*tt*freq) ); } );
        for(int vv = 0; vv < fInParam_.size(); vv++)
        {
            for(auto& dtc : fInParam_[vv].Ej_dtc_)
                dtc->fieldIn( fftFact_.data() );

            for(auto& dtc : fInParam_[vv].Ek_dtc_)
                dtc->fieldIn( fftFact_.data() );

            for(auto& dtc : fInParam_[vv].Hj_dtc_)
                dtc->fieldIn( fftFact_.data() );

            for(auto& dtc : fInParam_[vv].Hk_dtc_)
                dtc->fieldIn( fftFact_.data() );
        }
        ++t_step_;
    }

    void combineField(std::string fieldID, cplx_grid_ptr freqField, std::vector<std::shared_ptr<parallelStorageFreqDTC<T>>> dtcArr, std::vector<std::shared_ptr<masterImportDat>> combFieldVec, std::vector<int> fInAddInd)
    {
        for(int dd = 0; dd < dtcArr.size(); ++dd)
        {
            if(outProc_ == gridComm_.rank())
            {
                // Receive data from other processes and place it into the freq fields
                for(auto& getFields : combFieldVec)
                {
                    std::vector< std::vector< std::array<int,9> > > szProcOffset;
                    if(fieldID == "EJ")
                        szProcOffset = getFields->szProcOffsetEj_;
                    if(fieldID == "EK")
                        szProcOffset = getFields->szProcOffsetEk_;
                    if(fieldID == "HJ")
                        szProcOffset = getFields->szProcOffsetHj_;
                    if(fieldID == "HK")
                        szProcOffset = getFields->szProcOffsetHk_;
                    // Copy data from detectors to freq fields
                    if(szProcOffset[dd][0][0] == gridComm_.rank() )
                    {
                        for(auto& szOffProc : szProcOffset[dd])
                        {
                            for(int jj = 0; jj < szOffProc[2]-szOffProc[4]; ++jj)
                            {
                                for(int ii = 0; ii < szOffProc[1]-szOffProc[3]; ++ii)
                                {
                                    zaxpy_(nfreq_, 1.0/(static_cast<double>(dtcArr.size() * szProcOffset[dd].size() ) ), &dtcArr[dd]->outGrid()->point(0, ii+szOffProc[5], jj+szOffProc[6]), 1,   &freqField->point(0, ii + fInAddInd[0]+szOffProc[7], jj + fInAddInd[1]+szOffProc[8]), 1);
                                }
                            }
                        }
                    }
                    else if(szProcOffset[dd][0][0] != -1)
                    {
                        std::vector<cplx> temp_store(szProcOffset[dd][0][1] * szProcOffset[dd][0][2] * nfreq_,0.0);
                        gridComm_.recv(szProcOffset[dd][0][0], gridComm_.cantorTagGen(szProcOffset[dd][0][0], gridComm_.rank(), 3, 0), temp_store);
                        // Iterate over last index first since that 1 in 2D
                        for(auto& szOffProc : szProcOffset[dd])
                        {
                            for(int jj = 0; jj < szOffProc[2]-szOffProc[4]; ++jj)
                            {
                                for(int ii = 0; ii < szOffProc[1]-szOffProc[3]; ++ii)
                                {
                                    // Copy from vector that matches the way things are stored in the grid
                                    zaxpy_(nfreq_, 1.0/(static_cast<double>(dtcArr.size() * szProcOffset[dd].size() ) ), &temp_store[( (ii+szOffProc[5]) + (jj+szOffProc[6])*szOffProc[1])*nfreq_] , 1,   &freqField->point(0, ii + getFields->addIndex_[0]+szOffProc[7], jj + getFields->addIndex_[1]+szOffProc[8]), 1);
                                }
                            }
                        }
                    }
                }
            }
            else if(dtcArr[dd]->outGrid())
            {
                // Send data to master
                std::vector<cplx> sendVec(dtcArr[dd]->outGrid()->size(), 0.0);
                // Iterate over last index first since that 1 in 2D
                for(int jj = 0; jj < dtcArr[dd]->outGrid()->z(); ++jj)
                {
                    for(int ii = 0; ii < dtcArr[dd]->outGrid()->y(); ++ii)
                    {
                        // Copy into vector matching the way things are stored in the grid
                        zcopy_(nfreq_, &dtcArr[dd]->outGrid()->point(0,ii,jj), 1, &sendVec[(ii + jj*dtcArr[dd]->outGrid()->y()) * nfreq_], 1);
                    }
                }
                // Send the grid data to the output/collector process
                gridComm_.send(outProc_, gridComm_.cantorTagGen(gridComm_.rank(), outProc_, 3, 0), sendVec);
            }
        }
    }

    /**
     * @brief calculates the flux at the end of the calculation
     * @details uses the stored field infomration to fourier transform the fileds and cacluate the flux
     *
     * @param E_incd incident E field from TFSF surface
     * @param H_incd incident H field from the TFSF surface
     * @param off_incd An offset field for the incident flux so the E and H fields are calculated at the same point
     * @param TM true if 2D in the TM mode or 3D
    */
    void getFlux(std::vector<cplx> E_incd, std::vector<cplx> H_incd, std::vector<cplx> off_incd, bool TM)
    {
        for(int vv = 0; vv < fInParam_.size(); ++vv)
        {
            combineField("EJ", Ej_freq_[vv], fInParam_[vv].Ej_dtc_, fInParam_[vv].combineEjFields_, fInParam_[vv].addIndex_);
            combineField("EK", Ek_freq_[vv], fInParam_[vv].Ek_dtc_, fInParam_[vv].combineEkFields_, fInParam_[vv].addIndex_);
            combineField("HJ", Hj_freq_[vv], fInParam_[vv].Hj_dtc_, fInParam_[vv].combineHjFields_, fInParam_[vv].addIndex_);
            combineField("HK", Hk_freq_[vv], fInParam_[vv].Hk_dtc_, fInParam_[vv].combineHkFields_, fInParam_[vv].addIndex_);
        }
        if(outProc_ != gridComm_.rank())
            return;
        int nt = t_step_;
        std::ofstream f;
        f.open(fname_ + ".dat");
        cplx flux(0.0,0.0);
        cplx tempFlux(0.0,0.0);
        cplx flux_incd(0.0,0.0);
        // Storage for EjHk part of cross product
        std::vector<cplx_grid_ptr> flux_grid_ijk( fInParam_.size() );
        // Storage for EkHj part of cross product
        std::vector<cplx_grid_ptr> flux_grid_ikj( fInParam_.size() );

        // Construct the flux grids based on the sizes fo the compenent grids (either Ej/Hk exist, Ek/Hj exist, or both exist)
        for( int vv = 0; vv < flux_grid_ijk.size(); ++vv)
        {
            if(Ej_freq_[vv])
            {
                flux_grid_ijk[vv] = std::make_shared<Grid<cplx>>(std::array<int,3>( {{ Ej_freq_[vv]->z(), Ej_freq_[vv]->y(), 1 }} ) , std::array<double,3>( {{ Ej_freq_[vv]->dz(), Ej_freq_[vv]->dy(), 1 }} ) );
                // flux_grid_ikj[vv] = std::make_shared<Grid<cplx>>(std::array<int,3>( {{ Hk_freq_[vv]->z(), Hk_freq_[vv]->y(), 1 }} ) , std::array<double,3>( {{ Hk_freq_[vv]->dz(), Hk_freq_[vv]->dy(), 1 }} ) );
            }
            if(Ek_freq_[vv])
            {
                // flux_grid_ijk[vv] = std::make_shared<Grid<cplx>>(std::array<int,3>( {{ Ek_freq_[vv]->z(), Ek_freq_[vv]->y(), 1 }} ) , std::array<double,3>( {{ Ek_freq_[vv]->dz(), Ek_freq_[vv]->dy(), 1 }} ) );
                flux_grid_ikj[vv] = std::make_shared<Grid<cplx>>(std::array<int,3>( {{ Hj_freq_[vv]->z(), Hj_freq_[vv]->y(), 1 }} ) , std::array<double,3>( {{ Hj_freq_[vv]->dz(), Hj_freq_[vv]->dy(), 1 }} ) );
            }
        }
        // Calculate flux for every freqency
        for(int ff = 0; ff < nfreq_; ff++)
        {
            flux = 0.0;
            // Incdident flux from TFSF is 1D so no need to use E(H)j/k
            cplx Ej_inc(0.0,0.0), Hj_inc(0.0,0.0), Ek_inc(0.0,0.0), Hk_inc(0.0,0.0), offj_inc(0.0,0.0), offk_inc(0.0,0.0);
            // Calculate fourier transform of the incident field at current freq
            // TE mode has different offsets E fields not H fields
            if(TM)
            {
                for(int tt = 0; tt < E_incd.size(); ++tt)
                {
                    // Done to get correct incident feld information (correct phase factor for k fields)
                    Ej_inc   += getIncdField_( incd_prefactor_Ej_ *   E_incd[tt] ) * std::exp(cplx(0.0,freqList_.at(ff)*tt*(dt_ / static_cast<double>(timeInt_) ) ) );
                    Ek_inc   += getIncdField_( incd_prefactor_Ek_ *   E_incd[tt] ) * std::exp(cplx(0.0,freqList_.at(ff)*tt*(dt_ / static_cast<double>(timeInt_) ) ) );
                    Hj_inc   += getIncdField_( incd_prefactor_Hj_ *   H_incd[tt] ) * std::exp(cplx(0.0,freqList_.at(ff)*tt*(dt_ / static_cast<double>(timeInt_) ) ) );
                    Hk_inc   += getIncdField_( incd_prefactor_Hk_ *   H_incd[tt] ) * std::exp(cplx(0.0,freqList_.at(ff)*tt*(dt_ / static_cast<double>(timeInt_) ) ) );
                    offj_inc += getIncdField_( incd_prefactor_Hj_ * off_incd[tt] ) * std::exp(cplx(0.0,freqList_.at(ff)*tt*(dt_ / static_cast<double>(timeInt_) ) ) );
                    offk_inc += getIncdField_( incd_prefactor_Hk_ * off_incd[tt] ) * std::exp(cplx(0.0,freqList_.at(ff)*tt*(dt_ / static_cast<double>(timeInt_) ) ) );
                }

                flux_incd  = fluxIncdConv_ * ( Ej_inc * std::conj(Hk_inc) + Ej_inc * std::conj(offk_inc) ) / (2.0*std::pow(static_cast<double>(nt * timeInt_), 2.0) );
                flux_incd -= fluxIncdConv_ * ( Ek_inc * std::conj(Hj_inc) + Ek_inc * std::conj(offj_inc) ) / (2.0*std::pow(static_cast<double>(nt * timeInt_), 2.0) );
            }
            else
            {
                for(int tt = 0; tt < E_incd.size(); ++tt)
                {
                    // Done to get correct incident feld information (correct phase factor for k fields)
                    Ej_inc   += getIncdField_( incd_prefactor_Ej_ *   E_incd[tt] ) * std::exp(cplx(0.0,freqList_.at(ff)*tt*(dt_ / static_cast<double>(timeInt_) ) ) );
                    Ek_inc   += getIncdField_( incd_prefactor_Ek_ *   E_incd[tt] ) * std::exp(cplx(0.0,freqList_.at(ff)*tt*(dt_ / static_cast<double>(timeInt_) ) ) );
                    Hj_inc   += getIncdField_( incd_prefactor_Hj_ *   H_incd[tt] ) * std::exp(cplx(0.0,freqList_.at(ff)*tt*(dt_ / static_cast<double>(timeInt_) ) ) );
                    Hk_inc   += getIncdField_( incd_prefactor_Hk_ *   H_incd[tt] ) * std::exp(cplx(0.0,freqList_.at(ff)*tt*(dt_ / static_cast<double>(timeInt_) ) ) );
                    offj_inc += getIncdField_( incd_prefactor_Ej_ * off_incd[tt] ) * std::exp(cplx(0.0,freqList_.at(ff)*tt*(dt_ / static_cast<double>(timeInt_) ) ) );
                    offk_inc += getIncdField_( incd_prefactor_Ek_ * off_incd[tt] ) * std::exp(cplx(0.0,freqList_.at(ff)*tt*(dt_ / static_cast<double>(timeInt_) ) ) );
                }

                flux_incd  = fluxIncdConv_ * ( Ej_inc * std::conj(Hk_inc) + offj_inc * std::conj(Hk_inc) ) / (2.0*std::pow(static_cast<double>(nt * timeInt_), 2.0) );
                flux_incd -= fluxIncdConv_ * ( Ek_inc * std::conj(Hj_inc) + offk_inc * std::conj(Hj_inc) ) / (2.0*std::pow(static_cast<double>(nt * timeInt_), 2.0) );
            }

            for (int vv = 0; vv < fInParam_.size(); ++vv)
            {
                // temp fields to average over all points to get the value of each field in the center of the face
                cplx_grid_ptr tempE;
                cplx_grid_ptr tempH;
                // If Ej/Hk exist use them to build the temp fields, otherwise use Ek/Hj
                if(Ej_freq_[vv])
                {
                    tempE = std::make_shared<Grid<cplx>>(std::array<int,3>( {{ Ej_freq_[vv]->z(), Ej_freq_[vv]->y(), 1 }} ) , std::array<double,3>( {{ Ej_freq_[vv]->dz(), Ej_freq_[vv]->dy(), 1.0 }} ) );
                    tempH = std::make_shared<Grid<cplx>>(std::array<int,3>( {{ Hk_freq_[vv]->z(), Hk_freq_[vv]->y(), 1 }} ) , std::array<double,3>( {{ Hk_freq_[vv]->dz(), Hk_freq_[vv]->dy(), 1.0 }} ) );
                }
                else
                {
                    tempE = std::make_shared<Grid<cplx>>(std::array<int,3>( {{ Ek_freq_[vv]->z(), Ek_freq_[vv]->y(), 1 }} ) , std::array<double,3>( {{ Ek_freq_[vv]->dz(), Ek_freq_[vv]->dy(), 1.0 }} ) );
                    tempH = std::make_shared<Grid<cplx>>(std::array<int,3>( {{ Hj_freq_[vv]->z(), Hj_freq_[vv]->y(), 1 }} ) , std::array<double,3>( {{ Hj_freq_[vv]->dz(), Hj_freq_[vv]->dy(), 1.0 }} ) );
                }
                // Move freq field into continuous data structures
                if(Ej_freq_[vv])
                {
                    zcopy_(fInParam_[vv].sz_, &Ej_freq_[vv]->point(ff, 0, 0), nfreq_, tempE->data(), 1);
                    zcopy_(fInParam_[vv].sz_, &Hk_freq_[vv]->point(ff, 0, 0), nfreq_, tempH->data(), 1);
                }
                // Put Ej*conj(Hk) into flux_grid_ijk
                std::transform(tempE->data(), tempE->data()+tempE->size(), tempH->data(), flux_grid_ijk[vv]->data(), [&](cplx E, cplx H){return E*std::conj(H)/ std::pow(static_cast<double>(nt), 2.0); } );
                // fill the the E and H fields with zeros
                std::fill_n(tempE->data(), tempE->size(), 0.0);
                std::fill_n(tempH->data(), tempH->size(), 0.0);

                // Move freq field into continuous data structures
                if(Ek_freq_[vv])
                {
                    zcopy_(fInParam_[vv].sz_, &Ek_freq_[vv]->point(ff, 0, 0), nfreq_, tempE->data(), 1);
                    zcopy_(fInParam_[vv].sz_, &Hj_freq_[vv]->point(ff, 0, 0), nfreq_, tempH->data(), 1);
                }
                // Put Ek*conj(Hj) into flux_grid_ikj
                std::transform(tempE->data(), tempE->data()+tempE->size(), tempH->data(), flux_grid_ikj[vv]->data(), [&](cplx E, cplx H){return E*std::conj(H)/ std::pow(static_cast<double>(nt), 2.0); } );
                // Sum the flux_grid_ijk and flux_grid_ikj together
                zaxpy_(flux_grid_ijk[vv]->size(), -1.0, flux_grid_ikj[vv]->data(), 1, flux_grid_ijk[vv]->data(), 1);
            }
            // integrate flux over the spatial dimensions
            f << std::setw(18) << std::setprecision(15) << freqConv_ * freqList_[ff] << "\t" << std::setw(16) << std::setprecision(15) << std::abs(flux_incd) << "\t" << std::setw(16) << std::setprecision(15) << std::real(flux_incd) << "\t" << std::setw(16) << std::setprecision(15) << std::imag(flux_incd) << "\t";
            for(int vv = 0; vv < flux_grid_ijk.size(); vv++)
            {
                // Integrate over surface
                tempFlux = fluxConv_ * fInParam_[vv].weight_ * simps2D(flux_grid_ijk[vv]);
                // tempFlux = fInParam_[vv].weight_ * (flux_grid[vv][3], d_[0]);
                flux +=  tempFlux;
                f << std::setw(16) << std::setprecision(15) << std::abs(tempFlux) << "\t" << std::setw(16) << std::setprecision(15) << std::real(tempFlux) << "\t" << std::setw(16) << std::setprecision(15) << std::imag(tempFlux) << "\t";
            }
            if(flux_grid_ijk.size() > 1)
            {
                f << std::setw(16) << std::setprecision(15) << std::abs(flux) << "\t" << std::setw(16) << std::setprecision(15) << std::real(flux) << "\t" << std::setw(16) << std::setprecision(15) << std::imag(flux) << std::endl;
            }
            else
            {
                f << std::endl;
            }
        }
        if(save_)
            saveFields();
        return;
    }

    /**
     * @brief calculate the simpsion's rule approximation to an integral in 1D
     *
     * @param integrand array or vector of data to be integrated over
     * @param n number of steps to integrate over
     * @param d step size in direction of integration
     * @return numeric approx to the integral
     */
    // cplx simps(cplx* integrand, int n, double d)
    // {
    //     if(n == 1)
    //         return integrand[0]*d;
    //     cplx result(0.0,0.0);
    //     // How to handel an even number of points (Simpsons only defined for odd number so average over two sets)
    //     if(n % 2 == 0)
    //     {
    //         for(int ii = 0; ii < (n-2)/2; ii ++)
    //         {
    //             result += d/6.0 * (integrand[ii*2] + 4.0*integrand[ii*2+1] + integrand[(ii+1)*2]);
    //         }

    //         for(int ii = 1; ii < (n)/2; ii ++)
    //         {
    //             result += d/6.0 * (integrand[ii*2-1] + 4.0*integrand[ii*2] + integrand[ii*2+1]);
    //         }

    //         result += d/4.0*(integrand[0]+integrand[1] + integrand[n-1] + integrand[n-2]);
    //     }
    //     else
    //     {
    //         for(int ii = 0; ii < (n-1)/2; ii ++)
    //         {
    //             result += d/3.0 * (integrand[ii*2] + 4.0*integrand[ii*2+1] + integrand[(ii+1)*2]);
    //         }
    //     }
    //     return result;
    // }

    cplx simps(cplx* integrand, int n, double d)
    {
        if(n == 1)
            return integrand[0]*d;
        cplx result(0.0,0.0);
        if(n % 2 == 0)
        {
            for(int ii = 0; ii < (n-2)/2; ii ++)
                result += d/6.0 * (integrand[ii*2] + 4.0*integrand[ii*2+1] + integrand[(ii+1)*2]);

            for(int ii = 1; ii < (n)/2; ii ++)
                result += d/6.0 * (integrand[ii*2-1] + 4.0*integrand[ii*2] + integrand[ii*2+1]);

            result += d/4.0*(integrand[0]+integrand[1] + integrand[n-1] + integrand[n-2]);
        }
        else
        {
            for(int ii = 0; ii < (n-1)/2; ii ++)
                result += d/3.0 * (integrand[ii*2] + 4.0*integrand[ii*2+1] + integrand[(ii+1)*2]);
        }

        return result;
    }

    /**
     * @brief calculates the surface integral over a 2D grid using simpsions rule approximation
     *
     * @param grid pointer to grid that the surface integral is working over
     * @return the final integral value
     */
    cplx simps2D(cplx_grid_ptr grid)
    {
        if(grid->y() > 1 && grid->x() > 1)
        {
            std::vector<cplx>result(grid->y(),0.0);
            for(int jj = 0; jj < result.size(); ++jj)
                result[jj] = simps(&grid->point(0,jj), grid->x(), grid->dx());
            return simps(result.data(), result.size(), grid->dy() );
        }
        else if(grid->x() > 1)
            return simps(&grid->point(0,0), grid->x(), grid->dx());
        else if(grid->y() > 1)
            return simps(&grid->point(0,0), grid->y(), grid->dy());
    }

    /**
     * @brief calculate the surface integral for normalization purposes
     * @details calculates the area of the surface
     *
     * @param nx number of points of the surface in the x direction
     * @param ny number of points of the surface in the y direction
     * @param dx step size in the x direction
     * @param dy step size in the y direction
     * @return Area of the surface
     */
    cplx simps2DNorm(int nx, int ny, double dx, double dy)
    {

        std::vector<cplx>result(ny,0.0);
        for(int jj = 0; jj < result.size(); ++jj)
            result[jj] = simps(std::vector<cplx>( nx, cplx(1.0,0.0) ).data(), nx, dx);
        return simps(result.data(), result.size(), dy );
    }

    /**
     * @brief Saves te fields into a file to be used for other flux calculations
     */
    void saveFields()
    {
        // std::ofstream f;
        // f.open(fname_ + "_fields.dat", std::ios::binary | std::ios::out);
        // f.write( reinterpret_cast<char *>( &nfreq_ ), sizeof(nfreq_) );
        // f.write( reinterpret_cast<char *>( &sz_[0] ), sizeof(sz_[0]) );
        // f.write( reinterpret_cast<char *>( &sz_[1] ), sizeof(sz_[1]) );
        // f.write( reinterpret_cast<char *>( &sz_[2] ), sizeof(sz_[2]) );

        // int Ej_sz = Ej_freq_.size();
        // int Ek_sz = Ek_freq_.size();
        // int Hj_sz = Hj_freq_.size();
        // int Hk_sz = Hk_freq_.size();

        // f.write( reinterpret_cast<char *>( &Ej_sz ), sizeof( Ej_sz ) );
        // f.write( reinterpret_cast<char *>( &Ek_sz ), sizeof( Ek_sz ) );
        // f.write( reinterpret_cast<char *>( &Hj_sz ), sizeof( Hj_sz ) );
        // f.write( reinterpret_cast<char *>( &Hk_sz ), sizeof( Hk_sz ) );
        // cplx val(0,0);
        // for(int vv = 0; vv < fInParam_.size(); ++vv)
        // {
        //     // not just done with complex<double> since you can't convert it to char*
        //     for( auto& field : Ej_freq_[vv])
        //         f.write(reinterpret_cast<char *>( &field->point(0,0,0) ), sizeof(val)*field->x()*field->y()*field->z() );
        //     for( auto& field : Ek_freq_[vv])
        //         f.write(reinterpret_cast<char *>( &field->point(0,0,0) ), sizeof(val)*field->x()*field->y()*field->z() );
        //     for( auto& field : Hj_freq_[vv])
        //         f.write(reinterpret_cast<char *>( &field->point(0,0,0) ), sizeof(val)*field->x()*field->y()*field->z() );
        //     for( auto& field : Hk_freq_[vv])
        //         f.write(reinterpret_cast<char *>( &field->point(0,0,0) ), sizeof(val)*field->x()*field->y()*field->z() );
        // }
        // f.close();
        // return;
    }
    /**
     * @brief      loads incident fields into the calculations, weights to get negative
     *
     * @param  weight  typically -1.0, but how much to weight the field information here
     */
    void loadFields(double weight)
    {
        // std::ifstream f;
        // f.open(incd_fields_file_, std::ios::binary | std::ios::in);
        // int nx, ny, nz, nfreq;
        // f.read(reinterpret_cast<char *>( &nfreq ), sizeof( nfreq ) );
        // f.read(reinterpret_cast<char *>( &nx ), sizeof( nx ) );
        // f.read(reinterpret_cast<char *>( &ny ), sizeof( ny ) );
        // f.read(reinterpret_cast<char *>( &nz ), sizeof( nz ) );

        // int Ej_sz, Ek_sz, Hj_sz, Hk_sz;
        // f.read( reinterpret_cast<char *>( &Ej_sz), sizeof( Ej_sz ) );
        // f.read( reinterpret_cast<char *>( &Ek_sz), sizeof( Ek_sz ) );
        // f.read( reinterpret_cast<char *>( &Hj_sz), sizeof( Hj_sz ) );
        // f.read( reinterpret_cast<char *>( &Hk_sz), sizeof( Hk_sz ) );

        // if(nfreq_ != nfreq || sz_[0] != nx || sz_[1] != ny || sz_[2] != nz )
        //     throw std::logic_error("Given incident fields do not match size and frequency numbers for the current calculations");

        // if( ( Ej_sz != Ej_freq_.size() ) || ( Ek_sz != Ek_freq_.size() ) || ( Hj_sz != Hj_freq_.size() ) || ( Hk_sz != Hk_freq_.size() ) )
        //     throw std::logic_error("The size of the field vectors of the saved fields are not the same as those in the current calculation");

        // cplx val(0,0);
        // for(int vv = 0; vv < fInParam_.size(); vv++)
        // {
        //     for( auto& field : Ej_freq_[vv])
        //     {
        //         f.read(reinterpret_cast<char *>( &field->point(0,0,0) ), sizeof(val)*field->x()*field->y()*field->z() );
        //         zscal_(field->size(), weight, field->data(), 1);
        //     }
        //     for( auto& field : Ek_freq_[vv])
        //     {
        //         f.read(reinterpret_cast<char *>( &field->point(0,0,0) ), sizeof(val)*field->x()*field->y()*field->z());
        //         zscal_(field->size(), weight, field->data(), 1);
        //     }
        //     for( auto& field : Hj_freq_[vv])
        //     {
        //         f.read(reinterpret_cast<char *>( &field->point(0,0,0) ), sizeof(val)*field->x()*field->y()*field->z());
        //         zscal_(field->size(), weight, field->data(), 1);
        //     }
        //     for( auto& field : Hk_freq_[vv])
        //     {
        //         f.read(reinterpret_cast<char *>( &field->point(0,0,0) ), sizeof(val)*field->x()*field->y()*field->z());
        //         zscal_(field->size(), weight, field->data(), 1);
        //     }
        // }
        // f.close();
        // Make it minus to subtract out
    }
};

class parallelFluxDTCReal : public parallelFluxDTC<double>
{
    // typedef real_pgrid_ptr pgrid_ptr;

public:
    parallelFluxDTCReal(mpiInterface gridComm, std::string name, double weight, pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, std::array<int,3> loc, std::array<int,3> sz, bool cross_sec, bool save, bool load, int timeInt,   int nfreq, double fwidth, double fcen, DIRECTION propDir, std::array<double,3> d, double dt, double theta, double phi, double psi, double alpha, std::string incd_file, bool SI, double I0, double a);
    parallelFluxDTCReal(mpiInterface gridComm, std::string name, double weight, pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, std::array<int,3> loc, std::array<int,3> sz, bool cross_sec, bool save, bool load, int timeInt, double lamL,   double lamR,    int nLam, DIRECTION propDir, std::array<double,3> d, double dt, double theta, double phi, double psi, double alpha, std::string incd_file, bool SI, double I0, double a);
    // void fieldIn(double& tt);
    FieldInputParamsFlux makeParamIn(pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, DIRECTION dir, bool pl);
};

/**
 * Complex version of FluxDTC see base class for more descriptions
 */
class parallelFluxDTCCplx : public parallelFluxDTC<cplx>
{
    // typedef std::shared_ptr<parallelGrid<cplx> pgrid_ptr;
public:
    parallelFluxDTCCplx(mpiInterface gridComm, std::string name, double weight, pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, std::array<int,3> loc, std::array<int,3> sz, bool cross_sec, bool save, bool load, int timeInt,   int nfreq, double fwidth, double fcen, DIRECTION propDir, std::array<double,3> d, double dt, double theta, double phi, double psi, double alpha, std::string incd_file, bool SI, double I0, double a);
    parallelFluxDTCCplx(mpiInterface gridComm, std::string name, double weight, pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, std::array<int,3> loc, std::array<int,3> sz, bool cross_sec, bool save, bool load, int timeInt, double lamL,   double lamR,    int nLam, DIRECTION propDir, std::array<double,3> d, double dt, double theta, double phi, double psi, double alpha, std::string incd_file, bool SI, double I0, double a);
    // void fieldIn(double& tt);
    FieldInputParamsFlux makeParamIn(pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, DIRECTION dir, bool pl);
};
#endif