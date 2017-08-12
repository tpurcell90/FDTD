#ifndef PARALLEL_FDTD_FDTDFIELD
#define PARALLEL_FDTD_FDTDFIELD

#include <INPUTS/parallelInputs.hpp>
#include <DTC/parallelDTC_TXT.hpp>
#include <DTC/parallelDTC_COUT.hpp>
#include <DTC/parallelDTC_BIN.hpp>
#include <DTC/parallelDTC_FREQ.hpp>
#include <DTC/parallelFlux.hpp>
#include <DTC/toBitMap.hpp>
#include <SOURCE/parallelSourceNormal.hpp>
#include <SOURCE/parallelSourceOblique.hpp>
#include <SOURCE/parallelTFSF.hpp>
#include <UTIL/FDTD_up_eq.hpp>

/**
 * @brief The main FDTD propagator class
 * @details This class generates a propagator class that will update all electromagnetic fields
 */
template <typename T> class parallelFDTDFieldBase
{
protected:
    typedef std::shared_ptr<parallelGrid<T>> pgrid_ptr;
    typedef std::shared_ptr<parallelCPML<T>> pml_ptr;

    std::shared_ptr<mpiInterface> gridComm_; //!< mpi communicator for the propagator

    bool dielectricMatInPML_; //!< True if dielectric (constant or dispersive) material are in the PML's
    bool magMatInPML_; //!< True if magnetic (constant or dispersive) material are in the PML's

    int res_; //!< number of grid points per unit length
    int t_step_; //!< the number of time steps that happened

    int yExPBC_; //!< amount to extend the PBC for the Ex filed in the y direction
    int yEyPBC_; //!< amount to extend the PBC for the Ey filed in the y direction
    int yEzPBC_; //!< amount to extend the PBC for the Ez filed in the y direction

    int yHxPBC_; //!< amount to extend the PBC for the Hx filed in the y direction
    int yHyPBC_; //!< amount to extend the PBC for the Hy filed in the y direction
    int yHzPBC_; //!< amount to extend the PBC for the Hz filed in the y direction

    int pbcZMin_; //!< min value for z PBC
    int pbcZMax_; //!< min value for z PBC

    std::array<int,3> n_vec_; //!< the number of grid points in each direction
    std::array<int,3> ln_vec_; //!< the number of grid points in each direction for this process only
    std::array<int,3> pmlThickness_; //!< thickness of the PMLS in all directions
    std::array<double,3> d_; //!< the step size in all direction

    double dt_; //!< the time step  of the simulation
    double tcur_; //!< the current time of the simulation

    std::vector<std::shared_ptr<Obj>> objArr_; //!< vector containing all objects in the cell

    std::array<double,3> r_; //!< the vector descrbing the current location of the grid point for PBC's
    std::array<double,3> k_point_; //!< k-point vector for periodicity

    std::vector<std::shared_ptr<parallelDetectorBase<T> > > dtcArr_; //!< the vector of detectors in the cell
    std::vector<std::shared_ptr<parallelSourceBase<T> > > srcArr_; //!< the vector of all sources in the cell

    upLists upHx_; //!< the list of parameters used to update the Hx field containing : std::pair(std::array<int,8>( number of elements for the calculation, x start, y start, z start, x offset for j(y) spatial derivative, y offset for j(y) spatial derivative, z offset for j(y) spatial derivative, object array index ), std::array<double,2> ( scaling factor (assumed to be 1 right now since conductivity = 0.0), spatial derivative prefactor for j spatial derivative) )
    upLists upHy_; //!< the list of parameters used to update the Hy field containing : std::pair(std::array<int,8>( number of elements for the calculation, x start, y start, z start, x offset for j(z) spatial derivative, y offset for j(z) spatial derivative, z offset for j(z) spatial derivative, object array index ), std::array<double,2> ( scaling factor (assumed to be 1 right now since conductivity = 0.0), spatial derivative prefactor for j spatial derivative) )
    upLists upHz_; //!< the list of parameters used to update the Hz field containing : std::pair(std::array<int,8>( number of elements for the calculation, x start, y start, z start, x offset for j(x) spatial derivative, y offset for j(x) spatial derivative, z offset for j(x) spatial derivative, object array index ), std::array<double,2> ( scaling factor (assumed to be 1 right now since conductivity = 0.0), spatial derivative prefactor for j spatial derivative) )

    upLists upEx_; //!< the list of parameters used to update the Ex field containing : std::pair(std::array<int,8>( number of elements for the calculation, x start, y start, z start, x offset for j(y) spatial derivative, y offset for j(y) spatial derivative, z offset for j(y) spatial derivative, object array index ), std::array<double,2> ( scaling factor (assumed to be 1 right now since conductivity = 0.0), spatial derivative prefactor for j spatial derivative) )
    upLists upEy_; //!< the list of parameters used to update the Ey field containing : std::pair(std::array<int,8>( number of elements for the calculation, x start, y start, z start, x offset for j(z) spatial derivative, y offset for j(z) spatial derivative, z offset for j(z) spatial derivative, object array index ), std::array<double,2> ( scaling factor (assumed to be 1 right now since conductivity = 0.0), spatial derivative prefactor for j spatial derivative) )
    upLists upEz_; //!< the list of parameters used to update the Ez field containing : std::pair(std::array<int,8>( number of elements for the calculation, x start, y start, z start, x offset for j(x) spatial derivative, y offset for j(x) spatial derivative, z offset for j(x) spatial derivative, object array index ), std::array<double,2> ( scaling factor (assumed to be 1 right now since conductivity = 0.0), spatial derivative prefactor for j spatial derivative) )

    upLists upDx_; //!< the list of parameters used to update the Dx field and polarization fields containing : std::pair(std::array<int,8>( number of elements for the calculation, x start, y start, z start, x offset for j(y) spatial derivative, y offset for j(y) spatial derivative, z offset for j(y) spatial derivative, object array index ), std::array<double,2> ( scaling factor (assumed to be 1 right now since conductivity = 0.0), spatial derivative prefactor for j spatial derivative) )
    upLists upDy_; //!< the list of parameters used to update the Dy field and polarization fields containing : std::pair(std::array<int,8>( number of elements for the calculation, x start, y start, z start, x offset for j(z) spatial derivative, y offset for j(z) spatial derivative, z offset for j(z) spatial derivative, object array index ), std::array<double,2> ( scaling factor (assumed to be 1 right now since conductivity = 0.0), spatial derivative prefactor for j spatial derivative) )
    upLists upDz_; //!< the list of parameters used to update the Dz field and polarization fields containing : std::pair(std::array<int,8>( number of elements for the calculation, x start, y start, z start, x offset for j(x) spatial derivative, y offset for j(x) spatial derivative, z offset for j(x) spatial derivative, object array index ), std::array<double,2> ( scaling factor (assumed to be 1 right now since conductivity = 0.0), spatial derivative prefactor for j spatial derivative) )

    upLists upBx_; //!< the list of parameters used to update the Bx field and magnetization fields containing : std::pair(std::array<int,8>( number of elements for the calculation, x start, y start, z start, x offset for j(y) spatial derivative, y offset for j(y) spatial derivative, z offset for j(y) spatial derivative, object array index ), std::array<double,2> ( scaling factor (assumed to be 1 right now since conductivity = 0.0), spatial derivative prefactor for j spatial derivative) )
    upLists upBy_; //!< the list of parameters used to update the By field and magnetization fields containing : std::pair(std::array<int,8>( number of elements for the calculation, x start, y start, z start, x offset for j(z) spatial derivative, y offset for j(z) spatial derivative, z offset for j(z) spatial derivative, object array index ), std::array<double,2> ( scaling factor (assumed to be 1 right now since conductivity = 0.0), spatial derivative prefactor for j spatial derivative) )
    upLists upBz_; //!< the list of parameters used to update the Bz field and magnetization fields containing : std::pair(std::array<int,8>( number of elements for the calculation, x start, y start, z start, x offset for j(x) spatial derivative, y offset for j(x) spatial derivative, z offset for j(x) spatial derivative, object array index ), std::array<double,2> ( scaling factor (assumed to be 1 right now since conductivity = 0.0), spatial derivative prefactor for j spatial derivative) )

    std::function<void( std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr )> upHxFxn_; //!< function that will update the Hx field
    std::function<void( std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr )> upHyFxn_; //!< function that will update the Hy field
    std::function<void( std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr )> upHzFxn_; //!< function that will update the Hz field

    std::function<void( std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr )> upExFxn_; //!< function that will update the Ex field
    std::function<void( std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr )> upEyFxn_; //!< function that will update the Ey field
    std::function<void( std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr )> upEzFxn_; //!< function that will update the Ez field

    std::function<void(pml_ptr)> updateExPML_; //!< wrapper function to update Ex PMLs
    std::function<void(pml_ptr)> updateEyPML_; //!< wrapper function to update Ey PMLs
    std::function<void(pml_ptr)> updateEzPML_; //!< wrapper function to update Ez PMLs

    std::function<void(pml_ptr)> updateHxPML_; //!< wrapper function to update Hx PMLs
    std::function<void(pml_ptr)> updateHyPML_; //!< wrapper function to update Hy PMLs
    std::function<void(pml_ptr)> updateHzPML_; //!< wrapper function to update Hz PMLs

    std::function< void( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, T*, std::shared_ptr<Obj>) > upLorPxFxn_; //!< function that will update the Polarization fields
    std::function< void( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, T*, std::shared_ptr<Obj>) > upLorPyFxn_; //!< function that will update the Polarization fields
    std::function< void( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, T*, std::shared_ptr<Obj>) > upLorPzFxn_; //!< function that will update the Polarization fields

    std::function< void( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj>) > D2ExFxn_; //!< function that will convert the Dx field to the Ex fields
    std::function< void( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj>) > D2EyFxn_; //!< function that will convert the Dy field to the Ey fields
    std::function< void( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj>) > D2EzFxn_; //!< function that will convert the Dz field to the Ez fields

    std::function< void( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, T*, std::shared_ptr<Obj>) > upLorMxFxn_; //!< function that will update the Polarization fields
    std::function< void( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, T*, std::shared_ptr<Obj>) > upLorMyFxn_; //!< function that will update the Polarization fields
    std::function< void( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, T*, std::shared_ptr<Obj>) > upLorMzFxn_; //!< function that will update the Polarization fields

    std::function< void( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj>) > B2HxFxn_; //!< function that will convert the Dx field to the Ex fields
    std::function< void( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj>) > B2HyFxn_; //!< function that will convert the Dy field to the Ey fields
    std::function< void( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj>) > B2HzFxn_; //!< function that will convert the Dz field to the Ez fields

    std::function< void() > transferEx_; //!< function used to communicate the border values of the Ex field to each processor (either a simple return or Ex_->transferData() depending on circumstances)
    std::function< void() > transferEy_; //!< function used to communicate the border values of the Ey field to each processor (either a simple return or Ey_->transferData() depending on circumstances)
    std::function< void() > transferEz_; //!< function used to communicate the border values of the Ez field to each processor (either a simple return or Ez_->transferData() depending on circumstances)

    std::function< void() > transferHx_; //!< function used to communicate the border values of the Hx field to each processor (either a simple return or Hx_->transferData() depending on circumstances)
    std::function< void() > transferHy_; //!< function used to communicate the border values of the Hy field to each processor (either a simple return or Hy_->transferData() depending on circumstances)
    std::function< void() > transferHz_; //!< function used to communicate the border values of the Hz field to each processor (either a simple return or Hz_->transferData() depending on circumstances)

    std::function<void(pgrid_ptr, std::array<double, 3>&, int, int, int, int, int, int, int, double&, double&, double&)> pbcEx_; //!< function to apply PBC for the Ex field
    std::function<void(pgrid_ptr, std::array<double, 3>&, int, int, int, int, int, int, int, double&, double&, double&)> pbcEy_; //!< function to apply PBC for the Ey field
    std::function<void(pgrid_ptr, std::array<double, 3>&, int, int, int, int, int, int, int, double&, double&, double&)> pbcEz_; //!< function to apply PBC for the Ez field

    std::function<void(pgrid_ptr, std::array<double, 3>&, int, int, int, int, int, int, int, double&, double&, double&)> pbcHx_; //!< function to apply PBC for the Hx field
    std::function<void(pgrid_ptr, std::array<double, 3>&, int, int, int, int, int, int, int, double&, double&, double&)> pbcHy_; //!< function to apply PBC for the Hy field
    std::function<void(pgrid_ptr, std::array<double, 3>&, int, int, int, int, int, int, int, double&, double&, double&)> pbcHz_; //!< function to apply PBC for the Hz field

    std::vector<std::shared_ptr<parallelTFSFBase<T>>> tfsfArr_; //!< vector of all the TFSF objects

    std::vector<cplx> E_incd_; //!< vector of all the incident E field values
    std::vector<cplx> E_pl_incd_; //!< vector of all the incident E field values one point in front of the TFSF source start
    std::vector<cplx> H_incd_; //!< vector of all the incident H field values
    std::vector<cplx> H_mn_incd_; //!< vector of all the incident H field values one point behind the TFSF source start

    std::vector<T> scratch_; //!< vector for scratch operations

    std::array<int,8> axParams_; //!< Temp array to store all the axpy parameters for field updates
    std::array<double,2> prefactors_; //!< Temp array to store all the prefactor parameters for field updates

    std::vector<real_grid_ptr> weights_; //!< a map of the weights for each of the x y and z fields used to determine how to split up the grids for parallelization

    std::vector< std::shared_ptr< parallelDetectorFREQ_Base< T > > > dtcFreqArr_; //!< vector storing all dtcFREQ objects
    std::vector< std::shared_ptr< parallelFluxDTC< T > > > fluxArr_; //!< vector storing all flux objects

    std::shared_ptr<parallelGrid<int>> phys_Ex_; //!< Map of what objects are at each grid point for the Ex field.
    std::shared_ptr<parallelGrid<int>> phys_Ey_; //!< Map of what objects are at each grid point for the Ey field.
    std::shared_ptr<parallelGrid<int>> phys_Ez_; //!< Map of what objects are at each grid point for the Ez field.

    std::shared_ptr<parallelGrid<int>> phys_Hx_; //!< Map of what objects are at each grid point for the Hx field.
    std::shared_ptr<parallelGrid<int>> phys_Hy_; //!< Map of what objects are at each grid point for the Hy field.
    std::shared_ptr<parallelGrid<int>> phys_Hz_; //!< Map of what objects are at each grid point for the Hz field.

public:

    pgrid_ptr Hx_; //!< parallel grid corresponding to the Hx field
    pgrid_ptr Hy_; //!< parallel grid corresponding to the Hy field
    pgrid_ptr Hz_; //!< parallel grid corresponding to the Hz field

    pgrid_ptr Ex_; //!< parallel grid corresponding to the Ex field
    pgrid_ptr Ey_; //!< parallel grid corresponding to the Ey field
    pgrid_ptr Ez_; //!< parallel grid corresponding to the Ez field

    pgrid_ptr Bx_; //!< parallel grid corresponding to the Dx field
    pgrid_ptr By_; //!< parallel grid corresponding to the Dy field
    pgrid_ptr Bz_; //!< parallel grid corresponding to the Dz field

    pgrid_ptr Dx_; //!< parallel grid corresponding to the Dx field
    pgrid_ptr Dy_; //!< parallel grid corresponding to the Dy field
    pgrid_ptr Dz_; //!< parallel grid corresponding to the Dz field

    std::vector<pgrid_ptr> lorPx_; //!< a vector of Polarization fields for the x direction at the current time step
    std::vector<pgrid_ptr> lorPy_; //!< a vector of Polarization fields for the y direction at the current time step
    std::vector<pgrid_ptr> lorPz_; //!< a vector of Polarization fields for the z direction at the current time step

    std::vector<pgrid_ptr> prevLorPx_; //!< a vector of Polarization fields for the x direction at the previous time step
    std::vector<pgrid_ptr> prevLorPy_; //!< a vector of Polarization fields for the y direction at the previous time step
    std::vector<pgrid_ptr> prevLorPz_; //!< a vector of Polarization fields for the z direction at the previous time step

    std::vector<pgrid_ptr> lorMx_; //!< a vector of Polarization fields for the x direction at the current time step
    std::vector<pgrid_ptr> lorMy_; //!< a vector of Polarization fields for the y direction at the current time step
    std::vector<pgrid_ptr> lorMz_; //!< a vector of Polarization fields for the z direction at the current time step

    std::vector<pgrid_ptr> prevLorMx_; //!< a vector of Polarization fields for the x direction at the previous time step
    std::vector<pgrid_ptr> prevLorMy_; //!< a vector of Polarization fields for the y direction at the previous time step
    std::vector<pgrid_ptr> prevLorMz_; //!< a vector of Polarization fields for the z direction at the previous time step

    pml_ptr ExPML_; //!< PML for the Ex field
    pml_ptr EyPML_; //!< PML for the Ey field
    pml_ptr EzPML_; //!< PML for the Ez field
    pml_ptr HxPML_; //!< PML for the Hx field
    pml_ptr HyPML_; //!< PML for the Hy field
    pml_ptr HzPML_; //!< PML for the Hz field

    /**
     * @brief      Constructs a FDTD Propagator class
     *
     * @param[in]  IP        Input parameter object that read in values from a json input file
     * @param[in]  gridComm  A shared_ptr to the MPI interface for the calculation
     */
    parallelFDTDFieldBase(const parallelProgramInputs &IP, std::shared_ptr<mpiInterface> gridComm) :
        gridComm_(gridComm),
        dielectricMatInPML_(false),
        magMatInPML_(false),
        res_(IP.res_),
        t_step_(0),
        yExPBC_(0),
        yEyPBC_(0),
        yEzPBC_(0),
        yHxPBC_(0),
        yHyPBC_(0),
        yHzPBC_(0),
        pbcZMin_(0),
        pbcZMax_(1),
        pmlThickness_(IP.pmlThickness_),
        d_({{1.0/res_, 1.0/res_, 1.0/res_}}),
        dt_(IP.courant_ * d_[0]),
        tcur_(0),
        n_vec_( toN_vec(IP.size_ ) ),
        scratch_((n_vec_[0]+n_vec_[1]) * 2.0, 0.0),
        objArr_(IP.objArr_),
        k_point_(IP.k_point_),
        weights_()
    {
        // Reserve memory for all object vectors
        dtcArr_.reserve( IP.dtcLoc_.size() );
        srcArr_.reserve( IP.srcLoc_.size() );
        tfsfArr_.reserve( IP.tfsfLoc_.size() );
        dtcFreqArr_.reserve( IP.dtcLoc_.size() );
        fluxArr_.reserve( IP.fluxLoc_.size() );

        E_incd_   .reserve( ceil(IP.tMax_ / dt_ ) + 1 );
        E_pl_incd_.reserve( ceil(IP.tMax_ / dt_ ) + 1 );
        H_incd_   .reserve( ceil(IP.tMax_ / dt_ ) + 1 );
        H_mn_incd_.reserve( ceil(IP.tMax_ / dt_ ) + 1 );

        // Set up weights to scale where the process boundaries should be located (based on Instruction Calls for various objects)
        setupWeightsGrid(IP);

        // If only a 2D calculation set the number of points in the z direction to 1, otherwise like any other direction
        int nz = IP.size_[2] == 0 ? 1 : n_vec_[2]+2*gridComm_->npZ();

        // Construct and set up all the physical grids
        setupPhysFields(phys_Ex_, IP.periodic_, std::array<double,3>( {{ 0.5, 0.0, 0.0}} ), nz );
        setupPhysFields(phys_Ey_, IP.periodic_, std::array<double,3>( {{ 0.0, 0.5, 0.0}} ), nz );
        setupPhysFields(phys_Ez_, IP.periodic_, std::array<double,3>( {{ 0.0, 0.0, 0.5}} ), nz );

        setupPhysFields(phys_Hx_, IP.periodic_, std::array<double,3>( {{ 0.0, 0.5, 0.5}} ), nz );
        setupPhysFields(phys_Hy_, IP.periodic_, std::array<double,3>( {{ 0.5, 0.0, 0.5}} ), nz );
        setupPhysFields(phys_Hz_, IP.periodic_, std::array<double,3>( {{ 0.5, 0.5, 0.0}} ), nz );

        // Determine if magnetic or electric dielectric material are in the PMLs
        for(int pp = 0; pp < pmlThickness_.size(); ++pp)
        {
            // Determine what the i, j, k values for the PML are i is in the direction normal to the PML and j and k are inside the plane (i.e. if Left/Right ii = 0 (x) jj = 1(y) kk = 2(z))
            int cor_ii = pp;
            int cor_jj = (pp + 1) % 3;
            int cor_kk = (pp + 2) % 3;

            std::array<int,3> ptVec_ii = {0, 0, 0};
            std::array<int,3> ptVec_jj = {0, 0, 0};
            std::array<int,3> ptVec_kk = {0, 0, 0};

            // Constructs a unit vector to convert ii, jj, kk to xx, yy, zz
            ptVec_ii[cor_ii] = 1;
            ptVec_jj[cor_jj] = 1;
            ptVec_kk[cor_kk] = 1;
            int xx = 0, yy = 0 , zz = 0;
            int max_ii = pmlThickness_[cor_ii];
            // Max value for the left, bottom, or back PMLs
            // If the process starts outside the PML don't include it for dielectric/magnetic PML inclusion
            if(phys_Ex_->procLoc()[cor_ii] < pmlThickness_[cor_ii])
                max_ii = pmlThickness_[cor_ii] - phys_Ex_->procLoc()[cor_ii];
            else
                max_ii = 0;

            // If PML is inside the PML and extends outside of it then max is the max point inside the process
            if(max_ii >= phys_Ex_->ln_vec()[cor_ii])
                max_ii = phys_Ex_->ln_vec()[cor_ii]-1;

            // Min value for the right, top, or right PMLs
            int min_ii = n_vec_[cor_ii] - pmlThickness_[cor_ii] - 1;
            // If PML starts before the process starts set the min value to the start of the process, else if the PML starts inside the process set it to the local starting point, else don't include it
            if(phys_Ex_->procLoc()[cor_ii] >= n_vec_[cor_ii] - pmlThickness_[cor_ii])
                min_ii = 1;
            else if(phys_Ex_->procLoc()[cor_ii]+phys_Ex_->ln_vec()[cor_ii]-2 >= n_vec_[cor_ii] - pmlThickness_[cor_ii])
                min_ii = (n_vec_[cor_ii] - phys_Ex_->procLoc()[cor_ii]) - pmlThickness_[cor_ii];
            else
                min_ii = phys_Ex_->ln_vec()[cor_ii];

            // Loop over PML's area
            for(int jj = 1; jj < phys_Ex_->ln_vec()[cor_jj]-1; ++jj)
            {
                for(int kk = 1; kk < phys_Ex_->ln_vec()[cor_kk]-1; ++kk)
                {
                    // Bottom, Left, Back PML
                    for(int ii = 1; ii < max_ii; ++ii)
                    {
                        xx = ii*ptVec_ii[0]+jj*ptVec_jj[0]+kk*ptVec_kk[0];
                        yy = ii*ptVec_ii[1]+jj*ptVec_jj[1]+kk*ptVec_kk[1];
                        zz = ii*ptVec_ii[2]+jj*ptVec_jj[2]+kk*ptVec_kk[2];
                        if( ( objArr_[ phys_Ex_->point(xx,yy,zz) ]->epsInfty() > 1.0 ) || ( objArr_[ phys_Ex_->point(xx,yy,zz) ]->mat().size() > 1 ) )
                            dielectricMatInPML_ = true;
                        else if( ( objArr_[ phys_Ey_->point(xx,yy,zz) ]->epsInfty() > 1.0 ) || ( objArr_[ phys_Ey_->point(xx,yy,zz) ]->mat().size() > 1 ) )
                            dielectricMatInPML_ = true;
                        else if( ( objArr_[ phys_Ez_->point(xx,yy,zz) ]->epsInfty() > 1.0 ) || ( objArr_[ phys_Ez_->point(xx,yy,zz) ]->mat().size() > 1 ) )
                            dielectricMatInPML_ = true;

                        if( ( objArr_[ phys_Hx_->point(xx,yy,zz)]->muInfty() > 1.0) || ( objArr_[ phys_Hx_->point(xx,yy,zz)]->magMat().size() > 1) )
                            magMatInPML_ = true;
                        else if( ( objArr_[ phys_Hy_->point(xx,yy,zz)]->muInfty() > 1.0) || ( objArr_[ phys_Hy_->point(xx,yy,zz)]->magMat().size() > 1) )
                            magMatInPML_ = true;
                        else if( ( objArr_[ phys_Hz_->point(xx,yy,zz)]->muInfty() > 1.0) || ( objArr_[ phys_Hz_->point(xx,yy,zz)]->magMat().size() > 1) )
                            magMatInPML_ = true;
                    }
                    // Top, Right, Front PML
                    for(int ii = phys_Ex_->ln_vec()[cor_ii]-2; ii >= min_ii; --ii)
                    {
                        xx = ii*ptVec_ii[0] + jj*ptVec_jj[0] + kk*ptVec_kk[0];
                        yy = ii*ptVec_ii[1] + jj*ptVec_jj[1] + kk*ptVec_kk[1];
                        zz = ii*ptVec_ii[2] + jj*ptVec_jj[2] + kk*ptVec_kk[2];

                        if( ( objArr_[ phys_Ex_->point(xx,yy,zz) ]->epsInfty() > 1.0 ) || ( objArr_[ phys_Ex_->point(xx,yy,zz) ]->mat().size() > 1 ) )
                            dielectricMatInPML_ = true;
                        else if( ( objArr_[ phys_Ey_->point(xx,yy,zz) ]->epsInfty() > 1.0 ) || ( objArr_[ phys_Ey_->point(xx,yy,zz) ]->mat().size() > 1 ) )
                            dielectricMatInPML_ = true;
                        else if( ( objArr_[ phys_Ez_->point(xx,yy,zz) ]->epsInfty() > 1.0 ) || ( objArr_[ phys_Ez_->point(xx,yy,zz) ]->mat().size() > 1 ) )
                            dielectricMatInPML_ = true;

                        if( ( objArr_[ phys_Hx_->point(xx,yy,zz)]->muInfty() > 1.0) || ( objArr_[ phys_Hx_->point(xx,yy,zz)]->magMat().size() > 1) )
                            magMatInPML_ = true;
                        else if( ( objArr_[ phys_Hy_->point(xx,yy,zz)]->muInfty() > 1.0) || ( objArr_[ phys_Hy_->point(xx,yy,zz)]->magMat().size() > 1) )
                            magMatInPML_ = true;
                        else if( ( objArr_[ phys_Hz_->point(xx,yy,zz)]->muInfty() > 1.0) || ( objArr_[ phys_Hz_->point(xx,yy,zz)]->magMat().size() > 1) )
                            magMatInPML_ = true;
                    }
                }
            }
        }

        // See if any objects would require electric dispersive, or magnetic dispersive materials
        bool disp = false;
        bool magnetic = false;
        for(auto& obj : objArr_)
        {
            if(obj->magMat().size() > 1)
                magnetic = true;
            if(obj->mat().size() > 1 || obj->epsInfty() > 1.0)
                disp = true;
        }

        if(dielectricMatInPML_)
            disp = true;
        if(magMatInPML_)
            magnetic = true;
        // Initialize all girds
        if(IP.size_[2] != 0 || IP.pol_ == POLARIZATION::HZ || IP.pol_ == POLARIZATION::EX || IP.pol_ == POLARIZATION::EY)
        {
            // Defined TEz mode
            Ex_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_->npX(), n_vec_[1]+2*gridComm_->npY(), nz }}), d_, false);
            Hz_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_->npX(), n_vec_[1]+2*gridComm_->npY(), nz }}), d_, true);
            Ey_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_->npX(), n_vec_[1]+2*gridComm_->npY(), nz }}), d_, true);

            if(disp)
            {
                Dx_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_->npX(), n_vec_[1]+2*gridComm_->npY(), nz }}), d_, false);
                Dy_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_->npX(), n_vec_[1]+2*gridComm_->npY(), nz }}), d_, true);
            }
            if(magnetic)
            {
                Bz_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_->npX(), n_vec_[1]+2*gridComm_->npY(), nz }}), d_, true);
            }
            ln_vec_[0] = Hz_->local_x()-2;
            ln_vec_[1] = Hz_->local_y()-2;
            ln_vec_[2] = (nz == 1) ? 1 : Hz_->local_z()-2;
        }
        if(IP.size_[2] != 0 || IP.pol_ == POLARIZATION::EZ || IP.pol_ == POLARIZATION::HX || IP.pol_ == POLARIZATION::HY)
        {
            // Define TMz mode
            Hx_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }}),d_, true);
            Ez_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }}),d_,false);
            Hy_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }}),d_,false);

            if(disp)
            {
                Dz_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_->npX(), n_vec_[1]+2*gridComm_->npY(),  nz }}), d_, false);
            }
            if(magnetic)
            {
                Bx_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_->npX(), n_vec_[1]+2*gridComm_->npY(),  nz }}), d_, true);
                By_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_->npX(), n_vec_[1]+2*gridComm_->npY(),  nz }}), d_, false);
            }
            ln_vec_[0] = Ez_->local_x()-2;
            ln_vec_[1] = Ez_->local_y()-2;
            ln_vec_[2] = (nz == 1) ? 1 : Ez_->local_z()-2;
        }
        pbcZMin_ = (nz == 1) ? 0 : 1;
        pbcZMax_ = (nz == 1) ? 1 : ln_vec_[2];
        // Initialize object specific grids
        for(auto & obj :objArr_)
        {
            obj->setUpConsts(dt_);
            while((obj->mat().size() - 1) / 3.0 > lorPx_.size())
            {
                if(Dx_)
                {
                    prevLorPy_.push_back(std::make_shared<parallelGrid<T>>(gridComm_, false, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, true) );
                    prevLorPx_.push_back(std::make_shared<parallelGrid<T>>(gridComm_, false, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, false) );
                    lorPx_.push_back(std::make_shared<parallelGrid<T>>(gridComm_, false, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, false) );
                    lorPy_.push_back(std::make_shared<parallelGrid<T>>(gridComm_, false, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, true) );
                }
                if(Dz_)
                {
                    lorPz_.push_back(std::make_shared<parallelGrid<T>>(gridComm_, false, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, false) );
                    prevLorPz_.push_back(std::make_shared<parallelGrid<T>>(gridComm_, false, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, false) );
                }
            }
            while((obj->magMat().size() - 1) / 3.0 > lorMx_.size())
            {
                if(Bx_)
                {
                    prevLorMy_.push_back(std::make_shared<parallelGrid<T>>(gridComm_, false, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, true) );
                    prevLorMx_.push_back(std::make_shared<parallelGrid<T>>(gridComm_, false, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, false) );
                    lorMx_.push_back(std::make_shared<parallelGrid<T>>(gridComm_, false, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, false) );
                    lorMy_.push_back(std::make_shared<parallelGrid<T>>(gridComm_, false, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, true) );
                }
                if(Bz_)
                {
                    lorMz_.push_back(std::make_shared<parallelGrid<T>>(gridComm_, false, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, false) );
                    prevLorMz_.push_back(std::make_shared<parallelGrid<T>>(gridComm_, false, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, false) );
                }
            }
        }
        for(auto& xx : IP.inputMapSlicesX_)
            convertInputs2Map(IP, DIRECTION::X, xx);
        for(auto& yy : IP.inputMapSlicesY_)
            convertInputs2Map(IP, DIRECTION::Y, yy);
        for(auto& zz : IP.inputMapSlicesZ_)
            convertInputs2Map(IP, DIRECTION::Z, zz);
    }

    /**
     * @brief      Creates the update lists for a field
     *
     * @param[in]  physGrid  Map of all objects for the gird
     * @param[in]  pml       The pml for the grid
     * @param[in]  E         True if gird is an electric field
     * @param[in]  derivOff  An array describing the offset for the j spatial derivative
     * @param[in]  fieldEnd  The first grid point outside where the field is defined
     * @param[in]  d         grid spacing
     * @param      upU       upE/upH update list
     * @param      upD       The upD/upB update list
     */
    void initializeList(std::shared_ptr<parallelGrid<int>> physGrid, std::shared_ptr<parallelCPML<T>> pml, bool E, std::array<int,3> derivOff, std::array<int,3> fieldEnd, double d, upLists& upU, upLists& upD )
    {
        int zmin = physGrid->z() == 1 ? 0 : 1;
        int zmax = physGrid->z() == 1 ? 1 : physGrid->local_z()-1;
        std::vector<std::array<int, 5>> tempU, tempD;
        std::tie(tempU, tempD) = getBlasLists(  E, {{ 1, 1, zmin }}, {{ physGrid->ln_vec()[0]-1, physGrid->ln_vec()[1]-1, zmax }}, physGrid, pml);
        for(auto & up : tempU)
        {
            double ep_mu = E ? objArr_[up[4]]->epsInfty() : objArr_[up[4]]->muInfty();
            if( up[1] + physGrid->procLoc()[1] != fieldEnd[1] && up[2] + physGrid->procLoc()[2] != fieldEnd[2] )
            {
                if(up[3] + up[0] - 1 == fieldEnd[0])
                    upU.push_back(std::make_pair(std::array<int,8>({up[3]-1, up[0], up[1], up[2], derivOff[0], derivOff[1] , derivOff[2], up[4]}), std::array<double,2>({1.0, -1.0*dt_/(ep_mu*d)}) ) );
                else
                    upU.push_back(std::make_pair(std::array<int,8>({up[3]  , up[0], up[1], up[2], derivOff[0], derivOff[1] , derivOff[2], up[4]}), std::array<double,2>({1.0, -1.0*dt_/(ep_mu*d)}) ) );
            }
        }
        for(auto & up : tempD)
        {
            if( up[1] + physGrid->procLoc()[1] != fieldEnd[1] && up[2] + physGrid->procLoc()[2] != fieldEnd[2] )
            {
                if(up[3] + up[0] - 1 == fieldEnd[0])
                    upD.push_back(std::make_pair(std::array<int,8>({up[3]-1, up[0], up[1], up[2], derivOff[0], derivOff[1] , derivOff[2], up[4]}), std::array<double,2>({1.0, -1.0*dt_/d}) ) );
                else
                    upD.push_back(std::make_pair(std::array<int,8>({up[3]  , up[0], up[1], up[2], derivOff[0], derivOff[1] , derivOff[2], up[4]}), std::array<double,2>({1.0, -1.0*dt_/d}) ) );
            }
        }
    }

    /**
     * @brief      Constructs a DTC based off of the input parameters and puts it in the proper detector vector
     *
     * @param[in]  c             class type of the dtc (bin, bmp, cout, txt, freq)
     * @param[in]  grid          vector of the fields that need to be outputted
     * @param[in]  SI            true if outputting in SI units
     * @param[in]  loc           The location of the detectors lower left corner in grid points
     * @param[in]  sz            The size of the detector in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  fxn           Function used to modify base field data
     * @param[in]  txtType       if BMP what should be outputted to the text file
     * @param[in]  type          The type of the detector (Ex, Ey, Epow, etc)
     * @param[in]  freqList      The frequency list
     * @param[in]  timeInterval  The number of time steps per field output
     * @param[in]  a             unit length of the calculation
     * @param[in]  I0            unit current of the calculation
     * @param[in]  t_max         The time at the final time step
     */
    virtual void coustructDTC(DTCCLASS c, std::vector<pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, GRIDOUTFXN fxn, GRIDOUTTYPE txtType, DTCTYPE type, std::vector<double> freqList, double timeInterval, double a, double I0, double t_max) = 0;

    /**
     * @brief      Generates the FDTD update lists
     *
     * @param[in]  E         True if the field the list is being generated for is an electric field
     * @param[in]  min       An array containing the minimum grid points in all direction
     * @param[in]  max       An array containing the maximum grid points in all directions
     * @param[in]  physGrid  The map of all objects on that grid
     * @param[in]  pml       The cPML associated with the grid
     *
     * @return     A tuple containing all update lists for the grid containing (x start, y start, z start, number of elements to include, object parameters to use)
     */
    std::tuple< std::vector<std::array<int,5>>, std::vector<std::array<int,5>> > getBlasLists(bool E, std::array<int,3> min, std::array<int,3> max, std::shared_ptr<parallelGrid<int>> physGrid, std::shared_ptr<parallelCPML<T>> pml)
    {
        std::vector<std::array<int,5>> upULists;
        std::vector<std::array<int,5>> upDLists;

        // Include PMLs if material is inside them
        int PML_y_bot = ( (E && dielectricMatInPML_) || (!E && magMatInPML_) ) ? pml->lny_bot() : 0;
        int PML_y_top = ( (E && dielectricMatInPML_) || (!E && magMatInPML_) ) ? pml->lny_top() : 0;
        int PML_x_right = ( (E && dielectricMatInPML_) || (!E && magMatInPML_) ) ? pml->lnx_right() : 0;
        int PML_x_left  = ( (E && dielectricMatInPML_) || (!E && magMatInPML_) ) ? pml->lnx_left() : 0;
        int PML_z_back  = (Ez_ && Hz_ && ( (E && dielectricMatInPML_) || (!E && magMatInPML_) ) ) ? pml->lnz_back()  : 0;
        int PML_z_front = (Ez_ && Hz_ && ( (E && dielectricMatInPML_) || (!E && magMatInPML_) ) ) ? pml->lnz_front() : 0;

        std::shared_ptr<Obj> obj;
        // If inside the bottom PML everything should be updated with a D field
        for(int jj = min[1]; jj < min[1]+PML_y_bot; ++jj)
        {
            for(int kk = min[2]; kk < max[2]; ++kk )
            {
                int ii = min[0];
                while(ii < max[0])
                {
                    int iistore = ii;
                    // Check if points are in same object
                    while ( (ii < max[0]-1) && ( physGrid->point(ii,jj,kk) == physGrid->point(ii+1,jj,kk) ) )
                        ++ii;

                    obj = objArr_[ physGrid->point(iistore, jj, kk) ];
                    upDLists.push_back({{ iistore,jj,kk,ii-iistore+1,static_cast<int>(physGrid->point(iistore,jj,kk))}});

                    ++ii;
                }
            }
        }
        // Central region only the left, right, front and back PMLs should be 100% inside D lists
        for(int jj = min[1]+PML_y_bot; jj < max[1]-PML_y_top; ++jj)
        {
            for(int kk = min[2]; kk < min[2] + PML_z_back; ++kk )
            {
                int ii = min[0];
                while(ii < max[0])
                {
                    int iistore = ii;
                    // Check if points are in same object
                    while ( (ii < max[0]-1) && ( physGrid->point(ii,jj,kk) == physGrid->point(ii+1,jj,kk) ) )
                        ++ii;

                    objArr_[ physGrid->point(iistore, jj, kk) ];
                    upDLists.push_back({{ iistore,jj,kk,ii-iistore+1,static_cast<int>(physGrid->point(iistore,jj,kk))}});

                    ++ii;
                }
            }
            for(int kk = min[2]+PML_z_back; kk < max[2]-PML_z_front; ++kk)
            {
                int ii = min[0];
                while(ii < max[0])
                {
                    int iistore = ii;
                    // Check if points are in same object and whether or not at a PML boundary (always cut off at the boundary)
                    while ( (ii < max[0]-1) && ( physGrid->point(ii,jj,kk) == physGrid->point(ii+1,jj,kk) ) && (ii != PML_x_left ) && (ii != max[0]-PML_x_right ) )
                        ++ii;

                    obj = objArr_[ physGrid->point(iistore, jj, kk) ];
                    // If the point is not in a PML or dispersive material put it in update U field; else update D or
                    if( ( !E && obj->magMat().size() <= 1 ) || ( E && obj->mat().size() <= 1 && ii > PML_x_left && ii <=  max[0]-PML_x_right ) )
                        upULists.push_back({{ iistore,jj,kk,ii-iistore+1,static_cast<int>(physGrid->point(iistore,jj,kk))}});
                    else
                        upDLists.push_back({{ iistore,jj,kk,ii-iistore+1,static_cast<int>(physGrid->point(iistore,jj,kk))}});
                    ++ii;
                }
            }
            for(int kk = max[2]-PML_z_front; kk < max[2]; ++kk)
            {
                int ii = min[0];
                while(ii < max[0])
                {
                    int iistore = ii;
                    // Check if points are in same object
                    while ( (ii < max[0]-1) && ( physGrid->point(ii,jj,kk) == physGrid->point(ii+1,jj,kk) ) )
                        ++ii;
                    objArr_[ physGrid->point(iistore, jj, kk) ];
                    upDLists.push_back({{ iistore,jj,kk,ii-iistore+1,static_cast<int>(physGrid->point(iistore,jj,kk))}});
                    ++ii;
                }
            }
        }
        // If inside the top PML everything should be updated with a D field
        for(int jj = max[1]-PML_y_top; jj < max[1]; ++jj)
        {
            for(int kk = min[2]; kk < max[2]; ++kk )
            {
                int ii = min[0];
                while(ii < max[0])
                {
                    int iistore = ii;
                    // Check if points are in same object
                    while ( (ii < max[0]-1) && ( physGrid->point(ii,jj,kk) == physGrid->point(ii+1,jj,kk) ) )
                        ++ii;

                    objArr_[ physGrid->point(iistore, jj, kk) ];
                    upDLists.push_back({{ iistore,jj,kk,ii-iistore+1,static_cast<int>(physGrid->point(iistore,jj,kk))}});

                    ++ii;
                }
            }
        }
        return std::make_tuple(upULists, upDLists);
    }

    /**
     * @brief      Sets up the object map grids for each field
     *
     * @param      physGrid  The object map to be made
     * @param[in]  PBC       True if periodic boundary conditions used
     * @param[in]  offPt     values for the offset from the base grid point for each field
     * @param[in]  nz        number of grid points in the z direction
     */
    void setupPhysFields(std::shared_ptr<parallelGrid<int>>& physGrid, bool PBC, std::array<double,3> offPt, int nz)
    {
        physGrid = std::make_shared<parallelGrid<int> >(gridComm_, PBC, weights_, std::array<int,3>( {{ n_vec_[0]+2*gridComm_->npX(),n_vec_[1]+2*gridComm_->npY(), nz }} ), d_, false);
        // Test point used
        std::array<double,3> pt = {{ 0,0,0}};
        // Objects on the same points over write each other (last object made wins)
        int zmin = (nz == 1) ? 0 : 1;
        int zmax = (nz == 1) ? 1 : physGrid->ln_vec()[2]-1;
        for(int oo = 0; oo < objArr_.size(); ++oo)
        {
            // look at all local points only
            if(oo == 0 || objArr_[oo]->mat().size() > 1 || objArr_[oo]->epsInfty() != 1.0 || objArr_[oo]->magMat().size() > 1 || objArr_[oo]->muInfty() != 1.0 )
            {
                for(int ii = 1; ii < physGrid->ln_vec()[0]-1; ++ii)
                {
                    for(int jj = 1; jj < physGrid->ln_vec()[1]-1; ++jj)
                    {
                        for(int kk = zmin; kk < zmax; ++kk)
                        {
                            pt[0] = ( (ii-1) + offPt[0] + physGrid->procLoc()[0] - (n_vec_[0]-n_vec_[0] % 2)/2.0 )*d_[0];
                            pt[1] = ( (jj-1) + offPt[1] + physGrid->procLoc()[1] - (n_vec_[1]-n_vec_[1] % 2)/2.0 )*d_[1];
                            pt[2] = ( (kk-1) + offPt[2] + physGrid->procLoc()[2] - (n_vec_[2]-n_vec_[2] % 2)/2.0 )*d_[2];
                            if(objArr_[oo]->isObj(pt,d_[0])==true)
                            {
                                physGrid->point(ii,jj, kk) = oo;
                            }
                        }
                    }
                }
            }
        }
        // All borders of between the processors have a -1 to indicate they are borders
        zmin = 0;
        zmax = (nz == 1) ? 1 : physGrid->ln_vec()[2];
        if(nz != 1)
        {
            for(int ii = 0; ii < physGrid->ln_vec()[0]; ++ii)
            {
                for(int jj = 0; jj < physGrid->ln_vec()[1]; ++jj)
                {
                    physGrid->point(ii, jj, 0   ) = -1;
                    physGrid->point(ii, jj, zmax-1) = -1;
                }
            }
        }

        for(int ii = 0; ii < physGrid->ln_vec()[0]; ++ii)
        {
            for(int kk = zmin; kk < zmax; ++kk)
            {
                physGrid->point(ii, 0, kk) = -1;
                physGrid->point(ii, physGrid->ln_vec()[1]-1, kk) = -1;
            }
        }
        for(int kk = zmin; kk < zmax; ++kk)
        {
            for(int jj = 0; jj < physGrid->ln_vec()[1]; ++jj)
            {
                physGrid->point(0, jj, kk) = -1;
                physGrid->point(physGrid->ln_vec()[0]-1, jj, kk) = -1;
            }
        }
    }

    /**
     * @brief      Sets up the weight grid for dividing processes by taking average of H/E fields and for 3D calcs the flux region calcs
     *
     * @param[in]  IP    Input parameter object that is being used to construct the propagator
     */
    void setupWeightsGrid(const parallelProgramInputs &IP)
    {
        // test point will move across all grid points
        // Set weights based off of normal materials, number of axpy calls per time step (6 base + 3*n_lor_pol = 6 + 3*(objArr_[kk]->mat().size()-1)/3)*number of E fields
        std::array<double,3> pt = {{0,0,0}};
        weights_.push_back(std::make_shared<Grid<double>>(n_vec_, d_) );
        weights_.push_back(std::make_shared<Grid<double>>(n_vec_, d_) );
        weights_.push_back(std::make_shared<Grid<double>>(n_vec_, d_) );
        for(int oo = 0; oo < objArr_.size(); ++oo)
        {
            for(int ii = 0; ii < n_vec_[0]; ++ii)
            {
                for(int jj = 0; jj < n_vec_[1]; ++jj)
                {
                    for(int kk = 0; kk < n_vec_[2]; ++kk)
                    {
                        // Ex points located at ii+1/2, jj
                        pt[0] = (ii-(n_vec_[0]-1)/2.0+0.5)*d_[0];
                        pt[1] = (jj-(n_vec_[1]-1)/2.0    )*d_[1];
                        pt[2] = (kk-(n_vec_[2]-1)/2.0    )*d_[2];
                        if(objArr_[oo]->isObj(pt,d_[0]))
                            weights_[0]->point(ii,jj,kk) =  (2.0 + (2.0666666666666*static_cast<double>(objArr_[oo]->mat().size()-1) + 2.0666666666666*static_cast<double>(objArr_[oo]->magMat().size()-1) )/2.0 ) ;
                        // Ey points located ii, jj+1/2
                        pt[1] += 0.5*d_[1];
                        pt[0] -= 0.5*d_[0];
                        if(objArr_[oo]->isObj(pt,d_[0]))
                            weights_[1]->point(ii,jj,kk) =  (2.0 + (2.0666666666666*static_cast<double>(objArr_[oo]->mat().size()-1) + 2.0666666666666*static_cast<double>(objArr_[oo]->magMat().size()-1) )/2.0 ) ;
                        // Ez point is at ii, jj
                        pt[1] -= 0.5*d_[1];
                        pt[2] += 0.5*d_[2];
                        if(objArr_[oo]->isObj(pt,d_[0]))
                            weights_[2]->point(ii,jj,kk) =  (2.0 + (2.0666666666666*static_cast<double>(objArr_[oo]->mat().size()-1) + 2.0666666666666*static_cast<double>(objArr_[oo]->magMat().size()-1) )/2.0 ) ;
                    }
                }
            }
        }
        std::vector<double> copyZero(n_vec_[0]*n_vec_[2], 0.0);
        // Ey has no points on the ny-1 face
        std::copy_n( copyZero.data(), n_vec_[0]*n_vec_[2], &weights_[1]->point(0, n_vec_[1]-1, 0) );
        // Ex has no points on the nx-1 face
        for(int yy = 0; yy < n_vec_[1]; ++yy)
            dcopy_(n_vec_[2], copyZero.data(), 1, &weights_[0]->point(n_vec_[0]-1, yy, 0), weights_[0]->x() );
        // Ez has no points on the nz-1 face
        for(int yy = 0; yy < n_vec_[1]; ++yy)
            std::copy_n(copyZero.data(), n_vec_[0], &weights_[2]->point(0, yy, n_vec_[2]-1) );
        // Add PMLs
        for(int xx = 0; xx < IP.pmlThickness_[0]; ++xx)
        {
            // Ex field has no PML along the left and right
            for(int yy = 0; yy < weights_[0]->y(); ++ yy)
            {
                daxpy_(weights_[1]->z(), 1.0, std::vector<double>(n_vec_[2], 3.1).data(), 1, &weights_[1]->point(            xx, yy, 0), weights_[1]->x());
                daxpy_(weights_[2]->z(), 1.0, std::vector<double>(n_vec_[2], 3.1).data(), 1, &weights_[2]->point(            xx, yy, 0), weights_[2]->x());

                daxpy_(weights_[1]->z(), 1.0, std::vector<double>(n_vec_[2], 3.1).data(), 1, &weights_[1]->point(n_vec_[0]-1-xx, yy, 0), weights_[1]->x());
                daxpy_(weights_[2]->z(), 1.0, std::vector<double>(n_vec_[2], 3.1).data(), 1, &weights_[2]->point(n_vec_[0]-1-xx, yy, 0), weights_[2]->x());
            }
        }
        std::fill_n(copyZero.begin(), copyZero.size(), 3.1);
        for(int yy = 0; yy < IP.pmlThickness_[1]; ++yy)
        {
            // Ey field has no PML along the top and bottom
            std::transform(copyZero.begin(), copyZero.end(), &weights_[0]->point(0,             yy, 0), &weights_[0]->point(0,             yy, 0), [](double a, double b){return a+b;});
            std::transform(copyZero.begin(), copyZero.end(), &weights_[2]->point(0,             yy, 0), &weights_[2]->point(0,             yy, 0), [](double a, double b){return a+b;});

            std::transform(copyZero.begin(), copyZero.end(), &weights_[0]->point(0, n_vec_[1]-1-yy, 0), &weights_[0]->point(0, n_vec_[1]-1-yy, 0), [](double a, double b){return a+b;});
            std::transform(copyZero.begin(), copyZero.end(), &weights_[2]->point(0, n_vec_[1]-1-yy, 0), &weights_[2]->point(0, n_vec_[1]-1-yy, 0), [](double a, double b){return a+b;});
        }
        if(IP.size_[2] > 0)
        {
            for(int zz = 0; zz < IP.pmlThickness_[2]; ++zz)
            {
                // Ez field has no PML along the front/back
                for(int yy = 0; yy < weights_[0]->y(); ++ yy)
                {
                    std::transform(copyZero.begin(), copyZero.begin()+weights_[0]->x(), &weights_[0]->point(0, yy,             zz), &weights_[0]->point(0, yy,             zz), [](double a, double b){return a+b;});
                    std::transform(copyZero.begin(), copyZero.begin()+weights_[1]->x(), &weights_[1]->point(0, yy,             zz), &weights_[1]->point(0, yy,             zz), [](double a, double b){return a+b;});

                    std::transform(copyZero.begin(), copyZero.begin()+weights_[0]->x(), &weights_[0]->point(0, yy, n_vec_[2]-1-zz), &weights_[0]->point(0, yy, n_vec_[2]-1-zz), [](double a, double b){return a+b;});
                    std::transform(copyZero.begin(), copyZero.begin()+weights_[1]->x(), &weights_[1]->point(0, yy, n_vec_[2]-1-zz), &weights_[1]->point(0, yy, n_vec_[2]-1-zz), [](double a, double b){return a+b;});

                }
            }
        }
        // Computational cost of flux for 2D calcs is small so it is neglected.
        if(IP.size_[2] > 0)
        {
            // Include weights for flux regions
            for(int ff = 0; ff < IP.fluxLoc_.size(); ++ff)
            {
                std::vector<double> fluxWeight( std::max( std::max(n_vec_[0], n_vec_[1]), n_vec_[2] ), IP.fluxFreqList_[ff].size() );
                for(int yy = 0; yy < IP.fluxSz_[ff][1]; ++yy)
                {
                    daxpy_(IP.fluxSz_[ff][0], 1.0, fluxWeight.data(), 1, &weights_[0]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1]+yy, IP.fluxLoc_[ff][2]), 1);
                    daxpy_(IP.fluxSz_[ff][0], 1.0, fluxWeight.data(), 1, &weights_[1]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1]+yy, IP.fluxLoc_[ff][2]), 1);
                    daxpy_(IP.fluxSz_[ff][0], 1.0, fluxWeight.data(), 1, &weights_[2]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1]+yy, IP.fluxLoc_[ff][2]), 1);

                    daxpy_(IP.fluxSz_[ff][0], 1.0, fluxWeight.data(), 1, &weights_[0]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1]+yy, IP.fluxLoc_[ff][2]+IP.fluxSz_[ff][2]-1), 1);
                    daxpy_(IP.fluxSz_[ff][0], 1.0, fluxWeight.data(), 1, &weights_[1]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1]+yy, IP.fluxLoc_[ff][2]+IP.fluxSz_[ff][2]-1), 1);
                    daxpy_(IP.fluxSz_[ff][0], 1.0, fluxWeight.data(), 1, &weights_[2]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1]+yy, IP.fluxLoc_[ff][2]+IP.fluxSz_[ff][2]-1), 1);

                    daxpy_(IP.fluxSz_[ff][2], 1.0, fluxWeight.data(), 1, &weights_[0]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1]+yy, IP.fluxLoc_[ff][2]), 1);
                    daxpy_(IP.fluxSz_[ff][2], 1.0, fluxWeight.data(), 1, &weights_[1]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1]+yy, IP.fluxLoc_[ff][2]), 1);
                    daxpy_(IP.fluxSz_[ff][2], 1.0, fluxWeight.data(), 1, &weights_[2]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1]+yy, IP.fluxLoc_[ff][2]), 1);

                    daxpy_(IP.fluxSz_[ff][2], 1.0, fluxWeight.data(), 1, &weights_[0]->point(IP.fluxLoc_[ff][0]+IP.fluxSz_[ff][0]-1, IP.fluxLoc_[ff][1]+yy, IP.fluxLoc_[ff][2]), weights_[0]->x());
                    daxpy_(IP.fluxSz_[ff][2], 1.0, fluxWeight.data(), 1, &weights_[1]->point(IP.fluxLoc_[ff][0]+IP.fluxSz_[ff][0]-1, IP.fluxLoc_[ff][1]+yy, IP.fluxLoc_[ff][2]), weights_[1]->x());
                    daxpy_(IP.fluxSz_[ff][2], 1.0, fluxWeight.data(), 1, &weights_[2]->point(IP.fluxLoc_[ff][0]+IP.fluxSz_[ff][0]-1, IP.fluxLoc_[ff][1]+yy, IP.fluxLoc_[ff][2]), weights_[2]->x());
                }
                for(int zz = 0; zz < IP.fluxSz_[ff][2]; ++zz)
                {
                    daxpy_(IP.fluxSz_[ff][0], 1.0, fluxWeight.data(), 1, &weights_[0]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1], IP.fluxLoc_[ff][2]+zz), 1);
                    daxpy_(IP.fluxSz_[ff][0], 1.0, fluxWeight.data(), 1, &weights_[1]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1], IP.fluxLoc_[ff][2]+zz), 1);
                    daxpy_(IP.fluxSz_[ff][0], 1.0, fluxWeight.data(), 1, &weights_[2]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1], IP.fluxLoc_[ff][2]+zz), 1);

                    daxpy_(IP.fluxSz_[ff][0], 1.0, fluxWeight.data(), 1, &weights_[0]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1]+IP.fluxSz_[ff][1]-1, IP.fluxLoc_[ff][2]+zz), 1);
                    daxpy_(IP.fluxSz_[ff][0], 1.0, fluxWeight.data(), 1, &weights_[1]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1]+IP.fluxSz_[ff][1]-1, IP.fluxLoc_[ff][2]+zz), 1);
                    daxpy_(IP.fluxSz_[ff][0], 1.0, fluxWeight.data(), 1, &weights_[2]->point(IP.fluxLoc_[ff][0], IP.fluxLoc_[ff][1]+IP.fluxSz_[ff][1]-1, IP.fluxLoc_[ff][2]+zz), 1);
                }
            }
        }
    }
    /**
     * @brief return the current time
     * @return tcur_
     */
    inline double getTime(){return tcur_;}

    /**
     * @brief      steps the propagator forward one unit in time
     */
    void step()
    {
        // Update H/B fields
        updateH();
        updateB();

        // Include PML updates before transferring from B to H
        updateHxPML_(HxPML_);
        updateHyPML_(HyPML_);
        updateHzPML_(HzPML_);

        // Update all magnetization terms to the H field
        updateMagH();

        // Add the H incident field before stepping tfsf objects (like H updates) and then add the E incd (like normal updates)
        for(auto & tfsf : tfsfArr_)
        {
               H_incd_.push_back(tfsf->   H_incd());
            H_mn_incd_.push_back(tfsf->H_mn_incd());
            tfsf->updateFileds();
               E_incd_.push_back(tfsf->   E_incd());
            E_pl_incd_.push_back(tfsf->E_pl_incd());
        }
        for(auto & src :srcArr_)
            src->addPul(tcur_);

        // Transfer PBC and MPI related H field information
        pbcHx_(Hx_, k_point_, ln_vec_[0]+1, ln_vec_[1]+1, pbcZMax_+1, ln_vec_[0]+1, yHxPBC_, pbcZMin_, pbcZMax_  , d_[0], d_[1], d_[2] );
        pbcHy_(Hy_, k_point_, ln_vec_[0]+1, ln_vec_[1]+1, pbcZMax_+1, ln_vec_[0]  , yHyPBC_, pbcZMin_, pbcZMax_  , d_[0], d_[1], d_[2] );
        pbcHz_(Hz_, k_point_, ln_vec_[0]+1, ln_vec_[1]+1, pbcZMax_+1, ln_vec_[0]  , yHzPBC_, pbcZMin_, pbcZMax_+1, d_[0], d_[1], d_[2] );

        transferHx_();
        transferHy_();
        transferHz_();

        // Update E.D fields
        updateE();
        updateD();

        // Include PML updates before transferring from D to E
        updateExPML_(ExPML_);
        updateEyPML_(EyPML_);
        updateEzPML_(EzPML_);

        // Update all polarization terms to the E field
        updateDispE();

        // All E-field Updates should now be completed transfer border values for the E-fields
        pbcEx_(Ex_, k_point_, ln_vec_[0]+1, ln_vec_[1]+1, pbcZMax_+1, ln_vec_[0]  , yExPBC_, pbcZMin_, pbcZMax_+1, d_[0], d_[1], d_[2]);
        pbcEy_(Ey_, k_point_, ln_vec_[0]+1, ln_vec_[1]+1, pbcZMax_+1, ln_vec_[0]+1, yEyPBC_, pbcZMin_, pbcZMax_+1, d_[0], d_[1], d_[2]);
        pbcEz_(Ez_, k_point_, ln_vec_[0]+1, ln_vec_[1]+1, pbcZMax_+1, ln_vec_[0]+1, yEzPBC_, pbcZMin_, pbcZMax_  , d_[0], d_[1], d_[2]);

        transferEx_();
        transferEy_();
        transferEz_();

        // Increment time steps to before output as all fields should be updated to the next time step now
        tcur_ += dt_;
        ++t_step_;

        // Output all detector values
        for(auto & dtc : dtcArr_)
            if(t_step_ % dtc->timeInt() == 0)
                dtc->output(tcur_);
        for(auto & dtc : dtcFreqArr_)
            if(t_step_ % dtc->timeInt() == 0)
                dtc->output(tcur_);
        for(auto & flux : fluxArr_)
            if(t_step_ % flux->timeInt() == 0)
                flux->fieldIn(tcur_);
    }

    /**
     * @brief      Updates the H fields forward in time
     */
    void updateH()
    {
        for(auto& up : upHx_)
            upHxFxn_(std::get<0>(up), std::get<1>(up), Hx_, Ey_, Ez_);
        for(auto& up : upHy_)
            upHyFxn_(std::get<0>(up), std::get<1>(up), Hy_, Ez_, Ex_);
        for(auto& up : upHz_)
            upHzFxn_(std::get<0>(up), std::get<1>(up), Hz_, Ex_, Ey_);
    }

    /**
     * @brief      Updates the E fields forward in time
     */
    void updateE()
    {
        for(auto& up : upEx_)
            upExFxn_(std::get<0>(up), std::get<1>(up), Ex_, Hy_, Hz_);
        for(auto& up : upEy_)
            upEyFxn_(std::get<0>(up), std::get<1>(up), Ey_, Hz_, Hx_);
        for(auto& up : upEz_)
            upEzFxn_(std::get<0>(up), std::get<1>(up), Ez_, Hx_, Hy_);
    }

    /**
     * @brief      Updates the B fields forward in time
     */
    void updateB()
    {
        for(auto& up : upBx_)
            upHxFxn_(std::get<0>(up), std::get<1>(up), Bx_, Ey_, Ez_);
        for(auto& up : upBy_)
            upHyFxn_(std::get<0>(up), std::get<1>(up), By_, Ez_, Ex_);
        for(auto& up : upBz_)
            upHzFxn_(std::get<0>(up), std::get<1>(up), Bz_, Ex_, Ey_);
    }

    /**
     * @brief      Updates the D fields forward in time
     */
    void updateD()
    {
        for(auto& up : upDx_)
            upExFxn_(std::get<0>(up), std::get<1>(up), Dx_, Hy_, Hz_);
        for(auto& up : upDy_)
            upEyFxn_(std::get<0>(up), std::get<1>(up), Dy_, Hz_, Hx_);
        for(auto& up : upDz_)
            upEzFxn_(std::get<0>(up), std::get<1>(up), Dz_, Hx_, Hy_);
    }

    /**
     * @brief      Updates the polarization fields and adds them to the D field to get the E field for electrically dispersive materials
     */
    void updateDispE()
    {
        for(auto& up : upDx_)
        {
            upLorPxFxn_( std::get<0>(up), Ex_, lorPx_, prevLorPx_, scratch_.data(), objArr_[ std::get<0>(up)[7] ]);
            D2ExFxn_( std::get<0>(up), Dx_, Ex_, lorPx_, objArr_[ std::get<0>(up)[7] ]);
        }
        for(auto& up : upDy_)
        {
            upLorPyFxn_( std::get<0>(up), Ey_, lorPy_, prevLorPy_, scratch_.data(), objArr_[ std::get<0>(up)[7] ]);
            D2EyFxn_( std::get<0>(up), Dy_, Ey_, lorPy_, objArr_[ std::get<0>(up)[7] ]);
        }
        for(auto& up : upDz_)
        {
            upLorPzFxn_( std::get<0>(up), Ez_, lorPz_, prevLorPz_, scratch_.data(), objArr_[ std::get<0>(up)[7] ]);
            D2EzFxn_( std::get<0>(up), Dz_, Ez_, lorPz_, objArr_[ std::get<0>(up)[7] ]);
        }
    }

    /**
     * @brief      Updates the magnetization fields and adds them to the B field to get the H field for magnetically dispersive materials
     */
    void updateMagH()
    {
        for(auto& up : upBx_)
        {
            upLorMxFxn_( std::get<0>(up), Hx_, lorMx_, prevLorMx_, scratch_.data(), objArr_[ std::get<0>(up)[7] ]);
            B2HxFxn_( std::get<0>(up), Bx_, Hx_, lorMx_, objArr_[ std::get<0>(up)[7] ]);
        }
        for(auto& up : upBy_)
        {
            upLorMyFxn_( std::get<0>(up), Hy_, lorMy_, prevLorMy_, scratch_.data(), objArr_[ std::get<0>(up)[7] ]);
            B2HyFxn_( std::get<0>(up), By_, Hy_, lorMy_, objArr_[ std::get<0>(up)[7] ]);
        }
        for(auto& up : upBz_)
        {
            upLorMzFxn_( std::get<0>(up), Hz_, lorMz_, prevLorMz_, scratch_.data(), objArr_[ std::get<0>(up)[7] ]);
            B2HzFxn_( std::get<0>(up), Bz_, Hz_, lorMz_, objArr_[ std::get<0>(up)[7] ]);
        }
    }

    /**
     * @brief      Creates a BMP image for a slice of the input maps
     *
     * @param[in]  IP          Input parameters object used to create the FDTD propagator
     * @param[in]  sliceDir    Direction of the normal vector of the plane a slice is being taken of
     * @param[in]  sliceCoord  The value of the slice in real space along the normal direction to output the slice
     */
    void convertInputs2Map(const parallelProgramInputs &IP, DIRECTION sliceDir, double sliceCoord)
    {
        // Only do the output if on the first process
        if(gridComm_->rank() != 0)
            return;
        // iterator is what counts up to make each object have its value
        double iterator = 1;
        // will be set based on the direction if x normal ii=x, jj=y, kk=z
        int cor_ii = -1, cor_jj = -1, cor_kk = -1;
        std::string fname;
        // set the coordinate indexes needed to complete the operation
        if(sliceDir == DIRECTION::X)
        {
            cor_ii = 0;
            cor_jj = 1;
            cor_kk = 2;
            if(!Hz_ || !Ez_)
                throw std::logic_error("Slice in an YZ plane is not possible for a 2D calculation.");
            fname = "InputMap_YZ_plane_" + std::to_string(sliceCoord) + ".bmp";
        }
        else if(sliceDir == DIRECTION::Y)
        {
            cor_ii = 1;
            cor_jj = 0; // Easier to keep x as jj even though it should be kk for x
            cor_kk = 2;
            if(!Hz_ || !Ez_)
                throw std::logic_error("Slice in an XZ plane is not possible for a 2D calculation.");
            fname = "InputMap_XZ_plane_" + std::to_string(sliceCoord) + ".bmp";
        }
        else if(sliceDir == DIRECTION::Z)
        {
            cor_ii = 2;
            cor_jj = 0;
            cor_kk = 1;
            if( (!Hz_ || !Ez_) && sliceCoord != 0.0)
            {
                sliceCoord = 0.0;
            }
            fname = "InputMap_XY_plane_" + std::to_string(sliceCoord) + ".bmp";
        }
        else
            throw std::logic_error("Slice Direction must be X, Y, or Z");
        int sliceNum = (sliceCoord + IP.size_[cor_ii]/2.0 + 0.5*d_[cor_ii] )*IP.res_; // Do things in terms of grid points not actual values
        int map_nx = n_vec_[cor_jj];
        int map_ny = n_vec_[cor_kk];
        //construct map
        real_grid_ptr map = std::make_shared<Grid<double>>( std::array<int,3>({{map_nx, map_ny, 1}}) , std::array<double,3>({{ d_[cor_jj], d_[cor_kk], 1.0}}) );
        std::vector<double> ones(std::max(n_vec_[cor_jj], n_vec_[cor_kk]), 1);
        std::vector<double> includePML(ones.size(), 1.0);

        // PMLs regions initially set to 1
        for(int xx = 0; xx < IP.pmlThickness_[cor_jj]; ++xx)
        {
            dcopy_(map_ny, includePML.data(), 1, &map->point(xx             , 0 ), map->x() );
            dcopy_(map_ny, includePML.data(), 1, &map->point(map->x()-(1+xx), 0 ), map->x() );
        }
        for(int yy = 0; yy < IP.pmlThickness_[cor_kk]; ++yy)
        {
            dcopy_(map_nx, includePML.data(), 1, &map->point(0, yy              ), 1 );
            dcopy_(map_nx, includePML.data(), 1, &map->point(0, map->y()-(1+yy) ), 1 );
        }

        ++iterator;
        for(int dd = 0; dd < IP.dtcLoc_.size(); ++dd)
        {
            std::fill_n(ones.begin(), ones.size(), static_cast<double>(iterator) );
            if( sliceNum >= IP.dtcLoc_[dd][cor_ii] && sliceNum < IP.dtcSz_[dd][cor_ii] + IP.dtcLoc_[dd][cor_ii] )
            {
                if(IP.dtcSz_[dd][cor_jj] > 1)
                {
                    dcopy_(IP.dtcSz_[dd][cor_jj], ones.data(), 1, &map->point( IP.dtcLoc_[dd][cor_jj], IP.dtcLoc_[dd][cor_kk] ), 1);
                    if(IP.dtcSz_[dd][cor_kk] > 1)
                    {
                        dcopy_(IP.dtcSz_[dd][cor_jj]  , ones.data(), 1, &map->point( IP.dtcLoc_[dd][cor_jj]                          , IP.dtcLoc_[dd][cor_kk] + IP.dtcSz_[dd][cor_kk] - 1 ), 1);
                        dcopy_(IP.dtcSz_[dd][cor_kk]-2, ones.data(), 1, &map->point( IP.dtcLoc_[dd][cor_jj]                          , IP.dtcLoc_[dd][cor_kk] + 1)                         , map->x() );
                        dcopy_(IP.dtcSz_[dd][cor_kk]-2, ones.data(), 1, &map->point( IP.dtcLoc_[dd][cor_jj] + IP.dtcSz_[dd][cor_jj]-1, IP.dtcLoc_[dd][cor_kk] + 1)                         , map->x() );
                    }
                }
                else
                {
                    dcopy_(IP.dtcSz_[dd][cor_kk]  , ones.data(), 1, &map->point( IP.dtcLoc_[dd][cor_jj], IP.dtcLoc_[dd][cor_kk] ), map->x() );
                }
            }
            ++iterator;
        }
        // disregard the first object
        for(int tt = 0; tt < IP.tfsfLoc_.size(); ++tt)
        {
            std::fill_n(ones.begin(), ones.size(), static_cast<double>(iterator) );
            if( sliceNum >= IP.tfsfLoc_[tt][cor_ii] && sliceNum < IP.tfsfSize_[tt][cor_ii] + IP.tfsfLoc_[tt][cor_ii] )
            {
                if(IP.tfsfSize_[tt][cor_jj] > 1)
                {
                    dcopy_(IP.tfsfSize_[tt][cor_jj], ones.data(), 1, &map->point( IP.tfsfLoc_[tt][cor_jj], IP.tfsfLoc_[tt][cor_kk] ), 1);
                    if(IP.tfsfSize_[tt][cor_kk] > 1)
                    {
                        dcopy_(IP.tfsfSize_[tt][cor_jj]  , ones.data(), 1, &map->point( IP.tfsfLoc_[tt][cor_jj]                             , IP.tfsfLoc_[tt][cor_kk] + IP.tfsfSize_[tt][cor_kk] - 1 ), 1);
                        dcopy_(IP.tfsfSize_[tt][cor_kk]-2, ones.data(), 1, &map->point( IP.tfsfLoc_[tt][cor_jj]                             , IP.tfsfLoc_[tt][cor_kk] + 1)                            , map->x() );
                        dcopy_(IP.tfsfSize_[tt][cor_kk]-2, ones.data(), 1, &map->point( IP.tfsfLoc_[tt][cor_jj] + IP.tfsfSize_[tt][cor_jj]-1, IP.tfsfLoc_[tt][cor_kk] + 1)                            , map->x() );
                    }
                }
                else
                {
                    dcopy_(IP.tfsfSize_[tt][cor_kk]  , ones.data(), 1, &map->point( IP.tfsfLoc_[tt][cor_jj], IP.tfsfLoc_[tt][cor_kk] ), map->x() );
                }
            }
            ++iterator;
        }

        for(int ff = 0; ff < IP.fluxLoc_.size(); ++ff)
        {
            std::fill_n(ones.begin(), ones.size(), static_cast<double>(iterator) );
            if( sliceNum >= IP.fluxLoc_[ff][cor_ii] && sliceNum < IP.fluxSz_[ff][cor_ii] + IP.fluxLoc_[ff][cor_ii] )
            {
                if(IP.fluxSz_[ff][cor_jj] > 1)
                {
                    dcopy_(IP.fluxSz_[ff][cor_jj], ones.data(), 1, &map->point( IP.fluxLoc_[ff][cor_jj], IP.fluxLoc_[ff][cor_kk] ), 1);
                    if(IP.fluxSz_[ff][cor_kk] > 1)
                    {
                        dcopy_(IP.fluxSz_[ff][cor_jj]  , ones.data(), 1, &map->point( IP.fluxLoc_[ff][cor_jj]                          , IP.fluxLoc_[ff][cor_kk] + IP.fluxSz_[ff][cor_kk] - 1 ), 1);
                        dcopy_(IP.fluxSz_[ff][cor_kk]-2, ones.data(), 1, &map->point( IP.fluxLoc_[ff][cor_jj]                          , IP.fluxLoc_[ff][cor_kk] + 1)                         , map->x() );
                        dcopy_(IP.fluxSz_[ff][cor_kk]-2, ones.data(), 1, &map->point( IP.fluxLoc_[ff][cor_jj] + IP.fluxSz_[ff][cor_jj]-1, IP.fluxLoc_[ff][cor_kk] + 1)                         , map->x() );
                    }
                }
                else
                {
                    dcopy_(IP.fluxSz_[ff][cor_kk]  , ones.data(), 1, &map->point( IP.fluxLoc_[ff][cor_jj], IP.fluxLoc_[ff][cor_kk] ), map->x() );
                }
            }
            ++iterator;
        }

        for(int oo = 1; oo < objArr_.size(); ++oo)
        {
            // look at all local points only
            if(oo == 0 || objArr_[oo]->mat().size() > 1 || objArr_[oo]->epsInfty() != 1.0 || objArr_[oo]->magMat().size() > 1 || objArr_[oo]->muInfty() != 1.0 )
            {
                std::array<double,3>pt ={sliceCoord, sliceCoord, sliceCoord};
                for(int ii = 0; ii < map->x(); ++ii)
                {
                    for(int jj = 0; jj < map->y(); ++jj)
                    {
                        // split it up by component so you can see what goes where
                        pt[cor_jj] = ((ii)-(n_vec_[cor_jj]-1)/2.0)*d_[cor_jj];
                        pt[cor_kk] = ((jj)-(n_vec_[cor_kk]-1)/2.0)*d_[cor_kk];
                        if(objArr_[oo]->isObj(pt,d_[cor_jj])==true)
                            map->point(ii,jj) += static_cast<double>(iterator+oo)/3.0;

                        pt[cor_jj] = ((ii)+0.5-(n_vec_[cor_jj]-1)/2.0)*d_[cor_jj];
                        pt[cor_kk] = ((jj)-(n_vec_[cor_kk]-1)/2.0)*d_[cor_kk];
                        if(objArr_[oo]->isObj(pt,d_[cor_jj])==true)
                            map->point(ii,jj) += static_cast<double>(iterator+oo)/3.0;

                        pt[cor_jj] = ((ii)-(n_vec_[cor_jj]-1)/2.0)*d_[cor_jj];
                        pt[cor_kk] = ((jj)+0.5-(n_vec_[cor_kk]-1)/2.0)*d_[cor_kk];
                        if(objArr_[oo]->isObj(pt,d_[cor_jj])==true)
                            map->point(ii,jj) += static_cast<double>(iterator+oo)/3.0;
                    }
                }
            }
            else
            {
                for(int o1 = 1; o1 < objArr_.size(); ++o1)
                {
                    // If object is vacuum return to initial background values
                    std::array<double,3>pt ={sliceCoord, sliceCoord, sliceCoord};
                    for(int ii = 0; ii < map->x(); ++ii)
                    {
                        for(int jj = 0; jj < map->y(); ++jj)
                        {
                            pt[cor_jj] = ((ii)-(n_vec_[cor_jj]-1)/2.0)*d_[cor_jj];
                            pt[cor_kk] = ((jj)-(n_vec_[cor_kk]-1)/2.0)*d_[cor_kk];
                            if(oo != o1 && objArr_[oo]->isObj(pt,d_[cor_jj])==true && objArr_[o1]->isObj(pt,d_[cor_jj])==true)
                                map->point(ii,jj) -= static_cast<double>(iterator+o1)/3.0;

                            pt[cor_jj] = ((ii)+0.5-(n_vec_[cor_jj]-1)/2.0)*d_[cor_jj];
                            pt[cor_kk] = ((jj)-(n_vec_[cor_kk]-1)/2.0)*d_[cor_kk];
                            if(oo != o1 && objArr_[oo]->isObj(pt,d_[cor_jj])==true && objArr_[o1]->isObj(pt,d_[cor_jj])==true)
                                map->point(ii,jj) -= static_cast<double>(iterator+o1)/3.0;

                            pt[cor_jj] = ((ii)-(n_vec_[cor_jj]-1)/2.0)*d_[cor_jj];
                            pt[cor_kk] = ((jj)+0.5-(n_vec_[cor_kk]-1)/2.0)*d_[cor_kk];
                            if(oo != o1 && objArr_[oo]->isObj(pt,d_[cor_jj])==true && objArr_[o1]->isObj(pt,d_[cor_jj])==true)
                                map->point(ii,jj) -= static_cast<double>(iterator+o1)/3.0;
                        }
                    }
                }
            }
        }
        // OUtput to bmp
        GridToBitMap(map, fname, PLOTTYPE::MAG);
        return;
    }

    /**
     * @brief      Converts the double size vector into integers (number of grid points)
     *
     * @param[in]  size  The size of the field in real units
     *
     * @return     The size of the grid in number of grid points
     */
    inline std::array<int,3> toN_vec(std::array<double,3> size){ std::array<int,3> toRet; for(int ii = 0; ii < 3; ++ii) toRet[ii] = floor(res_ * size[ii] + 0.5) + 1; return toRet; }
    /**
     * @brief      Accessor function for dt_
     *
     * @return     dt_
     */
    inline double dt(){return dt_;}

    /**
     * @brief      Accessor function for E_incd_
     *
     * @return     E_incd_
     */
    inline std::vector<cplx>&    E_incd() {return    E_incd_;}
    /**
     * @brief      Accessor function for E_pl_incd_
     *
     * @return     E_pl_incd_
     */
    inline std::vector<cplx>& E_pl_incd() {return E_pl_incd_;}
    /**
     * @brief      Accessor function for H_incd_
     *
     * @return     H_incd_
     */
    inline std::vector<cplx>&    H_incd() {return    H_incd_;}
    /**
     * @brief      Accessor function for H_mn_incd_
     *
     * @return     H_mn_incd_
     */
    inline std::vector<cplx>& H_mn_incd() {return H_mn_incd_;}

    /**
     * @brief      Accessor function for fluxArr_
     *
     * @return     fluxArr_
     */
    inline std::vector<std::shared_ptr<parallelFluxDTC<T>>>& fluxArr(){return fluxArr_;}
    /**
     * @brief      Accessor function for dtcFreqArr_
     *
     * @return     dtcFreqArr_
     */
    inline std::vector<std::shared_ptr<parallelDetectorFREQ_Base<T>>>& dtcFreqArr() {return dtcFreqArr_;}
};

class parallelFDTDFieldReal : public parallelFDTDFieldBase<double>
{
public:
    /**
     * @brief      Constructs a FDTD Propagator class
     *
     * @param[in]  IP        Input parameter object that read in values from a json input file
     * @param[in]  gridComm  A shared_ptr to the MPI interface for the calculation
     */
    parallelFDTDFieldReal(parallelProgramInputs &IP, std::shared_ptr<mpiInterface> gridComm);


    /**
     * @brief      Constructs a DTC based off of the input parameters and puts it in the proper detector vector
     *
     * @param[in]  c             class type of the dtc (bin, bmp, cout, txt, freq)
     * @param[in]  grid          vector of the fields that need to be outputted
     * @param[in]  SI            true if outputting in SI units
     * @param[in]  loc           The location of the detectors lower left corner in grid points
     * @param[in]  sz            The size of the detector in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  fxn           Function used to modify base field data
     * @param[in]  txtType       if BMP what should be outputted to the text file
     * @param[in]  type          The type of the detector (Ex, Ey, Epow, etc)
     * @param[in]  freqList      The frequency list
     * @param[in]  timeInterval  The number of time steps per field output
     * @param[in]  a             unit length of the calculation
     * @param[in]  I0            unit current of the calculation
     * @param[in]  t_max         The time at the final time step
     */
    void coustructDTC(DTCCLASS c, std::vector<real_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, GRIDOUTFXN fxn, GRIDOUTTYPE txtType, DTCTYPE type, std::vector<double> freqList, double timeInterval, double a, double I0, double t_max);

};

class parallelFDTDFieldCplx : public parallelFDTDFieldBase<cplx>
{
public:
    /**
     * @brief      Constructs a FDTD Propagator class
     *
     * @param[in]  IP        Input parameter object that read in values from a json input file
     * @param[in]  gridComm  A shared_ptr to the MPI interface for the calculation
     */
    parallelFDTDFieldCplx(parallelProgramInputs &IP, std::shared_ptr<mpiInterface> gridComm);


    /**
     * @brief      Constructs a DTC based off of the input parameters and puts it in the proper detector vector
     *
     * @param[in]  c             class type of the dtc (bin, bmp, cout, txt, freq)
     * @param[in]  grid          vector of the fields that need to be outputted
     * @param[in]  SI            true if outputting in SI units
     * @param[in]  loc           The location of the detectors lower left corner in grid points
     * @param[in]  sz            The size of the detector in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  fxn           Function used to modify base field data
     * @param[in]  txtType       if BMP what should be outputted to the text file
     * @param[in]  type          The type of the detector (Ex, Ey, Epow, etc)
     * @param[in]  freqList      The frequency list
     * @param[in]  timeInterval  The number of time steps per field output
     * @param[in]  a             unit length of the calculation
     * @param[in]  I0            unit current of the calculation
     * @param[in]  t_max         The time at the final time step
     */
    void coustructDTC(DTCCLASS c, std::vector<cplx_pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, GRIDOUTFXN fxn, GRIDOUTTYPE txtType, DTCTYPE type, std::vector<double> freqList, double timeInterval, double a, double I0, double t_max);
};

#endif