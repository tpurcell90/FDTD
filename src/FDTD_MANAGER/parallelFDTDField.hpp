#ifndef PARALLEL_FDTD_FDTDFIELD
#define PARALLEL_FDTD_FDTDFIELD

#include <INPUTS/parallelInputs.hpp>
#include <DTC/parallelDTC_TXT.hpp>
#include <DTC/parallelDTC_COUT.hpp>
#include <DTC/parallelDTC_BIN.hpp>
#include <DTC/parallelDTC_FREQ.hpp>
#include <DTC/parallelFlux.hpp>
#include <SOURCE/parallelSourceNormal.hpp>
#include <SOURCE/parallelSourceOblique.hpp>
#include <SOURCE/parallelTFSF.hpp>
#include <PML/parallelPML.hpp>
#include <DTC/toBitMap.hpp>
#include <UTIL/FDTD_up_eq.hpp>

// #include <UTIL/FDTD_up_eq.hpp>
typedef cplx cplx;
typedef std::shared_ptr<std::vector<std::vector<std::array<double,5>>>> preConsts;

/**
 * @brief The main FDTD propagator class
 * @details This class generates a propagator class that will update all electromagnetic fields
 */
template <typename T> class parallelFDTDFieldBase
{
protected:
    typedef std::shared_ptr<parallelGrid<T>> pgrid_ptr;
    typedef std::shared_ptr<parallelCPML<T>> pml_ptr;
    mpiInterface & gridComm_; //!< mpi communicator for the propagator

    bool dielectricMatInPML_;

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

    upLists axHx_; //!< the list of parameters used to update the Hx field
    upLists axHy_; //!< the list of parameters used to update the Hy field
    upLists axHz_; //!< the list of parameters used to update the Hz field

    upLists axEx_; //!< the list of parameters used to update the Ex field
    upLists axEy_; //!< the list of parameters used to update the Ey field
    upLists axEz_; //!< the list of parameters used to update the Ez field

    upLists axDx_; //!< the list of parameters used to update the Dx field and polarization fields
    upLists axDy_; //!< the list of parameters used to update the Dy field and polarization fields
    upLists axDz_; //!< the list of parameters used to update the Dz field and polarization fields

    upLists axBx_; //!< the list of parameters used to update the Bx field and magnetization fields
    upLists axBy_; //!< the list of parameters used to update the By field and magnetization fields
    upLists axBz_; //!< the list of parameters used to update the Bz field and magnetization fields

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

    pgrid_ptr prevHx_; //!< parallel grid corresponding to the Hx field
    pgrid_ptr prevHy_; //!< parallel grid corresponding to the Hy field
    pgrid_ptr prevHz_; //!< parallel grid corresponding to the Hz field

    pgrid_ptr prevEx_; //!< parallel grid corresponding to the Ex field
    pgrid_ptr prevEy_; //!< parallel grid corresponding to the Ey field
    pgrid_ptr prevEz_; //!< parallel grid corresponding to the Ez field

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
     * @brief Constructor
     * @details Takes in the parallelProgramInputs structure and generates an FDTDField
     *
     * @param IP parallelProgramInputs
     * @param gridComm mpi communicator for the system
     */
    parallelFDTDFieldBase(parallelProgramInputs &IP, mpiInterface & gridComm) :
        gridComm_(gridComm),
        dielectricMatInPML_(false),
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
        dtcArr_.reserve( IP.dtcLoc_.size() );
        srcArr_.reserve( IP.srcLoc_.size() );
        tfsfArr_.reserve( IP.tfsfLoc_.size() );
        dtcFreqArr_.reserve( IP.dtcLoc_.size() );
        fluxArr_.reserve( IP.fluxLoc_.size() );

        E_incd_   .reserve( ceil(IP.tMax_ / dt_ ) + 1 );
        E_pl_incd_.reserve( ceil(IP.tMax_ / dt_ ) + 1 );
        H_incd_   .reserve( ceil(IP.tMax_ / dt_ ) + 1 );
        H_mn_incd_.reserve( ceil(IP.tMax_ / dt_ ) + 1 );
        setupWeightsGrid(IP);

        if(IP.size_[2] != 0)
        {
            Ex_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, false);
            Hz_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, true);
            Ey_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, true);
            Hx_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, true);
            Ez_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, false);
            Hy_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, false);

            bool disp = false;
            bool magnetic = false;
            for(auto& obj : objArr_)
            {
                if(obj->magMat().size() > 1)
                {
                    magnetic = true;
                }
                if(obj->mat().size() > 1)
                {
                    disp = true;
                }
            }
            if(disp)
            {
                Dx_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, false);
                Dy_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, true);
                Dz_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, false);
            }
            else
            {
                Dx_ = nullptr;
                Dy_ = nullptr;
                Dz_ = nullptr;
            }
            if(magnetic)
            {
                Bx_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, true);
                By_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, false);
                Bz_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, true);
            }
            else
            {
                Bx_ = nullptr;
                By_ = nullptr;
                Bz_ = nullptr;
            }

            ln_vec_[0] = Hz_->local_x()-2;
            ln_vec_[1] = Hz_->local_y()-2;
            ln_vec_[2] = Hz_->local_z()-2;
            pbcZMin_ = 1;
            pbcZMax_ = ln_vec_[2];
        }
        else if(IP.pol_ == POLARIZATION::HZ || IP.pol_ == POLARIZATION::EX || IP.pol_ == POLARIZATION::EY)
        {
            Ex_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(),n_vec_[1]+2*gridComm_.npY(), 1 }}),d_, false);
            Ey_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(),n_vec_[1]+2*gridComm_.npY(), 1 }}),d_, true);
            Ez_ = nullptr;


            Hx_ = nullptr;
            Hy_ = nullptr;
            Hz_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(),n_vec_[1]+2*gridComm_.npY(), 1 }}),d_, true);

            bool disp = false;
            bool magnetic = false;
            for(auto& obj : objArr_)
            {
                if(obj->magMat().size() > 1)
                {
                    magnetic = true;
                }
                if(obj->mat().size() > 1)
                {
                    disp = true;
                }
            }
            if(disp)
            {
                Dx_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, false);
                Dy_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, true);
            }
            else
            {
                Dx_ = nullptr;
                Dy_ = nullptr;
            }
            Dz_ = nullptr;
            Bx_ = nullptr;
            By_ = nullptr;
            if(magnetic)
            {
                Bz_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, true);
            }
            else
            {
                Bz_ = nullptr;
            }
            prevEx_ = nullptr;
            prevEy_ = nullptr;
            prevEz_ = nullptr;

            prevHx_ = nullptr;
            prevHy_ = nullptr;
            prevHz_ = nullptr;

            ln_vec_[0] = Hz_->local_x()-2;
            ln_vec_[1] = Hz_->local_y()-2;
            ln_vec_[2] = 1;
        }
        else
        {
            Ex_ = nullptr;
            Ey_ = nullptr;
            Ez_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(),n_vec_[1]+2*gridComm_.npY(), 1 }}),d_,false);

            Hx_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(),n_vec_[1]+2*gridComm_.npY(), 1 }}),d_, true);
            Hy_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(),n_vec_[1]+2*gridComm_.npY(), 1 }}),d_,false);
            Hz_ = nullptr;

            bool disp = false;
            bool magnetic = false;
            for(auto& obj : objArr_)
            {
                if(obj->magMat().size() > 1)
                {
                    magnetic = true;
                }
                if(obj->mat().size() > 1)
                {
                    disp = true;
                }
            }
            Dx_ = nullptr;
            Dy_ = nullptr;
            if(disp)
            {
                Dz_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, false);
            }
            else
            {
                Dz_ = nullptr;
            }
            if(magnetic)
            {
                Bx_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, true);
                By_ = std::make_shared<parallelGrid<T>> (gridComm_, IP.periodic_, weights_, std::array<int,3>({{ n_vec_[0]+2*gridComm_.npX(), n_vec_[1]+2*gridComm_.npY(),  n_vec_[2]+2*gridComm_.npZ() }}), d_, false);
            }
            else
            {
                Bx_ = nullptr;
                By_ = nullptr;
            }
            Bz_ = nullptr;

            prevEx_ = nullptr;
            prevEy_ = nullptr;
            prevEz_ = nullptr;

            prevHx_ = nullptr;
            prevHy_ = nullptr;
            prevHz_ = nullptr;

            ln_vec_[0] = Ez_->local_x()-2;
            ln_vec_[1] = Ez_->local_y()-2;
            ln_vec_[2] = 1;
        }
        for(auto& xx : IP.inputMapSlicesX_)
            convertInputs2Map(IP, DIRECTION::X, xx);
        for(auto& yy : IP.inputMapSlicesY_)
            convertInputs2Map(IP, DIRECTION::Y, yy);
        for(auto& zz : IP.inputMapSlicesZ_)
            convertInputs2Map(IP, DIRECTION::Z, zz);
    }
    /**
     * @brief Generates the the physical grids and then the parameter update lists
     * @details Uses the objArr_ information to determine which objects are at each grid point, and uses that to generate the parameter lists from there
     */
    void initializeLists(std::shared_ptr<parallelGridInt> phys_Ex, std::shared_ptr<parallelGridInt> phys_Ey, std::shared_ptr<parallelGridInt> phys_Ez, std::shared_ptr<parallelGridInt> phys_Hx, std::shared_ptr<parallelGridInt> phys_Hy, std::shared_ptr<parallelGridInt> phys_Hz)
    {
        std::tuple <int,int,int,int> BCs;
        if(Ex_)
            BCs = std::make_tuple(1, Ex_->local_x()-1, 1, Ex_->local_y()-1);
        else
            BCs = std::make_tuple(1, Ez_->local_x()-1, 1, Ez_->local_y()-1);
        int xmin, xmax, ymin, ymax, zmin, zmax;
        std::tie(xmin,xmax,ymin,ymax) = BCs;
        if(Ez_ && Hz_)
        {
            zmin = 1; zmax = Ex_->local_z()-1;
        }
        else
        {
            zmin = 0; zmax = 1;
        }
        if(Hz_)
        {
            // Get ax lists for the E and D fields via the functions
            std::vector<std::array<int, 5>> tempEx, tempDx;
            std::vector<std::array<int, 5>> tempEy, tempDy;
            std::vector<std::array<int, 5>> tempHz, tempBz;

            std::tie(tempEx, tempDx) = getAxLists(  true, {{ xmin, ymin, zmin }}, {{ xmax, ymax, zmax }}, phys_Ex);
            std::tie(tempEy, tempDy) = getAxLists(  true, {{ xmin, ymin, zmin }}, {{ xmax, ymax, zmax }}, phys_Ey);
            std::tie(tempHz, tempBz) = getAxLists( false, {{ xmin, ymin, zmin }}, {{ xmax, ymax, zmax }}, phys_Hz);

            axEx_.reserve( tempEx.size() );
            axEy_.reserve( tempEy.size() );
            axHz_.reserve( tempHz.size() );

            axDx_.reserve( tempDx.size() );
            axDy_.reserve( tempDy.size() );
            axBz_.reserve( tempBz.size() );

            for(auto & ax : tempEx)
            {
                double eps = objArr_[ax[4]]->epsInfty();
                // Ex does not update at right edge, if it contains that then cut it out, also calculate all pref prefactors needed
                if(ax[3] + ax[0] - 1 == n_vec_[0])
                    axEx_.push_back(std::make_pair(std::array<int,8>({ax[3]-1, ax[0], ax[1], ax[2], 0, -1, 0, ax[4]}), std::array<double,2>({1.0, -1.0*dt_/(eps*d_[1])}) ) );
                else
                    axEx_.push_back(std::make_pair(std::array<int,8>({ax[3]  , ax[0], ax[1], ax[2], 0, -1, 0, ax[4]}), std::array<double,2>({1.0, -1.0*dt_/(eps*d_[1])}) ) );
            }

            for(auto & ax : tempEy)
            {
                double eps = objArr_[ax[4]]->epsInfty();
                // Remove all top edge values for the y field/calculate all prefactors
                if( ax[1] + Ey_->procLoc()[1] != n_vec_[1] )
                    axEy_.push_back(std::make_pair(std::array<int,8>({ax[3], ax[0], ax[1], ax[2], 0, 0, -1, ax[4]}), std::array<double,2>({1.0, -1.0*dt_/(eps*d_[0])}) ) );
            }

            for(auto & ax : tempHz)
            {
                if( ax[1] + Hz_->procLoc()[1] != n_vec_[1] )
                {
                    double mu = objArr_[ax[4]]->muInfty();
                    if(ax[3] + ax[0] - 1 == n_vec_[0])
                        axHz_.push_back(std::make_pair(std::array<int,8>({ax[3]-1, ax[0], ax[1], ax[2], 1, 0, 0, ax[4]}) , std::array<double,2>({1.0, -1.0*dt_/(d_[0]*mu)}) ) );
                    else
                        axHz_.push_back(std::make_pair(std::array<int,8>({ax[3]  , ax[0], ax[1], ax[2], 1, 0, 0, ax[4]}) , std::array<double,2>({1.0, -1.0*dt_/(d_[0]*mu)}) ) );
                }
            }

            for(auto & ax : tempDx)
            {
                // Ex does not update at right edge, if it contains that then cut it out, also calculate all pref prefactors needed
                if(ax[3] + ax[0] - 1 == n_vec_[0])
                    axDx_.push_back(std::make_pair(std::array<int,8>({ax[3]-1, ax[0], ax[1], ax[2], 0, -1, 0, ax[4]}), std::array<double,2>({1.0, -1.0*dt_/d_[1]}) ) );
                else
                    axDx_.push_back(std::make_pair(std::array<int,8>({ax[3]  , ax[0], ax[1], ax[2], 0, -1, 0, ax[4]}), std::array<double,2>({1.0, -1.0*dt_/d_[1]}) ) );
            }
            for(auto & ax : tempDy)
            {
                // Remove all top edge values for the y field/calculate all prefactors
                if( ax[1] + Ey_->procLoc()[1] != n_vec_[1] )
                    axDy_.push_back(std::make_pair(std::array<int,8>({ax[3], ax[0], ax[1], ax[2], 0, 0, -1, ax[4]}), std::array<double,2>({1.0, -1.0*dt_/(d_[0])}) ) );
            }
            for(auto & ax : tempBz)
            {
                if( ax[1] + Hz_->procLoc()[1] != n_vec_[1] )
                {
                    if(ax[3] + ax[0] - 1 == n_vec_[0])
                        axBz_.push_back(std::make_pair(std::array<int,8>({ax[3]-1, ax[0], ax[1], ax[2], 1, 0, 0, ax[4]}) , std::array<double,2>({1.0, -1.0*dt_/(d_[0])}) ) );
                    else
                        axBz_.push_back(std::make_pair(std::array<int,8>({ax[3]  , ax[0], ax[1], ax[2], 1, 0, 0, ax[4]}) , std::array<double,2>({1.0, -1.0*dt_/(d_[0])}) ) );
                }
            }

        }
        if(Ez_)
        {
            std::vector<std::array<int, 5>> tempHx, tempBx;
            std::vector<std::array<int, 5>> tempHy, tempBy;
            std::vector<std::array<int, 5>> tempEz, tempDz;

            std::tie(tempHx, tempBx) = getAxLists( false, {{ xmin, ymin, zmin }}, {{ xmax, ymax, zmax }}, phys_Hx);
            std::tie(tempHy, tempBy) = getAxLists( false, {{ xmin, ymin, zmin }}, {{ xmax, ymax, zmax }}, phys_Hy);
            std::tie(tempEz, tempDz) = getAxLists(  true, {{ xmin, ymin, zmin }}, {{ xmax, ymax, zmax }}, phys_Ez);

            axHx_.reserve( tempHx.size() );
            axBx_.reserve( tempBx.size() );

            axHy_.reserve( tempHy.size() );
            axBy_.reserve( tempBy.size() );

            axEz_.reserve( tempEz.size() );
            axDz_.reserve( tempDz.size() );
            for(auto& ax : tempHx)
            {
                double mu = objArr_[ax[4]]->muInfty();
                if( ax[2] + Hy_->procLoc()[2] != n_vec_[2] && ax[1] + Ez_->procLoc()[1] != n_vec_[1] )
                    axHx_.push_back(std::make_pair(std::array<int,8>({ax[3], ax[0], ax[1], ax[2], 0, 1, 0, ax[4]} ), std::array<double,2>( { 1.0, -1.0*dt_/(mu*d_[0]) } ) ) );
            }
            // Hy has no right edge
            for(auto& ax : tempHy)
            {
                if( ax[2] + Hy_->procLoc()[2] != n_vec_[2] )
                {
                    double mu = objArr_[ax[4]]->muInfty();
                    if(ax[3] + ax[0] - 1 == n_vec_[0])
                        axHy_.push_back(std::make_pair(std::array<int,8>({ax[3]-1, ax[0], ax[1], ax[2], 0, 0, 1, ax[4]}) , std::array<double,2>( {1.0, -1.0*dt_/(mu*d_[0]) } ) ) );
                    else
                        axHy_.push_back(std::make_pair(std::array<int,8>({ax[3]  , ax[0], ax[1], ax[2], 0, 0, 1, ax[4]}) , std::array<double,2>( {1.0, -1.0*dt_/(mu*d_[0]) } ) ) );
                }
            }
            for(auto & ax : tempEz)
            {
                if( ax[2] + Ez_->procLoc()[2] != n_vec_[2] )
                {
                    double eps = objArr_[ax[4]]->epsInfty();
                    axEz_.push_back(std::make_pair(std::array<int,8>({ax[3], ax[0], ax[1], ax[2], -1, 0, 0, ax[4]}), std::array<double,2>({1.0, -1.0*dt_/(eps*d_[1])}) ) );
                }
            }

            for(auto& ax : tempBx)
            {
                if( ax[2] + Hx_->procLoc()[2] != n_vec_[2] && ax[1] + Hx_->procLoc()[1] != n_vec_[1] )
                    axBx_.push_back(std::make_pair(std::array<int,8>({ax[3], ax[0], ax[1], ax[2], 0, 1, 0, ax[4]}) , std::array<double,2>({1.0, -1.0*dt_/d_[0]}) ) );
            }
            // Hy has no right edge
            for(auto& ax : tempBy)
            {
                if( ax[2] + Hy_->procLoc()[2] != n_vec_[2] )
                {
                    if(ax[3] + ax[0] - 1 == n_vec_[0])
                        axBy_.push_back(std::make_pair(std::array<int,8>({ax[3]-1, ax[0], ax[1], ax[2], 0, 0, 1, ax[4]}) , std::array<double,2>({1.0, -1.0*dt_/d_[0]}) ) );
                    else
                        axBy_.push_back(std::make_pair(std::array<int,8>({ax[3]  , ax[0], ax[1], ax[2], 0, 0, 1, ax[4]}) , std::array<double,2>({1.0, -1.0*dt_/d_[0]}) ) );
                }
            }

            for(auto & ax : tempDz)
            {
                if( ax[2] + Ez_->procLoc()[2] != n_vec_[2] )
                    axDz_.push_back(std::make_pair(std::array<int,8>({ax[3], ax[0], ax[1], ax[2], -1, 0, 0, ax[4]}), std::array<double,2>({1.0, -1.0*dt_/(d_[1])}) ) );
            }
        }
    }

    /**
     * @brief      Constructs a DTC based off of the input parameters and puts it in the proper detector vector
     *
     * @param[in]  c             class type of the dtc (bin, bmp, cout, txt, freq)
     * @param[in]  fPow          True if outputting power
     * @param[in]  dtcNum        number of the detector
     * @param[in]  grid          vector of the fields that need to be outputted
     * @param[in]  SI            true if outputting in SI units
     * @param[in]  loc           The location of the detectors lower left corner in grid points
     * @param[in]  sz            The size of the detector in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  fxn           Function used to modify base field data
     * @param[in]  txtType       if BMP what should be outputted to the text file
     * @param[in]  type          The type of the detector (Ex, Ey, Epow, etc)
     * @param[in]  nfreq         The number of frequencies to keep track of if Freq detector
     * @param[in]  fcen          The center frequency of the freq range
     * @param[in]  fwidth        The frequency width of the detector
     * @param[in]  lamL          The left end point of the wavelength range
     * @param[in]  lamR          The right end point of the wavelength range
     * @param[in]  timeInterval  The number of time steps per field output
     * @param[in]  fristComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length of the calculation
     * @param[in]  I0            unit current of the calculation
     * @param[in]  t_max         The time at the final time step
     */
    virtual void coustructDTC(DTCCLASS c, bool fPow, int dtcNum, std::vector<pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, int nfreq, double fcen, double fwidth, double lamL, double lamR, double timeInterval, DIRECTION fristComp, double a, double I0, double t_max) = 0;

    /**
     * @brief      Constructs a frequency detector and returns a shared_ptr to it
     *
     * @param[in]  fPow          True if outputting power
     * @param[in]  dtcNum        number of the detector
     * @param      grid          vector of the fields that need to be outputted
     * @param[in]  SI            true if outputting in SI units
     * @param[in]  loc           The location of the detectors lower left corner in grid points
     * @param[in]  sz            The size of the detector in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  fxn           Function used to modify base field data
     * @param[in]  txtType       if BMP what should be outputted to the text file
     * @param[in]  type          The type of the detector (Ex, Ey, Epow, etc)
     * @param[in]  nfreq         The number of frequencies to keep track of if Freq detector
     * @param[in]  fcen          The center frequency of the freq range
     * @param[in]  fwidth        The frequency width of the detector
     * @param[in]  lamL          The left end point of the wavelength range
     * @param[in]  lamR          The right end point of the wavelength range
     * @param[in]  timeInterval  The number of time steps per field output
     * @param[in]  fristComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length of the calculation
     * @param[in]  I0            unit current of the calculation
     *
     * @return     A shared_ptr to the frequency detector
     */
    virtual std::shared_ptr<parallelDetectorFREQ_Base<T>> constructFreqDTC(bool fPow, int dtcNum, std::vector<pgrid_ptr> &grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, int nfreq, double fcen, double fwidth, double lamL, double lamR, double timeInterval, DIRECTION fristComp, double a, double I0) = 0;

    /**
     * @brief Creates geometry update lists for E and D fields
     * @details Creates the precurser lists that will be used to described the cell geometry  in the update parameter lists
     *
     * @param xmin minimum x grid point
     * @param xmax maximum x grid point
     * @param ymin minimum y grid point
     * @param ymax maximum y grid point
     * @param physGrid physical grid for the particular field
     * @return a pair of update preurser lists for the D and E fields
     */
    std::tuple< std::vector<std::array<int,5>>, std::vector<std::array<int,5>> > getAxLists(bool E, std::array<int,3> min, std::array<int,3> max, std::shared_ptr<parallelGridInt> physGrid)
    {
        std::vector<std::array<int,5>> axULists = {};
        std::vector<std::array<int,5>> axDLists = {};

        int PML_y_bot = 0;
        int PML_y_top = 0;
        int PML_x_right = 0;
        int PML_x_left = 0;
        int PML_z_back = 0;
        int PML_z_front = 0;
        if(dielectricMatInPML_ && Ez_ && Hz_)
        {
            PML_z_back  = EzPML_->lnz_back();
            PML_z_front = EzPML_->lnz_front();
            PML_y_bot   = ExPML_->lny_bot();
            PML_y_top   = ExPML_->lny_top();
            PML_x_right = EyPML_->lnx_right();
            PML_x_left  = EyPML_->lnx_left();
        }
        else if(dielectricMatInPML_ && Hz_)
        {
            PML_y_bot   = ExPML_->lny_bot();
            PML_y_top   = ExPML_->lny_top();
            PML_x_right = EyPML_->lnx_right();
            PML_x_left  = EyPML_->lnx_left();
        }
        else if(dielectricMatInPML_)
        {
            PML_y_bot   = EzPML_->lny_bot();
            PML_y_top   = EzPML_->lny_top();
            PML_x_right = EzPML_->lnx_right();
            PML_x_left  = EzPML_->lnx_left();
        }
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
                    while ( (ii < max[0]-1) && ( physGrid->point(ii,jj,kk) == physGrid->point(ii+1,jj,kk) ) ) ++ii;

                    axDLists.push_back({{ iistore,jj,kk,ii-iistore+1,static_cast<int>(physGrid->point(iistore,jj,kk))}});
                    ++ii;
                }
            }
        }
        // Centeral region only the left and right PMLs should be 100% inside D lists
        for(int jj = min[1]+PML_y_bot; jj < max[1]-PML_y_top; ++jj)
        {
            for(int kk = min[2]; kk < min[2] + PML_z_back; ++kk )
            {
                int ii = min[0];
                while(ii < max[0])
                {
                    int iistore = ii;
                    // Check if points are in same object
                    while ( (ii < max[0]-1) && ( physGrid->point(ii,jj,kk) == physGrid->point(ii+1,jj,kk) ) ) ++ii;
                    axDLists.push_back({{ iistore,jj,kk,ii-iistore+1,static_cast<int>(physGrid->point(iistore,jj,kk))}});
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
                    while ( (ii < max[0]-1) && ( physGrid->point(ii,jj,kk) == physGrid->point(ii+1,jj,kk) ) && (ii != PML_x_left ) && (ii != max[0]-PML_x_right ) ) ++ii;

                    // If the point is not in a PML or dispersive material put it in update E field; else update D
                    if( ( !E && objArr_[ physGrid->point(iistore,jj,kk) ]->magMat().size() <= 1 ) || (E && objArr_[ physGrid->point(iistore,jj,kk) ]->mat().size() <= 1 && ii > PML_x_left && ii <=  max[0]-PML_x_right ) )
                        axULists.push_back({{ iistore,jj,kk,ii-iistore+1,static_cast<int>(physGrid->point(iistore,jj,kk))}});
                    else
                        axDLists.push_back({{ iistore,jj,kk,ii-iistore+1,static_cast<int>(physGrid->point(iistore,jj,kk))}});
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
                    while ( (ii < max[0]-1) && ( physGrid->point(ii,jj,kk) == physGrid->point(ii+1,jj,kk) ) ) ++ii;
                    axDLists.push_back({{ iistore,jj,kk,ii-iistore+1,static_cast<int>(physGrid->point(iistore,jj,kk))}});
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
                    while ( (ii < max[0]-1) && ( physGrid->point(ii,jj,kk) == physGrid->point(ii+1,jj,kk) ) ) ++ii;
                    axDLists.push_back({{ iistore,jj,kk,ii-iistore+1,static_cast<int>(physGrid->point(iistore,jj,kk))}});
                    ++ii;
                }
            }
        }
        // std::cout << axULists.size() << '\t' << axDLists.size() << '\t' << axChiDLists.size() << std::endl;
        return std::make_tuple(axULists, axDLists);
    }
    /**
     * @brief Uses the objArr to generate the physical grids
     * @details goes through all the grid points and determines if there is an object is at the grid point and if so sets the point to that value
     *
     * @param physGrid physGrid to be set
     * @param offPt offset for all the coordinate
     * @param min minimum for all the coordinates
     * @param max maximum for all the coordinates
     */
    void setupPhysFields(std::shared_ptr<parallelGridInt> physGrid, std::array<double,3> offPt, std::array<int,3> min, std::array<int,3> max)
    {
        if( Ez_ && Hz_ )
        {
            // Test point used
            std::array<double,3> pt = {{ 0,0,0}};
            // Objects on the same points over write each other (last object made wins)
            for(int oo = 0; oo < objArr_.size(); ++oo)
            {
                // look at all local points only
                if(oo == 0 || objArr_[oo]->mat().size() > 2 || objArr_[oo]->epsInfty() != 1.0 || objArr_[oo]->magMat().size() > 2 || objArr_[oo]->muInfty() != 1.0 )
                {
                    for(int ii = min[0]; ii < max[0]; ++ii)
                    {
                        for(int jj = min[1]; jj < max[1]; ++jj)
                        {
                            for(int kk = min[2]; kk < max[2]; ++kk)
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
            for(int ii = min[0]-1; ii < max[0] + 1; ++ii)
            {
                for(int jj = min[1]-1; jj < max[1] + 1; ++jj)
                {
                    physGrid->point(ii, jj, 0     ) = -1;
                    physGrid->point(ii, jj, max[2]) = -1;
                }
            }
            for(int ii = min[0]-1; ii < max[0] + 1; ++ii)
            {
                for(int kk = min[2]-1; kk < max[2] + 1; ++kk)
                {
                    physGrid->point(ii, 0     , kk) = -1;
                    physGrid->point(ii, max[1], kk) = -1;
                }
            }
            for(int kk = min[2]-1; kk < max[2] + 1; ++kk)
            {
                for(int jj = min[1]-1; jj < max[1] + 1; ++jj)
                {
                    physGrid->point(0     , jj, kk) = -1;
                    physGrid->point(max[0], jj, kk) = -1;
                }
            }
        }
        else
        {
            // Test point used
            std::array<double,3> pt = {{(0.0, 0.0, 0.0)}};
            // Objects on the same points over write each other (last object made wins)
            for(int oo = 0; oo < objArr_.size(); ++oo)
            {
                // look at all local points only
                if(oo == 0 || objArr_[oo]->mat().size() > 2 || objArr_[oo]->epsInfty() != 1.0  || objArr_[oo]->magMat().size() > 2 || objArr_[oo]->muInfty() != 1.0)
                {
                    for(int ii = min[0]; ii < max[0]; ++ii)
                    {
                        for(int jj = min[1]; jj < max[1]; ++jj)
                        {
                            pt[0] = ((ii-1)+offPt[0]+physGrid->procLoc()[0]-(n_vec_[0]-1)/2.0)*d_[0];
                            pt[1] = ((jj-1)+offPt[1]+physGrid->procLoc()[1]-(n_vec_[1]-1)/2.0)*d_[1];
                            if(objArr_[oo]->isObj(pt,d_[0])==true)
                                physGrid->point(ii,jj) = oo;
                        }
                    }
                }
            }
            // All borders of between the processors have a -1 to indicate they are borders
            for(int ii = min[0]-1; ii < max[0] + 1; ii++)
            {
                physGrid->point(ii, 0)    = -1;
                physGrid->point(ii, max[1]) = -1;
            }
            for(int jj = min[1]-1; jj < max[1] + 1; ++jj)
            {
                physGrid->point(0, jj)    = -1;
                physGrid->point(max[0], jj) = -1;
            }
        }
    }

    /**
     * @brief      Is used to make the weights grid for dividing the FDTD fields evenly
     *
     * @param      IP    Input parameter from the constructor
     */
    void setupWeightsGrid(parallelProgramInputs &IP)
    {
        // test point will move across all grid points
        // Set weights based off of normal materials, number of axpy calls per time step (6 base + 3*n_lor_pol = 6 + 3*(objArr_[kk]->mat().size()-1)/3)*number of E fields
        // Separated by polarization since TE and TM grids look at different points in the same grid cell.
        if(n_vec_[2] > 1)
        {
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
            // Include weights for flux regions
            for(int ff = 0; ff < IP.fluxLoc_.size(); ++ff)
            {
                std::vector<double> fluxWeight( std::max( std::max(n_vec_[0], n_vec_[1]), n_vec_[2] ), IP.fluxNFreq_[ff] );
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
        else if(IP.pol_ == POLARIZATION::HZ || IP.pol_ == POLARIZATION::EX || IP.pol_ == POLARIZATION::EY)
        {

            std::array<double,3> pt = {{0.0, 0.0, 0.0}};
            weights_.push_back(std::make_shared<Grid<double>>(n_vec_, d_) );
            weights_.push_back(std::make_shared<Grid<double>>(n_vec_, d_) );
            // Objects on the same points over write each other (last object made wins)
            for(int oo = 0; oo < objArr_.size(); ++oo)
            {
                for(int ii = 0; ii < n_vec_[0]; ++ii)
                {
                    for(int jj = 0; jj < n_vec_[1]; ++jj)
                    {
                        // Ex points located at ii+1/2, jj
                        pt[0] = (ii-(n_vec_[0]-1)/2.0+0.5)*d_[0];
                        pt[1] = (jj-(n_vec_[1]-1)/2.0    )*d_[1];
                        if(objArr_[oo]->isObj(pt,d_[0])==true)
                            weights_[0]->point(ii,jj) =  (2.0 + (2.0666666666666*static_cast<double>(objArr_[oo]->mat().size()-1) + 2.0666666666666*static_cast<double>(objArr_[oo]->magMat().size()-1) ) / 2.0 ) ; // / objArr_[oo]->geo()[0] ;
                        // Ey points located ii, jj+1/2
                        pt[1] += 0.5*d_[1];
                        pt[0] -= 0.5*d_[0];
                        if(objArr_[oo]->isObj(pt,d_[0])==true)
                            weights_[1]->point(ii,jj) =  (2.0 + (2.0666666666666*static_cast<double>(objArr_[oo]->mat().size()-1) + 2.0666666666666*static_cast<double>(objArr_[oo]->magMat().size()-1) ) / 2.0 ) ; // / objArr_[oo]->geo()[0] ;
                    }
                }
            }
            // Ey has no points across the top edge
            dcopy_(n_vec_[0], std::vector<double>(n_vec_[0],0.0).data(), 1, &weights_[1]->point(0, n_vec_[1]-1), 1);
            // Ex has no points across right edge
            dcopy_(n_vec_[1], std::vector<double>(n_vec_[1],0.0).data(), 1, &weights_[0]->point(n_vec_[0]-1, 0), weights_[0]->x() );

            // Add PMLs
            for(int xx = 0; xx < IP.pmlThickness_[0]; ++xx)
            {
                // Ex field has no PML along the left and right
                daxpy_(weights_[1]->y(), 1.0, std::vector<double>(n_vec_[1], 3.1).data(), 1, &weights_[1]->point(      xx, 0), weights_[1]->x());
                daxpy_(weights_[1]->y(), 1.0, std::vector<double>(n_vec_[1], 3.1).data(), 1, &weights_[1]->point(n_vec_[0]-1-xx, 0), weights_[1]->x());
            }
            for(int yy = 0; yy < IP.pmlThickness_[1]; ++yy)
            {
                // Ey field has no PML along the top and bottom
                daxpy_(weights_[0]->x(), 1.0, std::vector<double>(n_vec_[0], 3.1).data(), 1, &weights_[0]->point(0,       yy), 1);
                daxpy_(weights_[0]->x(), 1.0, std::vector<double>(n_vec_[0], 3.1).data(), 1, &weights_[0]->point(0, n_vec_[1]-1-yy), 1);
            }
            // Hz has no points across the top and right edge
            // dcopy_(n_vec_[0], std::vector<double>(n_vec_[0],0.0).data(), 1, &weights_[2]->point(0, n_vec_[1]-1), 1);
            // dcopy_(n_vec_[1], std::vector<double>(n_vec_[1],0.0).data(), 1, &weights_[2]->point(n_vec_[0]-1, 0), weights_[1]->x());
        }
        else
        {
            std::array<double,3> pt = {{0.0, 0.0, 0.0}};
            weights_.push_back(std::make_shared<Grid<double>>(n_vec_, d_) );
            for(int oo = 0; oo < objArr_.size(); ++oo)
            {
                for(int ii = 0; ii < n_vec_[0]; ++ii)
                {
                    for(int jj = 0; jj < n_vec_[1]; ++jj)
                    {
                        // Ez point is at ii, jj
                        pt[0] = (ii-(n_vec_[0]-1)/2.0 )*d_[0];
                        pt[1] = (jj-(n_vec_[1]-1)/2.0 )*d_[1];
                        if(objArr_[oo]->isObj(pt,d_[0])==true)
                            weights_[0]->point(ii,jj) =  (2.0 + ( 2.0666666666666*static_cast<double>(objArr_[oo]->mat().size()-1) + 2.0666666666666*static_cast<double>(objArr_[oo]->magMat().size()-1) ) / 2.0 ) ; // / objArr_[oo]->geo()[0] ;
                    }
                }
            }
            // Add PMLs
            for(int xx = 0; xx < IP.pmlThickness_[0]; ++xx)
            {
                daxpy_(weights_[0]->y(), 1.0, std::vector<double>(n_vec_[1], 3.1).data(), 1, &weights_[0]->point(      xx, 0), weights_[0]->x());
                daxpy_(weights_[0]->y(), 1.0, std::vector<double>(n_vec_[1], 3.1).data(), 1, &weights_[0]->point(n_vec_[0]-1-xx, 0), weights_[0]->x());
            }
            for(int yy = 0; yy < IP.pmlThickness_[1]; ++yy)
            {
                daxpy_(weights_[0]->x(), 1.0, std::vector<double>(n_vec_[0], 3.1).data(), 1, &weights_[0]->point(0,       yy), 1);
                daxpy_(weights_[0]->x(), 1.0, std::vector<double>(n_vec_[0], 3.1).data(), 1, &weights_[0]->point(0, n_vec_[1]-1-yy), 1);
            }
            // Hx has no points across top edge
            // dcopy_(n_vec_[0], std::vector<double>(n_vec_[0],0.0).data(), 1, &weights_[0]->point(0, n_vec_[1]-1), 1);
            // Hy has no points across right edge
            // dcopy_(n_vec_[1], std::vector<double>(n_vec_[1],0.0).data(), 1, &weights_[1]->point(n_vec_[0]-1, 0), weights_[1]->x());
        }
        // Add time for flux regions
    }
    /**
     * @brief return the current time
     * @return tcur_
     */
    inline double & getTime(){return tcur_;}

    /**
     * @brief updates the grids one step in time
     * @details Runs each update function and propagates grids forward in time one step
     */
    void step()
    {
        // Update H fields
        updateH();
        updateMagH();

        updateHxPML_(HxPML_);
        updateHyPML_(HyPML_);
        updateHzPML_(HzPML_);
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

        pbcHx_(Hx_, k_point_, ln_vec_[0]+1, ln_vec_[1]+1, pbcZMax_+1, ln_vec_[0]+1, yHxPBC_, pbcZMin_, pbcZMax_  , d_[0], d_[1], d_[2] );
        pbcHy_(Hy_, k_point_, ln_vec_[0]+1, ln_vec_[1]+1, pbcZMax_+1, ln_vec_[0]  , yHyPBC_, pbcZMin_, pbcZMax_  , d_[0], d_[1], d_[2] );
        pbcHz_(Hz_, k_point_, ln_vec_[0]+1, ln_vec_[1]+1, pbcZMax_+1, ln_vec_[0]  , yHzPBC_, pbcZMin_, pbcZMax_+1, d_[0], d_[1], d_[2] );

        transferHx_();
        transferHy_();
        transferHz_();

        // Update E fields
        updateE();
        // Add dispersive materials and update E field PML's during upD function calls (done to include dispersive materials more easily)
        updateDispE();

        updateExPML_(ExPML_);
        updateEyPML_(EyPML_);
        updateEzPML_(EzPML_);

        // All E-field Updates should now be completed transfer border values for the E-fields
        pbcEx_(Ex_, k_point_, ln_vec_[0]+1, ln_vec_[1]+1, pbcZMax_+1, ln_vec_[0]  , yExPBC_, pbcZMin_, pbcZMax_+1, d_[0], d_[1], d_[2]);
        pbcEy_(Ey_, k_point_, ln_vec_[0]+1, ln_vec_[1]+1, pbcZMax_+1, ln_vec_[0]+1, yEyPBC_, pbcZMin_, pbcZMax_+1, d_[0], d_[1], d_[2]);
        pbcEz_(Ez_, k_point_, ln_vec_[0]+1, ln_vec_[1]+1, pbcZMax_+1, ln_vec_[0]+1, yEzPBC_, pbcZMin_, pbcZMax_  , d_[0], d_[1], d_[2]);

        transferEx_();
        transferEy_();
        transferEz_();

        // Increment time steps to before output as all fields should be updated to the next time step now
        tcur_ += dt_;
        t_step_ ++;

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
     * @brief Updates the magnetic fields
     * @details Propagates the magnetic fields forward in time using the magnetic field functors
     */
    void updateH()
    {
        for(auto& ax : axHx_)
            upHxFxn_(std::get<0>(ax), std::get<1>(ax), Hx_, Ey_, Ez_);
        for(auto& ax : axHy_)
            upHyFxn_(std::get<0>(ax), std::get<1>(ax), Hy_, Ez_, Ex_);
        for(auto& ax : axHz_)
            upHzFxn_(std::get<0>(ax), std::get<1>(ax), Hz_, Ex_, Ey_);
    }
    /**
     * @brief Updates the electric fields forward in time
     * @details Propagates the electric fields forward in time using the electric field functors
     */
    void updateE()
    {
        for(auto& ax : axEx_)
            upExFxn_(std::get<0>(ax), std::get<1>(ax), Ex_, Hy_, Hz_);
        for(auto& ax : axEy_)
            upEyFxn_(std::get<0>(ax), std::get<1>(ax), Ey_, Hz_, Hx_);
        for(auto& ax : axEz_)
            upEzFxn_(std::get<0>(ax), std::get<1>(ax), Ez_, Hx_, Hy_);
    }
    /**
     * @brief Updates the electric fields for dispersive materials
     * @details Updates he electric fields for the dispersive materials using the appropriate functors: Start off with updating the Lorentzian Polarization fields and the D fields (should be independent of each other), then combine in the E-fields; E = $\eps_\infty \vec{D} + \sum \vec{P}
     */
    void updateDispE()
    {
        for(auto& ax : axDx_)
        {
            upExFxn_(std::get<0>(ax), std::get<1>(ax), Dx_, Hy_, Hz_);
            upLorPxFxn_( std::get<0>(ax), Ex_, lorPx_, prevLorPx_, scratch_.data(), objArr_[ std::get<0>(ax)[7] ]);
            D2ExFxn_( std::get<0>(ax), Dx_, Ex_, lorPx_, objArr_[ std::get<0>(ax)[7] ]);
        }
        for(auto& ax : axDy_)
        {
            upEyFxn_(std::get<0>(ax), std::get<1>(ax), Dy_, Hz_, Hx_);
            upLorPyFxn_( std::get<0>(ax), Ey_, lorPy_, prevLorPy_, scratch_.data(), objArr_[ std::get<0>(ax)[7] ]);
            D2EyFxn_( std::get<0>(ax), Dy_, Ey_, lorPy_, objArr_[ std::get<0>(ax)[7] ]);
        }
        for(auto& ax : axDz_)
        {
            upEzFxn_(std::get<0>(ax), std::get<1>(ax), Dz_, Hx_, Hy_);
            upLorPzFxn_( std::get<0>(ax), Ez_, lorPz_, prevLorPz_, scratch_.data(), objArr_[ std::get<0>(ax)[7] ]);
            D2EzFxn_( std::get<0>(ax), Dz_, Ez_, lorPz_, objArr_[ std::get<0>(ax)[7] ]);
        }
    }

    /**
     * @brief Updates the magnetic fields for magnetically dispersive materials
     */
    void updateMagH()
    {
        for(auto& ax : axBx_)
        {
            upHxFxn_(std::get<0>(ax), std::get<1>(ax), Bx_, Ey_, Ez_);
            upLorMxFxn_( std::get<0>(ax), Hx_, lorMx_, prevLorMx_, scratch_.data(), objArr_[ std::get<0>(ax)[7] ]);
            B2HxFxn_( std::get<0>(ax), Bx_, Hx_, lorMx_, objArr_[ std::get<0>(ax)[7] ]);
        }
        for(auto& ax : axBy_)
        {
            upHyFxn_(std::get<0>(ax), std::get<1>(ax), By_, Ez_, Ex_);
            upLorMyFxn_( std::get<0>(ax), Hy_, lorMy_, prevLorMy_, scratch_.data(), objArr_[ std::get<0>(ax)[7] ]);
            B2HyFxn_( std::get<0>(ax), By_, Hy_, lorMy_, objArr_[ std::get<0>(ax)[7] ]);
        }
        for(auto& ax : axBz_)
        {
            upHzFxn_(std::get<0>(ax), std::get<1>(ax), Bz_, Ex_, Ey_);
            upLorMzFxn_( std::get<0>(ax), Hz_, lorMz_, prevLorMz_, scratch_.data(), objArr_[ std::get<0>(ax)[7] ]);
            B2HzFxn_( std::get<0>(ax), Bz_, Hz_, lorMz_, objArr_[ std::get<0>(ax)[7] ]);
        }
    }

    void convertInputs2Map(parallelProgramInputs &IP, DIRECTION sliceDir, double sliceCoord)
    {
        // Print out XY-Plane
        if(gridComm_.rank() != 0)
            return;
        double iterator = 1;
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
                std::cout << "WARNING: Setting sliceCoord to 0 since it is a 2D calculation." << std::endl;
                sliceCoord = 0.0;
            }
            fname = "InputMap_XY_plane_" + std::to_string(sliceCoord) + ".bmp";
        }
        else
            throw std::logic_error("Slice Direction must be X, Y, or Z");
        int sliceNum = (sliceCoord + IP.size_[cor_ii]/2.0 + 0.5*d_[cor_ii] )*IP.res_; // Do things in terms of grid points not actual values
        int map_nx = n_vec_[cor_jj];
        int map_ny = n_vec_[cor_kk];
        real_grid_ptr map = std::make_shared<Grid<double>>( std::array<int,3>({{map_nx, map_ny, 1}}) , std::array<double,3>({{ d_[cor_jj], d_[cor_kk], 1.0}}) );
        std::vector<double> ones(std::max(n_vec_[cor_jj], n_vec_[cor_kk]), 1);
        std::vector<double> includePML(ones.size(), 1.0);
        // PMLs are ones (basically background, but different to tell where they are)
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

        iterator++;
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
            if(oo == 0 || objArr_[oo]->mat().size() > 2 || objArr_[oo]->epsInfty() != 1.0 || objArr_[oo]->magMat().size() > 2 || objArr_[oo]->muInfty() != 1.0)
            {
                std::array<double,3>pt ={sliceCoord, sliceCoord, sliceCoord};
                for(int ii = 0; ii < map->x(); ++ii)
                {
                    for(int jj = 0; jj < map->y(); ++jj)
                    {
                        pt[cor_jj] = ((ii)-(n_vec_[cor_jj]-1)/2.0)*d_[cor_jj];
                        pt[cor_kk] = ((jj)-(n_vec_[cor_kk]-1)/2.0)*d_[cor_kk];
                        if(objArr_[oo]->isObj(pt,d_[cor_jj])==true)
                            map->point(ii,jj) += static_cast<double>(iterator)/3.0;

                        pt[cor_jj] = ((ii)+0.5-(n_vec_[cor_jj]-1)/2.0)*d_[cor_jj];
                        pt[cor_kk] = ((jj)-(n_vec_[cor_kk]-1)/2.0)*d_[cor_kk];
                        if(objArr_[oo]->isObj(pt,d_[cor_jj])==true)
                            map->point(ii,jj) += static_cast<double>(iterator)/3.0;

                        pt[cor_jj] = ((ii)-(n_vec_[cor_jj]-1)/2.0)*d_[cor_jj];
                        pt[cor_kk] = ((jj)+0.5-(n_vec_[cor_kk]-1)/2.0)*d_[cor_kk];
                        if(objArr_[oo]->isObj(pt,d_[cor_jj])==true)
                            map->point(ii,jj) += static_cast<double>(iterator)/3.0;
                    }
                }
            }
            ++iterator;
        }
        // to get to the next value do size of objArr +2

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
     * @brief Constructor
     * @details Takes in the parallelProgramInputs structure and generates an FDTDField
     *
     * @param IP parallelProgramInputs
     * @param gridComm mpi communicator for the system
     */
    parallelFDTDFieldReal(parallelProgramInputs &IP, mpiInterface & gridComm);

    /**
     * @brief      Constructs a frequency detector and returns a shared_ptr to it
     *
     * @param[in]  fPow          True if outputting power
     * @param[in]  dtcNum        number of the detector
     * @param      grid          vector of the fields that need to be outputted
     * @param[in]  SI            true if outputting in SI units
     * @param[in]  loc           The location of the detectors lower left corner in grid points
     * @param[in]  sz            The size of the detector in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  fxn           Function used to modify base field data
     * @param[in]  txtType       if BMP what should be outputted to the text file
     * @param[in]  type          The type of the detector (Ex, Ey, Epow, etc)
     * @param[in]  nfreq         The number of frequencies to keep track of if Freq detector
     * @param[in]  fcen          The center frequency of the freq range
     * @param[in]  fwidth        The frequency width of the detector
     * @param[in]  lamL          The left end point of the wavelength range
     * @param[in]  lamR          The right end point of the wavelength range
     * @param[in]  timeInterval  The number of time steps per field output
     * @param[in]  fristComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length of the calculation
     * @param[in]  I0            unit current of the calculation
     *
     * @return     A shared_ptr to the frequency detector
     */
    std::shared_ptr<parallelDetectorFREQ_Base<double>> constructFreqDTC(bool fPow, int dtcNum, std::vector<pgrid_ptr>& grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, int nfreq, double fcen, double fwidth, double lamL, double lamR, double timeInterval, DIRECTION fristComp, double a, double I0);

    /**
     * @brief      Constructs a DTC based off of the input parameters and puts it in the proper detector vector
     *
     * @param[in]  c             class type of the dtc (bin, bmp, cout, txt, freq)
     * @param[in]  fPow          True if outputting power
     * @param[in]  dtcNum        number of the detector
     * @param[in]  grid          vector of the fields that need to be outputted
     * @param[in]  SI            true if outputting in SI units
     * @param[in]  loc           The location of the detectors lower left corner in grid points
     * @param[in]  sz            The size of the detector in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  fxn           Function used to modify base field data
     * @param[in]  txtType       if BMP what should be outputted to the text file
     * @param[in]  type          The type of the detector (Ex, Ey, Epow, etc)
     * @param[in]  nfreq         The number of frequencies to keep track of if Freq detector
     * @param[in]  fcen          The center frequency of the freq range
     * @param[in]  fwidth        The frequency width of the detector
     * @param[in]  lamL          The left end point of the wavelength range
     * @param[in]  lamR          The right end point of the wavelength range
     * @param[in]  timeInterval  The number of time steps per field output
     * @param[in]  fristComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length of the calculation
     * @param[in]  I0            unit current of the calculation
     * @param[in]  t_max         The time at the final time step
     */
    void coustructDTC(DTCCLASS c, bool fPow, int dtcNum, std::vector<pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, int nfreq, double fcen, double fwidth, double lamL, double lamR, double timeInterval, DIRECTION fristComp, double a, double I0, double t_max);
};

class parallelFDTDFieldCplx : public parallelFDTDFieldBase<cplx>
{
public:
    /**
     * @brief Constructor
     * @details Takes in the parallelProgramInputs structure and generates an FDTDField
     *
     * @param IP parallelProgramInputs
     * @param gridComm mpi communicator for the system
     */
    parallelFDTDFieldCplx(parallelProgramInputs &IP, mpiInterface & gridComm);

    /**
     * @brief      Constructs a frequency detector and returns a shared_ptr to it
     *
     * @param[in]  fPow          True if outputting power
     * @param[in]  dtcNum        number of the detector
     * @param      grid          vector of the fields that need to be outputted
     * @param[in]  SI            true if outputting in SI units
     * @param[in]  loc           The location of the detectors lower left corner in grid points
     * @param[in]  sz            The size of the detector in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  fxn           Function used to modify base field data
     * @param[in]  txtType       if BMP what should be outputted to the text file
     * @param[in]  type          The type of the detector (Ex, Ey, Epow, etc)
     * @param[in]  nfreq         The number of frequencies to keep track of if Freq detector
     * @param[in]  fcen          The center frequency of the freq range
     * @param[in]  fwidth        The frequency width of the detector
     * @param[in]  lamL          The left end point of the wavelength range
     * @param[in]  lamR          The right end point of the wavelength range
     * @param[in]  timeInterval  The number of time steps per field output
     * @param[in]  fristComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length of the calculation
     * @param[in]  I0            unit current of the calculation
     *
     * @return     A shared_ptr to the frequency detector
     */
    std::shared_ptr<parallelDetectorFREQ_Base<cplx>> constructFreqDTC(bool fPow, int dtcNum, std::vector<pgrid_ptr>& grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, int nfreq, double fcen, double fwidth, double lamL, double lamR, double timeInterval, DIRECTION fristComp, double a, double I0);

    /**
     * @brief      Constructs a DTC based off of the input parameters and puts it in the proper detector vector
     *
     * @param[in]  c             class type of the dtc (bin, bmp, cout, txt, freq)
     * @param[in]  fPow          True if outputting power
     * @param[in]  dtcNum        number of the detector
     * @param[in]  grid          vector of the fields that need to be outputted
     * @param[in]  SI            true if outputting in SI units
     * @param[in]  loc           The location of the detectors lower left corner in grid points
     * @param[in]  sz            The size of the detector in grid points
     * @param[in]  out_name      The output file name
     * @param[in]  fxn           Function used to modify base field data
     * @param[in]  txtType       if BMP what should be outputted to the text file
     * @param[in]  type          The type of the detector (Ex, Ey, Epow, etc)
     * @param[in]  nfreq         The number of frequencies to keep track of if Freq detector
     * @param[in]  fcen          The center frequency of the freq range
     * @param[in]  fwidth        The frequency width of the detector
     * @param[in]  lamL          The left end point of the wavelength range
     * @param[in]  lamR          The right end point of the wavelength range
     * @param[in]  timeInterval  The number of time steps per field output
     * @param[in]  fristComp     The first component: is loc[0] in the x y or z direction
     * @param[in]  a             unit length of the calculation
     * @param[in]  I0            unit current of the calculation
     * @param[in]  t_max         The time at the final time step
     */
    void coustructDTC(DTCCLASS c, bool fPow, int dtcNum, std::vector<pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, DTCTYPE type, int nfreq, double fcen, double fwidth, double lamL, double lamR, double timeInterval, DIRECTION fristComp, double a, double I0, double t_max);
};

#endif