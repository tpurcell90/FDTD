#ifndef FDTD_TFSF
#define FDTD_TFSF

#include <SOURCE/Pulse.hpp>
#include <unordered_map>

struct paramStoreTFSF
{
    int incdStart_j_; //!< incident start for the j grid
    int incdStart_k_; //!< incident start for the k grid
    int strideField_; //1< stride for the real grids
    int strideIncd_; //!< stride for the incident grids
    int addIncdProp_; //!< value added to incd start as looping over transCor2
    cplx prefactor_j_; //!< prefactor for the j grid
    cplx prefactor_k_; //!< prefactor for the k grid
    std::array<int,2> szTrans_j_; //!< size of the transverse directions the j grid
    std::array<int,2> szTrans_k_; //!< size of the transverse directions the k grid
    std::array<int,3> loc_; //!< location of the lower left corner of the tfsf region
    std::array<int,3> addVec_; //!< what to add when looping over transCor2
};

/**
 * @brief TFSF region for introducing plane waves into the system
 */
template <typename T> class parallelTFSFBase
{
protected:
    std::shared_ptr<mpiInterface> gridComm_; //!< Communicator for MPI calls in the grid
    int gridLen_; //!< Length of the 1D auxiliary grid for the incident fields
    double phi_; //!< angle of plane wave's k-vector in the xy plane
    double theta_; //!< polar angle of the plane wave's k-vector
    double psi_; //!< angle describing the polarization of the plane wave from the vector \vec{k}\times e_{z} (if k is along z psi_ = phi_)
    double alpha_; //!< phase difference for circularly polarized light
    double phiPrefactCalc_; //!< phi used in the calculation of the prefactors (may change if 3D calculation)
    double psiPrefactCalc_; //!< psi used in the calculation of the prefactors (may change if chiral light)
    double dt_; //!< time step
    double dx_; //!< step size of 1D auxillary grid
    double t_step_; //!< current integer time step
    std::array<int, 3> sz_; //!< size of the total field region
    std::array<int, 3> loc_; //!< location of the tfsf origin (a corner)
    std::vector<double> cdd_; //!< PML prefactor D-field's effect on itself
    std::vector<double> cdh_; //!< PML prefactor  H-field's effect on D
    std::vector<double> chh_; //!< PML prefactor H-field's effect on itself
    std::vector<double> chb0_; //!< PML prefactor prevB field's effect on H
    std::vector<double> chb1_; //!< PML prefactor B field's effect on H
    std::vector<std::shared_ptr<PulseBase>> pul_; //!< PulseBase of incident wave
    std::vector<cplx> D_old_; //!< old D field values
    std::vector<cplx> B_old_; //1< old B field values

    std::shared_ptr<parallelGrid<T>> Ex_; //!< grid_ptr to the Ex field
    std::shared_ptr<parallelGrid<T>> Ey_; //!< grid_ptr to the Ey field
    std::shared_ptr<parallelGrid<T>> Ez_; //!< grid_ptr to the Ez field
    std::shared_ptr<parallelGrid<T>> Hx_; //!< grid_ptr to the Hx field
    std::shared_ptr<parallelGrid<T>> Hy_; //!< grid_ptr to the Hy field
    std::shared_ptr<parallelGrid<T>> Hz_; //!< grid_ptr to the Hz field

    std::shared_ptr<Grid<cplx>> E_incd_; //!< TFSF incident E field
    std::shared_ptr<Grid<cplx>> H_incd_; //!< TFSF incident H field
    std::shared_ptr<Grid<cplx>> D_incd_; //!< TFSF incident D field
    std::shared_ptr<Grid<cplx>> B_incd_; //!< TFSF incident B field

    std::vector<cplx> incdTransfer_; //!< scartch space to transfer a cplx pulse with cplx prefactor to a real grid
    int originQuadrent_; //!< which quadrant the origin is in 1 is bot left increases in a counter clockwise direction

    std::shared_ptr<paramStoreTFSF> botSurE_; //!< parameter structure describing parameters needed to add the TFSF incident E field to the bot surface of the TFSF region
    std::shared_ptr<paramStoreTFSF> botSurH_; //!< parameter structure describing parameters needed to add the TFSF incident H field to the bot surface of the TFSF region

    std::shared_ptr<paramStoreTFSF> topSurE_; //!< parameter structure describing parameters needed to add the TFSF incident E field to the top surface of the TFSF region
    std::shared_ptr<paramStoreTFSF> topSurH_; //!< parameter structure describing parameters needed to add the TFSF incident H field to the top surface of the TFSF region

    std::shared_ptr<paramStoreTFSF> leftSurE_; //!< parameter structure describing parameters needed to add the TFSF incident E field to the left surface of the TFSF region
    std::shared_ptr<paramStoreTFSF> leftSurH_; //!< parameter structure describing parameters needed to add the TFSF incident H field to the left surface of the TFSF region

    std::shared_ptr<paramStoreTFSF> rightSurE_; //!< parameter structure describing parameters needed to add the TFSF incident E field to the right surface of the TFSF region
    std::shared_ptr<paramStoreTFSF> rightSurH_; //!< parameter structure describing parameters needed to add the TFSF incident H field to the right surface of the TFSF region

    std::shared_ptr<paramStoreTFSF> frontSurE_; //!< parameter structure describing parameters needed to add the TFSF incident E field to the front surface of the TFSF region
    std::shared_ptr<paramStoreTFSF> frontSurH_; //!< parameter structure describing parameters needed to add the TFSF incident H field to the front surface of the TFSF region

    std::shared_ptr<paramStoreTFSF> backSurE_; //!< parameter structure describing parameters needed to add the TFSF incident E field to the back surface of the TFSF region
    std::shared_ptr<paramStoreTFSF> backSurH_; //!< parameter structure describing parameters needed to add the TFSF incident H field to the back surface of the TFSF region

    std::function< void(std::shared_ptr<parallelGrid<T>>, std::shared_ptr<parallelGrid<T>>, std::shared_ptr<Grid<cplx>>, cplx* incdTransfer_, std::shared_ptr<paramStoreTFSF>) > addEBot_; //!< Function to add the TFSF incident E field to the bot surface of the TFSF region
    std::function< void(std::shared_ptr<parallelGrid<T>>, std::shared_ptr<parallelGrid<T>>, std::shared_ptr<Grid<cplx>>, cplx* incdTransfer_, std::shared_ptr<paramStoreTFSF>) > addETop_; //!< Function to add the TFSF incident E field to the top surface of the TFSF region
    std::function< void(std::shared_ptr<parallelGrid<T>>, std::shared_ptr<parallelGrid<T>>, std::shared_ptr<Grid<cplx>>, cplx* incdTransfer_, std::shared_ptr<paramStoreTFSF>) > addELeft_; //!< Function to add the TFSF incident E field to the left surface of the TFSF region
    std::function< void(std::shared_ptr<parallelGrid<T>>, std::shared_ptr<parallelGrid<T>>, std::shared_ptr<Grid<cplx>>, cplx* incdTransfer_, std::shared_ptr<paramStoreTFSF>) > addERight_; //!< Function to add the TFSF incident E field to the right surface of the TFSF region
    std::function< void(std::shared_ptr<parallelGrid<T>>, std::shared_ptr<parallelGrid<T>>, std::shared_ptr<Grid<cplx>>, cplx* incdTransfer_, std::shared_ptr<paramStoreTFSF>) > addEBack_; //!< Function to add the TFSF incident E field to the front surface of the TFSF region
    std::function< void(std::shared_ptr<parallelGrid<T>>, std::shared_ptr<parallelGrid<T>>, std::shared_ptr<Grid<cplx>>, cplx* incdTransfer_, std::shared_ptr<paramStoreTFSF>) > addEFront_; //!< Function to add the TFSF incident E field to the back surface of the TFSF region

    std::function< void(std::shared_ptr<parallelGrid<T>>, std::shared_ptr<parallelGrid<T>>, std::shared_ptr<Grid<cplx>>, cplx* incdTransfer_, std::shared_ptr<paramStoreTFSF>) > addHBot_; //!< Function to add the TFSF incident H field to the bot surface of the TFSF region
    std::function< void(std::shared_ptr<parallelGrid<T>>, std::shared_ptr<parallelGrid<T>>, std::shared_ptr<Grid<cplx>>, cplx* incdTransfer_, std::shared_ptr<paramStoreTFSF>) > addHTop_; //!< Function to add the TFSF incident H field to the top surface of the TFSF region
    std::function< void(std::shared_ptr<parallelGrid<T>>, std::shared_ptr<parallelGrid<T>>, std::shared_ptr<Grid<cplx>>, cplx* incdTransfer_, std::shared_ptr<paramStoreTFSF>) > addHLeft_; //!< Function to add the TFSF incident H field to the left surface of the TFSF region
    std::function< void(std::shared_ptr<parallelGrid<T>>, std::shared_ptr<parallelGrid<T>>, std::shared_ptr<Grid<cplx>>, cplx* incdTransfer_, std::shared_ptr<paramStoreTFSF>) > addHRight_; //!< Function to add the TFSF incident H field to the right surface of the TFSF region
    std::function< void(std::shared_ptr<parallelGrid<T>>, std::shared_ptr<parallelGrid<T>>, std::shared_ptr<Grid<cplx>>, cplx* incdTransfer_, std::shared_ptr<paramStoreTFSF>) > addHBack_; //!< Function to add the TFSF incident H field to the front surface of the TFSF region
    std::function< void(std::shared_ptr<parallelGrid<T>>, std::shared_ptr<parallelGrid<T>>, std::shared_ptr<Grid<cplx>>, cplx* incdTransfer_, std::shared_ptr<paramStoreTFSF>) > addHFront_; //!< Function to add the TFSF incident H field to the back surface of the TFSF region


public:
    /**
     * @brief      Construct a TFSF surface
     *
     * @param[in]  gridComm  The mpiInterface of the grids
     * @param[in]  locO      location of TFSF origin
     * @param[in]  sz        size of the total field region
     * @param[in]  theta     polar angle of the plane wave's k-vector
     * @param[in]  phi       angle of plane wave's k-vector in the xy plane
     * @param[in]  psi       angle describing the polarization of the plane wave from the vector \vec{k}\times e_{z} (if k is along z psi_ = phi_)
     * @param[in]  circPol   POLARIZATION::R if R polarized, L if L polarized, linear if anything else
     * @param[in]  kLenRelJ  ratio between the size of the axis oriented along psi to that perpendicular to it for elliptically polarized light
     * @param[in]  dx        step size of incident fields
     * @param[in]  dt        time step
     * @param[in]  pul       incident pulseBase
     * @param[in]  Ex        grid_ptr to the Ex field
     * @param[in]  Ey        grid_ptr to the Ey field
     * @param[in]  Ez        grid_ptr to the Ez field
     * @param[in]  Hx        grid_ptr to the Hx field
     * @param[in]  Hy        grid_ptr to the Hy field
     * @param[in]  Hz        grid_ptr to the Hz field
     */
    parallelTFSFBase(std::shared_ptr<mpiInterface> gridComm, std::array<int,3> locO, std::array<int,3> sz, double theta, double phi, double psi, POLARIZATION circPol, double kLenRelJ, double dx, double dt, std::vector<std::shared_ptr<PulseBase>> pul, std::shared_ptr<parallelGrid<T>> Ex, std::shared_ptr<parallelGrid<T>> Ey, std::shared_ptr<parallelGrid<T>> Ez, std::shared_ptr<parallelGrid<T>> Hx, std::shared_ptr<parallelGrid<T>> Hy, std::shared_ptr<parallelGrid<T>> Hz) :
        gridComm_(gridComm),
        gridLen_(0),
        phi_(phi),
        theta_(theta),
        psi_(psi),
        alpha_(0.0),
        phiPrefactCalc_(phi),
        psiPrefactCalc_(psi),
        dt_(dt),
        dx_(dx),
        t_step_(0),
        sz_(sz),
        loc_(locO),
        cdd_(20,0.0),
        cdh_(20,0.0),
        chh_(20,0.0),
        chb0_(20,0.0),
        chb1_(20,0.0),
        pul_(pul),
        D_old_(20,0.0),
        B_old_(20,0.0),
        Ex_(Ex),
        Ey_(Ey),
        Ez_(Ez),
        Hx_(Hx),
        Hy_(Hy),
        Hz_(Hz),
        incdTransfer_( sz[isamax_(sz.size(), sz.data(), 1)-1] )
    {
        if(std::any_of(sz_.begin(), sz_.end(), [](int a){return a == 1; }))
            gridLen_ = 2.0 * std::accumulate(sz_.begin(), sz_.end(), 0);
        else
            gridLen_ = std::accumulate(sz_.begin(), sz_.end(), 0);

        E_incd_ = std::make_shared<Grid<cplx>>( std::array<int,3>({{gridLen_+25,1,1}}), std::array<double,3>({{dx_, dx_, 1}}) );
        H_incd_ = std::make_shared<Grid<cplx>>( std::array<int,3>({{gridLen_+24,1,1}}), std::array<double,3>({{dx_, dx_, 1}}) );
        D_incd_ = std::make_shared<Grid<cplx>>( std::array<int,3>({{20,1,1}}), std::array<double,3>({{dx_, dx_, 1}}) );
        B_incd_ = std::make_shared<Grid<cplx>>( std::array<int,3>({{20,1,1}}), std::array<double,3>({{dx_, dx_, 1}}) );


        for(int ii = 0; ii < 20; ii++)
        {
            double sigE = 0.8*(4.0+1)/dx_ * pow( static_cast<double>(ii)      / 19.0 , 4.0);
            double sigH = 0.8*(4.0+1)/dx_ * pow((static_cast<double>(ii) +0.5)/ 19.0 , 4.0);
            cdd_[ii]  = (2*1.0 - sigE*dt_) / (2*1.0 + sigE*dt_);
            cdh_[ii]  = (2* dt_) / (dx_ * (2*1.0 + sigE*dt_));
            chh_[ii]  = (2*1.0 - sigH*dt_) / (2*1.0 + sigH*dt_);
            chb0_[ii] = (2*1.0 + 0.0*dt_) / (2*1.0 + sigH*dt_);
            chb1_[ii] = (2*1.0 - 0.0*dt_) / (2*1.0 + sigH*dt_);
        }

        // Determine the origin quadrant
        if(phi_ <= M_PI/2.0 && phi_ > 0)
            originQuadrent_ = 1;
        else if(phi_ <= M_PI && phi_ > M_PI/2.0)
            originQuadrent_ = 2;
        else if(phi_ <= 3*M_PI/2.0 && phi_ > M_PI)
            originQuadrent_ = 3;
        else if((phi_ <= 3*M_PI/2.0 && phi_ >= 2.0*M_PI) || phi_ == 0)
            originQuadrent_ = 4;

        // Flip which side the corner is in
        if(theta > M_PI/2.0 || theta_ < 0.0)
            originQuadrent_ *= -1;

        // Phi does not matter when along z-axis set it to make it easier to calculate the field
        if(theta_ == 0 || theta_ == M_PI)
            phiPrefactCalc_ = M_PI/2.0;
        // Determine the psi_/phi_ and alpha_ to match the elliptical shape of the light
        if(circPol == POLARIZATION::R || circPol == POLARIZATION::L )
        {
            // Get alpha parameters
            double c = pow(kLenRelJ, 2.0);
            // phi/psi control the light polarization angle
            psiPrefactCalc_ = 0.5 * asin( sqrt( ( pow(cos(2*psi_),2)*4*c + pow( (1+c)*sin(2*psi_), 2.0) ) / pow(1.0+c, 2.0) ) );
            alpha_ = acos( ( (c - 1.0)*sin(2*psi_) ) / sqrt( pow(cos(2.0*psi_),2.0)*4.0*c + pow( (1.0+c)*sin(2.0*psi_), 2.0) ) );
            if(std::abs( std::tan(psi_) ) > 1)
                psiPrefactCalc_ = M_PI/2.0 - psiPrefactCalc_;
            if(circPol == POLARIZATION::R)
                alpha_ *= -1.0;
        }
        std::shared_ptr<parallelGrid<T>> zField = Ez_ ? Ez_ : Hz_;
        bool twoDim = (Ez_ && Hz_) ? false : true;
        // If TFSF is only a point in 2D and the calc is not a 1D replacement or a line in 3D throw an error since TFSF acts on a plane
        if( ( !twoDim && (sz_[0] <= 1 && sz_[1] <= 1) || (sz_[0] <= 1 && sz_[2] <= 1) || (sz_[1] <= 1 && sz_[2] <= 1) ) || ( twoDim && sz_[0] <=1 && sz_[1] <=1 && zField->local_x() != 3 && zField->local_y() != 3) )
            throw std::logic_error("TFSF sources are not point sources");
        // If an xy-plane, and not 2D (would be a box) make z-normal planes
        if(sz_[0] > 1 && sz_[1] > 1 && !twoDim)
        {
            // Make sure the cell is not a box and then check direction of propagation
            if( sz_[2] <= 1 && ( theta_ != 0.0 && theta_ != M_PI && Ez && Hz) )
            {
                throw std::logic_error("Direction of propagation of TFSF surface is parallel to the main surface.");
            }
            // If a box or a wave traveling in the +z direction
            if( theta_ == 0.0 || sz_[2] > 1)
            {
                backSurE_   = genSurface(zField, DIRECTION::Z, false, true, circPol, kLenRelJ);
                backSurH_   = genSurface(zField, DIRECTION::Z, false, false, circPol, kLenRelJ);
            }
            // If a box or a wave traveling in the -z direction
            if( theta_ == M_PI || sz_[2] > 1)
            {
                frontSurE_   = genSurface(zField, DIRECTION::Z, true, true, circPol, kLenRelJ);
                frontSurH_   = genSurface(zField, DIRECTION::Z, true, false, circPol, kLenRelJ);
            }
        }
        // If in 3D both sz[1] and sz[2] need to be greater than 1; 2D calcs need to have only sz[1] > 1 or be 1D along y
        if( ( sz_[1] > 1 && (sz_[2] > 1 || twoDim ) ) || (twoDim && (zField->local_y() == 3) ) )
        {
            // Make sure the cell is not a box and then check direction of propagation
            if( sz_[1] <= 1 && ( std::abs(originQuadrent_) == 1 ||  std::abs(originQuadrent_) == 3 || theta_ == 0.0 || theta_ == M_PI) )
            {
                throw std::logic_error("Direction of propagation of TFSF surface is parallel to the main surface.");
            }
            // If a box or a wave traveling in the +x direction
            if( std::abs(originQuadrent_) == 2 || sz_[0] > 1)
            {
                rightSurE_ = genSurface(zField, DIRECTION::X, true, true, circPol, kLenRelJ);
                rightSurH_ = genSurface(zField, DIRECTION::X, true, false, circPol, kLenRelJ);
            }
            // If a box or a wave traveling in the -x direction
            if( std::abs(originQuadrent_) == 4 || sz_[0] > 1)
            {
                leftSurE_  = genSurface(zField, DIRECTION::X, false, true, circPol, kLenRelJ);
                leftSurH_  = genSurface(zField, DIRECTION::X, false, false, circPol, kLenRelJ);
            }
        }
        // If in 3D both sz[0] and sz[2] need to be greater than 1; 2D calcs need to have only sz[0] > 1 or be 1D along x
        if( (sz_[0] > 1 && ( sz_[2] > 1 || twoDim ) ) || (twoDim && (zField->local_x() == 3) ) )
        {
            // Make sure the cell is not a box and then check direction of propagation
            if( sz_[0] <= 1 && ( std::abs(originQuadrent_) == 2 || std::abs(originQuadrent_) == 4  || theta_ == 0.0 || theta_ == M_PI) )
            {
                throw std::logic_error("Direction of propagation of TFSF surface is parallel to the main surface.");
            }
            // If a box or a wave traveling in the +y direction
            if(std::abs( originQuadrent_ ) == 1 || sz_[1] > 1 )
            {
                botSurE_ = genSurface(zField, DIRECTION::Y, false, true, circPol, kLenRelJ);
                botSurH_ = genSurface(zField, DIRECTION::Y, false, false, circPol, kLenRelJ);
            }
            // If a box or a wave traveling in the -y direction
            if(std::abs( originQuadrent_ ) == 3 || sz_[1] > 1 )
            {
                topSurE_  = genSurface(zField, DIRECTION::Y, true, true, circPol, kLenRelJ);
                topSurH_  = genSurface(zField, DIRECTION::Y, true, false, circPol, kLenRelJ);
            }
        }
    }
    /**
     * @brief      generates a tfsf surface data structures
     *
     * @param[in]  zField    either the Ez or Hz field for TE or TM mode
     * @param[in]  dir       direction of the surface, left and right X, top and bottom Y
     * @param[in]  pl        if the surface is top or right then true, else false
     * @param[in]  E         true if the source is for the E field
     * @param[in]  circPol   POLARIZATION::R if R polarized, L if L polarized, linear if anything else
     * @param[in]  kLenRelJ  ratio between the size of the axis oriented along psi to that perpendicular to it for elliptically polarized light
     *
     * @return     a TFSF Surface data structure
     */
    std::shared_ptr<paramStoreTFSF> genSurface(std::shared_ptr<parallelGrid<T>> zField, DIRECTION dir, bool pl, bool E, POLARIZATION circPol, double kLenRelJ)
    {
        paramStoreTFSF sur;
        sur.loc_ = {{-1, -1, -1}};
        sur.szTrans_j_ = {{-1, -1}};
        sur.szTrans_k_ = {{-1, -1}};
        int planeStartCor =  -1; // Coordinate of the plane (perpendicular to surface)
        int cor = -1; // the index in the loc/sz vectors corresponding to the i direction
        int transCor1 = -1; // the index in the loc/sz vectors corresponding to the j or k direction (not z)
        int transCor2 = -1; // the index in the loc/sz vectors corresponding to the j or k direction (z if possible)
        int local_add = -1; // local size of the field in the i direction
        int trans1_local_add = -1; // local size of the gird in the trans1 direction
        int trans2_local_add = -1; // local size of the gird in the trans2 direction
        // Set prefactors for updating H field (used on E_inc) use of XOR since sign flips for E/H and pl/mn
        sur.prefactor_j_ = ( E ^ pl ) ?      dt_/dx_ : -1.0*dt_/dx_;
        sur.prefactor_k_ = ( E ^ pl ) ? -1.0*dt_/dx_ :      dt_/dx_;
        if(dir ==DIRECTION::X)
        {
            sur.addVec_ = {{ 0, 0, 1 }};
            cor = 0;
            transCor1 = 1;
            transCor2 = 2;
            sur.strideField_ = zField->local_x()*zField->local_z();
        }
        else if( dir == DIRECTION::Y)
        {
            sur.addVec_ = {{ 0, 0, 1 }};
            cor = 1;
            transCor1 = 0;
            transCor2 = 2;
            sur.strideField_ = 1;
        }
        else if( dir == DIRECTION::Z)
        {
            sur.addVec_ = {{ 0, 1, 0 }};
            cor = 2;
            transCor1 = 0;
            transCor2 = 1;
            sur.strideField_ = 1;
        }
        local_add = zField->ln_vec()[cor] - 2;
        trans1_local_add = zField->ln_vec()[transCor1] - 2;
        trans2_local_add = zField->ln_vec()[transCor2] - 2;
        // Determine the psi_/phi_ and alpha_ to match the elliptical shape of the light
        if(circPol == POLARIZATION::R || circPol == POLARIZATION::L )
        {
            // If circularly polarized do checks to see which fields get the complex $e^{i\alpha}$ prefactor multiplied to it
            if(theta_ == 0 || theta_ == M_PI)
            {
                if(dir == DIRECTION::X && !E)
                    sur.prefactor_k_ *= exp( cplx( 0, alpha_ ) );
                else if(dir == DIRECTION::Y && E)
                    sur.prefactor_j_ *= exp( cplx( 0, alpha_ ) );
                else if(dir == DIRECTION::Z)
                    E ? sur.prefactor_k_ *= exp( cplx( 0, alpha_ ) ) : sur.prefactor_j_ *= exp( cplx( 0, alpha_ ) ) ;
            }
            else if(std::abs(originQuadrent_) == 1 || std::abs(originQuadrent_) == 3)
            {
                if(dir == DIRECTION::X && E)
                    sur.prefactor_j_ *= exp( cplx( 0, alpha_ ) );
                else if(dir == DIRECTION::Y)
                    E ? sur.prefactor_k_ *= exp( cplx( 0, alpha_ ) ) : sur.prefactor_j_ *= exp( cplx( 0, alpha_ ) ) ;
                else if(dir == DIRECTION::Z && !E)
                    sur.prefactor_k_ *= exp( cplx( 0, alpha_ ) );
            }
            else
            {
                if(dir == DIRECTION::X)
                    E ? sur.prefactor_k_ *= exp( cplx( 0, alpha_ ) ) : sur.prefactor_j_ *= exp( cplx( 0, alpha_ ) ) ;
                else if(dir == DIRECTION::Y && !E)
                    sur.prefactor_k_ *= exp( cplx( 0, alpha_ ) );
                else if(dir == DIRECTION::Z && E)
                    sur.prefactor_j_ *= exp( cplx( 0, alpha_ ) );
            }
        }
        // Angle of incidence modifications
        if(Hz_ && Ez_)
        {
            if( dir == DIRECTION::X)
            {
                // if E modify Hz_inc, else Ez_inc
                sur.prefactor_j_ *= E ?  -1.0 * cos(psiPrefactCalc_)*sin(theta_) : sin(psiPrefactCalc_)*sin(theta_);
                // if E modify Hy_inc, else Ey_inc
                sur.prefactor_k_ *= E ?  -1.0 * ( sin(psiPrefactCalc_)*cos(phiPrefactCalc_) - cos(psiPrefactCalc_)*cos(theta_)*sin(phiPrefactCalc_) ) : -1.0 * ( cos(psiPrefactCalc_)*cos(phiPrefactCalc_) + sin(psiPrefactCalc_)*cos(theta_)*sin(phiPrefactCalc_) );
            }
            else if( dir == DIRECTION::Y)
            {
                // if E modify Hx_inc, else Ex_inc
                sur.prefactor_j_ *= E ?  ( sin(psiPrefactCalc_)*sin(phiPrefactCalc_) + cos(psiPrefactCalc_)*cos(theta_)*cos(phiPrefactCalc_) ) : (cos(psiPrefactCalc_)*sin(phiPrefactCalc_) - sin(psiPrefactCalc_)*cos(theta_)*cos(phiPrefactCalc_) );
                // if E modify Hz_inc, else Ez_inc
                sur.prefactor_k_ *= E ?  -1.0 * cos(psiPrefactCalc_)*sin(theta_) : sin(psiPrefactCalc_)*sin(theta_);
            }
            else if( dir == DIRECTION::Z)
            {
                // if E modify Hy_inc, else Ey_inc
                sur.prefactor_j_ *= E ?  -1.0 * ( sin(psiPrefactCalc_)*cos(phiPrefactCalc_) - cos(psiPrefactCalc_)*cos(theta_)*sin(phiPrefactCalc_) ) : -1.0 * ( cos(psiPrefactCalc_)*cos(phiPrefactCalc_) + sin(psiPrefactCalc_)*cos(theta_)*sin(phiPrefactCalc_) );
                // if E modify Hx_inc, else Ex_inc
                sur.prefactor_k_ *= E ?  ( sin(psiPrefactCalc_)*sin(phiPrefactCalc_) + cos(psiPrefactCalc_)*cos(theta_)*cos(phiPrefactCalc_) ) : (cos(psiPrefactCalc_)*sin(phiPrefactCalc_) - sin(psiPrefactCalc_)*cos(theta_)*cos(phiPrefactCalc_) );
            }
        }
        else if(Hz_)
        {
            // Z inc field is unaffected only do for x and y fields
            if(!E)
                dir==DIRECTION::X ? sur.prefactor_k_ *= std::cos(phi_) : sur.prefactor_j_ *= -1.0*std::sin(phi_);
        }
        else if(Ez_)
        {
            // Z inc field is unaffected only do for x and y fields
            if(E)
                dir==DIRECTION::X ? sur.prefactor_k_ *= -1.0*std::cos(phi_) : sur.prefactor_j_ *= std::sin(phi_);
        }

        // Where does the plane start?
        if(!pl)
            E ? planeStartCor = loc_[cor] : planeStartCor = loc_[cor] - 1;
        else
            planeStartCor = loc_[cor] + sz_[cor]-1;

        // Does the process contain the plane?
        if(planeStartCor >= zField->procLoc()[cor] && planeStartCor < zField->procLoc()[cor] + local_add)
        {
            sur.loc_[cor] = planeStartCor - zField->procLoc()[cor]+1;
            // Find the location in the transverse directions
            // Does the location of the plane start within the process or is it a continuation into this process?
            if( loc_[transCor1] >= zField->procLoc()[transCor1] && loc_[transCor1] < zField->procLoc()[transCor1] + trans1_local_add )
                sur.loc_[transCor1] = loc_[transCor1] - zField->procLoc()[transCor1] + 1;
            else if(loc_[transCor1] < zField->procLoc()[transCor1] && loc_[transCor1] + sz_[transCor1] > zField->procLoc()[transCor1])
                sur.loc_[transCor1] = 1;

            // Only a transCor2 if 3D calc
            if(Hz_ && Ez_)
            {
                // Does the location of the plane start within the process or is it a continuation into this process?
                if( loc_[transCor2] >= zField->procLoc()[transCor2] && loc_[transCor2] < zField->procLoc()[transCor2] + trans2_local_add )
                    sur.loc_[transCor2] = loc_[transCor2] - zField->procLoc()[transCor2] + 1;
                else if(loc_[transCor2] < zField->procLoc()[transCor2] && loc_[transCor2] + sz_[transCor2] > zField->procLoc()[transCor2])
                    sur.loc_[transCor2] = 1;
            }
            else
            {
                sur.loc_[2] = 0;
            }
            // Does the process actually have the plane?
            if(sur.loc_[transCor1] != -1 && (!Hz_ || !Ez_ || sur.loc_[transCor2] != -1 ) )
            {
                // Does the plane go out of the process?
                if(sz_[transCor1] + loc_[transCor1] > zField->procLoc()[transCor1] + trans1_local_add)
                {
                    // Yes? then set size to be size of process - loc
                    sur.szTrans_j_[0] = trans1_local_add - sur.loc_[transCor1] + 1;
                    sur.szTrans_k_[0] = trans1_local_add - sur.loc_[transCor1] + 1;
                }
                else
                {
                    // No? Set based on the size inside the process

                    // For the Ey field update from j0+1/2 to j1-1/2; Hy j0 to j1; Ex field update from i0+1/2 to i1-1/2; Hz i0 to i1; Ez field update from i0 to i1; Hz i0+1/2 to i1-1/2
                    sur.szTrans_j_[0] = ( !E ^ DIRECTION::Y == dir ) ? loc_[transCor1] + sz_[transCor1] - (zField->procLoc()[transCor1] + sur.loc_[transCor1] - 1) : loc_[transCor1] + sz_[transCor1] - (zField->procLoc()[transCor1] + sur.loc_[transCor1] - 1) - 1;
                    // For the Hz field update from j0+1/2 to j1-1/2; Ez j0 to j1; Hy field update from i0+1/2 to i1-1/2; Ey i0 to i1; Hx field update from i0 to i1; Ex i0+1/2 to i1-1/2
                    sur.szTrans_k_[0] = (  E ^ DIRECTION::Y == dir ) ? loc_[transCor1] + sz_[transCor1] - (zField->procLoc()[transCor1] + sur.loc_[transCor1] - 1) : loc_[transCor1] + sz_[transCor1] - (zField->procLoc()[transCor1] + sur.loc_[transCor1] - 1) - 1;
                }
                if(Hz_ && Ez_)
                {
                    // Does the plane go out of the process in the trans2 Direction?
                    if(sz_[transCor2] + loc_[transCor2] > zField->procLoc()[transCor2] + trans2_local_add)
                    {
                        // Yes? then set size to be size of process - loc
                        sur.szTrans_j_[1] = trans2_local_add - sur.loc_[transCor2] + 1;
                        sur.szTrans_k_[1] = trans2_local_add - sur.loc_[transCor2] + 1;
                    }
                    else
                    {
                        // No? Set based on the size inside the process

                        // For the Hy field update from k0+1/2 to k1-1/2; Ey k0 to k1; Hz field update from k0 to i1; Ez k0+1/2 to k1-1/2; the Hx field update from j0+1/2 to j1-1/2; Ez i0 to j1
                        sur.szTrans_j_[1] = (  E ^ DIRECTION::Y == dir )? loc_[transCor2] + sz_[transCor2] - (zField->procLoc()[transCor2] + sur.loc_[transCor2] - 1) : loc_[transCor2] + sz_[transCor2] - (zField->procLoc()[transCor2] + sur.loc_[transCor2] - 1) - 1;
                        // For the Ez field update from k0+1/2 to k1-1/2; Hz k0 to k1; Ex field update from k0 to i1; Hx k0+1/2 to k1-1/2; the Ey field update from j0+1/2 to j1-1/2; Hy i0 to j1
                        sur.szTrans_k_[1] = ( !E ^ DIRECTION::Y == dir )? loc_[transCor2] + sz_[transCor2] - (zField->procLoc()[transCor2] + sur.loc_[transCor2] - 1) : loc_[transCor2] + sz_[transCor2] - (zField->procLoc()[transCor2] + sur.loc_[transCor2] - 1) - 1;
                    }
                }
                else
                {
                    // If 2D szTrans[1] = 1 by default
                    sur.szTrans_j_[1] = 1;
                    sur.szTrans_k_[1] = 1;
                }
                bool pl_prop = (theta_ == 0 || (std::fmod(theta_, M_PI) != 0 && abs(originQuadrent_) == 1) || (std::fmod(theta_, M_PI) != 0 && abs(originQuadrent_) == 4) ) ? true : false;
                // Default to assume that the dir == direction of propagation
                if( pl ^ pl_prop )
                {
                    !E ? sur.incdStart_j_ = 2 : sur.incdStart_j_ = 1;
                    !E ? sur.incdStart_k_ = 2 : sur.incdStart_k_ = 1;
                }
                else
                {
                    sur.incdStart_j_ = sz_[cor] + 1;
                    sur.incdStart_k_ = sz_[cor] + 1;
                }
                // Assume that the direction used is in a direction that does not change
                sur.strideIncd_ = 0;
                sur.addIncdProp_ = 0;

                // Correct the assumptions here
                if(theta_ == 0 || theta_ == M_PI)
                {
                    if(dir != DIRECTION::Z)
                    {
                        sur.addIncdProp_ = (theta_ == 0) ? 1 : -1;
                        if(theta_ == 0)
                        {
                            // The start for the Ez, Hz, Ey, Hy, Hx, and Ex incident fields will be the same
                            sur.incdStart_j_  = 2 + (zField->procLoc()[2] - loc_[2]) + (sur.loc_[2] - 1);
                            sur.incdStart_k_  = 2 + (zField->procLoc()[2] - loc_[2]) + (sur.loc_[2] - 1);
                        }
                        else
                        {
                            // The start for both the Ez, Hy, and Ex are the same, and starts for the Hz, Hx, and Ey incident fields will be the same
                            sur.incdStart_j_  = ( E ^ (dir == DIRECTION::Y) ) ? 2 + sz_[2] - (zField->procLoc()[2] - loc_[2] + (sur.loc_[2] ) ) : 2 + sz_[2] - 1 - (zField->procLoc()[2] - loc_[2] + (sur.loc_[2] ) );
                            sur.incdStart_k_  = (!E ^ (dir == DIRECTION::Y) ) ? 2 + sz_[2] - (zField->procLoc()[2] - loc_[2] + (sur.loc_[2] ) ) : 2 + sz_[2] - 1 - (zField->procLoc()[2] - loc_[2] + (sur.loc_[2] ) );
                        }
                    }
                }
                else if(abs(originQuadrent_) == 1 || abs(originQuadrent_) == 3)
                {
                    if(dir != DIRECTION::Y)
                    {
                        if(dir == DIRECTION::X)
                            sur.strideIncd_ = (abs(originQuadrent_) == 1) ? 1 : -1;
                        else if(dir == DIRECTION::Z)
                            sur.addIncdProp_ = (abs(originQuadrent_) == 1) ? 1 : -1;

                        if( abs(originQuadrent_) == 1 )
                        {
                            // The start location for incd fields is the same for Ex, Ey, Ez, Hx, Hy, and Hz
                            sur.incdStart_j_  = 2 + (zField->procLoc()[1] - loc_[1]) + (sur.loc_[1] - 1);
                            sur.incdStart_k_  = 2 + (zField->procLoc()[1] - loc_[1]) + (sur.loc_[1] - 1);
                        }
                        else
                        {
                            if(dir == DIRECTION::X )
                            {
                                // The Ez, Hx and Hy incident fields have the same start but the Hz, Ex, Ey fields have the same but different value
                                sur.incdStart_j_  = !E ? 2 + sz_[1] - (zField->procLoc()[1] - loc_[1] + (sur.loc_[1] - 1) + sur.szTrans_j_[0] ) : 2 + sz_[1] - 1 - (zField->procLoc()[1] - loc_[1] + (sur.loc_[1] - 1) + sur.szTrans_j_[0] );
                                sur.incdStart_k_  =  E ? 2 + sz_[1] - (zField->procLoc()[1] - loc_[1] + (sur.loc_[1] - 1) + sur.szTrans_k_[0] ) : 2 + sz_[1] - 1 - (zField->procLoc()[1] - loc_[1] + (sur.loc_[1] - 1) + sur.szTrans_k_[0] );
                            }
                            else if(dir == DIRECTION::Z)
                            {
                                // The Ez, Hx and Hy incident fields have the same start but the Hz, Ex, Ey fields have the same but different value
                                sur.incdStart_j_  =  E ? 2 + sz_[1] - (zField->procLoc()[1] - loc_[1] + (sur.loc_[1] ) ) : 2 + sz_[1] - 1 - (zField->procLoc()[1] - loc_[1] + (sur.loc_[1] ) );
                                sur.incdStart_k_  = !E ? 2 + sz_[1] - (zField->procLoc()[1] - loc_[1] + (sur.loc_[1] ) ) : 2 + sz_[1] - 1 - (zField->procLoc()[1] - loc_[1] + (sur.loc_[1] ) );
                            }
                        }
                    }
                }
                else if(abs(originQuadrent_) == 4 || abs(originQuadrent_) == 2)
                {
                    if(dir != DIRECTION::X)
                    {
                        sur.strideIncd_ = (std::abs(originQuadrent_) == 4) ? 1 : -1;
                        if(abs(originQuadrent_) == 4 )
                        {
                            // The Ex, Hx, Ey, Hy, Ez, and Hz incident field starts are all the same
                            sur.incdStart_j_  = 2 + (zField->procLoc()[0] - loc_[0]) + (sur.loc_[0] - 1);
                            sur.incdStart_k_  = 2 + (zField->procLoc()[0] - loc_[0]) + (sur.loc_[0] - 1);
                        }
                        else
                        {
                            // The incident Ex, Ey, and Hz incident fields have the same start, while Hx, Hy and Ez incident fields have the same start point
                            sur.incdStart_j_  = ( E ^ ( dir == DIRECTION::Z) ) ? 2 + sz_[0] - (zField->procLoc()[0] - loc_[0] + (sur.loc_[0] - 1) + sur.szTrans_j_[0] ) : 2 + sz_[0] - 1 - (zField->procLoc()[0] - loc_[0] + (sur.loc_[0] - 1) + sur.szTrans_j_[0] );
                            sur.incdStart_k_  = (!E ^ ( dir == DIRECTION::Z) ) ? 2 + sz_[0] - (zField->procLoc()[0] - loc_[0] + (sur.loc_[0] - 1) + sur.szTrans_k_[0] ) : 2 + sz_[0] - 1 - (zField->procLoc()[0] - loc_[0] + (sur.loc_[0] - 1) + sur.szTrans_k_[0] );
                        }
                    }
                }
            }
            else
                return nullptr;
            return std::make_shared<paramStoreTFSF>(sur);
        }
        else
        {
            return nullptr;
        }
    }

    /**
     * @brief      updates the the fields using update functors
     */
    void updateFileds()
    {
        addHBot_  (Hz_, Hx_, E_incd_, incdTransfer_.data(), botSurH_);
        addHTop_  (Hz_, Hx_, E_incd_, incdTransfer_.data(), topSurH_);

        addHLeft_ (Hy_, Hz_, E_incd_, incdTransfer_.data(), leftSurH_);
        addHRight_(Hy_, Hz_, E_incd_, incdTransfer_.data(), rightSurH_);

        addHBack_ (Hx_, Hy_, E_incd_, incdTransfer_.data(), backSurH_);
        addHFront_(Hx_, Hy_, E_incd_, incdTransfer_.data(), frontSurH_);

        step();

        addEBot_  (Ez_, Ex_, H_incd_, incdTransfer_.data(), botSurE_);
        addETop_  (Ez_, Ex_, H_incd_, incdTransfer_.data(), topSurE_);

        addELeft_ (Ey_, Ez_, H_incd_, incdTransfer_.data(), leftSurE_);
        addERight_(Ey_, Ez_, H_incd_, incdTransfer_.data(), rightSurE_);

        addEBack_ (Ex_, Ey_, H_incd_, incdTransfer_.data(), backSurE_);
        addEFront_(Ex_, Ey_, H_incd_, incdTransfer_.data(), frontSurE_);
    }
    /**
     * @return origin Location
     */
    inline std::vector<int> & loc()  {return loc_;}
    /**
     * @return size of total field region
     */
    inline std::vector<int> & size() {return sz_;}
    /**
     * @return angle of incidence of the light
     */
    inline double phi(){return phi_;}

    /**
     * @brief      returns prefactor phi factor
     *
     * @return     phiPrefactCalc_
     */
    inline double phiPreFact(){return phiPrefactCalc_;}

    /**
     * @return angle of incidence of the light
     */
    inline double theta(){return theta_;}

    /**
     * @return angle of polarization of the light
     */
    inline double psi(){return psi_;}

    /**
     * @brief      returns prefactor psi factor
     *
     * @return     psiPrefactCalc_
     */
    inline double psiPreFact(){return psiPrefactCalc_;}

    /**
     * @return phase angle of incidence of the light if circularly polarized
     */
    inline double alpha(){return alpha_;}
    /**
     * @return origin quadrant
     */
    inline int quadrant(){return originQuadrent_;}
    /**
     * @return length of 1D auxiliary fields
     */
    inline int gridLen(){return gridLen_;}

    /**
     * @brief      Accessor function for the incident E field
     */
    inline cplx E_incd() {return E_incd_->point(2,0);}
    /**
     * @brief      Accessor function for the incident E field one grid point further from the start in the direction of propagation
     */
    inline cplx E_pl_incd() {return E_incd_->point(3,0);}
    /**
     * @brief      Accessor function for the incident H field
     */
    inline cplx H_incd() {return H_incd_->point(2,0);}
    /**
     * @brief      Accessor function for the incident H field one grid point further from the start in the direction  oppisite of propagation
     */
    inline cplx H_mn_incd() {return H_incd_->point(1,0);}

    /**
     * @brief Moves the incident fields forward one time step
     * @details Uses the 1D FDTD equations to propagate the fields in time
     */
    void step()
    {
        // Update incident H field
        zaxpy_(gridLen_+4,      dt_/dx_, &E_incd_->point(0,0), 1, &H_incd_->point(0,0), 1);
        zaxpy_(gridLen_+4, -1.0*dt_/dx_, &E_incd_->point(1,0), 1, &H_incd_->point(0,0), 1);

        // Copy the auxiliary fields to incident fields
        zcopy_(20, &B_incd_->point(0,0), 1, B_old_.data(), 1);
        zcopy_(20, &D_incd_->point(0,0), 1, D_old_.data(), 1);

        // Do the H PMLs
        zaxpy_(20,      dt_/dx_, &E_incd_->point(gridLen_+4,0), 1, &B_incd_->point(0,0), 1);
        zaxpy_(20, -1.0*dt_/dx_, &E_incd_->point(gridLen_+5,0), 1, &B_incd_->point(0,0), 1);
        for(int ii = 0; ii < 20; ii ++)
            H_incd_->point(gridLen_+4+ii,0) = chh_[ii] * H_incd_->point(gridLen_+4+ii,0) + chb0_[ii] * B_incd_->point(ii,0) - chb1_[ii]*B_old_[ii];

        // Add the pulses
        E_incd_ -> point(0,0) = 0.0;
        for(auto& pul : pul_)
            E_incd_ -> point(0,0) += pul->pulse(static_cast<double>(t_step_)*dt_);

        // Update H incident fields
        zaxpy_(gridLen_+3,      dt_/dx_, &H_incd_->point(0,0), 1, &E_incd_->point(1,0), 1);
        zaxpy_(gridLen_+3, -1.0*dt_/dx_, &H_incd_->point(1,0), 1, &E_incd_->point(1,0), 1);

        // PMLs for E
        for(int ii = 0; ii < 20; ii++)
            D_incd_->point(ii,0) = cdd_[ii] * D_incd_->point(ii,0) + cdh_[ii] * (H_incd_->point(gridLen_+3+ii,0) - H_incd_->point(gridLen_+4+ii,0));
        zaxpy_(20,  1.0, &D_incd_->point(0,0), 1, &E_incd_->point(gridLen_+4,0), 1);
        zaxpy_(20, -1.0, D_old_.data()       , 1, &E_incd_->point(gridLen_+4,0), 1);
        t_step_++;
    }
};

namespace tfsfUpdateFxnReal
{
    /**
     * @brief      Adds a tfsf fields to two components of the real grids.
     *
     * @param[in]  grid_j        grid pointers to the grid j
     * @param[in]  grid_k        grid pointers to the grid k
     * @param[in]  incd          The incd field value
     * @param[in]  incdTransfer  pointer to the start of scratch space vector fro transfer
     * @param[in]  sur           The surface parameter struct
     */
    void addTFSFTwoComp (real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur);

    /**
     * @brief      Adds a tfsf fields to the j component of the real grids.
     *
     * @param[in]  grid_j        grid pointers to the grid j
     * @param[in]  grid_k        grid pointers to the grid k
     * @param[in]  incd          The incd field value
     * @param[in]  incdTransfer  pointer to the start of scratch space vector fro transfer
     * @param[in]  sur           The surface parameter struct
     */
    void addTFSFOneCompJ(real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur);

    /**
     * @brief      Adds a tfsf fields to k component of the real grids.
     *
     * @param[in]  grid_j        grid pointers to the grid j
     * @param[in]  grid_k        grid pointers to the grid k
     * @param[in]  incd          The incd field value
     * @param[in]  incdTransfer  pointer to the start of scratch space vector fro transfer
     * @param[in]  sur           The surface parameter struct
     */
    void addTFSFOneCompK(real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur);
    /**
     * @brief      does grid->gridTransfer() for grid pointers
     *
     * @param[in]  grid grid that needs grid->gridTransfer to it
     */
    void transferDat(real_pgrid_ptr grid);
}

namespace tfsfUpdateFxnCplx
{
    /**
     * @brief      Adds a tfsf fields to two components of the real grids.
     *
     * @param[in]  grid_j        grid pointers to the grid j
     * @param[in]  grid_k        grid pointers to the grid k
     * @param[in]  incd          The incd field value
     * @param[in]  incdTransfer  pointer to the start of scratch space vector fro transfer
     * @param[in]  sur           The surface parameter struct
     */
    void addTFSFTwoComp (cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur);
    /**
     * @brief      Adds a tfsf fields to the j component of the real grids.
     *
     * @param[in]  grid_j        grid pointers to the grid j
     * @param[in]  grid_k        grid pointers to the grid k
     * @param[in]  incd          The incd field value
     * @param[in]  incdTransfer  pointer to the start of scratch space vector fro transfer
     * @param[in]  sur           The surface parameter struct
     */
    void addTFSFOneCompJ(cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur);
    /**
     * @brief      Adds a tfsf fields to k component of the real grids.
     *
     * @param[in]  grid_j        grid pointers to the grid j
     * @param[in]  grid_k        grid pointers to the grid k
     * @param[in]  incd          The incd field value
     * @param[in]  incdTransfer  pointer to the start of scratch space vector fro transfer
     * @param[in]  sur           The surface parameter struct
     */
    void addTFSFOneCompK(cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur);
    /**
     * @brief      does grid->gridTransfer() for grid pointers
     *
     * @param[in]  grid grid that needs grid->gridTransfer to it
     */
    void transferDat(cplx_pgrid_ptr);
}

class parallelTFSFReal : public parallelTFSFBase<double>
{
public:

    /**
     * @brief      Construct a TFSF surface
     *
     * @param[in]  gridComm  The mpiInterface of the grids
     * @param[in]  locO      location of TFSF origin
     * @param[in]  sz        size of the total field region
     * @param[in]  theta     polar angle of the plane wave's k-vector
     * @param[in]  phi       angle of plane wave's k-vector in the xy plane
     * @param[in]  psi       angle describing the polarization of the plane wave from the vector \vec{k}\times e_{z} (if k is along z psi_ = phi_)
     * @param[in]  circPol   POLARIZATION::R if R polarized, L if L polarized, linear if anything else
     * @param[in]  kLenRelJ  ratio between the size of the axis oriented along psi to that perpendicular to it for elliptically polarized light
     * @param[in]  dx        step size of incident fields
     * @param[in]  dt        time step
     * @param[in]  pul       incident pulseBase
     * @param[in]  Ex        grid_ptr to the Ex field
     * @param[in]  Ey        grid_ptr to the Ey field
     * @param[in]  Ez        grid_ptr to the Ez field
     * @param[in]  Hx        grid_ptr to the Hx field
     * @param[in]  Hy        grid_ptr to the Hy field
     * @param[in]  Hz        grid_ptr to the Hz field
     */
    parallelTFSFReal(std::shared_ptr<mpiInterface> gridComm, std::array<int,3> locO, std::array<int,3> sz, double theta, double phi, double psi, POLARIZATION circPol, double kLenRelJ, double dx, double dt,  std::vector<std::shared_ptr<PulseBase>> pul, real_pgrid_ptr Ex, real_pgrid_ptr Ey, real_pgrid_ptr Ez, real_pgrid_ptr Hx, real_pgrid_ptr Hy, real_pgrid_ptr Hz);
};
class parallelTFSFCplx : public parallelTFSFBase<cplx>
{
public:

    /**
     * @brief      Construct a TFSF surface
     *
     * @param[in]  gridComm  The mpiInterface of the grids
     * @param[in]  locO      location of TFSF origin
     * @param[in]  sz        size of the total field region
     * @param[in]  theta     polar angle of the plane wave's k-vector
     * @param[in]  phi       angle of plane wave's k-vector in the xy plane
     * @param[in]  psi       angle describing the polarization of the plane wave from the vector \vec{k}\times e_{z} (if k is along z psi_ = phi_)
     * @param[in]  circPol   POLARIZATION::R if R polarized, L if L polarized, linear if anything else
     * @param[in]  kLenRelJ  ratio between the size of the axis oriented along psi to that perpendicular to it for elliptically polarized light
     * @param[in]  dx        step size of incident fields
     * @param[in]  dt        time step
     * @param[in]  pul       incident pulseBase
     * @param[in]  Ex        grid_ptr to the Ex field
     * @param[in]  Ey        grid_ptr to the Ey field
     * @param[in]  Ez        grid_ptr to the Ez field
     * @param[in]  Hx        grid_ptr to the Hx field
     * @param[in]  Hy        grid_ptr to the Hy field
     * @param[in]  Hz        grid_ptr to the Hz field
     */
    parallelTFSFCplx(std::shared_ptr<mpiInterface> gridComm, std::array<int,3> locO, std::array<int,3> sz, double theta, double phi, double psi, POLARIZATION circPl, double kLenRelJ, double dx, double dt,  std::vector<std::shared_ptr<PulseBase>> pul, cplx_pgrid_ptr Ex, cplx_pgrid_ptr Ey, cplx_pgrid_ptr Ez, cplx_pgrid_ptr Hx, cplx_pgrid_ptr Hy, cplx_pgrid_ptr Hz);
};

#endif