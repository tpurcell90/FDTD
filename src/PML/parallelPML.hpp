#ifndef PARALLEL_FDTD_PML
#define PARALLEL_FDTD_PML

#include <src/GRID/parallelGrid.hpp>
#include <src/OBJECTS/Obj.hpp>


struct updatePsiParams
{
    int transSz_; //!< size in the transvers diretion (how many points to do mkl operations on)
    int stride_; //!< stride for the mkl operations
    double b_; //!< the b parameter as defined in chapter 7 of Taflove
    double c_; //!< the c parameter as defined in chapter 7 of Taflove
    double cOff_; //!< the c parameter as defined in chapter 7 of Taflove, but with the oppisite sing of c
    std::array<int,3> loc_; //!< the starting point of the mkl operations
    std::array<int,3> locOff_; //!< the starting point of teh mkl operations offset by the corret value
};
/**
 * @brief a stoarge struct for adding ths psi fields into the the EM fields
 * @details Stores all necessary inforation to add teh psi fileds to the EM grids

 */
struct updateGridParams
{
    int nAx_; //!< size of mkl operator
    int stride_; //!< tride for mkl operator
    double Db_; //!< psi factor add on prefactor
    std::array<int,3> loc_; //!< starting point of MKL operations
};

/**
 * @brief Implimentation of the CPMLs for the parallel FDTD fields
 * @details An implimientation of CPML for absorbing the EM fields at the cell boundaries
 *
 */
template <typename T> class parallelCPML
{
    typedef std::shared_ptr<parallelGrid<T>> pgrid_ptr;
protected:
    mpiInterface & gridComm_; //!< MPI Communicator
    POLARIZATION pol_i_; //!< Polarization of the Field that the CPML is acting on
    DIRECTION i_; //!< Direction corresponding to polarization of the field the CPML is acting on
    DIRECTION j_; //!< Direction next in the cycle from i_; i.e. if i_ is DIRECTION::Y then j_ is DIRECTION::Z
    DIRECTION k_; //!< Direction last one in the cycle from i_; i.e. if i_ is DIRECTION::Y then j_ is DIRECTION::X
    double m_; //!< Exponent for calculating CPML constants sigma and kappa
    double ma_; //!< Exponent for calculating CPML constant a
    double sigmaMax_; //!< Maximum value of sigma
    double kappaMax_; //!< Maximum value of kappa
    double aMax_; //!< Maximum value of a
    double dt_; //!< time step
    std::array<int,3> n_vec_; //!< vector storing the thickness of the PML in all directions
    std::array<int,3> ln_vec_pl_; //!< vector storing the local thickness of the PMLs in the positive directions (top, right, and front)
    std::array<int,3> ln_vec_mn_; //!< vector storing the local thickness of the PMLs in the positive directions (bottom, left, and back)
    std::array<double,3> d_; //!< vector storing the step size in all directions
    std::vector<double> eta_eff_top_; //!< eta effective along top PML
    std::vector<double> eta_eff_bot_; //!< eta effective along bottom PML
    std::vector<double> eta_eff_left_; //!< eta effective along left PML
    std::vector<double> eta_eff_right_; //!< eta effective along right PML
    std::vector<double> eta_eff_front_; //!< eta effective along right PML
    std::vector<double> eta_eff_back_; //!< eta effective along right PML
    std::function<void(std::vector<updateGridParams>&, std::vector<updatePsiParams>&, pgrid_ptr, pgrid_ptr, pgrid_ptr)> upPsi_j_; //!< update function for the $\\psi_j$ field
    std::function<void(std::vector<updateGridParams>&, std::vector<updatePsiParams>&, pgrid_ptr, pgrid_ptr, pgrid_ptr)> upPsi_k_; //!< update function for the $\\psi_k$ field

public:
    pgrid_ptr grid_i_; //!< FDTD field polarized in the i direction (E/H) is the same as it is in pol_i
    pgrid_ptr grid_j_; //!< FDTD field polarized in the j direction (E/H) is the opposite as it is in pol_i
    pgrid_ptr grid_k_; //!< FDTD field polarized in the k direction (E/H) is the opposite as it is in pol_i
    pgrid_ptr psi_j_; //!< CPML helper field $\\psi$ polarized in the j direction
    pgrid_ptr psi_k_; //!< CPML helper field $\\psi$ polarized in the k direction
    std::vector<updatePsiParams>  updateListPsi_j_; //!< update parameters for $\\psi_j$ field
    std::vector<updatePsiParams>  updateListPsi_k_; //!< update parameters for $\\psi_k$ field
    std::vector<updateGridParams> updateListGrid_j_; //!< update parameters for the updates in the j direction
    std::vector<updateGridParams> updateListGrid_k_; //!< update parameters for the updates in the k direction

    /**
     * @brief      { function_description }
     *
     * @param[in]  gridComm  mpi communicator
     * @param[in]  weights   The weights
     * @param[in]  grid_i    shared pointer to grid_i (field PML is being applied to)
     * @param[in]  grid_j    shared pointer to grid_j (field polarized in the j direction of the PML; nullptr if none)
     * @param[in]  grid_k    shared pointer to grid_k (field polarized in the k direction of the PML; nullptr if none)
     * @param[in]  pol_i     polarization of grid_o
     * @param[in]  n_vec     Vector storing thickness of the PMLs in all directions
     * @param[in]  m         scalling factor for sigma
     * @param[in]  ma        scaling factor for a
     * @param[in]  aMax      max a value
     * @param[in]  d         vector storing the step sizes in all directions
     * @param[in]  dt        time step
     * @param[in]  phys_Ex   The physical ex grid
     * @param[in]  phys_Ey   The physical ey grid
     * @param[in]  objArr    The object arr
     */
    parallelCPML(mpiInterface gridComm, std::vector<real_grid_ptr> weights, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k, POLARIZATION pol_i, std::array<int,3> n_vec, double m, double ma, double aMax, std::array<double,3> d, double dt, int_pgrid_ptr phys_Ex, int_pgrid_ptr phys_Ey, std::vector<std::shared_ptr<Obj>> objArr) :
        gridComm_(gridComm),
        pol_i_(pol_i),
        n_vec_(n_vec),
        ln_vec_mn_({{ 0,0,0,}}),
        ln_vec_pl_({{ 0,0,0,}}),
        m_(m),
        ma_(ma),
        aMax_(aMax),
        d_(d),
        dt_(dt),
        grid_i_(grid_i),
        grid_j_(grid_j),
        grid_k_(grid_k)
    {
        if(!grid_i)
            throw std::logic_error("PMLs have to have a real field to be associated with, grid_i_ can't be a nullptr.");
        kappaMax_ = 1.0;
        sigmaMax_ = 0.8*(m_+1)/d_[0];
        if(pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EX)
        {
            i_ = DIRECTION::X;
            j_ = DIRECTION::Y;
            k_ = DIRECTION::Z;
        }
        else if(pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EY)
        {
            i_ = DIRECTION::Y;
            j_ = DIRECTION::Z;
            k_ = DIRECTION::X;
        }
        else
        {
            i_ = DIRECTION::Z;
            j_ = DIRECTION::X;
            k_ = DIRECTION::Y;
        }
        // In 2D z direction is assumed to be isotropic
        if(grid_k_ && (grid_i_->local_z() != 1 || j_ != DIRECTION::Z) )
        {
            psi_j_ = std::make_shared<parallelGrid<T>>(gridComm, grid_k_->PBC(), weights, grid_k_->n_vec(), grid_k_->d() );
        }
        else
            psi_j_ = nullptr;

        // In 2D z direction is assumed to be isotropic
        if(grid_j_ && (grid_i_->local_z() != 1 || k_ != DIRECTION::Z) )
        {
            psi_k_ = std::make_shared<parallelGrid<T>>(gridComm, grid_j_->PBC(), weights, grid_j_->n_vec(), grid_j_->d() );
        }
        else
            psi_k_ = nullptr;

        if(!psi_k_ && ! psi_j_)
            throw std::logic_error("Well now this is awkward, you apperently are telling the parallelCPML that you only have one field, light needs at least 2 so FAIL!!!!");
        genDatStruct();
        genEtaEff(phys_Ex,phys_Ey, objArr);
        initalizeGrid(phys_Ex, phys_Ey, objArr);
        gridComm_.barrier();
    }

    /**
     * @brief      finds teh eta_eff values for all PMLs
     * @details    no physEz since it is not needed in averaging
     *
     * @param[in]  phys_Ex  physical grid for the y direction
     * @param[in]  phys_Ey  physical grid for the x direction
     * @param[in]  objArr   vector of all objects in the cell
     */
    void genEtaEff(int_pgrid_ptr phys_Ex, int_pgrid_ptr phys_Ey, std::vector<std::shared_ptr<Obj>> objArr)
    {
        gridComm_.barrier();
        for(int ii = 0; ii < n_vec_[0] + gridComm_.npX()*2; ii ++)
        {
            std::vector<int> physEyVals = phys_Ey->getPlaneYZ(ii);

            double eps_sum = 0.0;
            // Sum up all the dielectric constants in the plane, if it is a boundary region don't add it
            for(int vv = 0; vv < physEyVals.size(); vv++)
                if(physEyVals[vv] != -1)
                    eps_sum += objArr[physEyVals[vv]]->epsInfty();
            // Find the average dielectric constant if outside the boundary regions
            double normEpsSum = phys_Ex->local_z() != 1 ? physEyVals.size() - 2*(gridComm_.npY()*phys_Ex->z() + gridComm_.npZ()*phys_Ex->y() ) + 2*(gridComm_.npY() + 1) : physEyVals.size() - 2*gridComm_.npY();
            if(eta_eff_left_.size() < n_vec_[0] && eps_sum != 0.0)
            {
                eta_eff_left_.push_back( eps_sum/normEpsSum );
            }
        }
        for(int ii = phys_Ey->x() - 2; ii >= phys_Ey->x() - (n_vec_[0] + gridComm_.npX()*2); ii --)
        {
            std::vector<int> physEyVals = phys_Ey->getPlaneYZ(ii);

            double eps_sum = 0.0;
            // Sum up all the dielectric constants in the plane, if it is a boundary region don't add it
            for(int vv = 0; vv < physEyVals.size(); vv++)
                if(physEyVals[vv] != -1)
                    eps_sum += objArr[physEyVals[vv]]->epsInfty();
            // Find the average dielectric constant if outside the boundary regions
            double normEpsSum = phys_Ex->local_z() != 1 ? physEyVals.size() - 2*(gridComm_.npY()*phys_Ex->z() + gridComm_.npZ()*phys_Ex->y() ) + 2*(gridComm_.npY() + 1) : physEyVals.size() - 2*gridComm_.npY();
            if(eta_eff_right_.size() < n_vec_[0] && eps_sum != 0.0)
            {
                eta_eff_right_.push_back( eps_sum/normEpsSum );
            }
        }
        for(int ii = 0; ii < n_vec_[1] + gridComm_.npY()*2; ii ++)
        {
            std::vector<int> physExVals = phys_Ex->getPlaneXZ(ii);
            double eps_sum = 0.0;
            // Sum up all the dielectric constants in the plane, if it is a boundary region don't add it
            for(int vv = 0; vv < physExVals.size(); vv++)
                if(physExVals[vv] != -1)
                    eps_sum += objArr[physExVals[vv]]->epsInfty();

            // Find the average dielectric constant if outside the boundary regions
            double normEpsSum = phys_Ex->local_z() != 1 ? (physExVals.size() - 2*(gridComm_.npZ()*phys_Ex->x() + gridComm_.npX()*phys_Ex->z() ) + 4) : (physExVals.size() - 2*( gridComm_.npX() ) );
            if(eta_eff_bot_.size() < n_vec_[1] && eps_sum != 0)
            {
                eta_eff_bot_.push_back(eps_sum/normEpsSum ); //Not based on numProcs since in xz planes only 1 proc just corners
            }
        }
        for(int ii = phys_Ex->y() - 2; ii >= phys_Ex->y() - (n_vec_[1] + gridComm_.npY()*2); ii --)
        {
            std::vector<int> physExVals = phys_Ex->getPlaneXZ(ii);
            double eps_sum = 0.0;
            // Sum up all the dielectric constants in the plane, if it is a boundary region don't add it
            for(int vv = 0; vv < physExVals.size(); vv++)
                if(physExVals[vv] != -1)
                    eps_sum += objArr[physExVals[vv]]->epsInfty();

            // Find the average dielectric constant if outside the boundary regions
            double normEpsSum = phys_Ex->local_z() != 1 ? (physExVals.size() - 2*(gridComm_.npZ()*phys_Ex->x() + gridComm_.npX()*phys_Ex->z() ) + 4) : (physExVals.size() - 2*( gridComm_.npX() ) );
            if(eta_eff_top_.size() < n_vec_[1] && eps_sum != 0)
            {
                eta_eff_top_.push_back( eps_sum/normEpsSum ); //Not based on numProcs since in xz planes only 1 proc just corners
            }
        }

        if(grid_i_->local_z() != 1)
        {
            for(int ii = 0; ii < n_vec_[2] + gridComm_.npZ()*2; ii ++)
            {
                std::vector<int> physExVals = phys_Ex->getPlaneXY(ii);
                double eps_sum = 0.0;
                // Sum up all the dielectric constants in the plane, if it is a boundary region don't add it
                for(int vv = 0; vv < physExVals.size(); vv++)
                    if(physExVals[vv] != -1)
                        eps_sum += objArr[physExVals[vv]]->epsInfty();

                // Find the average dielectric constant if outside the boundary regions
                if(eta_eff_back_.size() < n_vec_[2] && eps_sum != 0)
                    eta_eff_back_.push_back(eps_sum/(physExVals.size() - 2*(gridComm_.npY()*phys_Ex->x() + gridComm_.npX()*phys_Ex->y() ) + 2*(gridComm_.npY() + 1) ) );
            }
            for(int ii = phys_Ex->z() - 2; ii >= phys_Ex->z() - (n_vec_[2] + gridComm_.npZ()*2); ii --)
            {
                std::vector<int> physExVals = phys_Ex->getPlaneXY(ii);
                double eps_sum = 0.0;
                // Sum up all the dielectric constants in the plane, if it is a boundary region don't add it
                for(int vv = 0; vv < physExVals.size(); vv++)
                    if(physExVals[vv] != -1)
                        eps_sum += objArr[physExVals[vv]]->epsInfty();

                // Find the average dielectric constant if outside the boundary regions
                if(eta_eff_front_.size() < n_vec_[2] && eps_sum != 0)
                    eta_eff_front_.push_back(eps_sum/(physExVals.size() - 2*(gridComm_.npY()*phys_Ex->x() + gridComm_.npX()*phys_Ex->y() ) + 2*(gridComm_.npY() + 1) ) );
            }
        }
    }

    /**
     * @brief determine the local nx and ny values for all directions on the local processor
     */
    void genDatStruct()
    {
        if(grid_i_->procLoc()[0] < n_vec_[0])
            ln_vec_mn_[0] = grid_i_->procLoc()[0] + grid_i_->local_x()-2 < n_vec_[0] ? grid_i_->local_x()-2 : n_vec_[0] - grid_i_->procLoc()[0];
        if(grid_i_->procLoc()[0] + grid_i_->local_x()-2 > grid_i_->x() - 2*gridComm_.npX() - n_vec_[0])
            ln_vec_pl_[0] = grid_i_->procLoc()[0] > grid_i_->x() - 2*gridComm_.npX() - n_vec_[0] ? grid_i_->local_x()-2 : grid_i_->procLoc()[0] + grid_i_->local_x() - 2 - (grid_i_->x() - 2*gridComm_.npX() - n_vec_[0]);;

        if(grid_i_->procLoc()[1] < n_vec_[1])
            ln_vec_mn_[1] = grid_i_->procLoc()[1] + grid_i_->local_y()-2 < n_vec_[1] ? grid_i_->local_y()-2 : n_vec_[1] - grid_i_->procLoc()[1];
        if(grid_i_->procLoc()[1] + grid_i_->local_y()-2 > grid_i_->y() - gridComm_.npY()*2 - n_vec_[1])
            ln_vec_pl_[1] = grid_i_->procLoc()[1] > grid_i_->y() - gridComm_.npY()*2 - n_vec_[1] ? grid_i_->local_y()-2 : grid_i_->procLoc()[1] + grid_i_->local_y() - 2 - (grid_i_->y() - 2*gridComm_.npY() - n_vec_[1]);

        if(grid_i_->local_z() != 1)
        {
            if(grid_i_->procLoc()[2] < n_vec_[2])
                ln_vec_mn_[2] = grid_i_->procLoc()[2] + grid_i_->local_z()-2 < n_vec_[2] ? grid_i_->local_z()-2 : n_vec_[2] - grid_i_->procLoc()[2];

            if(grid_i_->procLoc()[2] + grid_i_->local_z()-2 > grid_i_->z() - gridComm_.npZ()*2 - n_vec_[2])
                ln_vec_pl_[2] = grid_i_->procLoc()[2] > grid_i_->z() - gridComm_.npZ()*2 - n_vec_[2] ? grid_i_->local_z()-2 : grid_i_->procLoc()[2] + grid_i_->local_z() - 2 - (grid_i_->z() - 2*gridComm_.npZ() - n_vec_[2]);
        }
    }

    /**
     * @brief      Generate the psi field update lists.
     *
     * @param[in]  dir      Direction of the psi field.
     * @param[in]  pl       True if top, right or front.
     * @param[in]  startPt  Where the PMLs starts.
     * @param[in]  pmlEdge  Where the PML ends
     * @param[in]  nDir     Thickness in that direction.
     * @param[in]  dirMax   How far to iterate over.
     *
     * @return     The psi up list.
     */
    std::vector<updatePsiParams> getPsiUpList(DIRECTION dir, bool pl, int startPt,  int pmlEdge, int nDir, int dirMax)
    {
        std::vector<updatePsiParams> paramList;
        updatePsiParams param;

        bool negC = true;
        int transSz2 = 1;
        std::vector<double> *eta_eff;
        if(dir == DIRECTION::X)
        {
            if(grid_i_->local_z() != 1)
            {
                if(pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EZ)
                    param.transSz_ -= 1;
                param.transSz_= grid_i_->local_z()-2;
                transSz2 = grid_i_->local_y()-2;
                if( (pol_i_ == POLARIZATION::EY || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HZ ) && gridComm_.rank() == gridComm_.size() - 1)
                    transSz2 -= 1;
            }
            else
            {
                param.transSz_= grid_i_->local_y()-2;
                if( (pol_i_ == POLARIZATION::EY || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HZ ) && gridComm_.rank() == gridComm_.size() - 1)
                    param.transSz_ -= 1;
            }
            param.stride_ = grid_i_->local_x();


            if(pol_i_ == POLARIZATION::EZ || pol_i_ == POLARIZATION::EY)
                negC = false;

            if(pl)
                eta_eff = &eta_eff_right_;
            else
                eta_eff = &eta_eff_left_;
        }
        else if(dir == DIRECTION::Y)
        {
            param.transSz_= grid_i_->local_x()-2;
            if(pol_i_ == POLARIZATION::EX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HZ )
                param.transSz_ -= 1;

            if(grid_i_->local_z() != 1)
            {
                transSz2 = grid_i_->local_z()-2;
                if( (pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EZ ) )
                    transSz2 -= 1;
            }
            param.stride_ = 1;

            if(pol_i_ == POLARIZATION::EZ || pol_i_ == POLARIZATION::EX)
                negC = false;

            if(pl)
                eta_eff = &eta_eff_top_;
            else
                eta_eff = &eta_eff_bot_;
        }
        else if(dir == DIRECTION::Z)
        {
            param.transSz_= grid_i_->local_x()-2;
            if(pol_i_ == POLARIZATION::EX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HZ )
                param.transSz_ -= 1;

            transSz2 = grid_i_->local_y()-2;
            if( (pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::EY ) && gridComm_.rank() == gridComm_.size() - 1)
                transSz2 -= 1;

            param.stride_ = 1;

            if(pol_i_ == POLARIZATION::EY || pol_i_ == POLARIZATION::EX)
                negC = false;

            if(pl)
                eta_eff = &eta_eff_front_;
            else
                eta_eff = &eta_eff_back_;
        }
        else
            throw std::logic_error("PMLs need to be in X,Y, or Z direction it can't be NONE");

        double distOff = 0.0;
        if(pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HZ)
            distOff = pl ? -0.5 : 0.5;

        int dist = 0;
        int ccStart = 0;
        if( (dir == DIRECTION::Z) && (pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EZ) )
            pl ? ccStart = 1 : dirMax -= 1;
        else if( dir == DIRECTION::Y && (gridComm_.rank() == gridComm_.size()-1) && (pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EY) )
            pl ? ccStart = 1 : dirMax -= 1;
        else if( dir == DIRECTION::X && (pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EX) )
            pl ? ccStart = 1 : dirMax -= 1;

        int zStart = 0;
        if(grid_i_->local_z() != 1)
            zStart = 1;

        for(int cc = ccStart; cc < dirMax; cc++)
        {
            if(dir == DIRECTION::X)
                dist = pmlEdge > nDir ? (grid_i_->x()-2*gridComm_.npX()) - (startPt - cc) : cc+startPt;
            else if(dir == DIRECTION::Y)
                dist = pmlEdge > nDir ? (grid_i_->y()-2*gridComm_.npY()) - (startPt - cc) : cc+startPt;
            else if(dir == DIRECTION::Z)
                dist = pmlEdge > nDir ? (grid_i_->z()-2*gridComm_.npZ()) - (startPt - cc) : cc+startPt;

            double sig = sigma( (static_cast<double>(dist) + distOff), static_cast<double>(nDir-1), eta_eff->at(dist));
            double kap = kappa( (static_cast<double>(dist) + distOff), static_cast<double>(nDir-1) );
            double a   = aVal ( (static_cast<double>(dist) + distOff), static_cast<double>(nDir-1) );

            param.b_    = b( sig, a, kap );
            param.c_    = c( sig, a, kap );
            param.cOff_ = c( sig, a, kap );
            negC ? param.c_ *= -1.0 : param.cOff_ *= -1.0;

            for(int jj = 0; jj < transSz2 ; ++jj)
            {
                if(dir == DIRECTION::X)
                {
                    if(grid_i_->local_z() != 1)
                    {
                        pl ? param.loc_ = {grid_i_->local_x()-2 - cc, jj+1, zStart} : param.loc_ = {cc+1, jj+1, zStart};
                        if(pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HZ)
                            pl ? param.locOff_ = {grid_i_->local_x()-1 - cc, 1+jj, zStart} : param.locOff_ = {cc+2, 1+jj, zStart};
                        else
                            pl ? param.locOff_ = {grid_i_->local_x()-3 - cc, 1+jj, zStart} : param.locOff_ = {cc  , 1+jj, zStart};
                    }
                    else
                    {
                        pl ? param.loc_ = {grid_i_->local_x()-2 - cc, 1, zStart+jj} : param.loc_ = {cc+1, 1, zStart+jj};
                        if(pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HZ)
                            pl ? param.locOff_ = {grid_i_->local_x()-1 - cc, 1, zStart+jj} : param.locOff_ = {cc+2, 1, zStart+jj};
                        else
                            pl ? param.locOff_ = {grid_i_->local_x()-3 - cc, 1, zStart+jj} : param.locOff_ = {cc  , 1, zStart+jj};
                    }
                }
                else if(dir == DIRECTION::Y)
                {
                    pl ? param.loc_ = {1, grid_i_->local_y()-2 - cc, zStart+jj} : param.loc_ = {1, cc+1, zStart+jj};
                    if(pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HZ)
                        pl ? param.locOff_ = {1, grid_i_->local_y()-1 - cc, zStart+jj} : param.locOff_ = {1, cc+2, zStart+jj};
                    else
                        pl ? param.locOff_ = {1, grid_i_->local_y()-3 - cc, zStart+jj} : param.locOff_ = {1, cc  , zStart+jj};
                }
                else if(dir == DIRECTION::Z)
                {
                    pl ? param.loc_ = {1, jj+1, grid_i_->local_z()-2 - cc} : param.loc_ = {1, jj+1, cc+1};
                    if(pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HZ)
                        pl ? param.locOff_ = {1, jj+1, grid_i_->local_z()-1 - cc} : param.locOff_ = {1, jj+1, cc+2};
                    else
                        pl ? param.locOff_ = {1, jj+1, grid_i_->local_z()-3 - cc} : param.locOff_ = {1, jj+1, cc  };
                }

                paramList.push_back(param);
            }
        }
        return paramList;
    }

    /**
     * @brief generate the lists that add psi to grid
     *
     * @param[in] xmin minimum in the x direction
     * @param[in] xmax maximum in the x direction
     * @param[in] ymin minimum in the y direction
     * @param[in] ymax maximum in the y direction
     * @param[in] zmin minimum in the z direction
     * @param[in] zmax maximum in the z direction
     * @param[in] stride stride for mkl operations
     * @param[in] physGrid physical grid for direction
     * @return ax lists to generate the lists that add psi to grid
     */
    std::vector<std::array<int,5>> getAxLists(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int stride, int_pgrid_ptr physGrid)
    {
        std::vector<std::array<int,5>> axLists;
        if(stride == 1)
        {
            for(int kk = zmin; kk < zmax; ++kk)
            {
                for(int jj = ymin; jj < ymax; ++jj)
                {
                    int ii = xmin;
                    while(ii < xmax)
                    {
                        int iistore = ii;
                        while (ii < xmax-1 && physGrid->point(ii,jj,kk) == physGrid->point(ii+1,jj,kk) )
                            ++ii;
                        axLists.push_back( {{ iistore, jj, kk, ii-iistore+1, static_cast<int>(physGrid->point(iistore,jj,kk)) }} );
                        ++ii;
                    }
                }
            }
        }
        else
        {
            for(int kk = zmin; kk < zmax; ++kk)
            {
                for(int ii = xmin; ii < xmax; ++ii)
                {
                    int jj = ymin;
                    while(jj < ymax)
                    {
                        int jjstore = jj;
                        while (jj < ymax-1 && physGrid->point(ii,jj,kk) == physGrid->point(ii,jj+1,kk) )
                            ++jj;
                        axLists.push_back( {{ ii, jjstore, kk, jj-jjstore+1, static_cast<int>(physGrid->point(ii,jjstore,kk)) }} );
                        ++jj;
                    }
                }
            }
        }
        return axLists;
    }

    /**
     * @brief      Generates the lists that add psi fields to grid_i
     * @details    no phys_Ez since there will always be a break at Ez
     *
     * @param[in]  dir      direction of pml
     * @param[in]  pl       true if top or right
     * @param[in]  phys_Ex  The physical ex grids
     * @param[in]  phys_Ey  The physical ey grids
     * @param[in]  objArr   list of objects in the cell
     *
     * @return     The grid up list.
     */
    std::vector<updateGridParams> getGridUpList(DIRECTION dir, bool pl,  int_pgrid_ptr phys_Ex, int_pgrid_ptr phys_Ey, std::vector<std::shared_ptr<Obj>> objArr)
    {
        std::vector<updateGridParams> paramList;
        updateGridParams param;

        bool negD = true;

        if(dir == DIRECTION::X)
        {
            param.stride_ = grid_i_->local_x()*grid_i_->local_z();
            int xmin, xmax, ymin, ymax, zmin, zmax;
            ymin = 1;
            ymax = grid_i_->local_y()-1;
            if( ( gridComm_.rank() == gridComm_.size()-1 ) && (pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EY) )
                ymax -= 1;
            if(pl)
            {
                xmin = grid_i_->local_x() - ln_vec_pl_[0] - 1;
                xmax = grid_i_->local_x() -1;
                if( pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::EX || pol_i_ == POLARIZATION::HY)
                    xmax -= 1;
            }
            else
            {
                xmin = 1;
                xmax = ln_vec_mn_[0];
                if( pol_i_ == POLARIZATION::EZ || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EY)
                    xmax += 1;
            }
            if(grid_i_->local_z() == 1)
            {
                zmin = 0;
                zmax = 1;
            }
            else
            {
                zmin = 1;
                zmax = grid_i_->local_z() - 1;
                if(pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EZ)
                    zmax -= 1;
            }
            std::vector<std::array<int,5>> axParams = getAxLists(xmin, xmax, ymin, ymax, zmin, zmax, param.stride_, phys_Ey);
            for(auto & axList : axParams)
            {
                param.loc_ = {axList[0], axList[1], axList[2]};
                param.nAx_ = axList[3];
                param.Db_    = dt_ / d_[0];
                if(pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::EY)
                    param.Db_*=-1.0;
                paramList.push_back(param);
            }
        }
        else if(dir == DIRECTION::Y)
        {
            param.stride_ = 1;
            int xmin, xmax, ymin, ymax, zmin, zmax;
            xmin = 1;
            xmax = grid_i_->local_x()-1;
            if( pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EX)
                xmax -= 1;
            if(pl)
            {
                ymin = grid_i_->local_y() - ln_vec_pl_[1] - 1;
                ymax = grid_i_->local_y() -1;
                if( pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EY)
                    ymax -= 1;
            }
            else
            {
                ymin = 1;
                ymax = ln_vec_mn_[1];
                if( pol_i_ == POLARIZATION::EZ || pol_i_ == POLARIZATION::EX || pol_i_ == POLARIZATION::HY)
                    ymax += 1;
            }
            if(grid_i_->local_z() == 1)
            {
                zmin = 0;
                zmax = 1;
            }
            else
            {
                zmin = 1;
                zmax = grid_i_->local_z() - 1;
                if(pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EZ)
                    zmax -= 1;
            }
            std::vector<std::array<int,5>> axParams = getAxLists(xmin, xmax, ymin, ymax, zmin , zmax, param.stride_, phys_Ex);
            for(auto & axList : axParams)
            {
                param.loc_ = {axList[0], axList[1], axList[2]};
                param.nAx_ = axList[3];
                param.Db_    = dt_ / d_[1];
                if(pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EZ)
                    param.Db_*=-1.0;
                paramList.push_back(param);
            }
        }
        else if(dir == DIRECTION::Z)
        {
            param.stride_ = 1;
            int xmin, xmax, ymin, ymax, zmin, zmax;
            ymin = 1;
            ymax = grid_i_->local_y()-1;
            if( ( gridComm_.rank() == gridComm_.size()-1 ) && (pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EY) )
                ymax -= 1;
            if(pl)
            {
                zmin = grid_i_->local_z() - ln_vec_pl_[2] - 1;
                zmax = grid_i_->local_z() -1;
                if( pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EZ)
                    zmax -= 1;
            }
            else
            {
                zmin = 1;
                zmax = ln_vec_mn_[2];
                if( pol_i_ == POLARIZATION::EX || pol_i_ == POLARIZATION::EY || pol_i_ == POLARIZATION::HZ)
                    zmax += 1;
            }
            xmin = 1;
            xmax = grid_i_->local_x() - 1;
            if(pol_i_ == POLARIZATION::EX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HZ)
                xmax -= 1;
            std::vector<std::array<int,5>> axParams = getAxLists(xmin, xmax, ymin, ymax, zmin, zmax, param.stride_, phys_Ey);
            for(auto & axList : axParams)
            {
                param.loc_ = {axList[0], axList[1], axList[2]};
                param.nAx_ = axList[3];
                param.Db_    = dt_ / d_[2];
                if(pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EX)
                    param.Db_*=-1.0;
                paramList.push_back(param);
            }
        }
        else
            throw std::logic_error("Direction for grid updates has to be X, Y or Z not NONE");

        return paramList;
    }

    /**
     * @brief      Creates the psi update lists and the addition lists
     * @details    No phys_Ez since there will always be a break at the Ez grid
     *
     * @param[in]  phys_Ex  The physical ex grids
     * @param[in]  phys_Ey  The physical ey grids
     * @param[in]  objArr   The lists of objects in the grid
     */
    void initalizeGrid(int_pgrid_ptr phys_Ex, int_pgrid_ptr phys_Ey, std::vector<std::shared_ptr<Obj>> objArr)
    {
        std::vector<updatePsiParams> psiTemp;
        std::vector<updateGridParams> gridTemp;
        int dirMax;
        if(ln_vec_pl_[1] > 0 && i_ != DIRECTION::Y)
        {
            if(grid_k_ && j_ == DIRECTION::Y)
            {
                psiTemp  = getPsiUpList (DIRECTION::Y, true, grid_i_->procLoc()[1] + grid_i_->local_y()-2,  grid_i_->y()-2*gridComm_.npY(), n_vec_[1], ln_vec_pl_[1]);
                gridTemp = getGridUpList(DIRECTION::Y, true, phys_Ex, phys_Ey, objArr);
                updateListPsi_j_.reserve(updateListPsi_j_.size() + psiTemp.size());
                updateListPsi_j_.insert(updateListPsi_j_.end(), psiTemp.begin(), psiTemp.end());
                updateListGrid_j_.reserve(updateListGrid_j_.size() + gridTemp.size());
                updateListGrid_j_.insert(updateListGrid_j_.end(), gridTemp.begin(), gridTemp.end());
            }
            else if(grid_j_ && k_ == DIRECTION::Y)
            {
                psiTemp  = getPsiUpList (DIRECTION::Y, true, grid_i_->procLoc()[1] + grid_i_->local_y()-2,  grid_i_->y()-2*gridComm_.npY(), n_vec_[1], ln_vec_pl_[1]);
                gridTemp = getGridUpList(DIRECTION::Y, true, phys_Ex, phys_Ey, objArr);
                updateListPsi_k_.reserve(updateListPsi_k_.size() + psiTemp.size());
                updateListPsi_k_.insert(updateListPsi_k_.end(), psiTemp.begin(), psiTemp.end());
                updateListGrid_k_.reserve(updateListGrid_k_.size() + gridTemp.size());
                updateListGrid_k_.insert(updateListGrid_k_.end(), gridTemp.begin(), gridTemp.end());
            }
        }
        if(ln_vec_mn_[1] > 0 && i_ != DIRECTION::Y)
        {
            if(grid_k_ && j_ == DIRECTION::Y)
            {
                psiTemp  = getPsiUpList (DIRECTION::Y, false, grid_i_->procLoc()[1], n_vec_[1], n_vec_[1], ln_vec_mn_[1]);
                gridTemp = getGridUpList(DIRECTION::Y, false, phys_Ex, phys_Ey, objArr);

                updateListPsi_j_.reserve(updateListPsi_j_.size() + psiTemp.size());
                updateListPsi_j_.insert(updateListPsi_j_.end(), psiTemp.begin(), psiTemp.end());
                updateListGrid_j_.reserve(updateListGrid_j_.size() + gridTemp.size());
                updateListGrid_j_.insert(updateListGrid_j_.end(), gridTemp.begin(), gridTemp.end());
            }
            else if(grid_j_ && k_ == DIRECTION::Y)
            {
                psiTemp  = getPsiUpList (DIRECTION::Y, false, grid_i_->procLoc()[1], n_vec_[1], n_vec_[1], ln_vec_mn_[1]);
                gridTemp = getGridUpList(DIRECTION::Y, false, phys_Ex, phys_Ey, objArr);

                updateListPsi_k_.reserve(updateListPsi_k_.size() + psiTemp.size());
                updateListPsi_k_.insert(updateListPsi_k_.end(), psiTemp.begin(), psiTemp.end());
                updateListGrid_k_.reserve(updateListGrid_k_.size() + gridTemp.size());
                updateListGrid_k_.insert(updateListGrid_k_.end(), gridTemp.begin(), gridTemp.end());

            }
        }
        if(ln_vec_mn_[0] > 0 && i_ != DIRECTION::X)
        {
            if(grid_k_ && j_ == DIRECTION::X)
            {
                psiTemp  = getPsiUpList (DIRECTION::X, false, grid_i_->procLoc()[0], n_vec_[0], n_vec_[0], ln_vec_mn_[0]);
                gridTemp = getGridUpList(DIRECTION::X, false, phys_Ex, phys_Ey, objArr);

                updateListPsi_j_.reserve(updateListPsi_j_.size() + psiTemp.size());
                updateListPsi_j_.insert(updateListPsi_j_.end(), psiTemp.begin(), psiTemp.end());
                updateListGrid_j_.reserve(updateListGrid_j_.size() + gridTemp.size());
                updateListGrid_j_.insert(updateListGrid_j_.end(), gridTemp.begin(), gridTemp.end());
            }
            else if(grid_j_ && k_ == DIRECTION::X)
            {
                psiTemp  = getPsiUpList (DIRECTION::X, false, grid_i_->procLoc()[0], n_vec_[0], n_vec_[0], ln_vec_mn_[0]);
                gridTemp = getGridUpList(DIRECTION::X, false, phys_Ex, phys_Ey, objArr);

                updateListPsi_k_.reserve(updateListPsi_k_.size() + psiTemp.size());
                updateListPsi_k_.insert(updateListPsi_k_.end(), psiTemp.begin(), psiTemp.end());
                updateListGrid_k_.reserve(updateListGrid_k_.size() + gridTemp.size());
                updateListGrid_k_.insert(updateListGrid_k_.end(), gridTemp.begin(), gridTemp.end());
            }
        }
        if(ln_vec_pl_[0] > 0 && i_ != DIRECTION::X)
        {
            if(grid_k_ && j_ == DIRECTION::X)
            {
                psiTemp  = getPsiUpList (DIRECTION::X, true, grid_i_->procLoc()[0] + grid_i_->local_x()-2,  grid_i_->x()-2*gridComm_.npX(), n_vec_[0], ln_vec_pl_[0]);
                gridTemp = getGridUpList(DIRECTION::X, true, phys_Ex, phys_Ey, objArr);

                updateListPsi_j_.reserve(updateListPsi_j_.size() + psiTemp.size());
                updateListPsi_j_.insert(updateListPsi_j_.end(), psiTemp.begin(), psiTemp.end());
                updateListGrid_j_.reserve(updateListGrid_j_.size() + gridTemp.size());
                updateListGrid_j_.insert(updateListGrid_j_.end(), gridTemp.begin(), gridTemp.end());
            }
            else if(grid_j_ && k_ == DIRECTION::X)
            {
                psiTemp  = getPsiUpList (DIRECTION::X, true, grid_i_->procLoc()[0] + grid_i_->local_x()-2,  grid_i_->x()-2*gridComm_.npX(), n_vec_[0], ln_vec_pl_[0]);
                gridTemp = getGridUpList(DIRECTION::X, true, phys_Ex, phys_Ey, objArr);

                updateListPsi_k_.reserve(updateListPsi_k_.size() + psiTemp.size());
                updateListPsi_k_.insert(updateListPsi_k_.end(), psiTemp.begin(), psiTemp.end());
                updateListGrid_k_.reserve(updateListGrid_k_.size() + gridTemp.size());
                updateListGrid_k_.insert(updateListGrid_k_.end(), gridTemp.begin(), gridTemp.end());
            }
        }


        if(ln_vec_mn_[2] > 0 && i_ != DIRECTION::Z)
        {
            if(grid_k_ && j_ == DIRECTION::Z)
            {
                psiTemp  = getPsiUpList (DIRECTION::Z, false, grid_i_->procLoc()[2], n_vec_[2], n_vec_[2], ln_vec_mn_[2]);
                gridTemp = getGridUpList(DIRECTION::Z, false, phys_Ex, phys_Ey, objArr);

                updateListPsi_j_.reserve(updateListPsi_j_.size() + psiTemp.size());
                updateListPsi_j_.insert(updateListPsi_j_.end(), psiTemp.begin(), psiTemp.end());
                updateListGrid_j_.reserve(updateListGrid_j_.size() + gridTemp.size());
                updateListGrid_j_.insert(updateListGrid_j_.end(), gridTemp.begin(), gridTemp.end());
            }
            else if(grid_j_ && k_ == DIRECTION::Z)
            {
                psiTemp  = getPsiUpList (DIRECTION::Z, false, grid_i_->procLoc()[2], n_vec_[2], n_vec_[2], ln_vec_mn_[2]);
                gridTemp = getGridUpList(DIRECTION::Z, false, phys_Ex, phys_Ey, objArr);

                updateListPsi_k_.reserve(updateListPsi_k_.size() + psiTemp.size());
                updateListPsi_k_.insert(updateListPsi_k_.end(), psiTemp.begin(), psiTemp.end());
                updateListGrid_k_.reserve(updateListGrid_k_.size() + gridTemp.size());
                updateListGrid_k_.insert(updateListGrid_k_.end(), gridTemp.begin(), gridTemp.end());
            }
        }
        if(ln_vec_pl_[2] > 0 && i_ != DIRECTION::Z)
        {
            if(grid_k_ && j_ == DIRECTION::Z)
            {
                psiTemp  = getPsiUpList (DIRECTION::Z, true, grid_i_->procLoc()[2] + grid_i_->local_z()-2,  grid_i_->z()-2*gridComm_.npZ(), n_vec_[2], ln_vec_pl_[2]);
                gridTemp = getGridUpList(DIRECTION::Z, true, phys_Ex, phys_Ey, objArr);

                updateListPsi_j_.reserve(updateListPsi_j_.size() + psiTemp.size());
                updateListPsi_j_.insert(updateListPsi_j_.end(), psiTemp.begin(), psiTemp.end());
                updateListGrid_j_.reserve(updateListGrid_j_.size() + gridTemp.size());
                updateListGrid_j_.insert(updateListGrid_j_.end(), gridTemp.begin(), gridTemp.end());
            }
            else if(grid_j_ && k_ == DIRECTION::Z)
            {
                psiTemp  = getPsiUpList (DIRECTION::Z, true, grid_i_->procLoc()[0] + grid_i_->local_z()-2,  grid_i_->z()-2*gridComm_.npZ(), n_vec_[2], ln_vec_pl_[2]);
                gridTemp = getGridUpList(DIRECTION::Z, true, phys_Ex, phys_Ey, objArr);

                updateListPsi_k_.reserve(updateListPsi_k_.size() + psiTemp.size());
                updateListPsi_k_.insert(updateListPsi_k_.end(), psiTemp.begin(), psiTemp.end());
                updateListGrid_k_.reserve(updateListGrid_k_.size() + gridTemp.size());
                updateListGrid_k_.insert(updateListGrid_k_.end(), gridTemp.begin(), gridTemp.end());
            }
        }
    }
    /**
     * @brief updates the girds
     */
    void updateGrid()
    {
        upPsi_j_(updateListGrid_j_, updateListPsi_j_, grid_i_, psi_j_, grid_k_);
        upPsi_k_(updateListGrid_k_, updateListPsi_k_, grid_i_, psi_k_, grid_j_);
    }

    /**
     * @brief      calculates b from Taflove chapter 7
     *
     * @param[in]  sig   sigma value at that point
     * @param[in]  a     value at that point
     * @param[in]  kap   kappa value at that point
     *
     * @return     b
     */
    inline double b(double sig, double a, double kap)
    {
        return std::exp(-1.0 * (sig/kap + a) * dt_);
    }

    /**
     * @brief      calculates c from Taflove chapter 7
     *
     * @param[in]  sigma value at taht point
     * @param[in]  value at that point
     * @param[in]  kappa value at that point
     *
     * @return     c
     */
    inline double c(double sig, double a, double kap)
    {
        return (sig == 0 && a == 0) ? 0 : sig / (sig*kap + std::pow(kap,2.0)*a) * (std::exp(-1.0 * (sig/kap + a) * dt_) - 1.0);
    }

    /**
     * @brief      calculates a from Taflove chapter 7
     *
     * @param[in]  ii     current point location
     * @param[in]  iiMax  maximum ii value
     *
     * @return     a
     */
    inline double aVal(double ii, double iiMax)
    {
        return (0.0 <= ii && ii <= iiMax) ? aMax_ * std::pow( (iiMax - ii) / iiMax, ma_) : 0.0;
    }

    /**
     * @brief      calculates kappa from Taflove chapter 7
     *
     * @param[in]  ii current point location
     * @param[in]  iiMax maximum ii value
     *
     * @return     kappa
     */
    inline double kappa(double ii, double iiMax)
    {
        return (0.0 <= ii && ii <= iiMax) ? 1.0 + (kappaMax_ - 1.0) * std::pow((iiMax - ii) / iiMax , m_) : 1.0;
    }

    /**
     * @brief      {calculates sigma from Taflove chapter 7
     *
     * @param[in]  ii       ii current point location
     * @param[in]  iiMax    iiMax maximum ii value
     * @param[in]  eta_eff  The effective $\\eta$ for the point
     *
     * @return     sigma
     */
    inline double sigma(double ii, double iiMax, double eta_eff)
    {
        return (0.0 <= ii && ii <= iiMax) ? sigmaMax_ / eta_eff * pow((iiMax - ii) / iiMax , m_) : 0.0;
    }
    /**
     * @brief      ln_vec_right accessor function
     *
     * @return     ln_vec_right
     */
    inline int lnx_right() {return ln_vec_pl_[0];}
    /**
     * @brief      ln_vec_left accessor function
     *
     * @return     ln_vec_left
     */
    inline int lnx_left()  {return ln_vec_mn_[0] ;}
    /**
     * @brief      ln_vec_top accessor function
     *
     * @return     ln_vec_top
     */
    inline int lny_top()   {return ln_vec_pl_[1]  ;}
    /**
     * @brief      ln_vec_bot accessor function
     *
     * @return     ln_vec_bot
     */
    inline int lny_bot()   {return ln_vec_mn_[1]  ;}

     /**
     * @brief      ln_vec_top accessor function
     *
     * @return     ln_vec_top
     */
    inline int lnz_back()   {return ln_vec_pl_[2]  ;}
    /**
     * @brief      ln_vec_bot accessor function
     *
     * @return     ln_vec_bot
     */
    inline int lnz_front()   {return ln_vec_mn_[2]  ;}
};
/**
 * @brief functions for psi updates
 *
 */
namespace pmlUpdateFxnReal
{
    /**
     * @brief      Adds psi to grid_i_
     *
     * @param      gridParamList  list of parameters to update grid_i_
     * @param      psiParamList   ist of parameters to update Psi field
     * @param[in]  grid_i         grid_i_
     * @param[in]  psi            psi field to update
     * @param[in]  grid           grid used to update psi
     */
    void addPsi(std::vector<updateGridParams> &gridParamList, std::vector<updatePsiParams> &psiParamList, real_pgrid_ptr grid_i, real_pgrid_ptr psi, real_pgrid_ptr grid);

    /**
     * @brief      update psi field
     *
     * @param      list to loop over to update psi
     * @param[in]  psi field to update
     * @param[in]  grid used to update psi
     */
    void updatePsiField(std::vector<updatePsiParams>&paramList, real_pgrid_ptr psi , real_pgrid_ptr grid);
}
namespace pmlUpdateFxnCplx
{
    /**
     * @brief      Adds psi to grid_i_
     *
     * @param      gridParamList  list of parameters to update grid_i_
     * @param      psiParamList   ist of parameters to update Psi field
     * @param[in]  grid_i         grid_i_
     * @param[in]  psi            psi field to update
     * @param[in]  grid           grid used to update psi
     */
    void addPsi(std::vector<updateGridParams> &gridParamList, std::vector<updatePsiParams> &psiParamList, cplx_pgrid_ptr grid_i, cplx_pgrid_ptr psi, cplx_pgrid_ptr grid);

    /**
     * @brief      update psi field
     *
     * @param      list to loop over to update psi
     * @param[in]  psi field to update
     * @param[in]  grid used to update psi
     */
    void updatePsiField(std::vector<updatePsiParams>&paramList , cplx_pgrid_ptr psi , cplx_pgrid_ptr grid);
}

class parallelCPMLReal : public parallelCPML<double>
{
public:
    /**
     * @brief      { function_description }
     *
     * @param[in]  gridComm  mpi communicator
     * @param[in]  weights   The weights
     * @param[in]  grid_i    shared pointer to grid_i (field PML is being applied to)
     * @param[in]  grid_j    shared pointer to grid_j (field polarized in the j direction of the PML; nullptr if none)
     * @param[in]  grid_k    shared pointer to grid_k (field polarized in the k direction of the PML; nullptr if none)
     * @param[in]  pol_i     polarization of grid_o
     * @param[in]  n_vec     Vector storing thickness of the PMLs in all directions
     * @param[in]  m         scalling factor for sigma
     * @param[in]  ma        scaling factor for a
     * @param[in]  aMax      max a value
     * @param[in]  d         vector storing the step sizes in all directions
     * @param[in]  dt        time step
     * @param[in]  phys_Ex   The physical ex grid
     * @param[in]  phys_Ey   The physical ey grid
     * @param[in]  objArr    The object arr
     */
    parallelCPMLReal(mpiInterface gridComm, std::vector<real_grid_ptr> weights, real_pgrid_ptr grid_i, real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, POLARIZATION pol_i, std::array<int,3> n_vec, double m, double ma, double aMax, std::array<double,3> d, double dt, int_pgrid_ptr phys_Ex, int_pgrid_ptr phys_Ey, std::vector<std::shared_ptr<Obj>> objArr);
};

class parallelCPMLCplx : public parallelCPML<cplx>
{
public:
    /**
     * @brief      { function_description }
     *
     * @param[in]  gridComm  mpi communicator
     * @param[in]  weights   The weights
     * @param[in]  grid_i    shared pointer to grid_i (field PML is being applied to)
     * @param[in]  grid_j    shared pointer to grid_j (field polarized in the j direction of the PML; nullptr if none)
     * @param[in]  grid_k    shared pointer to grid_k (field polarized in the k direction of the PML; nullptr if none)
     * @param[in]  pol_i     polarization of grid_o
     * @param[in]  n_vec     Vector storing thickness of the PMLs in all directions
     * @param[in]  m         scalling factor for sigma
     * @param[in]  ma        scaling factor for a
     * @param[in]  aMax      max a value
     * @param[in]  d         vector storing the step sizes in all directions
     * @param[in]  dt        time step
     * @param[in]  phys_Ex   The physical ex grid
     * @param[in]  phys_Ey   The physical ey grid
     * @param[in]  objArr    The object arr
     */
    parallelCPMLCplx(mpiInterface gridComm, std::vector<real_grid_ptr> weights, std::shared_ptr<parallelGrid<cplx > > grid_i, std::shared_ptr<parallelGrid<cplx > > grid_j, std::shared_ptr<parallelGrid<cplx > > grid_k, POLARIZATION pol_i, std::array<int,3> n_vec, double m, double ma, double aMax, std::array<double,3> d, double dt, int_pgrid_ptr phys_Ex, int_pgrid_ptr phys_Ey, std::vector<std::shared_ptr<Obj>> objArr);
};
#endif