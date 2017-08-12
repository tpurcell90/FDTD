#ifndef PARALLEL_FDTD_PML
#define PARALLEL_FDTD_PML

#include <src/OBJECTS/Obj.hpp>

/**
 * @brief      parameters to update the $\psi$ fields for CPMLs
 */
struct updatePsiParams
{
    int transSz_; //!< size in the transverse direction (how many points to do blas operations on)
    int stride_; //!< stride for the blas operations
    double b_; //!< the b parameter as defined in chapter 7 of Taflove
    double c_; //!< the c parameter as defined in chapter 7 of Taflove
    double cOff_; //!< the c parameter as defined in chapter 7 of Taflove, but with the opposite sign of c
    std::array<int,3> loc_; //!< the starting point of the blas operations
    std::array<int,3> locOff_; //!< the starting point of the blas operations offset by the correct value
};
/**
 * @brief      Parameters to update the fields using the $\psi$ fields
 */
struct updateGridParams
{
    int nAx_; //!< size of blas operator
    int stride_; //!< stride for blas operator
    double Db_; //!< psi factor add on prefactor
    std::array<int,3> loc_; //!< starting point of MKL operations
};


template <typename T> class parallelCPML
{
    typedef std::shared_ptr<parallelGrid<T>> pgrid_ptr;
protected:
    std::shared_ptr<mpiInterface> gridComm_; //!< MPI Communicator
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
    parallelCPML(std::shared_ptr<mpiInterface> gridComm, std::vector<real_grid_ptr> weights, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k, POLARIZATION pol_i, std::array<int,3> n_vec, double m, double ma, double aMax, std::array<double,3> d, double dt, int_pgrid_ptr physGrid, std::vector<std::shared_ptr<Obj>> objArr) :
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
        // Determine the PML's i,j,k directions based on the polarization of the field, direction i shares the direction of the field polarization
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
        // In 2D z direction is assumed to be isotropic, but 3D it is not
        if(grid_k_ && (grid_i_->local_z() != 1 || j_ != DIRECTION::Z) )
        {
            psi_j_ = std::make_shared<parallelGrid<T>>(gridComm, grid_k_->PBC(), weights, grid_k_->n_vec(), grid_k_->d() );
        }
        else
            psi_j_ = nullptr;

        // In 2D z direction is assumed to be isotropic, but in 3D it is not
        if(grid_j_ && (grid_i_->local_z() != 1 || k_ != DIRECTION::Z) )
        {
            psi_k_ = std::make_shared<parallelGrid<T>>(gridComm, grid_j_->PBC(), weights, grid_j_->n_vec(), grid_j_->d() );
        }
        else
            psi_k_ = nullptr;

        if(!psi_k_ && ! psi_j_)
            throw std::logic_error("Both PML auxiliary fields are undefined this likely means that there is an issue with grid_i or grid_k");

        findLnVecs();

        // Calculate the mean $\eta_r$ values for each PML
        genEtaEff(physGrid, DIRECTION::X, false, objArr, eta_eff_left_);
        genEtaEff(physGrid, DIRECTION::X,  true, objArr, eta_eff_right_);

        genEtaEff(physGrid, DIRECTION::Y, false, objArr, eta_eff_bot_);
        genEtaEff(physGrid, DIRECTION::Y,  true, objArr, eta_eff_top_);


        if(grid_i_->z() != 1)
        {
            genEtaEff(physGrid, DIRECTION::Z, false, objArr, eta_eff_back_);
            genEtaEff(physGrid, DIRECTION::Z,  true, objArr, eta_eff_front_);
        }
        // Initialize all update lists
        initalizeLists(ln_vec_mn_[0], DIRECTION::X, false, grid_i_->procLoc()[0]                       , n_vec_[0]                      , n_vec_[0], ln_vec_mn_[0]);
        initalizeLists(ln_vec_pl_[0], DIRECTION::X,  true, grid_i_->procLoc()[0] + grid_i_->local_x()-2, grid_i_->x()-2*gridComm_->npX(), n_vec_[0], ln_vec_pl_[0]);
        initalizeLists(ln_vec_mn_[1], DIRECTION::Y, false, grid_i_->procLoc()[1]                       , n_vec_[1]                      , n_vec_[1], ln_vec_mn_[1]);
        initalizeLists(ln_vec_pl_[1], DIRECTION::Y,  true, grid_i_->procLoc()[1] + grid_i_->local_y()-2, grid_i_->y()-2*gridComm_->npY(), n_vec_[1], ln_vec_pl_[1]);
        initalizeLists(ln_vec_mn_[2], DIRECTION::Z, false, grid_i_->procLoc()[2]                       , n_vec_[2]                      , n_vec_[2], ln_vec_mn_[2]);
        initalizeLists(ln_vec_pl_[2], DIRECTION::Z,  true, grid_i_->procLoc()[2] + grid_i_->local_z()-2, grid_i_->z()-2*gridComm_->npZ(), n_vec_[2], ln_vec_pl_[2]);
    }

    /**
     * @brief      calculates the effective mean value of $\eta_r$ for a plane in the FDTD cell
     *
     * @param[in]  physVals        Slice from the object map grid
     * @param[in]  objArr          The object arr
     * @param[in]  normEtaSumDiff  The value to subtract from physVals.size() due to the -1 border
     *
     * @return     value of $\sqrt{\varepsilon_{r,eff} \mu_{r,eff}}$
     */
    double genEtaEffVal(const std::vector<int>& physVals, std::vector<std::shared_ptr<Obj>> objArr, double normEtaSumDiff)
    {
        double eps_sum = 0.0;
        double mu_sum = 0.0;
        // Sum up all the dielectric constants in the plane, if it is a boundary region don't add it
        for(auto& val : physVals)
        {
            if(val != -1)
            {
                eps_sum += objArr[val]->epsInfty() / (physVals.size() - normEtaSumDiff);
                mu_sum += objArr[val]->muInfty() / (physVals.size() - normEtaSumDiff);
            }
        }
        // Find the average dielectric constant if outside the boundary regions
        return eps_sum * mu_sum;
    }

    /**
     * @brief      For the entire thickness of the PML find the values of the effective wave impedence
     *
     * @param[in]  physGrid      The object map grid
     * @param[in]  planeNormDir  Direction of the normal vector to the PML
     * @param[in]  pl            True if a pl PML
     * @param[in]  objArr        The object arr
     * @param      eta_eff       Vector storing the effective wave impedance for the PML
     */
    void genEtaEff(int_pgrid_ptr physGrid, DIRECTION planeNormDir, bool pl, std::vector<std::shared_ptr<Obj>> objArr, std::vector<double>& eta_eff)
    {
        // Find cor_ii, jj, kk based on the the normal vector of the PML planes
        // npII is the number of processes in the ii direction
        // normEtaSumDiff removes the -1 boundary from the object
        int cor_ii = 0, cor_jj = 0,  cor_kk = 0, npII = 0;
        double normEtaSumDiff = 0.0;
        if(planeNormDir == DIRECTION::X)
        {
            cor_ii = 0, cor_jj = 1, cor_kk = 2;
            npII = gridComm_->npX();
            //For 3D: total elements in boundary = 2*(number processors in J dir * size in z + number of processors in K dir * size in j - 4 * num procs j * num procs k) last term to remove corners
            normEtaSumDiff = physGrid->local_z() != 1 ? 2*(gridComm_->npY()*physGrid->z() + gridComm_->npZ()*physGrid->y() ) - 4*(gridComm_->npY() * gridComm_->npZ()) : 2*gridComm_->npY();
        }
        else if(planeNormDir == DIRECTION::Y)
        {
            cor_ii = 1, cor_jj = 2, cor_kk = 0;
            npII = gridComm_->npY();
            //For 3D: total elements in boundary = 2*(number processors in J dir * size in z + number of processors in K dir * size in j - 4 * num procs j * num procs k) last term to remove corners
            normEtaSumDiff = physGrid->local_z() != 1 ? 2*(gridComm_->npZ()*physGrid->x() + gridComm_->npX()*physGrid->z() ) - 4*(gridComm_->npZ() * gridComm_->npX()) : (2*( gridComm_->npX() ) );
        }
        else if(planeNormDir == DIRECTION::Z)
        {
            cor_ii = 2, cor_jj = 0, cor_kk = 1;
            npII = gridComm_->npZ();
            //For 3D: total elements in boundary = 2*(number processors in J dir * size in z + number of processors in K dir * size in j - 4 * num procs j * num procs k) last term to remove corners
            normEtaSumDiff = 2*(gridComm_->npY()*physGrid->x() + gridComm_->npX()*physGrid->y() ) - 4*(gridComm_->npY() * gridComm_->npX());
        }
        // find values to loop over
        int min = pl ? physGrid->n_vec()[cor_ii] - (n_vec_[cor_ii] + npII*2) : 0;
        int max = pl ? physGrid->n_vec()[cor_ii] - 1 : n_vec_[cor_ii] + npII*2;

        // loop over those values
        for(int ii = min; ii < max; ++ii)
        {
            double etaVal = 0.0;
            if(planeNormDir == DIRECTION::X)
                etaVal = genEtaEffVal( physGrid->getPlaneYZ(ii), objArr, normEtaSumDiff);
            else if(planeNormDir == DIRECTION::Y)
                etaVal = genEtaEffVal( physGrid->getPlaneXZ(ii), objArr, normEtaSumDiff);
            else if(planeNormDir == DIRECTION::Z)
            {
                etaVal = genEtaEffVal( physGrid->getPlaneXY(ii), objArr, normEtaSumDiff);
            }
            // Find the average dielectric constant if outside the boundary regions
            if(eta_eff.size() < n_vec_[cor_ii] && etaVal != 0.0)
                eta_eff.push_back( etaVal );
        }
    }

    /**
     * @brief      Finds each process's local PML thicknesses
     */
    void findLnVecs()
    {
        // 1) Does the process inside the left PML? yes-> 2) Does the PML extend outside the process's region yes-> ln_vec_mn_[0] = local size of process, no-> ln_vec_mn_[0] = PML thickness-where the process started
        if(grid_i_->procLoc()[0] < n_vec_[0])
            ln_vec_mn_[0] = grid_i_->procLoc()[0] + grid_i_->local_x()-2 < n_vec_[0] ? grid_i_->local_x()-2 : n_vec_[0] - grid_i_->procLoc()[0];

        // 1) Does the right PML begin before the process's region end? yes-> 2) Does the process begin after the PML begins? yes->ln_vec_pl_[0] = local size of process, no-> ln_vec_pl_[0] = the process's size - PML thickness
        if(grid_i_->procLoc()[0] + grid_i_->local_x()-2 > grid_i_->x() - 2*gridComm_->npX() - n_vec_[0])
            ln_vec_pl_[0] = grid_i_->procLoc()[0] > grid_i_->x() - 2*gridComm_->npX() - n_vec_[0] ? grid_i_->local_x()-2 : grid_i_->procLoc()[0] + grid_i_->local_x() - 2 - (grid_i_->x() - 2*gridComm_->npX() - n_vec_[0]);;

        // 1) Does the process inside the bottom PML? yes-> 2) Does the PML extend outside the process's region yes-> ln_vec_mn_[1] = local size of process, no-> ln_vec_mn_[1] = PML thickness-where the process started
        if(grid_i_->procLoc()[1] < n_vec_[1])
            ln_vec_mn_[1] = grid_i_->procLoc()[1] + grid_i_->local_y()-2 < n_vec_[1] ? grid_i_->local_y()-2 : n_vec_[1] - grid_i_->procLoc()[1];

        // 1) Does the top PML begin before the process's region end? yes-> 2) Does the process begin after the PML begins? yes->ln_vec_pl_[1] = local size of process, no-> ln_vec_pl_[1] = the process's size - PML thickness
        if(grid_i_->procLoc()[1] + grid_i_->local_y()-2 > grid_i_->y() - gridComm_->npY()*2 - n_vec_[1])
            ln_vec_pl_[1] = grid_i_->procLoc()[1] > grid_i_->y() - gridComm_->npY()*2 - n_vec_[1] ? grid_i_->local_y()-2 : grid_i_->procLoc()[1] + grid_i_->local_y() - 2 - (grid_i_->y() - 2*gridComm_->npY() - n_vec_[1]);

        if(grid_i_->local_z() != 1)
        {
            // 1) Does the process inside the back PML? yes-> 2) Does the PML extend outside the process's region yes-> ln_vec_mn_[2] = local size of process, no-> ln_vec_mn_[2] = PML thickness-where the process started
            if(grid_i_->procLoc()[2] < n_vec_[2])
                ln_vec_mn_[2] = grid_i_->procLoc()[2] + grid_i_->local_z()-2 < n_vec_[2] ? grid_i_->local_z()-2 : n_vec_[2] - grid_i_->procLoc()[2];

            // 1) Does the front PML begin before the process's region end? yes-> 2) Does the process begin after the PML begins? yes->ln_vec_pl_[2] = local size of process, no-> ln_vec_pl_[2] = the process's size - PML thickness
            if(grid_i_->procLoc()[2] + grid_i_->local_z()-2 > grid_i_->z() - gridComm_->npZ()*2 - n_vec_[2])
                ln_vec_pl_[2] = grid_i_->procLoc()[2] > grid_i_->z() - gridComm_->npZ()*2 - n_vec_[2] ? grid_i_->local_z()-2 : grid_i_->procLoc()[2] + grid_i_->local_z() - 2 - (grid_i_->z() - 2*gridComm_->npZ() - n_vec_[2]);
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

        int cor_trans1 = -1, cor_trans2 = -1, cor_norm = -1;
        int trans1FieldOff = 0, trans2FieldOff = 0;
        int ccStart = 0; // where to start the loop over pmls
        // Determine the coordinate system for the PML update lists compared to x,y,z
        if(dir == DIRECTION::X)
        {
            // 3D calculation do blas operations over z coordinate; for 2D do blas operations over y coordinate
            cor_norm = 0;
            cor_trans1 = (grid_i_->local_z() == 1) ? 1 : 2;
            cor_trans2 = (grid_i_->local_z() == 1) ? 2 : 1;

            param.stride_ = grid_i_->local_x();
            // What eta effective vector should be used
            eta_eff = !pl ? &eta_eff_left_ : &eta_eff_right_;
            // Determines sign of the spatial derivative
            negC = (pol_i_ == POLARIZATION::EZ || pol_i_ == POLARIZATION::EY) ? false : true;

            // Ey, Hx, and Hz fields have one less point in the y direction
            if( (pol_i_ == POLARIZATION::EY || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HZ ) && gridComm_->rank() == gridComm_->size() - 1)
                (grid_i_->local_z() == 1) ? trans1FieldOff = 1 : trans2FieldOff = 1;

            // Ez, Hx, and Hy fields have one less point in the z direction
            if( (pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EZ ) )
                (grid_i_->local_z() == 1) ? trans2FieldOff = 1 : trans1FieldOff = 1;

            // PML start / Max conditions depending on where in the map they are
            if( pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EX)
                pl ? ccStart = 1 : dirMax -= 1;
        }
        else if(dir == DIRECTION::Y)
        {
            cor_norm = 1;
            cor_trans1 = 0;
            cor_trans2 = 2;

            param.stride_ = 1;
            // What eta effective vector should be used
            eta_eff = !pl ? &eta_eff_bot_ : &eta_eff_top_;
            // Determines sign of the spatial derivative
            negC = (pol_i_ == POLARIZATION::EZ || pol_i_ == POLARIZATION::EX) ? false : true;
            // Ex, Hy, and Hz fields have one less point in the x direction
            if(pol_i_ == POLARIZATION::EX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HZ )
                trans1FieldOff = 1;
            // Ez, Hx, and Hy fields have one less point in the z direction
            if( (pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EZ ) )
                trans2FieldOff = 1;
            // Hx, Ey, and Hz fields all have one less point in the z direction
            if( ( gridComm_->rank() == gridComm_->size()-1) && (pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EY) )
                pl ? ccStart = 1 : dirMax -= 1;
        }
        else if(dir == DIRECTION::Z)
        {
            cor_norm = 2;
            cor_trans1 = 0;
            cor_trans2 = 1;
            param.stride_ = 1;
            // What eta effective vector should be used
            eta_eff = !pl ? &eta_eff_back_ : &eta_eff_front_;
            // Determines sign of the spatial derivative
            negC = (pol_i_ == POLARIZATION::EY || pol_i_ == POLARIZATION::EX) ? false : true;
            // Ex, Hy, and Hz fields have one less point in the x direction
            if(pol_i_ == POLARIZATION::EX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HZ )
                trans1FieldOff = 1;
            // Ey, Hx, and Hz fields have one less point in the y direction
            if( (pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::EY ) && gridComm_->rank() == gridComm_->size() - 1)
                trans2FieldOff = 1;
            // Hx, Hy, and Ez fields all have one less point in the z direction
            if( pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EZ )
                pl ? ccStart = 1 : dirMax -= 1;
        }
        else
        {
            throw std::logic_error("PMLs need to be in X,Y, or Z direction it can't be NONE");
        }
        param.transSz_ = grid_i_->ln_vec()[cor_trans1] - 2 - trans1FieldOff;
        transSz2       = grid_i_->ln_vec()[cor_trans2] - 2 - trans2FieldOff;
        if(grid_i_->local_z() == 1 && cor_trans2 == 2)
            transSz2 = 1;

        // Magnetic fields are all 0.5 units off the main gird points in a direction not along its polarization (which are not included.)
        int dist = 0;
        double distOff = 0.0;
        if(pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HZ)
            distOff = pl ? -0.5 : 0.5;

        // Get the constants for all directions by looping over the PML thickness
        for(int cc = ccStart; cc < dirMax; ++cc)
        {
            // calculate distance from edge (facing into the cell)
            dist = pmlEdge > nDir ? (grid_i_->n_vec()[cor_norm]-2*gridComm_->npArr()[cor_norm]) - (startPt - cc) : cc+startPt;

            double sig = sigma( (static_cast<double>(dist) + distOff), static_cast<double>(nDir-1), eta_eff->at(dist));
            double kap = kappa( (static_cast<double>(dist) + distOff), static_cast<double>(nDir-1) );
            double a   = aVal ( (static_cast<double>(dist) + distOff), static_cast<double>(nDir-1) );

            param.b_    = b( sig, a, kap );
            param.c_    = c( sig, a, kap );
            param.cOff_ = c( sig, a, kap );
            negC ? param.c_ *= -1.0 : param.cOff_ *= -1.0;

            // initialize the loc and locOff to -1, -1, -1
            param.loc_    = {-1, -1, -1};
            param.locOff_ = {-1, -1, -1};

            // blas operations will always start at 1 (transSz)
            param.loc_[cor_trans1] = 1;
            param.locOff_[cor_trans1] = 1;
            // along the transSz2 direction  get the offset location from the update and then put the param set at the end of the vector
            for(int jj = 0; jj < transSz2 ; ++jj)
            {
                // param.loc_[cor_norm] is based on the distance and the locOff_[cor_norm] would be + 1 if an E field polarization, and - 1 if an H field polarization
                param.loc_[cor_norm] = pl ? grid_i_->ln_vec()[cor_norm] - 2 - cc : cc+1;
                if( pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::HZ )
                    param.locOff_[cor_norm] = pl ? grid_i_->ln_vec()[cor_norm]-1-cc : cc + 2;
                else
                    param.locOff_[cor_norm] = pl ? grid_i_->ln_vec()[cor_norm]-3-cc : cc;

                // looping parameters start at 1+jj
                param.loc_[cor_trans2] = 1+jj;
                param.locOff_[cor_trans2] = 1+jj;

                // if 2D z coordinate will be 0 not one so subtract 1
                if(grid_i_->local_z() == 1)
                {
                    param.loc_[2] -= 1;
                    param.locOff_[2] -= 1;
                }
                paramList.push_back(param);
            }
        }
        return paramList;
    }

    /**
     * @brief generate the lists that add psi to grid
     *
     * @param[in] min an array describing minimum values in each direction
     * @param[in] max an array describing maximum values in each direction
     * @param[in] stride stride for blas operations
     * @param[in] physGrid physical grid for direction
     * @return ax lists to generate the lists that add psi to grid
     */
    std::vector<std::array<int,5>> getAxLists(std::array<int,3> min, std::array<int,3> max, int stride)
    {
        // Loop over the PMLS and add the PML updates are independent of the object maps since they act on the D field if not in vacuum, the code is commented out in case this needs to change.
        std::vector<std::array<int,5>> axLists;
        if(stride == 1)
        {
            for(int kk = min[2]; kk < max[2]; ++kk)
            {
                for(int jj = min[1]; jj < max[1]; ++jj)
                {
                    axLists.push_back( {{ min[0], jj, kk, max[0]-min[0]+1, 0 }} );
                    // int ii = min[0];
                    // while(ii < max[0])
                    // {
                    //     int iistore = ii;
                    //     while (ii < max[0]-1 && physGrid->point(ii,jj,kk) == physGrid->point(ii+1,jj,kk) )
                    //         ++ii;
                    //     axLists.push_back( {{ iistore, jj, kk, ii-iistore+1, static_cast<int>(physGrid->point(iistore,jj,kk)) }} );
                    //     ++ii;
                    // }
                }
            }
        }
        else
        {
            // If 2D do the blas operations on the Y direction, if 3D do blas operations in the Z direction
            if(grid_i_->local_z()==1)
            {
                for(int kk = min[2]; kk < max[2]; ++kk)
                {
                    for(int ii = min[0]; ii < max[0]; ++ii)
                    {
                        axLists.push_back( {{ ii, min[1], kk, max[1]-min[1]+1, 0 }} );
                        // int jj = min[1];
                        // while(jj < max[1])
                        // {
                        //     int jjstore = jj;
                        //     while (jj < max[1]-1 && physGrid->point(ii,jj,kk) == physGrid->point(ii,jj+1,kk) )
                        //         ++jj;
                        //     axLists.push_back( {{ ii, jjstore, kk, jj-jjstore+1, static_cast<int>(physGrid->point(ii,jjstore,kk)) }} );
                        //     ++jj;
                        // }
                    }
                }
            }
            else
            {
                for(int jj = min[1]; jj < max[1]; ++jj)
                {
                    for(int ii = min[0]; ii < max[0]; ++ii)
                    {
                        axLists.push_back( {{ ii, jj, min[2], max[2]-min[2]+1, 0 }} );
                        // int kk = min[2];
                        // while(kk < max[2])
                        // {
                        //     int kkstore = kk;
                        //     while (kk < max[2]-1 && physGrid->point(ii,jj,kk) == physGrid->point(ii,jj,kk+1) )
                        //         ++kk;
                        //     axLists.push_back( {{ ii, jj, kkstore, kk-kkstore+1, static_cast<int>(physGrid->point(ii,jj,kkstore)) }} );
                        //     ++kk;
                        // }
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
    std::vector<updateGridParams> getGridUpList(DIRECTION dir, bool pl)
    {
        std::vector<updateGridParams> paramList;
        updateGridParams param;

        std::array<int,3> min = {0,0,0};
        std::array<int,3> max = {0,0,0};
        int cor_ii = 0, cor_jj = 0, cor_kk = 0, off_ii = 0;
        // Determine the x, y, and z min/max for each direction. For transverse directions it is the field size limitations and the normal direction is based on the PML thickness
        if(dir == DIRECTION::X)
        {
            param.stride_ = grid_i_->local_x();
            cor_ii = 0; cor_jj = 1; cor_kk = 2;

            // Prefactor sign can be determined by looking at Taflove Chapter 7
            param.Db_    = dt_ / d_[0];
            if(pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::EY)
                param.Db_*=-1.0;

            // offset in the ii direction is based on taking the mirror image of the point in the opposite PML
            if(pl && ( pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::EX || pol_i_ == POLARIZATION::HY ) )
                off_ii = -1;
            else if( !pl && ( pol_i_ == POLARIZATION::EZ || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EY ) )
                off_ii = 1;
        }
        else if(dir == DIRECTION::Y)
        {
            param.stride_ = 1;
            cor_ii = 1; cor_jj = 2; cor_kk = 0;

            // Prefactor sign can be determined by looking at Taflove Chapter 7
            param.Db_    = dt_ / d_[1];
            if(pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EZ)
                param.Db_*=-1.0;

            // offset in the ii direction is based on taking the mirror image of the point in the opposite PML
            if(pl && ( pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EY ) )
                off_ii = -1;
            else if( !pl && ( pol_i_ == POLARIZATION::EZ || pol_i_ == POLARIZATION::EX || pol_i_ == POLARIZATION::HY ) )
                off_ii = 1;
        }
        else if(dir == DIRECTION::Z)
        {
            param.stride_ = 1;
            cor_ii = 2; cor_jj = 0; cor_kk = 1;

            // Prefactor sign can be determined by looking at Taflove Chapter 7
            param.Db_    = dt_ / d_[2];
            if(pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EX)
                param.Db_*=-1.0;

            // offset in the ii direction is based on taking the mirror image of the point in the opposite PML
            if(pl && ( pol_i_ == POLARIZATION::EZ || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY ) )
                off_ii = -1;
            else if( !pl && ( pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::EX || pol_i_ == POLARIZATION::EY ) )
                off_ii = 1;
        }
        else
            throw std::logic_error("Direction for grid updates has to be X, Y or Z not NONE");

        min[cor_ii] += pl ? grid_i_->ln_vec()[cor_ii] - ln_vec_pl_[cor_ii] - 1  : 1 ;
        max[cor_ii] += pl ? grid_i_->ln_vec()[cor_ii] -1  + off_ii: ln_vec_mn_[cor_ii] + off_ii;

        min[cor_jj] = 1;
        max[cor_jj] = grid_i_->ln_vec()[cor_jj]-1;

        min[cor_kk] = 1;
        max[cor_kk] = grid_i_->ln_vec()[cor_kk]-1;

        // Ex, Hy, Hz fields have one less point in the x direction
        if( cor_ii != 0 && ( pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::EX || pol_i_ == POLARIZATION::HY ) )
            max[0] -= 1;
        // Hx, Ey, Hz fields have one less point in the y direction
        if( cor_ii != 1 && ( gridComm_->rank() == gridComm_->size()-1 ) && (pol_i_ == POLARIZATION::HZ || pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::EY) )
            max[1] -= 1;
        // Hx, Hy, Ez fields have one less point in the z direction
        if( cor_ii != 2 && grid_i_->local_z() != 1 && ( pol_i_ == POLARIZATION::HX || pol_i_ == POLARIZATION::HY || pol_i_ == POLARIZATION::EZ ) )
            max[2] -= 1;


        std::vector<std::array<int,5>> axParams = getAxLists(min, max, param.stride_);
        for(auto & axList : axParams)
        {
            param.loc_ = {axList[0], axList[1], axList[2]};
            param.nAx_ = axList[3];
            paramList.push_back(param);
        }

        return paramList;
    }

    /**
     * @brief      Fills the psi and grid update lists with the correct values
     *
     * @param[in]  ln_pml   local size of the PML thickness
     * @param[in]  dir      The direction of polarization of the $\psi$ field
     * @param[in]  pl       True if top, right or front.
     * @param[in]  startPt  Where the PMLs starts.
     * @param[in]  pmlEdge  Where the PML ends
     * @param[in]  nDir     Thickness in that direction.
     * @param[in]  dirMax   How far to iterate over.
     * @param      psiList   The $\psi$ update list
     * @param      gridList  The grid_i update list
     */
    void fillLists(DIRECTION dir, bool pl, int startPt, int pmlEdge, int nDir, int dirMax, std::vector<updatePsiParams>& psiList, std::vector<updateGridParams>& gridList)
    {
        std::vector<updatePsiParams> psiTemp;
        std::vector<updateGridParams> gridTemp;

        psiTemp  = getPsiUpList (dir, pl, startPt, pmlEdge, nDir, dirMax);
        gridTemp = getGridUpList(dir, pl);

        psiList.reserve(psiList.size() + psiTemp.size());
        psiList.insert(psiList.end(), psiTemp.begin(), psiTemp.end());
        gridList.reserve(gridList.size() + gridTemp.size());
        gridList.insert(gridList.end(), gridTemp.begin(), gridTemp.end());
    }

    /**
     * @brief      Initializes the update list with the correct parameters
     *
     * @param[in]  ln_pml   local size of the PML thickness
     * @param[in]  dir      The direction of polarization of the $\psi$ field
     * @param[in]  pl       True if top, right or front.
     * @param[in]  startPt  Where the PMLs starts.
     * @param[in]  pmlEdge  Where the PML ends
     * @param[in]  nDir     Thickness in that direction.
     * @param[in]  dirMax   How far to iterate over.
     */
    void initalizeLists(int ln_pml, DIRECTION dir, bool pl, int startPt, int pmlEdge, int nDir, int dirMax)
    {
        if(ln_pml > 0 && i_ != dir)
        {
            if(grid_k_ && j_ == dir)
                fillLists( dir, pl, startPt, pmlEdge, nDir, dirMax, updateListPsi_j_, updateListGrid_j_);
            else if(grid_j_ && k_ == dir)
                fillLists( dir, pl, startPt, pmlEdge, nDir, dirMax, updateListPsi_k_, updateListGrid_k_);
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
     * @brief      Constructor class
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
    parallelCPMLReal(std::shared_ptr<mpiInterface> gridComm, std::vector<real_grid_ptr> weights, real_pgrid_ptr grid_i, real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, POLARIZATION pol_i, std::array<int,3> n_vec, double m, double ma, double aMax, std::array<double,3> d, double dt, int_pgrid_ptr physGrid, std::vector<std::shared_ptr<Obj>> objArr);
};

class parallelCPMLCplx : public parallelCPML<cplx>
{
public:
    /**
     * @brief      Constructor class
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
    parallelCPMLCplx(std::shared_ptr<mpiInterface> gridComm, std::vector<real_grid_ptr> weights, std::shared_ptr<parallelGrid<cplx > > grid_i, std::shared_ptr<parallelGrid<cplx > > grid_j, std::shared_ptr<parallelGrid<cplx > > grid_k, POLARIZATION pol_i, std::array<int,3> n_vec, double m, double ma, double aMax, std::array<double,3> d, double dt, int_pgrid_ptr physGrid, std::vector<std::shared_ptr<Obj>> objArr);
};
#endif