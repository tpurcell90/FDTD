#include <DTC/parallelFlux.hpp>

parallelFluxDTCReal::parallelFluxDTCReal(std::shared_ptr<mpiInterface> gridComm, std::string name, double weight, pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, std::array<int,3> loc, std::array<int,3> sz, bool cross_sec, bool save, bool load, int timeInt, std::vector<double> freqList, DIRECTION propDir, std::array<double,3> d, double dt, double theta, double phi, double psi, double alpha, std::string incd_file, bool SI, double I0, double a) :
    parallelFluxDTC<double>(gridComm, name, weight, Ex,Ey,Ez, Hx, Hy, Hz, loc, sz, cross_sec, save, load, timeInt, freqList, propDir, d, dt, theta, phi, psi, alpha, incd_file, SI, I0, a)
{
    getIncdField_ = [](cplx a){return cplx(std::real(a), 0.0); };
    if( !(Ez && Hz ) )
    {
        if(sz_[0] > 1 && sz_[1] > 1)
        {
            // Surface at each end of the box
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Y, false) );
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Y, true ) );

            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::X, false) );
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::X, true ) );
        }
        else if(sz_[0] == 1 && ( (Hz && Hz->local_x() != 3) || (Ez && Ez->local_x() !=3) ) )
        {
            // Single surface assume positive weight
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::X, true) );
        }
        else
        {
            // Single surface assume positive weight
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Y, true) );
        }
    }
    else
    {
        if(sz_[0] > 1 && sz_[1] > 1 && sz_[2] > 1)
        {
            // A plane for each surface of the box
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::X, true) );
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::X, false) );

            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Y, true) );
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Y, false) );

            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Z, true) );
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Z, false) );
        }
        else if(sz_[0] == 1)
        {
            // Single surface assume positive weight
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::X, true) );
        }
        else if(sz_[1] == 1)
        {
            // Single surface assume positive weight
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Y, true) );
        }
        else if(sz_[2] == 1)
        {
            // Single surface assume positive weight
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Z, true) );
        }
    }
}

parallelFluxDTCCplx::parallelFluxDTCCplx(std::shared_ptr<mpiInterface> gridComm, std::string name, double weight, pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, std::array<int,3> loc, std::array<int,3> sz, bool cross_sec, bool save, bool load, int timeInt, std::vector<double> freqList, DIRECTION propDir, std::array<double,3> d, double dt, double theta, double phi, double psi, double alpha, std::string incd_file, bool SI, double I0, double a) :
    parallelFluxDTC<cplx>(gridComm, name, weight, Ex,Ey,Ez, Hx, Hy, Hz, loc, sz, cross_sec, save, load, timeInt, freqList, propDir, d, dt, theta, phi, psi, alpha, incd_file, SI, I0, a)
{
    getIncdField_ = [](cplx a){return a; };
    if(sz_.size() == 2)
    {
        if(sz_[0] > 1 && sz_[1] > 1)
        {
            // Surface at each end of the box
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::X, true ) );
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::X, false) );

            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Y, true ) );
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Y, false) );
        }
        else if(sz_[0] == 1 && ( (Hz && Hz->local_x() != 3) || (Ez && Ez->local_x() !=3) ) )
        {
            // Single surface assume positive weight
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::X, true) );
        }
        else
        {
            // Single surface assume positive weight
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Y, true) );
        }
    }
    else
    {
        if(sz_[0] > 1 && sz_[1] > 1 && sz_[2] > 1)
        {
            // A plane for each surface of the box
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::X, true) );
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::X, false) );

            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Y, true) );
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Y, false) );

            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Z, true) );
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Z, false) );
        }
        else if(sz_[0] == 1)
        {
            // Single surface assume positive weight
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::X, true) );
        }
        else if(sz_[1] == 1)
        {
            // Single surface assume positive weight
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Y, true) );
        }
        // Single surface assume positive weight
        else if(sz_[2] == 1)
        {
            fInParam_.push_back(makeParamIn(Ex,Ey,Ez, Hx, Hy, Hz, DIRECTION::Z, true) );
        }
    }
}

parallelFluxDTCCplx::FieldInputParamsFlux parallelFluxDTCCplx::makeParamIn(pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, DIRECTION dir, bool pl)
{
    FieldInputParamsFlux to_return;
    masterImportDat toMaster;

    int cor = -1; //!< Coordinate of the direction of copying
    int transCor1 = -1; //!< Coordinate of the direction of the main loop for transferring data
    int transCor2 = -1; //!< Coordinate of the direction of the second loop for transferring data
    int corI = -1; //!< coordinate of the I direction
    int corJ = -1; //!< coordinate of the J direction
    int corK = -1; //!< coordinate of the K direction
    std::array<int,3> sz = {{ 0, 0, 0 }};
    pgrid_ptr Ej; //!< grid pointer to the Ej field
    pgrid_ptr Ek; //!< grid pointer to the Ek field
    pgrid_ptr Hj; //!< grid pointer to the Hj field
    pgrid_ptr Hk; //!< grid pointer to the Hk field

    // Set up the values of coordinates that are needed, and the size of the surfaces

    if(dir == DIRECTION::X)
    {
        Ej = Ey; Hj = Hy; Ek = Ez; Hk = Hz;
        // Coordinates corresponding to how it is in terms of adding from main girds to detector grids
        if( (Ej && Ej->local_z() == 1) || (Ek && Ek->local_z() == 1) )
        {
            cor = 1; transCor1 = 2; transCor2 = 0;
        }
        else
        {
            cor = 2; transCor1 = 1; transCor2 = 0;
        }

        // ijk coordinates
        corI = 0; corJ = 1; corK = 2;
        sz = { { 1, sz_[1], sz_[2] } };
    }
    else if(dir == DIRECTION::Y)
    {
        Ej = Ez; Hj = Hz; Ek = Ex; Hk = Hx;
        // Coordinates corresponding to how it is in terms of adding from main girds to detector grids
        cor = 0; transCor1 = 2; transCor2 = 1;
        // ijk coordinates
        corI = 1; corJ = 2; corK = 0;
        sz = { { sz_[0], 1, sz_[2] } };
    }
    else
    {
        Ej = Ex; Hj = Hx; Ek = Ey; Hk = Hy;
        // Coordinates corresponding to how it is in terms of adding from main girds to detector grids
        cor = 0; transCor1 = 1; transCor2 = 2;
        // ijk coordinates
        corI = 2; corJ = 0; corK = 1;
        sz = { { sz_[0], sz_[1], 1 } };
    }

    std::array<int,3> loc(loc_);
    if(pl)
    {
        loc[transCor2] = loc_[transCor2] + sz_[transCor2] - 1;
    }

    std::array<int,3> locOff(loc);
        // total size of the surface is the area of the sizes in both direction
    to_return.sz_ = sz[cor]*sz[transCor1];
    // addIndex set to see where in the surface components to place this processes data
    if(Ez)
        toMaster.addIndex_ = {{ Ez->procLoc()[transCor1] - loc_[transCor1], Ez->procLoc()[cor] - loc_[cor] }};
    else
        toMaster.addIndex_ = {{ Hz->procLoc()[transCor1] - loc_[transCor1], Hz->procLoc()[cor] - loc_[cor] }};

    // If addIndex < 0 then the object starts in this process
    if(toMaster.addIndex_[0] < 0)
        toMaster.addIndex_[0] = 0;
    if(toMaster.addIndex_[1] < 0)
        toMaster.addIndex_[1] = 0;

    if(Ez && Hz)
    {
        std::array<int,3> szOff(sz);
        // Make the fields all centered around the center point of the face that we are looking at
        locOff[corK] -= 1;
         szOff[corK] += 1;
        to_return.Ek_dtc_.push_back( std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Ek, dir, locOff, szOff, freqList_) );

        to_return.Hk_dtc_.push_back( std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Hk, dir, locOff, szOff, freqList_) );
        locOff[corI] -= 1;
        to_return.Hk_dtc_.push_back( std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Hk, dir, locOff, szOff, freqList_) );

        locOff[corK] += 1;
         szOff[corK] -= 1;
        locOff[corJ] -= 1;
         szOff[corJ] += 1;

        to_return.Hj_dtc_.push_back( std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Hj, dir, locOff, szOff, freqList_) );
        locOff[corI] += 1;
        to_return.Hj_dtc_.push_back( std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Hj, dir, locOff, szOff, freqList_) );

        to_return.Ej_dtc_.push_back( std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Ej, dir, locOff, szOff, freqList_) );
    }
    else if(Ez)
    {
        if(Ek)
        {
            locOff[corI] += -1;
            // Ez
            to_return.Ek_dtc_.push_back( std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Ek, dir, loc   , sz, freqList_) );
            // Hy center around Ez points
            to_return.Hj_dtc_.push_back( std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Hj, dir, loc   , sz, freqList_) );
            to_return.Hj_dtc_.push_back( std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Hj, dir, locOff, sz, freqList_) );
        }
        else
        {
            locOff[corI] += -1;
            // Ez
            to_return.Ej_dtc_.push_back( std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Ej, dir, loc   , sz, freqList_) );
            // Hy center around Ez points
            to_return.Hk_dtc_.push_back( std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Hk, dir, loc   , sz, freqList_) );
            to_return.Hk_dtc_.push_back( std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Hk, dir, locOff, sz, freqList_) );
        }
    }
    else
    {
        if(Ej)
        {
            locOff[corI] += 1;
            // Ey
            to_return.Ej_dtc_.push_back(std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Ej, dir, loc   , sz, freqList_) );
            // Hz
            to_return.Hk_dtc_.push_back(std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Hk, dir, loc   , sz, freqList_) );
            // Center around Hz points
            to_return.Ej_dtc_.push_back(std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Ej, dir, locOff, sz, freqList_) );
        }
        else
        {
            locOff[corI] += 1;
            // Ey
            to_return.Ek_dtc_.push_back(std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Ek, dir, loc   , sz, freqList_) );
            // Hz
            to_return.Hj_dtc_.push_back(std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Hj, dir, loc   , sz, freqList_) );
            // Center around Hz points
            to_return.Ek_dtc_.push_back(std::make_shared<parallelStorageFreqDTCCplx>(outProc_, Ek, dir, locOff, sz, freqList_) );
        }
    }
    // For consistency with boxes everything point outward is positive
    pl  ? to_return.weight_ = 1.0 : to_return.weight_ = -1.0;

    if( gridComm_->rank() == outProc_ )
    {
        // Make all Ej, Ek, Hj, and Hk freq grids for collection and output on the process with rank == outProc
        if(to_return.Ej_dtc_.size() > 0)
            Ej_freq_.push_back( std::make_shared<Grid<cplx>>( std::array<int,3>( {{ nfreq_, sz[transCor1], sz[cor] }} ), std::array<double,3>( {{dLam_, d_[transCor1], d_[cor] }} ) ) ) ;
        else
            Ej_freq_.push_back(nullptr);
        if(to_return.Ek_dtc_.size() > 0)
            Ek_freq_.push_back( std::make_shared<Grid<cplx>>( std::array<int,3>( {{ nfreq_, sz[transCor1], sz[cor] }} ), std::array<double,3>( {{dLam_, d_[transCor1], d_[cor] }} ) ) ) ;
        else
            Ek_freq_.push_back(nullptr);
        if(to_return.Hj_dtc_.size() > 0)
            Hj_freq_.push_back( std::make_shared<Grid<cplx>>( std::array<int,3>( {{ nfreq_, sz[transCor1], sz[cor] }} ), std::array<double,3>( {{dLam_, d_[transCor1], d_[cor] }} ) ) ) ;
        else
            Hj_freq_.push_back(nullptr);
        if(to_return.Hk_dtc_.size() > 0)
            Hk_freq_.push_back( std::make_shared<Grid<cplx>>( std::array<int,3>( {{ nfreq_, sz[transCor1], sz[cor] }} ), std::array<double,3>( {{dLam_, d_[transCor1], d_[cor] }} ) ) ) ;
        else
            Hk_freq_.push_back(nullptr);
    }
    else
    {
        Ej_freq_.push_back(nullptr);
        Ek_freq_.push_back(nullptr);
        Hj_freq_.push_back(nullptr);
        Hk_freq_.push_back(nullptr);
    }

    // Get the final importing field parameters constructed
    toMaster.szProcOffsetEj_ = constructSzProcOffsetLists( to_return.Ej_dtc_, toMaster.addIndex_, Ej, Ek, corJ, cor, transCor1 );
    toMaster.szProcOffsetEk_ = constructSzProcOffsetLists( to_return.Ek_dtc_, toMaster.addIndex_, Ek, Ej, corK, cor, transCor1 );
    toMaster.szProcOffsetHj_ = constructSzProcOffsetLists( to_return.Hj_dtc_, toMaster.addIndex_, Hj, Hk, corJ, cor, transCor1 );
    toMaster.szProcOffsetHk_ = constructSzProcOffsetLists( to_return.Hk_dtc_, toMaster.addIndex_, Hk, Hj, corK, cor, transCor1 );

    // have the master process collect all necessary data for importing fields at the end of the calculation
    if(gridComm_->rank() == outProc_)
    {
        // Gather all the toMaster on the collection/output process and add them to the appropriate vectors
        std::vector<std::shared_ptr<masterImportDat>> masterProc;
        std::vector<masterImportDat> allProcs;
        mpi::gather(*gridComm_, toMaster, allProcs, outProc_);
        for(auto & proc : allProcs)
        {
            if(proc.szProcOffsetEj_.size() > 0)
                to_return.combineEjFields_.push_back(std::make_shared<masterImportDat>(proc) );
            if(proc.szProcOffsetEk_.size() > 0)
                to_return.combineEkFields_.push_back(std::make_shared<masterImportDat>(proc) );
            if(proc.szProcOffsetHj_.size() > 0)
                to_return.combineHjFields_.push_back(std::make_shared<masterImportDat>(proc) );
            if(proc.szProcOffsetHk_.size() > 0)
                to_return.combineHkFields_.push_back(std::make_shared<masterImportDat>(proc) );
        }
        if(to_return.combineHkFields_.empty() && to_return.combineHjFields_.empty() && to_return.combineEkFields_.empty() && to_return.combineEjFields_.empty())
            throw std::logic_error("One of the flux surfaces is outside the FDTD cell.");
    }
    else
    {
        // If not output/collection then send it and make the combine field vectors empty
        mpi::gather(*gridComm_, toMaster,  outProc_);
        to_return.combineEjFields_ = {};
        to_return.combineEkFields_ = {};
        to_return.combineHjFields_ = {};
        to_return.combineHkFields_ = {};
    }
    to_return.addIndex_ = toMaster.addIndex_;
    return to_return;
}

parallelFluxDTCReal::FieldInputParamsFlux parallelFluxDTCReal::makeParamIn(pgrid_ptr Ex, pgrid_ptr Ey, pgrid_ptr Ez, pgrid_ptr Hx, pgrid_ptr Hy, pgrid_ptr Hz, DIRECTION dir, bool pl)
{
    FieldInputParamsFlux to_return;
    masterImportDat toMaster;

    int cor = -1; //!< Coordinate of the direction of copying
    int transCor1 = -1; //!< Coordinate of the direction of the main loop for transferring data
    int transCor2 = -1; //!< Coordinate of the direction of the second loop for transferring data
    int corI = -1; //!< coordinate of the I direction
    int corJ = -1; //!< coordinate of the J direction
    int corK = -1; //!< coordinate of the K direction
    std::array<int,3> sz = {{ 0, 0, 0 }};
    pgrid_ptr Ej; //!< grid pointer to the Ej field
    pgrid_ptr Ek; //!< grid pointer to the Ek field
    pgrid_ptr Hj; //!< grid pointer to the Hj field
    pgrid_ptr Hk; //!< grid pointer to the Hk field

    // Set up the values of coordinates that are needed, and the size of the surfaces

    if(dir == DIRECTION::X)
    {
        Ej = Ey; Hj = Hy; Ek = Ez; Hk = Hz;
        // Coordinates corresponding to how it is in terms of adding from main girds to detector grids
        if( (Ej && Ej->local_z() == 1) || (Ek && Ek->local_z() == 1) )
        {
            cor = 1; transCor1 = 2; transCor2 = 0;
        }
        else
        {
            cor = 2; transCor1 = 1; transCor2 = 0;
        }

        // ijk coordinates
        corI = 0; corJ = 1; corK = 2;
        sz = { { 1, sz_[1], sz_[2] } };
    }
    else if(dir == DIRECTION::Y)
    {
        Ej = Ez; Hj = Hz; Ek = Ex; Hk = Hx;
        // Coordinates corresponding to how it is in terms of adding from main girds to detector grids
        cor = 0; transCor1 = 2; transCor2 = 1;
        // ijk coordinates
        corI = 1; corJ = 2; corK = 0;
        sz = { { sz_[0], 1, sz_[2] } };
    }
    else
    {
        Ej = Ex; Hj = Hx; Ek = Ey; Hk = Hy;
        // Coordinates corresponding to how it is in terms of adding from main girds to detector grids
        cor = 0; transCor1 = 1; transCor2 = 2;
        // ijk coordinates
        corI = 2; corJ = 0; corK = 1;
        sz = { { sz_[0], sz_[1], 1 } };
    }

    std::array<int,3> loc(loc_);
    if(pl)
    {
        loc[transCor2] = loc_[transCor2] + sz_[transCor2] - 1;
    }

    std::array<int,3> locOff(loc);
        // total size of the surface is the area of the sizes in both direction
    to_return.sz_ = sz[cor]*sz[transCor1];
    // addIndex set to see where in the surface components to place this processes data
    if(Ez)
        toMaster.addIndex_ = {{ Ez->procLoc()[transCor1] - loc_[transCor1], Ez->procLoc()[cor] - loc_[cor] }};
    else
        toMaster.addIndex_ = {{ Hz->procLoc()[transCor1] - loc_[transCor1], Hz->procLoc()[cor] - loc_[cor] }};

    // If addIndex < 0 object starts in this process, so set it to 0
    if(toMaster.addIndex_[0] < 0)
        toMaster.addIndex_[0] = 0;

    if(toMaster.addIndex_[1] < 0)
        toMaster.addIndex_[1] = 0;

    if(Ez && Hz)
    {
        std::array<int,3> szOff(sz);
        // Make the fields all centered around the center point of the face that we are looking at
        locOff[corK] -= 1;
        szOff[corK] += 1;

        to_return.Ek_dtc_.push_back( std::make_shared<parallelStorageFreqDTCReal>(outProc_, Ek, dir, locOff, szOff, freqList_) );

        to_return.Hk_dtc_.push_back( std::make_shared<parallelStorageFreqDTCReal>(outProc_, Hk, dir, locOff, szOff, freqList_) );
        locOff[corI] -= 1;
        to_return.Hk_dtc_.push_back( std::make_shared<parallelStorageFreqDTCReal>(outProc_, Hk, dir, locOff, szOff, freqList_) );

        locOff[corK] += 1;
        szOff[corK] -= 1;

        locOff[corJ] -= 1;
        szOff[corJ] += 1;

        to_return.Hj_dtc_.push_back( std::make_shared<parallelStorageFreqDTCReal>(outProc_, Hj, dir, locOff, szOff, freqList_) );
        locOff[corI] += 1;
        to_return.Hj_dtc_.push_back( std::make_shared<parallelStorageFreqDTCReal>(outProc_, Hj, dir, locOff, szOff, freqList_) );

        to_return.Ej_dtc_.push_back( std::make_shared<parallelStorageFreqDTCReal>(outProc_, Ej, dir, locOff, szOff, freqList_) );
    }
    else if(Ez)
    {
        if(Ek)
        {
            locOff[corI] += -1;
            // Ez
            to_return.Ek_dtc_.push_back( std::make_shared<parallelStorageFreqDTCReal>(outProc_, Ek, dir, loc   , sz, freqList_) );
            // Hy center around Ez points
            to_return.Hj_dtc_.push_back( std::make_shared<parallelStorageFreqDTCReal>(outProc_, Hj, dir, loc   , sz, freqList_) );
            to_return.Hj_dtc_.push_back( std::make_shared<parallelStorageFreqDTCReal>(outProc_, Hj, dir, locOff, sz, freqList_) );
        }
        else
        {
            locOff[corI] += -1;
            // Ez
            to_return.Ej_dtc_.push_back( std::make_shared<parallelStorageFreqDTCReal>(outProc_, Ej, dir, loc   , sz, freqList_) );
            // Hy center around Ez points
            to_return.Hk_dtc_.push_back( std::make_shared<parallelStorageFreqDTCReal>(outProc_, Hk, dir, loc   , sz, freqList_) );
            to_return.Hk_dtc_.push_back( std::make_shared<parallelStorageFreqDTCReal>(outProc_, Hk, dir, locOff, sz, freqList_) );
        }
    }
    else
    {
        if(Ej)
        {
            locOff[corI] += 1;
            // Ey
            to_return.Ej_dtc_.push_back(std::make_shared<parallelStorageFreqDTCReal>(outProc_, Ej, dir, loc   , sz, freqList_) );
            // Hz
            to_return.Hk_dtc_.push_back(std::make_shared<parallelStorageFreqDTCReal>(outProc_, Hk, dir, loc   , sz, freqList_) );
            // Center around Hz points
            to_return.Ej_dtc_.push_back(std::make_shared<parallelStorageFreqDTCReal>(outProc_, Ej, dir, locOff, sz, freqList_) );
        }
        else
        {
            locOff[corI] += 1;
            // Ey
            to_return.Ek_dtc_.push_back(std::make_shared<parallelStorageFreqDTCReal>(outProc_, Ek, dir, loc   , sz, freqList_) );
            // Hz
            to_return.Hj_dtc_.push_back(std::make_shared<parallelStorageFreqDTCReal>(outProc_, Hj, dir, loc   , sz, freqList_) );
            // Center around Hz points
            to_return.Ek_dtc_.push_back(std::make_shared<parallelStorageFreqDTCReal>(outProc_, Ek, dir, locOff, sz, freqList_) );
        }
    }
    // For consistency with boxes everything point outward is positive
    pl  ? to_return.weight_ = 1.0 : to_return.weight_ = -1.0;

    if( gridComm_->rank() == outProc_ )
    {
        // Make all Ej, Ek, Hj, and Hk freq grids for collection and output on the process with rank == outProc
        if(to_return.Ej_dtc_.size() > 0)
            Ej_freq_.push_back( std::make_shared<Grid<cplx>>( std::array<int,3>( {{ nfreq_, sz[transCor1], sz[cor] }} ), std::array<double,3>( {{dLam_, d_[transCor1], d_[cor] }} ) ) ) ;
        else
            Ej_freq_.push_back(nullptr);
        if(to_return.Ek_dtc_.size() > 0)
            Ek_freq_.push_back( std::make_shared<Grid<cplx>>( std::array<int,3>( {{ nfreq_, sz[transCor1], sz[cor] }} ), std::array<double,3>( {{dLam_, d_[transCor1], d_[cor] }} ) ) ) ;
        else
            Ek_freq_.push_back(nullptr);
        if(to_return.Hj_dtc_.size() > 0)
            Hj_freq_.push_back( std::make_shared<Grid<cplx>>( std::array<int,3>( {{ nfreq_, sz[transCor1], sz[cor] }} ), std::array<double,3>( {{dLam_, d_[transCor1], d_[cor] }} ) ) ) ;
        else
            Hj_freq_.push_back(nullptr);
        if(to_return.Hk_dtc_.size() > 0)
            Hk_freq_.push_back( std::make_shared<Grid<cplx>>( std::array<int,3>( {{ nfreq_, sz[transCor1], sz[cor] }} ), std::array<double,3>( {{dLam_, d_[transCor1], d_[cor] }} ) ) ) ;
        else
            Hk_freq_.push_back(nullptr);
    }
    else
    {
        Ej_freq_.push_back(nullptr);
        Ek_freq_.push_back(nullptr);
        Hj_freq_.push_back(nullptr);
        Hk_freq_.push_back(nullptr);
    }
    // Construct parameters needed to import fields at the end of the calculation
    toMaster.szProcOffsetEj_ = constructSzProcOffsetLists( to_return.Ej_dtc_, toMaster.addIndex_, Ej, Ek, corJ, cor, transCor1 );
    toMaster.szProcOffsetEk_ = constructSzProcOffsetLists( to_return.Ek_dtc_, toMaster.addIndex_, Ek, Ej, corK, cor, transCor1 );
    toMaster.szProcOffsetHj_ = constructSzProcOffsetLists( to_return.Hj_dtc_, toMaster.addIndex_, Hj, Hk, corJ, cor, transCor1 );
    toMaster.szProcOffsetHk_ = constructSzProcOffsetLists( to_return.Hk_dtc_, toMaster.addIndex_, Hk, Hj, corK, cor, transCor1 );

    // If master process collect all parameters
    if(gridComm_->rank() == outProc_)
    {
        // Gather all the toMaster on the collection/output process and add them to the appropriate vectors
        std::vector<std::shared_ptr<masterImportDat>> masterProc;
        std::vector<masterImportDat> allProcs;
        mpi::gather(*gridComm_, toMaster, allProcs, outProc_);
        for(auto & proc : allProcs)
        {
            if(proc.szProcOffsetEj_.size() > 0)
                to_return.combineEjFields_.push_back(std::make_shared<masterImportDat>(proc) );
            if(proc.szProcOffsetEk_.size() > 0)
                to_return.combineEkFields_.push_back(std::make_shared<masterImportDat>(proc) );
            if(proc.szProcOffsetHj_.size() > 0)
                to_return.combineHjFields_.push_back(std::make_shared<masterImportDat>(proc) );
            if(proc.szProcOffsetHk_.size() > 0)
                to_return.combineHkFields_.push_back(std::make_shared<masterImportDat>(proc) );
        }
        if(to_return.combineHkFields_.empty() && to_return.combineHjFields_.empty() && to_return.combineEkFields_.empty() && to_return.combineEjFields_.empty())
            throw std::logic_error("One of the flux surfaces is outside the FDTD cell.");
    }
    else
    {
        // If not output/collection then send it and make the combine field vectors empty
        mpi::gather(*gridComm_, toMaster,  outProc_);
        to_return.combineEjFields_ = {};
        to_return.combineEkFields_ = {};
        to_return.combineHjFields_ = {};
        to_return.combineHkFields_ = {};
    }
    to_return.addIndex_ = toMaster.addIndex_;
    return to_return;
}