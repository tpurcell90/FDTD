#include <UTIL/FDTD_up_eq.hpp>

void FDTDCompUpdateFxnReal::OneCompCurlJ (std::array<int,8>& axParams, std::array<double,2>& prefactors, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k)
{
    // Finite difference of the j derivative components in the curl
    daxpy_(axParams[0],      prefactors[1], &grid_j->point(axParams[1]            ,axParams[2]            ,axParams[3]            ), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    daxpy_(axParams[0], -1.0*prefactors[1], &grid_j->point(axParams[1]+axParams[6],axParams[2]+axParams[4],axParams[3]+axParams[5]), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    return;
}

void FDTDCompUpdateFxnReal::OneCompCurlK (std::array<int,8>& axParams, std::array<double,2>& prefactors, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k)
{
    // Finite differnce of the k derivatives components in the curl
    daxpy_(axParams[0], -1.0*prefactors[1], &grid_k->point(axParams[1]            ,axParams[2]            ,axParams[3]            ), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    daxpy_(axParams[0],      prefactors[1], &grid_k->point(axParams[1]+axParams[4],axParams[2]+axParams[5],axParams[3]+axParams[6]), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    return;
}

void FDTDCompUpdateFxnReal::TwoCompCurl (std::array<int,8>& axParams, std::array<double,2>& prefactors, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k)
{
    // Finite difference of the j derivative components in the curl
    daxpy_(axParams[0],      prefactors[1], &grid_j->point(axParams[1]            ,axParams[2]            ,axParams[3]            ), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    daxpy_(axParams[0], -1.0*prefactors[1], &grid_j->point(axParams[1]+axParams[6],axParams[2]+axParams[4],axParams[3]+axParams[5]), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    // Finite difference of the k derivative components in the curl
    daxpy_(axParams[0], -1.0*prefactors[1], &grid_k->point(axParams[1]            ,axParams[2]            ,axParams[3]            ), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    daxpy_(axParams[0],      prefactors[1], &grid_k->point(axParams[1]+axParams[4],axParams[2]+axParams[5],axParams[3]+axParams[6]), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    return;
}

void FDTDCompUpdateFxnCplx::OneCompCurlJ (std::array<int,8>& axParams, std::array<double,2>& prefactors, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k)
{
    // Finite difference of the j derivative components in the curl
    zaxpy_(axParams[0],      prefactors[1], &grid_j->point(axParams[1]            ,axParams[2]            ,axParams[3]            ), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    zaxpy_(axParams[0], -1.0*prefactors[1], &grid_j->point(axParams[1]+axParams[6],axParams[2]+axParams[4],axParams[3]+axParams[5]), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    return;
}

void FDTDCompUpdateFxnCplx::OneCompCurlK (std::array<int,8>& axParams, std::array<double,2>& prefactors, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k)
{
    // Finite differnce of the k derivatives components in the curl
    zaxpy_(axParams[0], -1.0*prefactors[1], &grid_k->point(axParams[1]            ,axParams[2]            ,axParams[3]            ), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    zaxpy_(axParams[0],      prefactors[1], &grid_k->point(axParams[1]+axParams[4],axParams[2]+axParams[5],axParams[3]+axParams[6]), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    return;
}

void FDTDCompUpdateFxnCplx::TwoCompCurl (std::array<int,8>& axParams, std::array<double,2>& prefactors, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k)
{
    // Finite difference of the j derivative components in the curl
    zaxpy_(axParams[0],      prefactors[1], &grid_j->point(axParams[1]            ,axParams[2]            ,axParams[3]            ), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    zaxpy_(axParams[0], -1.0*prefactors[1], &grid_j->point(axParams[1]+axParams[6],axParams[2]+axParams[4],axParams[3]+axParams[5]), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    // Finite difference of the k derivative components in the curl
    zaxpy_(axParams[0], -1.0*prefactors[1], &grid_k->point(axParams[1]            ,axParams[2]            ,axParams[3]            ), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    zaxpy_(axParams[0],      prefactors[1], &grid_k->point(axParams[1]+axParams[4],axParams[2]+axParams[5],axParams[3]+axParams[6]), 1, &grid_i->point(axParams[1],axParams[2],axParams[3]), 1);
    return;
}

void FDTDCompUpdateFxnReal::UpdateLorPol(std::array<int,8>& axParams, pgrid_ptr grid_i, std::vector<pgrid_ptr> & lorPi, std::vector<pgrid_ptr> & prevLorPi, double* jstore, std::shared_ptr<Obj> obj)
{
    for(int pp = 0; pp < obj->alpha().size(); ++pp)
    {
        // Store current values of the Lorentzian polarizations into jstore
        dcopy_(axParams[0], &lorPi[pp]->point(axParams[1], axParams[2], axParams[3]),  1, jstore, 1);

        // Update the polarizations as done in Taflove Ch. 7
        dscal_(axParams[0], obj->alpha()[pp],     &lorPi[pp]->point(axParams[1], axParams[2], axParams[3]),  1);
        daxpy_(axParams[0], obj->   xi()[pp], &prevLorPi[pp]->point(axParams[1], axParams[2], axParams[3]),  1, &lorPi[pp] ->point(axParams[1], axParams[2], axParams[3]), 1);
        daxpy_(axParams[0], obj->gamma()[pp],        &grid_i->point(axParams[1], axParams[2], axParams[3]),  1, &lorPi[pp] ->point(axParams[1], axParams[2], axParams[3]), 1);
        // reset prevLorP with the previously stored values in jstore
        dcopy_(axParams[0], jstore,1, &prevLorPi[pp]->point(axParams[1], axParams[2], axParams[3]), 1);
    }
    return;
}

void FDTDCompUpdateFxnCplx::UpdateLorPol(std::array<int,8>& axParams, pgrid_ptr grid_i, std::vector<pgrid_ptr> & lorPi, std::vector<pgrid_ptr> & prevLorPi, cplx* jstore, std::shared_ptr<Obj> obj)
{
    for(int pp = 0; pp < obj->alpha().size(); ++pp)
    {
        // Store current values of the Lorentzian polarizations into jstore
        zcopy_(axParams[0], &lorPi[pp]->point(axParams[1], axParams[2], axParams[3]),  1, jstore, 1);

        // Update the polarizations as done in Taflove Ch. 7
        zscal_(axParams[0], obj->alpha()[pp],     &lorPi[pp]->point(axParams[1], axParams[2], axParams[3]),  1);
        zaxpy_(axParams[0], obj->   xi()[pp], &prevLorPi[pp]->point(axParams[1], axParams[2], axParams[3]),  1, &lorPi[pp] ->point(axParams[1], axParams[2], axParams[3]), 1);
        zaxpy_(axParams[0], obj->gamma()[pp],        &grid_i->point(axParams[1], axParams[2], axParams[3]),  1, &lorPi[pp] ->point(axParams[1], axParams[2], axParams[3]), 1);
        // reset prevLorP with the previously stored values in jstore
        zcopy_(axParams[0], jstore,1, &prevLorPi[pp]->point(axParams[1], axParams[2], axParams[3]), 1);
    }
    return;
}

void FDTDCompUpdateFxnReal::UpdateLorMag(std::array<int,8>& axParams, pgrid_ptr grid_i, std::vector<pgrid_ptr> & lorMi, std::vector<pgrid_ptr> & prevLorMi, double* jstore, std::shared_ptr<Obj> obj)
{
    for(int mm = 0; mm < obj->magAlpha().size(); ++mm)
    {
        // Store current values of the Lorentzian polarizations into jstore
        dcopy_(axParams[0], &lorMi[mm]->point(axParams[1], axParams[2], axParams[3]),  1, jstore, 1);

        // Update the polarizations as done in Taflove Ch. 7
        dscal_(axParams[0], obj->magAlpha()[mm],     &lorMi[mm]->point(axParams[1], axParams[2], axParams[3]),  1);
        daxpy_(axParams[0], obj->   magXi()[mm], &prevLorMi[mm]->point(axParams[1], axParams[2], axParams[3]),  1, &lorMi[mm] ->point(axParams[1], axParams[2], axParams[3]), 1);
        daxpy_(axParams[0], obj->magGamma()[mm],        &grid_i->point(axParams[1], axParams[2], axParams[3]),  1, &lorMi[mm] ->point(axParams[1], axParams[2], axParams[3]), 1);
        // reset prevLorP with the previously stored values in jstore
        dcopy_(axParams[0], jstore,1, &prevLorMi[mm]->point(axParams[1], axParams[2], axParams[3]), 1);
    }
    return;
}

void FDTDCompUpdateFxnCplx::UpdateLorMag(std::array<int,8>& axParams, pgrid_ptr grid_i, std::vector<pgrid_ptr> & lorMi, std::vector<pgrid_ptr> & prevLorMi, cplx* jstore, std::shared_ptr<Obj> obj)
{
    for(int mm = 0; mm < obj->magAlpha().size(); ++mm)
    {
        // Store current values of the Lorentzian polarizations into jstore
        zcopy_(axParams[0], &lorMi[mm]->point(axParams[1], axParams[2], axParams[3]),  1, jstore, 1);

        // Update the polarizations as done in Taflove Ch. 7
        zscal_(axParams[0], obj->magAlpha()[mm],     &lorMi[mm]->point(axParams[1], axParams[2], axParams[3]),  1);
        zaxpy_(axParams[0], obj->   magXi()[mm], &prevLorMi[mm]->point(axParams[1], axParams[2], axParams[3]),  1, &lorMi[mm] ->point(axParams[1], axParams[2], axParams[3]), 1);
        zaxpy_(axParams[0], obj->magGamma()[mm],        &grid_i->point(axParams[1], axParams[2], axParams[3]),  1, &lorMi[mm] ->point(axParams[1], axParams[2], axParams[3]), 1);
        // reset prevLorP with the previously stored values in jstore
        zcopy_(axParams[0], jstore,1, &prevLorMi[mm]->point(axParams[1], axParams[2], axParams[3]), 1);
    }
    return;
}

void FDTDCompUpdateFxnReal::DtoE(std::array<int,8>& axParams, pgrid_ptr Di, pgrid_ptr Ei, std::vector<pgrid_ptr> & lorPi, std::shared_ptr<Obj> obj)
{
    double eps = obj->epsInfty();
    // Set the E field to the D field
    dcopy_(axParams[0], &Di->point(axParams[1],axParams[2], axParams[3]), 1,&Ei->point(axParams[1],axParams[2], axParams[3]), 1);
    // Scale E field by 1/eps (Done because E = 1/eps D)
    dscal_(axParams[0], 1.0/eps, &Ei->point(axParams[1],axParams[2], axParams[3]), 1);
    // Add all Polarizations
    for(int pp = 0; pp < obj->alpha().size(); ++pp)
        daxpy_(axParams[0], -1.0/eps, &lorPi[pp]->point(axParams[1],axParams[2], axParams[3]), 1,&Ei->point(axParams[1],axParams[2], axParams[3]), 1);
    return;

}

void FDTDCompUpdateFxnCplx::DtoE(std::array<int,8>& axParams, pgrid_ptr Di, pgrid_ptr Ei, std::vector<pgrid_ptr> & lorPi, std::shared_ptr<Obj> obj)
{
    cplx eps( obj->epsInfty(), 0.0 );
    // Set the E field to the D field
    zcopy_(axParams[0], &Di->point(axParams[1],axParams[2], axParams[3]), 1,&Ei->point(axParams[1],axParams[2], axParams[3]), 1);
    // Scale E field by 1/eps (Done because E = 1/eps D)
    zscal_(axParams[0], 1.0/eps, &Ei->point(axParams[1],axParams[2], axParams[3]), 1);
    // Add all Polarizations
    for(int pp = 0; pp < obj->alpha().size(); ++pp)
        zaxpy_(axParams[0], -1.0/eps, &lorPi[pp]->point(axParams[1],axParams[2], axParams[3]), 1,&Ei->point(axParams[1],axParams[2], axParams[3]), 1);
    return;
}

void FDTDCompUpdateFxnReal::BtoH(std::array<int,8>& axParams, pgrid_ptr Bi, pgrid_ptr Hi, std::vector<pgrid_ptr> & lorMi, std::shared_ptr<Obj> obj)
{
    double mu = obj->muInfty();
    // Set the E field to the D field
    dcopy_(axParams[0], &Bi->point(axParams[1],axParams[2], axParams[3]), 1, &Hi->point(axParams[1],axParams[2], axParams[3]), 1);
    // Scale E field by 1/eps (Done because E = 1/eps D)
    dscal_(axParams[0], 1.0/mu, &Hi->point(axParams[1],axParams[2], axParams[3]), 1);
    // Add all Polarizations
    for(int mm = 0; mm < obj->magAlpha().size(); ++mm)
        daxpy_(axParams[0], -1.0/mu, &lorMi[mm]->point(axParams[1],axParams[2], axParams[3]), 1,&Hi->point(axParams[1],axParams[2], axParams[3]), 1);
    return;
}

void FDTDCompUpdateFxnCplx::BtoH(std::array<int,8>& axParams, pgrid_ptr Bi, pgrid_ptr Hi, std::vector<pgrid_ptr> & lorMi, std::shared_ptr<Obj> obj)
{
    cplx mu( obj->muInfty(), 0.0 );
    // Set the E field to the D field
    zcopy_(axParams[0], &Bi->point(axParams[1],axParams[2], axParams[3]), 1, &Hi->point(axParams[1],axParams[2], axParams[3]), 1);
    // Scale E field by 1/eps (Done because E = 1/eps D)
    zscal_(axParams[0], 1.0/mu, &Hi->point(axParams[1],axParams[2], axParams[3]), 1);
    // Add all Polarizations
    for(int mm = 0; mm < obj->magAlpha().size(); ++mm)
        zaxpy_(axParams[0], -1.0/mu, &lorMi[mm]->point(axParams[1],axParams[2], axParams[3]), 1,&Hi->point(axParams[1],axParams[2], axParams[3]), 1);
    return;
}

void FDTDCompUpdateFxnReal::applyPBC(pgrid_ptr fUp, std::array<double,3> & k_point, int nx, int ny, int nz, int xmax, int ymax, int zmin, int zmax, double & dx, double & dy, double & dz)
{
    // if real fields then the PBC is just copying from one side to the other
    if(zmin != 0)
    {
        for(int jj = 1; jj < ny; ++jj)
        {
            dcopy_(nz-1, &fUp->point(xmax-1, jj, 1     ), fUp->local_x(), &fUp->point(0   , jj,      1), fUp->local_x() );
            dcopy_(nz-1, &fUp->point(1     , jj, 1     ), fUp->local_x(), &fUp->point(xmax, jj,      1), fUp->local_x() );

            dcopy_(nx-1, &fUp->point(1     , jj, zmax-1), 1             , &fUp->point(1   , jj, zmin-1), 1 );
            dcopy_(nx-1, &fUp->point(1     , jj, zmin  ), 1             , &fUp->point(1   , jj, zmax  ), 1 );
        }

        // // X edges
        dcopy_(nx-1, &fUp->point(1,      1, zmin  ), 1, &fUp->point(1, ymax  , zmax  ), 1);
        dcopy_(nx-1, &fUp->point(1, ymax-1, zmin  ), 1, &fUp->point(1, 0     , zmax  ), 1);
        dcopy_(nx-1, &fUp->point(1,      1, zmax-1), 1, &fUp->point(1, ymax  , zmin-1), 1);
        dcopy_(nx-1, &fUp->point(1, ymax-1, zmax-1), 1, &fUp->point(1, 0     , zmin-1), 1);

        // Y edges
        dcopy_(ny-1, &fUp->point(     1, 1, zmin  ), fUp->local_x()*fUp->local_z(), &fUp->point(xmax  , 1, zmax  ), fUp->local_x()*fUp->local_z());
        dcopy_(ny-1, &fUp->point(xmax-1, 1, zmin  ), fUp->local_x()*fUp->local_z(), &fUp->point(0     , 1, zmax  ), fUp->local_x()*fUp->local_z());
        dcopy_(ny-1, &fUp->point(     1, 1, zmax-1), fUp->local_x()*fUp->local_z(), &fUp->point(xmax  , 1, zmin-1), fUp->local_x()*fUp->local_z());
        dcopy_(ny-1, &fUp->point(xmax-1, 1, zmax-1), fUp->local_x()*fUp->local_z(), &fUp->point(0     , 1, zmin-1), fUp->local_x()*fUp->local_z());

        // Z edges
        dcopy_(nz-1, &fUp->point(1     , 1     , 1), fUp->local_x(), &fUp->point(xmax, ymax, 1), fUp->local_x());
        dcopy_(nz-1, &fUp->point(xmax-1, 1     , 1), fUp->local_x(), &fUp->point(0   , ymax, 1), fUp->local_x());
        dcopy_(nz-1, &fUp->point(1     , ymax-1, 1), fUp->local_x(), &fUp->point(xmax, 0   , 1), fUp->local_x());
        dcopy_(nz-1, &fUp->point(xmax-1, ymax-1, 1), fUp->local_x(), &fUp->point(0   , 0   , 1), fUp->local_x());

        // //Corners
        fUp->point(xmax, ymax, zmax  ) = fUp->point(1     , 1     , zmin  );
        fUp->point(0   , ymax, zmax  ) = fUp->point(xmax-1, 1     , zmin  );
        fUp->point(xmax, 0   , zmax  ) = fUp->point(1     , ymax-1, zmin  );
        fUp->point(0   , 0   , zmax  ) = fUp->point(xmax-1, ymax-1, zmin  );

        fUp->point(xmax, ymax, zmin-1) = fUp->point(1     , 1     , zmax-1);
        fUp->point(0   , ymax, zmin-1) = fUp->point(xmax-1, 1     , zmax-1);
        fUp->point(xmax, 0   , zmin-1) = fUp->point(1     , ymax-1, zmax-1);
        fUp->point(0   , 0   , zmin-1) = fUp->point(xmax-1, ymax-1, zmax-1);
    }
    else
    {
        dcopy_(fUp->local_y(), &fUp->point(xmax-1, 0), fUp->local_x(), &fUp->point(0   , 0), fUp->local_x() );
        dcopy_(fUp->local_y(), &fUp->point(1     , 0), fUp->local_x(), &fUp->point(xmax, 0), fUp->local_x() );
    }
}

void FDTDCompUpdateFxnCplx::applyPBC(pgrid_ptr fUp, std::array<double,3> & k_point, int nx, int ny, int nz, int xmax, int ymax, int zmin, int zmax, double & dx, double & dy, double & dz)
{
    if(zmin != 0)
    {
        for(int jj = 1; jj < ny; ++jj)
        {
            zcopy_(nz-1, &fUp->point(xmax-1, jj, 1     ), fUp->local_x(), &fUp->point(0   , jj,      1), fUp->local_x() );
            zcopy_(nz-1, &fUp->point(1     , jj, 1     ), fUp->local_x(), &fUp->point(xmax, jj,      1), fUp->local_x() );

            zcopy_(nx-1, &fUp->point(1     , jj, zmax-1), 1             , &fUp->point(1   , jj, zmin-1), 1 );
            zcopy_(nx-1, &fUp->point(1     , jj, zmin  ), 1             , &fUp->point(1   , jj, zmax  ), 1 );
        }

        // // X edges
        zcopy_(nx-1, &fUp->point(1,      1, zmin  ), 1, &fUp->point(1, ymax  , zmax  ), 1);
        zcopy_(nx-1, &fUp->point(1, ymax-1, zmin  ), 1, &fUp->point(1, 0     , zmax  ), 1);
        zcopy_(nx-1, &fUp->point(1,      1, zmax-1), 1, &fUp->point(1, ymax  , zmin-1), 1);
        zcopy_(nx-1, &fUp->point(1, ymax-1, zmax-1), 1, &fUp->point(1, 0     , zmin-1), 1);

        // // Y edges
        zcopy_(ny-1, &fUp->point(     1, 1, zmin  ), fUp->local_x()*fUp->local_z(), &fUp->point(xmax  , 1, zmax  ), fUp->local_x()*fUp->local_z());
        zcopy_(ny-1, &fUp->point(xmax-1, 1, zmin  ), fUp->local_x()*fUp->local_z(), &fUp->point(0     , 1, zmax  ), fUp->local_x()*fUp->local_z());
        zcopy_(ny-1, &fUp->point(     1, 1, zmax-1), fUp->local_x()*fUp->local_z(), &fUp->point(xmax  , 1, zmin-1), fUp->local_x()*fUp->local_z());
        zcopy_(ny-1, &fUp->point(xmax-1, 1, zmax-1), fUp->local_x()*fUp->local_z(), &fUp->point(0     , 1, zmin-1), fUp->local_x()*fUp->local_z());

        // // Z edges
        zcopy_(nz-1, &fUp->point(1     , 1     , 1), fUp->local_x(), &fUp->point(xmax, ymax, 1), fUp->local_x());
        zcopy_(nz-1, &fUp->point(xmax-1, 1     , 1), fUp->local_x(), &fUp->point(0   , ymax, 1), fUp->local_x());
        zcopy_(nz-1, &fUp->point(1     , ymax-1, 1), fUp->local_x(), &fUp->point(xmax, 0   , 1), fUp->local_x());
        zcopy_(nz-1, &fUp->point(xmax-1, ymax-1, 1), fUp->local_x(), &fUp->point(0   , 0   , 1), fUp->local_x());

        // //Corners
        fUp->point(xmax, ymax, zmax  ) = fUp->point(1     , 1     , zmin  );
        fUp->point(0   , ymax, zmax  ) = fUp->point(xmax-1, 1     , zmin  );
        fUp->point(xmax, 0   , zmax  ) = fUp->point(1     , ymax-1, zmin  );
        fUp->point(0   , 0   , zmax  ) = fUp->point(xmax-1, ymax-1, zmin  );

        fUp->point(xmax, ymax, zmin-1) = fUp->point(1     , 1     , zmax-1);
        fUp->point(0   , ymax, zmin-1) = fUp->point(xmax-1, 1     , zmax-1);
        fUp->point(xmax, 0   , zmin-1) = fUp->point(1     , ymax-1, zmax-1);
        fUp->point(0   , 0   , zmin-1) = fUp->point(xmax-1, ymax-1, zmax-1);
    }
    else
    {
        zcopy_(fUp->local_y(), &fUp->point(xmax-1, 0), fUp->local_x(), &fUp->point(0   , 0), fUp->local_x() );
        zcopy_(fUp->local_y(), &fUp->point(1     , 0), fUp->local_x(), &fUp->point(xmax, 0), fUp->local_x() );
    }
    // std::array<double, 3> r;
    // for(int kk = zmin; kk < nz; ++kk)
    // {
    //     for(int jj = 1; jj < ny; ++jj)
    //     {
    //         r = {{(xmax-1)*dx,(jj-1)*dy, kk > 0 ? kk*dz : (kk-1)*dz}};
    //         fUp->point(0,jj,kk) = (std::exp(cplx(0.0,-1.0) * cplx(std::inner_product(r.begin(),r.end(), k_point.begin(),0) ) ) ) * fUp->point(xmax-1,jj,kk);
    //         r = {{0*dx,(jj-1)*dy, kk > 0 ? kk*dz : (kk-1)*dz}};
    //         fUp->point(xmax,jj,kk) = (std::exp(cplx(0.0,-1.0) * cplx(std::inner_product(r.begin(),r.end(), k_point.begin(),0) ) ) ) * fUp->point(1,jj,kk);
    //     }
    // }
    // if(zmin != 0)
    // {
    //     for(int ii = 1; ii < nx; ++ii)
    //     {
    //         for(int jj = 1; jj < ny; ++jj)
    //         {
    //             r = {{(ii-1)*dx,(jj-1)*dy, (zmax-1)*dz }};
    //             fUp->point(ii,jj,0) = (std::exp(cplx(0.0,-1.0) * cplx(std::inner_product(r.begin(),r.end(), k_point.begin(),0) ) ) ) * fUp->point(ii,jj,zmax-1);
    //             r = {{(ii-1)*dx,(jj-1)*dy, 0*dz}};
    //             fUp->point(ii,jj,zmax) = (std::exp(cplx(0.0,-1.0) * cplx(std::inner_product(r.begin(),r.end(), k_point.begin(),0) ) ) ) * fUp->point(ii,jj,zmin);
    //         }
    //     }
    // }
}

void FDTDCompUpdateFxnReal::applyPBC1Proc(pgrid_ptr fUp, std::array<double,3> & k_point, int nx, int ny, int nz, int xmax, int ymax, int zmin, int zmax, double & dx, double & dy, double & dz)
{
    // if real fields then the PBC is just copying from one side to the other
    if(zmin != 0)
    {
        for(int kk = zmin; kk < nz; ++kk)
        {
            dcopy_(nx-1, &fUp->point(1   , ymax-1, kk),  1, &fUp->point(1 , 0   , kk), 1 );
            dcopy_(nx-1, &fUp->point(1   , 1     , kk),  1, &fUp->point(1 , ymax, kk), 1 );
        }
        for(int jj = 1; jj < ny; ++jj)
        {
            dcopy_(nz-1, &fUp->point(xmax-1, jj, 1), fUp->local_x(), &fUp->point(0   , jj, 1), fUp->local_x() );
            dcopy_(nz-1, &fUp->point(1     , jj, 1), fUp->local_x(), &fUp->point(xmax, jj, 1), fUp->local_x() );

            dcopy_(nx-1, &fUp->point(1, jj, zmax-1), 1, &fUp->point(1, jj, zmin-1), 1 );
            dcopy_(nx-1, &fUp->point(1, jj, zmin  ), 1, &fUp->point(1, jj, zmax  ), 1 );
        }

        // X edges
        dcopy_(nx-1, &fUp->point(1,      1, zmin  ), 1, &fUp->point(1, ymax  , zmax  ), 1);
        dcopy_(nx-1, &fUp->point(1, ymax-1, zmin  ), 1, &fUp->point(1, 0     , zmax  ), 1);
        dcopy_(nx-1, &fUp->point(1,      1, zmax-1), 1, &fUp->point(1, ymax  , zmin-1), 1);
        dcopy_(nx-1, &fUp->point(1, ymax-1, zmax-1), 1, &fUp->point(1, 0     , zmin-1), 1);

        // Y edges
        dcopy_(ny-1, &fUp->point(     1, 1, zmin  ), fUp->local_x()*fUp->local_z(), &fUp->point(xmax  , 1, zmax  ), fUp->local_x()*fUp->local_z());
        dcopy_(ny-1, &fUp->point(xmax-1, 1, zmin  ), fUp->local_x()*fUp->local_z(), &fUp->point(0     , 1, zmax  ), fUp->local_x()*fUp->local_z());
        dcopy_(ny-1, &fUp->point(     1, 1, zmax-1), fUp->local_x()*fUp->local_z(), &fUp->point(xmax  , 1, zmin-1), fUp->local_x()*fUp->local_z());
        dcopy_(ny-1, &fUp->point(xmax-1, 1, zmax-1), fUp->local_x()*fUp->local_z(), &fUp->point(0     , 1, zmin-1), fUp->local_x()*fUp->local_z());

        // Z edges
        dcopy_(nz-1, &fUp->point(1     , 1     , 1), fUp->local_x(), &fUp->point(xmax, ymax, 1), fUp->local_x());
        dcopy_(nz-1, &fUp->point(xmax-1, 1     , 1), fUp->local_x(), &fUp->point(0   , ymax, 1), fUp->local_x());
        dcopy_(nz-1, &fUp->point(1     , ymax-1, 1), fUp->local_x(), &fUp->point(xmax, 0   , 1), fUp->local_x());
        dcopy_(nz-1, &fUp->point(xmax-1, ymax-1, 1), fUp->local_x(), &fUp->point(0   , 0   , 1), fUp->local_x());

        //Corners
        fUp->point(xmax, ymax, zmax  ) = fUp->point(1     , 1     , zmin  );
        fUp->point(0   , ymax, zmax  ) = fUp->point(xmax-1, 1     , zmin  );
        fUp->point(xmax, 0   , zmax  ) = fUp->point(1     , ymax-1, zmin  );
        fUp->point(0   , 0   , zmax  ) = fUp->point(xmax-1, ymax-1, zmin  );

        fUp->point(xmax, ymax, zmin-1) = fUp->point(1     , 1     , zmax-1);
        fUp->point(0   , ymax, zmin-1) = fUp->point(xmax-1, 1     , zmax-1);
        fUp->point(xmax, 0   , zmin-1) = fUp->point(1     , ymax-1, zmax-1);
        fUp->point(0   , 0   , zmin-1) = fUp->point(xmax-1, ymax-1, zmax-1);
    }
    else
    {
        dcopy_(nx-1, &fUp->point(1     , ymax-1),  1            , &fUp->point(1   , 0   ),              1 );
        dcopy_(nx-1, &fUp->point(1     , 1     ),  1            , &fUp->point(1   , ymax),              1 );
        dcopy_(ny-1, &fUp->point(xmax-1, 1     ), fUp->local_x(), &fUp->point(0   , 1   ), fUp->local_x() );
        dcopy_(ny-1, &fUp->point(1     , 1     ), fUp->local_x(), &fUp->point(xmax, 1   ), fUp->local_x() );
    }
}

void FDTDCompUpdateFxnCplx::applyPBC1Proc(pgrid_ptr fUp, std::array<double,3> & k_point, int nx, int ny, int nz, int xmax, int ymax, int zmin, int zmax, double & dx, double & dy, double & dz)
{
    if(zmin != 0)
    {
        for(int kk = zmin; kk < nz; ++kk)
        {
            zcopy_(nx-1, &fUp->point(1   , ymax-1, kk),  1, &fUp->point(1 , 0   , kk), 1 );
            zcopy_(nx-1, &fUp->point(1   , 1     , kk),  1, &fUp->point(1 , ymax, kk), 1 );
        }
        for(int jj = 1; jj < ny; ++jj)
        {
            zcopy_(nz-1, &fUp->point(xmax-1, jj, 1), fUp->local_x(), &fUp->point(0   , jj, 1), fUp->local_x() );
            zcopy_(nz-1, &fUp->point(1     , jj, 1), fUp->local_x(), &fUp->point(xmax, jj, 1), fUp->local_x() );

            zcopy_(nx-1, &fUp->point(1, jj, zmax-1), 1, &fUp->point(1, jj, zmin-1), 1 );
            zcopy_(nx-1, &fUp->point(1, jj, zmin  ), 1, &fUp->point(1, jj, zmax  ), 1 );
        }

        // X edges
        zcopy_(nx-1, &fUp->point(1,      1, zmin  ), 1, &fUp->point(1, ymax  , zmax  ), 1);
        zcopy_(nx-1, &fUp->point(1, ymax-1, zmin  ), 1, &fUp->point(1, 0     , zmax  ), 1);
        zcopy_(nx-1, &fUp->point(1,      1, zmax-1), 1, &fUp->point(1, ymax  , zmin-1), 1);
        zcopy_(nx-1, &fUp->point(1, ymax-1, zmax-1), 1, &fUp->point(1, 0     , zmin-1), 1);

        // Y edges
        zcopy_(ny-1, &fUp->point(     1, 1, zmin  ), fUp->local_x()*fUp->local_z(), &fUp->point(xmax  , 1, zmax  ), fUp->local_x()*fUp->local_z());
        zcopy_(ny-1, &fUp->point(xmax-1, 1, zmin  ), fUp->local_x()*fUp->local_z(), &fUp->point(0     , 1, zmax  ), fUp->local_x()*fUp->local_z());
        zcopy_(ny-1, &fUp->point(     1, 1, zmax-1), fUp->local_x()*fUp->local_z(), &fUp->point(xmax  , 1, zmin-1), fUp->local_x()*fUp->local_z());
        zcopy_(ny-1, &fUp->point(xmax-1, 1, zmax-1), fUp->local_x()*fUp->local_z(), &fUp->point(0     , 1, zmin-1), fUp->local_x()*fUp->local_z());

        // Z edges
        zcopy_(nz-1, &fUp->point(1     , 1     , 1), fUp->local_x(), &fUp->point(xmax, ymax, 1), fUp->local_x());
        zcopy_(nz-1, &fUp->point(xmax-1, 1     , 1), fUp->local_x(), &fUp->point(0   , ymax, 1), fUp->local_x());
        zcopy_(nz-1, &fUp->point(1     , ymax-1, 1), fUp->local_x(), &fUp->point(xmax, 0   , 1), fUp->local_x());
        zcopy_(nz-1, &fUp->point(xmax-1, ymax-1, 1), fUp->local_x(), &fUp->point(0   , 0   , 1), fUp->local_x());

        //Corners
        fUp->point(xmax, ymax, zmax  ) = fUp->point(1     , 1     , zmin  );
        fUp->point(0   , ymax, zmax  ) = fUp->point(xmax-1, 1     , zmin  );
        fUp->point(xmax, 0   , zmax  ) = fUp->point(1     , ymax-1, zmin  );
        fUp->point(0   , 0   , zmax  ) = fUp->point(xmax-1, ymax-1, zmin  );

        fUp->point(xmax, ymax, zmin-1) = fUp->point(1     , 1     , zmax-1);
        fUp->point(0   , ymax, zmin-1) = fUp->point(xmax-1, 1     , zmax-1);
        fUp->point(xmax, 0   , zmin-1) = fUp->point(1     , ymax-1, zmax-1);
        fUp->point(0   , 0   , zmin-1) = fUp->point(xmax-1, ymax-1, zmax-1);
    }
    else
    {
        zcopy_(nx-1, &fUp->point(1     , ymax-1),  1            , &fUp->point(1   , 0   ),              1 );
        zcopy_(nx-1, &fUp->point(1     , 1     ),  1            , &fUp->point(1   , ymax),              1 );
        zcopy_(ny-1, &fUp->point(xmax-1, 1     ), fUp->local_x(), &fUp->point(0   , 1   ), fUp->local_x() );
        zcopy_(ny-1, &fUp->point(1     , 1     ), fUp->local_x(), &fUp->point(xmax, 1   ), fUp->local_x() );
    }
    // std::array<double, 3> r;
    // for(int kk = zmin; kk < zmax; ++kk)
    // {
    //     for(int ii = 1; ii < nx; ++ii)
    //     {
    //         // Calculate the r factor and then update the fields with the complex values
    //         r = {{(ii-1)*dx,(ymax-1)*dy, kk > 0 ? kk*dz : (kk-1)*dz }};
    //         fUp->point(ii,0,kk) = (std::exp(cplx(0.0,-1.0) * cplx(std::inner_product(r.begin(),r.end(), k_point.begin(),0) ) ) ) * fUp->point(ii,ymax-1,kk);
    //         r = {{(ii-1)*dx,0.0*dy, kk > 0 ? kk*dz : (kk-1)*dz }};
    //         fUp->point(ii,ymax,kk) = (std::exp(cplx(0.0,-1.0) * cplx(std::inner_product(r.begin(),r.end(), k_point.begin(),0) ) ) ) * fUp->point(ii,1,kk);
    //     }
    // }
    // for(int kk = zmin; kk < zmax; ++kk)
    // {
    //     for(int jj = 1; jj < ny; ++jj)
    //     {
    //         r = {{(xmax-1)*dx,(jj-1)*dy, kk > 0 ? kk*dz : (kk-1)*dz }};
    //         fUp->point(0,jj,kk) = (std::exp(cplx(0.0,-1.0) * cplx(std::inner_product(r.begin(),r.end(), k_point.begin(),0) ) ) ) * fUp->point(xmax-1,jj,kk);
    //         r = {{0*dx,(jj-1)*dy, kk > 0 ? kk*dz : (kk-1)*dz }};
    //         fUp->point(xmax,jj,kk) = (std::exp(cplx(0.0,-1.0) * cplx(std::inner_product(r.begin(),r.end(), k_point.begin(),0) ) ) ) * fUp->point(1,jj,kk);
    //     }
    // }
    // if(zmin != 0)
    // {
    //     for(int ii = 1; ii < nx; ++ii)
    //     {
    //         for(int jj = 1; jj < ny; ++jj)
    //         {
    //             r = {{(ii-1)*dx,(jj-1)*dy, (zmax-1)*dz }};
    //             fUp->point(ii,jj,0) = (std::exp(cplx(0.0,-1.0) * cplx(std::inner_product(r.begin(),r.end(), k_point.begin(),0) ) ) ) * fUp->point(ii,jj,zmax-1);
    //             r = {{(ii-1)*dx,(jj-1)*dy, 0*dz}};
    //             fUp->point(ii,jj,zmax) = (std::exp(cplx(0.0,-1.0) * cplx(std::inner_product(r.begin(),r.end(), k_point.begin(),0) ) ) ) * fUp->point(ii,jj,zmin);
    //         }
    //     }
    // }
}