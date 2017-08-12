#include <SOURCE/parallelTFSF.hpp>
void tfsfUpdateFxnReal::addTFSFTwoComp(real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur)
{
    for(int ll = 0; ll < sur->szTrans_j_[1]; ++ll)
    {
        zcopy_(sur->szTrans_j_[0], &incd->point(sur->incdStart_j_+ll*sur->addIncdProp_,0), sur->strideIncd_, incdTransfer, 1);
        zscal_(sur->szTrans_j_[0], sur->prefactor_j_, incdTransfer, 1);
        daxpy_(sur->szTrans_j_[0], 1.0, reinterpret_cast<double*>(incdTransfer ), 2, &grid_j->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1], sur->loc_[2]+ll*sur->addVec_[2]), sur->strideField_);
    }
    for(int ll = 0; ll < sur->szTrans_k_[1]; ++ll)
    {
        zcopy_(sur->szTrans_k_[0], &incd->point(sur->incdStart_k_+ll*sur->addIncdProp_,0), sur->strideIncd_, incdTransfer, 1);
        zscal_(sur->szTrans_k_[0], sur->prefactor_k_, incdTransfer, 1);
        daxpy_(sur->szTrans_k_[0], 1.0, reinterpret_cast<double*>(incdTransfer ), 2, &grid_k->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1], sur->loc_[2]+ll*sur->addVec_[2]), sur->strideField_);
    }
}

void tfsfUpdateFxnReal::addTFSFOneCompJ(real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur)
{
    for(int ll = 0; ll < sur->szTrans_j_[1]; ++ll)
    {
        // std::transform(&incd->point(sur->incdStart_j_+ll*sur->addIncdProp_,0), &incd->point(sur->incdStart_j_+ll*sur->addIncdProp_,0) + sur->szTrans_j_[0], &grid_j->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1]), &grid_j->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1]), [&](cplx incd, double field) { return field + std::real(incd * sur->prefactor_j_); } );
        // daxpy_(sur->szTrans_j_[0], sur->prefactor_j_, &incd->point(sur->incdStart_j_+ll*sur->addIncdProp_,0), sur->strideIncd_, &grid_j->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1]), sur->strideField_);
        zcopy_(sur->szTrans_j_[0], &incd->point(sur->incdStart_j_+ll*sur->addIncdProp_,0), sur->strideIncd_, incdTransfer, 1);
        zscal_(sur->szTrans_j_[0], sur->prefactor_j_, incdTransfer, 1);
        daxpy_(sur->szTrans_j_[0], 1.0, reinterpret_cast<double*>(incdTransfer ), 2, &grid_j->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1], sur->loc_[2]+ll*sur->addVec_[2]), sur->strideField_);

    }
}

void tfsfUpdateFxnReal::addTFSFOneCompK(real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur)
{
    for(int ll = 0; ll < sur->szTrans_k_[1]; ++ll)
    {
        // std::transform(&incd->point(sur->incdStart_k_+ll*sur->addIncdProp_,0), &incd->point(sur->incdStart_k_+ll*sur->addIncdProp_,0) + sur->szTrans_k_[0], &grid_k->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1]), &grid_k->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1]), [&](cplx incd, double field) { return field + std::real(incd * sur->prefactor_k_); } );
        // daxpy_(sur->szTrans_k_[0], sur->prefactor_k_, &incd->point(sur->incdStart_k_+ll*sur->addIncdProp_,0), sur->strideIncd_, &grid_k->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1]), sur->strideField_);
        zcopy_(sur->szTrans_k_[0], &incd->point(sur->incdStart_k_+ll*sur->addIncdProp_,0), sur->strideIncd_, incdTransfer, 1);
        zscal_(sur->szTrans_k_[0], sur->prefactor_k_, incdTransfer, 1);
        daxpy_(sur->szTrans_k_[0], 1.0, reinterpret_cast<double*>(incdTransfer ), 2, &grid_k->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1], sur->loc_[2]+ll*sur->addVec_[2]), sur->strideField_);
    }
}

void tfsfUpdateFxnReal::transferDat(real_pgrid_ptr grid)
{
    grid->transferDat();
}

void tfsfUpdateFxnCplx::addTFSFTwoComp(cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur)
{
    for(int ll = 0; ll < sur->szTrans_j_[1]; ++ll)
        zaxpy_(sur->szTrans_j_[0], sur->prefactor_j_, &incd->point(sur->incdStart_j_+ll*sur->addIncdProp_,0), sur->strideIncd_, &grid_j->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1], sur->loc_[2]+ll*sur->addVec_[2]), sur->strideField_);
    for(int ll = 0; ll < sur->szTrans_k_[1]; ++ll)
        zaxpy_(sur->szTrans_k_[0], sur->prefactor_k_, &incd->point(sur->incdStart_k_+ll*sur->addIncdProp_,0), sur->strideIncd_, &grid_k->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1], sur->loc_[2]+ll*sur->addVec_[2]), sur->strideField_);
}

void tfsfUpdateFxnCplx::addTFSFOneCompJ(cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur)
{
    for(int ll = 0; ll < sur->szTrans_j_[1]; ++ll)
        zaxpy_(sur->szTrans_j_[0], sur->prefactor_j_, &incd->point(sur->incdStart_j_+ll*sur->addIncdProp_,0), sur->strideIncd_, &grid_j->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1]), sur->strideField_);
}

void tfsfUpdateFxnCplx::addTFSFOneCompK(cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur)
{
    for(int ll = 0; ll < sur->szTrans_k_[1]; ++ll)
        zaxpy_(sur->szTrans_k_[0], sur->prefactor_k_, &incd->point(sur->incdStart_k_+ll*sur->addIncdProp_,0), sur->strideIncd_, &grid_k->point(sur->loc_[0]+ll*sur->addVec_[0], sur->loc_[1]+ll*sur->addVec_[1]), sur->strideField_);
}
void tfsfUpdateFxnCplx::transferDat(cplx_pgrid_ptr grid)
{
    grid->transferDat();
}


parallelTFSFReal::parallelTFSFReal(std::shared_ptr<mpiInterface> gridComm, std::array<int,3> loc, std::array<int,3> sz, double theta, double phi, double psi, POLARIZATION circPol, double kLenRelJ, double dx, double dt, std::vector<std::shared_ptr<PulseBase>> pul, real_pgrid_ptr Ex, real_pgrid_ptr Ey, real_pgrid_ptr Ez, real_pgrid_ptr Hx, real_pgrid_ptr Hy, real_pgrid_ptr Hz) :
    parallelTFSFBase<double>(gridComm, loc, sz, theta, phi, psi, circPol, kLenRelJ, dx, dt, pul, Ex, Ey, Ez, Hx, Hy,  Hz)
{
    if(Ez_ && Hz_)
    {
        if(botSurE_)
            addEBot_ = tfsfUpdateFxnReal::addTFSFTwoComp;
        else
            addEBot_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        if(botSurH_)
            addHBot_ = tfsfUpdateFxnReal::addTFSFTwoComp;
        else
            addHBot_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(topSurE_)
            addETop_ = tfsfUpdateFxnReal::addTFSFTwoComp;
        else
            addETop_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(topSurH_)
            addHTop_ = tfsfUpdateFxnReal::addTFSFTwoComp;
        else
            addHTop_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(leftSurE_)
            addELeft_ = tfsfUpdateFxnReal::addTFSFTwoComp;
        else
            addELeft_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(leftSurH_)
            addHLeft_ = tfsfUpdateFxnReal::addTFSFTwoComp;
        else
            addHLeft_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(rightSurE_)
            addERight_ = tfsfUpdateFxnReal::addTFSFTwoComp;
        else
            addERight_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(rightSurH_)
            addHRight_ = tfsfUpdateFxnReal::addTFSFTwoComp;
        else
            addHRight_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(backSurE_)
            addEBack_ = tfsfUpdateFxnReal::addTFSFTwoComp;
        else
            addEBack_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(backSurH_)
            addHBack_ = tfsfUpdateFxnReal::addTFSFTwoComp;
        else
            addHBack_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(frontSurE_)
            addEFront_ = tfsfUpdateFxnReal::addTFSFTwoComp;
        else
            addEFront_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(frontSurH_)
            addHFront_ = tfsfUpdateFxnReal::addTFSFTwoComp;
        else
            addHFront_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
    }
    else if(Hz)
    {
        if(botSurH_)
        {
            addHBot_ = tfsfUpdateFxnReal::addTFSFOneCompJ;
        }
        else
        {
            addHBot_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(botSurE_)
        {
            addEBot_ = tfsfUpdateFxnReal::addTFSFOneCompK;
        }
        else
        {
            addEBot_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(topSurH_)
        {
            addHTop_ = tfsfUpdateFxnReal::addTFSFOneCompJ;
        }
        else
        {
            addHTop_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(topSurE_)
        {
            addETop_ = tfsfUpdateFxnReal::addTFSFOneCompK;
        }
        else
        {
            addETop_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(leftSurH_)
        {
            addHLeft_ = tfsfUpdateFxnReal::addTFSFOneCompK;
        }
        else
        {
            addHLeft_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(leftSurE_)
        {
            addELeft_ = tfsfUpdateFxnReal::addTFSFOneCompJ;
        }
        else
        {
            addELeft_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(rightSurH_)
        {
            addHRight_ = tfsfUpdateFxnReal::addTFSFOneCompK;
        }
        else
        {
            addHRight_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(rightSurE_)
        {
            addERight_ = tfsfUpdateFxnReal::addTFSFOneCompJ;
        }
        else
        {
            addERight_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        addEBack_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        addHBack_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        addEFront_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        addHFront_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
    }
    else if(Ez)
    {
        if(botSurH_)
        {
            addHBot_ = tfsfUpdateFxnReal::addTFSFOneCompK;
        }
        else
        {
            addHBot_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(botSurE_)
        {
            addEBot_ = tfsfUpdateFxnReal::addTFSFOneCompJ;
        }
        else
        {
            addEBot_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(topSurH_)
        {
            addHTop_ = tfsfUpdateFxnReal::addTFSFOneCompK;
        }
        else
        {
            addHTop_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(topSurE_)
        {
            addETop_ = tfsfUpdateFxnReal::addTFSFOneCompJ;
        }
        else
        {
            addETop_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(leftSurH_)
        {
            addHLeft_ = tfsfUpdateFxnReal::addTFSFOneCompJ;
        }
        else
        {
            addHLeft_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(leftSurE_)
        {
            addELeft_ = tfsfUpdateFxnReal::addTFSFOneCompK;
        }
        else
        {
            addELeft_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(rightSurH_)
        {
            addHRight_ = tfsfUpdateFxnReal::addTFSFOneCompJ;
        }
        else
        {
            addHRight_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(rightSurE_)
        {
            addERight_ = tfsfUpdateFxnReal::addTFSFOneCompK;
        }
        else
        {
            addERight_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        addEBack_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        addHBack_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        addEFront_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        addHFront_ = [](real_pgrid_ptr grid_j, real_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
    }
}
parallelTFSFCplx::parallelTFSFCplx(std::shared_ptr<mpiInterface> gridComm, std::array<int,3> loc, std::array<int,3> sz, double theta, double phi, double psi, POLARIZATION circPol, double kLenRelJ, double dx, double dt, std::vector<std::shared_ptr<PulseBase>> pul, cplx_pgrid_ptr Ex, cplx_pgrid_ptr Ey, cplx_pgrid_ptr Ez, cplx_pgrid_ptr Hx, cplx_pgrid_ptr Hy, cplx_pgrid_ptr Hz) :
    parallelTFSFBase<cplx>(gridComm, loc, sz, theta, phi, psi, circPol, kLenRelJ, dx, dt, pul, Ex, Ey, Ez, Hx, Hy, Hz)
{
    if(Ez_ && Hz_)
    {
        if(botSurE_)
            addEBot_ = tfsfUpdateFxnCplx::addTFSFTwoComp;
        else
            addEBot_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        if(botSurH_)
            addHBot_ = tfsfUpdateFxnCplx::addTFSFTwoComp;
        else
            addHBot_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(topSurE_)
            addETop_ = tfsfUpdateFxnCplx::addTFSFTwoComp;
        else
            addETop_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(topSurH_)
            addHTop_ = tfsfUpdateFxnCplx::addTFSFTwoComp;
        else
            addHTop_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(leftSurE_)
            addELeft_ = tfsfUpdateFxnCplx::addTFSFTwoComp;
        else
            addELeft_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(leftSurH_)
            addHLeft_ = tfsfUpdateFxnCplx::addTFSFTwoComp;
        else
            addHLeft_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(rightSurE_)
            addERight_ = tfsfUpdateFxnCplx::addTFSFTwoComp;
        else
            addERight_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(rightSurH_)
            addHRight_ = tfsfUpdateFxnCplx::addTFSFTwoComp;
        else
            addHRight_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(backSurE_)
            addEBack_ = tfsfUpdateFxnCplx::addTFSFTwoComp;
        else
            addEBack_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(backSurH_)
            addHBack_ = tfsfUpdateFxnCplx::addTFSFTwoComp;
        else
            addHBack_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(frontSurE_)
            addEFront_ = tfsfUpdateFxnCplx::addTFSFTwoComp;
        else
            addEFront_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};

        if(frontSurH_)
            addHFront_ = tfsfUpdateFxnCplx::addTFSFTwoComp;
        else
            addHFront_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
    }
    else if(Hz)
    {
        if(botSurH_)
        {
            addHBot_ = tfsfUpdateFxnCplx::addTFSFOneCompJ;
        }
        else
        {
            addHBot_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(botSurE_)
        {
            addEBot_ = tfsfUpdateFxnCplx::addTFSFOneCompK;
        }
        else
        {
            addEBot_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(topSurH_)
        {
            addHTop_ = tfsfUpdateFxnCplx::addTFSFOneCompJ;
        }
        else
        {
            addHTop_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(topSurE_)
        {
            addETop_ = tfsfUpdateFxnCplx::addTFSFOneCompK;
        }
        else
        {
            addETop_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(leftSurH_)
        {
            addHLeft_ = tfsfUpdateFxnCplx::addTFSFOneCompK;
        }
        else
        {
            addHLeft_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(leftSurE_)
        {
            addELeft_ = tfsfUpdateFxnCplx::addTFSFOneCompJ;
        }
        else
        {
            addELeft_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(rightSurH_)
        {
            addHRight_ = tfsfUpdateFxnCplx::addTFSFOneCompK;
        }
        else
        {
            addHRight_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(rightSurE_)
        {
            addERight_ = tfsfUpdateFxnCplx::addTFSFOneCompJ;
        }
        else
        {
            addERight_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        addEBack_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        addHBack_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        addEFront_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        addHFront_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
    }
    else if(Ez)
    {
        if(botSurH_)
        {
            addHBot_ = tfsfUpdateFxnCplx::addTFSFOneCompK;
        }
        else
        {
            addHBot_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(botSurE_)
        {
            addEBot_ = tfsfUpdateFxnCplx::addTFSFOneCompJ;
        }
        else
        {
            addEBot_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(topSurH_)
        {
            addHTop_ = tfsfUpdateFxnCplx::addTFSFOneCompK;
        }
        else
        {
            addHTop_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(topSurE_)
        {
            addETop_ = tfsfUpdateFxnCplx::addTFSFOneCompJ;
        }
        else
        {
            addETop_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(leftSurH_)
        {
            addHLeft_ = tfsfUpdateFxnCplx::addTFSFOneCompJ;
        }
        else
        {
            addHLeft_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(leftSurE_)
        {
            addELeft_ = tfsfUpdateFxnCplx::addTFSFOneCompK;
        }
        else
        {
            addELeft_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(rightSurH_)
        {
            addHRight_ = tfsfUpdateFxnCplx::addTFSFOneCompJ;
        }
        else
        {
            addHRight_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        if(rightSurE_)
        {
            addERight_ = tfsfUpdateFxnCplx::addTFSFOneCompK;
        }
        else
        {
            addERight_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        }

        addEBack_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        addHBack_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        addEFront_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
        addHFront_ = [](cplx_pgrid_ptr grid_j, cplx_pgrid_ptr grid_k, cplx_grid_ptr incd, cplx* incdTransfer, std::shared_ptr<paramStoreTFSF> sur){return;};
    }
}