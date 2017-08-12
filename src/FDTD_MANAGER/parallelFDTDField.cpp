#include <FDTD_MANAGER/parallelFDTDField.hpp>
#include <iomanip>

parallelFDTDFieldReal::parallelFDTDFieldReal(parallelProgramInputs &IP, std::shared_ptr<mpiInterface> gridComm) :
    parallelFDTDFieldBase<double>(IP, gridComm)
{
    // If there is an Hz field set up all TE mode functions otherwise set them to do nothing
    if(Hz_)
    {
        // Initialize the PMLs
        if(magMatInPML_)
        {
            HzPML_   = std::make_shared<parallelCPMLReal>(gridComm_, weights_, Bz_, Ex_, Ey_, POLARIZATION::HZ, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Hz_, objArr_);
        }
        else
        {
            HzPML_   = std::make_shared<parallelCPMLReal>(gridComm_, weights_, Hz_, Ex_, Ey_, POLARIZATION::HZ, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Hz_, objArr_);
        }
        if(dielectricMatInPML_)
        {
            ExPML_   = std::make_shared<parallelCPMLReal>(gridComm_, weights_, Dx_, Hy_, Hz_, POLARIZATION::EX, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Ex_, objArr_);
            EyPML_   = std::make_shared<parallelCPMLReal>(gridComm_, weights_, Dy_, Hz_, Hx_, POLARIZATION::EY, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Ey_, objArr_);
        }
        else
        {
            ExPML_   = std::make_shared<parallelCPMLReal>(gridComm_, weights_, Ex_, Hy_, Hz_, POLARIZATION::EX, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Ex_, objArr_);
            EyPML_   = std::make_shared<parallelCPMLReal>(gridComm_, weights_, Ey_, Hz_, Hx_, POLARIZATION::EY, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Ey_, objArr_);
        }
        // Fill the update Lists
        initializeList(phys_Hz_, HzPML_, false, std::array<int,3>( {{ 1,  0,  0 }} ), std::array<int,3>( {{n_vec_[0]  , n_vec_[1]  , n_vec_[2]+1 }} ), d_[0], upHz_, upBz_);
        initializeList(phys_Ex_, ExPML_,  true, std::array<int,3>( {{ 0, -1,  0 }} ), std::array<int,3>( {{n_vec_[0]  , n_vec_[1]+1, n_vec_[2]+1 }} ), d_[1], upEx_, upDx_);
        initializeList(phys_Ey_, EyPML_,  true, std::array<int,3>( {{ 0,  0, -1 }} ), std::array<int,3>( {{n_vec_[0]+1, n_vec_[1]  , n_vec_[2]+1 }} ), d_[0], upEy_, upDy_);

        if(!Ez_)
        {
            upExFxn_ = &FDTDCompUpdateFxnReal::OneCompCurlK;
            upEyFxn_ = &FDTDCompUpdateFxnReal::OneCompCurlJ;
        }
        else
        {
            upExFxn_ = &FDTDCompUpdateFxnReal::TwoCompCurl;
            upEyFxn_ = &FDTDCompUpdateFxnReal::TwoCompCurl;
        }

        upHzFxn_ = &FDTDCompUpdateFxnReal::TwoCompCurl;

        updateExPML_ = [](pml_ptr pml){pml->updateGrid();};
        updateEyPML_ = [](pml_ptr pml){pml->updateGrid();};
        updateHzPML_ = [](pml_ptr pml){pml->updateGrid();};

        upLorPxFxn_ = &FDTDCompUpdateFxnReal::UpdateLorPol;
        upLorPyFxn_ = &FDTDCompUpdateFxnReal::UpdateLorPol;
        upLorMzFxn_ = &FDTDCompUpdateFxnReal::UpdateLorMag;

        D2ExFxn_ = &FDTDCompUpdateFxnReal::DtoE;
        D2EyFxn_ = &FDTDCompUpdateFxnReal::DtoE;

        B2HzFxn_ = &FDTDCompUpdateFxnReal::BtoH;

        if(IP.periodic_ && gridComm_->size() > 1)
        {
            pbcEx_ = &FDTDCompUpdateFxnReal::applyPBC;
            pbcEy_ = &FDTDCompUpdateFxnReal::applyPBC;
            pbcHz_ = &FDTDCompUpdateFxnReal::applyPBC;

            yExPBC_ = ln_vec_[1]+1;
            yEyPBC_ = ln_vec_[1]+1;
            yHzPBC_ = ln_vec_[1]+1;
            // Field definition buffer
            if(gridComm_->size()-1 == gridComm_->rank())
            {
                yHzPBC_ -= 1;
                yEyPBC_ -= 1;
            }
        }
        else if(IP.periodic_)
        {
            pbcEx_ = &FDTDCompUpdateFxnReal::applyPBC1Proc;
            pbcEy_ = &FDTDCompUpdateFxnReal::applyPBC1Proc;
            pbcHz_ = &FDTDCompUpdateFxnReal::applyPBC1Proc;

            yHzPBC_ = ln_vec_[1];
            yExPBC_ = ln_vec_[1]+1;
            yEyPBC_ = ln_vec_[1];
        }
        else
        {
            pbcEx_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
            pbcEy_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
            pbcHz_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
        }

        // Transfer fields should do nothing if only on one processor
        if(gridComm_->size() > 1)
        {
            transferEx_ = [=](){Ex_->transferDat();};
            transferEy_ = [=](){Ey_->transferDat();};
            transferHz_ = [=](){Hz_->transferDat();};
        }
        else
        {
            transferEx_ = [](){return;};
            transferEy_ = [](){return;};
            transferHz_ = [](){return;};
        }
    }
    else
    {
        upHzFxn_ = [](std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr){return;};

        upExFxn_ = [](std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr){return;};;
        upEyFxn_ = [](std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr){return;};;

        updateExPML_ = [](pml_ptr pml){return;};
        updateEyPML_ = [](pml_ptr pml){return;};
        updateHzPML_ = [](pml_ptr pml){return;};

        upLorPxFxn_ = []( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, double*, std::shared_ptr<Obj> ){return;};
        upLorPyFxn_ = []( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, double*, std::shared_ptr<Obj> ){return;};
        upLorMzFxn_ = []( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, double*, std::shared_ptr<Obj> ){return;};

        D2ExFxn_ = []( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj> ){return;};
        D2EyFxn_ = []( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj> ){return;};
        B2HzFxn_ = []( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj> ){return;};

        pbcEx_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
        pbcEy_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
        pbcHz_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};

        transferEx_ = [](){return;};
        transferEy_ = [](){return;};
        transferHz_ = [](){return;};
    }
    // If there is an Ez field set up all TM functions otherwise set them to do nothing
    if(Ez_)
    {
        //initialize the PMLs
        if(magMatInPML_)
        {
            HxPML_   = std::make_shared<parallelCPMLReal>(gridComm_, weights_, Bx_, Ey_, Ez_, POLARIZATION::HX, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Hx_, objArr_);
            HyPML_   = std::make_shared<parallelCPMLReal>(gridComm_, weights_, By_, Ez_, Ex_, POLARIZATION::HY, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Hy_, objArr_);
        }
        else
        {
            HxPML_   = std::make_shared<parallelCPMLReal>(gridComm_, weights_, Hx_, Ey_, Ez_, POLARIZATION::HX, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Hx_, objArr_);
            HyPML_   = std::make_shared<parallelCPMLReal>(gridComm_, weights_, Hy_, Ez_, Ex_, POLARIZATION::HY, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Hy_, objArr_);
        }
        if(dielectricMatInPML_)
        {
            EzPML_   = std::make_shared<parallelCPMLReal>(gridComm_, weights_, Dz_, Hx_, Hy_, POLARIZATION::EZ, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Ez_, objArr_);
        }
        else
        {
            EzPML_   = std::make_shared<parallelCPMLReal>(gridComm_, weights_, Ez_, Hx_, Hy_, POLARIZATION::EZ, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Ez_, objArr_);
        }
        // Fill update lists
        initializeList(phys_Ez_, EzPML_,  true, std::array<int,3>( {{-1,  0,  0 }} ), std::array<int,3>( {{n_vec_[0]+1, n_vec_[1]+1, n_vec_[2]   }} ), d_[0], upEz_, upDz_);
        initializeList(phys_Hx_, HxPML_, false, std::array<int,3>( {{ 0,  1,  0 }} ), std::array<int,3>( {{n_vec_[0]+1, n_vec_[1]  , n_vec_[2]   }} ), d_[1], upHx_, upBx_);
        initializeList(phys_Hy_, HyPML_, false, std::array<int,3>( {{ 0,  0,  1 }} ), std::array<int,3>( {{n_vec_[0]  , n_vec_[1]+1, n_vec_[2]   }} ), d_[0], upHy_, upBy_);

        if(!Hz_)
        {
            upHxFxn_ = &FDTDCompUpdateFxnReal::OneCompCurlK;
            upHyFxn_ = &FDTDCompUpdateFxnReal::OneCompCurlJ;
        }
        else
        {
            upHxFxn_ = &FDTDCompUpdateFxnReal::TwoCompCurl;
            upHyFxn_ = &FDTDCompUpdateFxnReal::TwoCompCurl;
        }

        upEzFxn_ = &FDTDCompUpdateFxnReal::TwoCompCurl;

        upLorMxFxn_ = &FDTDCompUpdateFxnReal::UpdateLorMag;
        upLorMyFxn_ = &FDTDCompUpdateFxnReal::UpdateLorMag;

        upLorPzFxn_ = &FDTDCompUpdateFxnReal::UpdateLorPol;

        B2HxFxn_ = &FDTDCompUpdateFxnReal::BtoH;
        B2HyFxn_ = &FDTDCompUpdateFxnReal::BtoH;

        D2EzFxn_ = &FDTDCompUpdateFxnReal::DtoE;

        updateHxPML_ = [](pml_ptr pml){pml->updateGrid();};
        updateHyPML_ = [](pml_ptr pml){pml->updateGrid();};
        updateEzPML_ = [](pml_ptr pml){pml->updateGrid();};

        if(IP.periodic_ && gridComm_->size() > 1)
        {
            pbcHx_ = &FDTDCompUpdateFxnReal::applyPBC;
            pbcHy_ = &FDTDCompUpdateFxnReal::applyPBC;
            pbcEz_ = &FDTDCompUpdateFxnReal::applyPBC;

            yHxPBC_ = ln_vec_[1]+1;
            yHyPBC_ = ln_vec_[1]+1;
            yEzPBC_ = ln_vec_[1]+1;
            //Field definition buffer
            if(gridComm_->size()-1 == gridComm_->rank())
                yHxPBC_ = ln_vec_[1];
        }
        else if(IP.periodic_)
        {
            pbcHx_ = &FDTDCompUpdateFxnReal::applyPBC1Proc;
            pbcHy_ = &FDTDCompUpdateFxnReal::applyPBC1Proc;
            pbcEz_ = &FDTDCompUpdateFxnReal::applyPBC1Proc;
            yHxPBC_ = ln_vec_[1];
            yHyPBC_ = ln_vec_[1]+1;
            yEzPBC_ = ln_vec_[1]+1;
        }
        else
        {
            pbcHx_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
            pbcHy_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
            pbcEz_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
        }

        // Transfer fields should do nothing if only on one processor
        if(gridComm_->size() > 1)
        {
            transferHx_ = [=](){Hx_->transferDat();};
            transferHy_ = [=](){Hy_->transferDat();};
            transferEz_ = [=](){Ez_->transferDat();};

        }
        else
        {
            transferEz_ = [](){return;};
            transferHx_ = [](){return;};
            transferHy_ = [](){return;};
        }
    }
    else
    {
        upHxFxn_ = [](std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr){return;};;
        upHyFxn_ = [](std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr){return;};

        upEzFxn_ = [](std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr){return;};

        upLorMxFxn_ = []( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, double*, std::shared_ptr<Obj> ){return;};
        upLorMyFxn_ = []( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, double*, std::shared_ptr<Obj> ){return;};
        upLorPzFxn_ = []( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, double*, std::shared_ptr<Obj> ){return;};

        B2HxFxn_ = []( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj> ){return;};
        B2HyFxn_ = []( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj> ){return;};
        D2EzFxn_ = []( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj> ){return;};

        updateHxPML_ = [](pml_ptr pml){return;};
        updateHyPML_ = [](pml_ptr pml){return;};
        updateEzPML_ = [](pml_ptr pml){return;};

        pbcHx_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
        pbcHy_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
        pbcEz_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};

        transferEz_ = [](){return;};
        transferHx_ = [](){return;};
        transferHy_ = [](){return;};
    }
    // Construct all soft sources
    for(int ss = 0; ss < IP.srcPol_.size(); ss++)
    {
        // Make the pulse (including all pulses to be used)
        std::vector<std::shared_ptr<PulseBase>> pul;
        for(int pp = 0; pp < IP.srcPulShape_[ss].size(); pp ++)
        {
            if(IP.srcPulShape_[ss][pp] == PLSSHAPE::CONTINUOUS)
                pul.push_back(std::make_shared<PulseCont>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp], dt_));
            else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::BH)
                pul.push_back(std::make_shared<PulseBH>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp], dt_));
            else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::RECT)
                pul.push_back(std::make_shared<PulseRect>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp], dt_));
            else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::GAUSSIAN)
                pul.push_back(std::make_shared<PulseGauss>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp], dt_));
            else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::RICKER)
                pul.push_back(std::make_shared<PulseRicker>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp], dt_));
            else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::RAMP_CONT)
                pul.push_back(std::make_shared<PulseRampCont>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp], dt_));
        }
        //  Make the source act on any of the fields
        if(IP.srcPol_[ss] == POLARIZATION::EX)
        {
            if(int( round(IP.srcPhi_[ss]) ) % 90 == 0)
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalReal>(gridComm_, pul, Ex_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            else
            {
                srcArr_.push_back( std::make_shared<parallelSourceObliqueReal >(gridComm_, pul, Ex_, POLARIZATION::EX, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceObliqueReal >(gridComm_, pul, Ey_, POLARIZATION::EY, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
            }
        }
        else if(IP.srcPol_[ss] == POLARIZATION::EY)
        {
            if(int( round(IP.srcPhi_[ss]) ) % 90 == 0)
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalReal>(gridComm_, pul, Ey_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            else
            {
                srcArr_.push_back( std::make_shared<parallelSourceObliqueReal >(gridComm_, pul, Ex_, POLARIZATION::EX, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceObliqueReal >(gridComm_, pul, Ey_, POLARIZATION::EY, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
            }
        }
        else if(IP.srcPol_[ss] == POLARIZATION::EZ)
        {
            if(int( round(IP.srcPhi_[ss]) ) % 90 == 0)
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalReal>(gridComm_, pul, Ez_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            else
            {
                srcArr_.push_back( std::make_shared<parallelSourceObliqueReal >(gridComm_, pul, Ez_, POLARIZATION::EZ, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
            }
        }
        else if(IP.srcPol_[ss] == POLARIZATION::HX)
        {
            if(int( round(IP.srcPhi_[ss]) ) % 90 == 0)
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalReal>(gridComm_, pul, Hx_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            else
            {
                srcArr_.push_back( std::make_shared<parallelSourceObliqueReal >(gridComm_, pul, Hy_, POLARIZATION::HY, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceObliqueReal >(gridComm_, pul, Hx_, POLARIZATION::HX, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
            }
        }
        else if(IP.srcPol_[ss] == POLARIZATION::HY)
        {
            if(int( round(IP.srcPhi_[ss]) ) % 90 == 0)
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalReal>(gridComm_, pul, Hy_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            else
            {
                srcArr_.push_back( std::make_shared<parallelSourceObliqueReal >(gridComm_, pul, Hx_, POLARIZATION::HX, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceObliqueReal >(gridComm_, pul, Hy_, POLARIZATION::HY, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
            }
        }
        else if(IP.srcPol_[ss] == POLARIZATION::HZ)
        {
            if(int( round(IP.srcPhi_[ss]) ) % 90 == 0)
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalReal>(gridComm_, pul, Hz_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            else
            {
                srcArr_.push_back( std::make_shared<parallelSourceObliqueReal >(gridComm_, pul, Hz_, POLARIZATION::HZ, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
            }
        }
        else
        {
            double axRat = IP.srcEllipticalKratio_[ss];
            double psi = IP.srcPsi_[ss];
            double psiPrefactCalc = psi;
            double alphaOff = 0.0;
            double prefactor_k_ = 1.0;
            double prefactor_j_ = 1.0;
            double c = pow(axRat, 2.0);

            // phi/psi control the light polarization angle
            psiPrefactCalc = 0.5 * asin( sqrt( ( pow(cos(2.0*psi),2.0)*4.0*c + pow( (1.0+c)*sin(2.0*psi), 2.0) ) / pow(1.0+c, 2.0) ) );
            alphaOff = acos( ( (c - 1.0)*sin(2.0*psi) ) / sqrt( pow(cos(2.0*psi),2.0)*4.0*c + pow( (1.0+c)*sin(2.0*psi), 2.0) ) );
            if(std::abs( std::tan(psi) ) > 1)
                psiPrefactCalc = M_PI/2.0 - psiPrefactCalc;
            if(IP.srcPol_[ss] == POLARIZATION::R)
                alphaOff *= -1.0;

            if( isamin_(IP.srcSz_[ss].size(), IP.srcSz_[ss].data(), 1)-1 == 0 )
            {
                prefactor_j_ *= -1.0 * cos(psiPrefactCalc);
                prefactor_k_ *= -1.0 * sin(psiPrefactCalc);
            }
            else if( isamin_(IP.srcSz_[ss].size(), IP.srcSz_[ss].data(), 1)-1 == 1 )
            {
                prefactor_j_ *=  -1.0 * cos(psiPrefactCalc);
                prefactor_k_ *=         sin(psiPrefactCalc);
            }
            else if( isamin_(IP.srcSz_[ss].size(), IP.srcSz_[ss].data(), 1)-1 == 2 )
            {
                prefactor_j_ *= -1.0*sin(psiPrefactCalc);
                prefactor_k_ *=      cos(psiPrefactCalc);
            }

            for(auto pulse : pul)
                pulse->modE0(prefactor_j_);
            std::vector<std::shared_ptr<PulseBase>> phaseOffPul;
            for(int pp = 0; pp < IP.srcPulShape_[ss].size(); pp ++)
            {
                if(IP.srcPulShape_[ss][pp] == PLSSHAPE::CONTINUOUS)
                    phaseOffPul.push_back(std::make_shared<PulseCont>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp]*prefactor_k_, dt_, alphaOff));
                else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::BH)
                    phaseOffPul.push_back(std::make_shared<PulseBH>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp]*prefactor_k_, dt_, alphaOff));
                else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::RECT)
                    phaseOffPul.push_back(std::make_shared<PulseRect>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp]*prefactor_k_, dt_, alphaOff));
                else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::GAUSSIAN)
                    phaseOffPul.push_back(std::make_shared<PulseGauss>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp]*prefactor_k_, dt_, alphaOff));
                else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::RICKER)
                    phaseOffPul.push_back(std::make_shared<PulseRicker>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp]*prefactor_k_, dt_, alphaOff));
                else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::RAMP_CONT)
                    phaseOffPul.push_back(std::make_shared<PulseRampCont>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp]*prefactor_k_, dt_, alphaOff));
            }
            if(isamin_(IP.srcSz_[ss].size(), IP.srcSz_[ss].data(), 1)-1 == 0 )
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalReal>(gridComm_, pul, Ey_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceNormalReal>(gridComm_, phaseOffPul, Ez_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            if(isamin_(IP.srcSz_[ss].size(), IP.srcSz_[ss].data(), 1)-1 == 1 )
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalReal>(gridComm_, pul, Ex_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceNormalReal>(gridComm_, phaseOffPul, Ez_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            if(isamin_(IP.srcSz_[ss].size(), IP.srcSz_[ss].data(), 1)-1 == 2 )
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalReal>(gridComm_, pul, Ex_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceNormalReal>(gridComm_, phaseOffPul, Ey_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
        }
    }

    // Construct all TFSF surfaces
    for(int tt = 0; tt < IP.tfsfSize_.size(); tt++)
    {
        if(IP.tfsfSize_[tt][0] != 0.0 || IP.tfsfSize_[tt][1] != 0.0)
        {
        // Make the pulse (including all pulses to be used)
            std::vector<std::shared_ptr<PulseBase>> pul;
            for(int pp = 0; pp < IP.tfsfPulShape_[tt].size(); pp ++)
            {
                if(IP.tfsfPulShape_[tt][pp] == PLSSHAPE::CONTINUOUS)
                    pul.push_back(std::make_shared<PulseCont>(IP.tfsfPulFxn_[tt][pp], IP.tfsfEmax_[tt][pp], dt_));
                else if(IP.tfsfPulShape_[tt][pp] == PLSSHAPE::BH)
                    pul.push_back(std::make_shared<PulseBH>(IP.tfsfPulFxn_[tt][pp], IP.tfsfEmax_[tt][pp], dt_));
                else if(IP.tfsfPulShape_[tt][pp] == PLSSHAPE::RECT)
                    pul.push_back(std::make_shared<PulseRect>(IP.tfsfPulFxn_[tt][pp], IP.tfsfEmax_[tt][pp], dt_));
                else if(IP.tfsfPulShape_[tt][pp] == PLSSHAPE::GAUSSIAN)
                    pul.push_back(std::make_shared<PulseGauss>(IP.tfsfPulFxn_[tt][pp], IP.tfsfEmax_[tt][pp], dt_));
                else if(IP.tfsfPulShape_[tt][pp] == PLSSHAPE::RICKER)
                    pul.push_back(std::make_shared<PulseRicker>(IP.tfsfPulFxn_[tt][pp], IP.tfsfEmax_[tt][pp], dt_));
                else if(IP.tfsfPulShape_[tt][pp] == PLSSHAPE::RAMP_CONT)
                    pul.push_back(std::make_shared<PulseRampCont>(IP.tfsfPulFxn_[tt][pp], IP.tfsfEmax_[tt][pp], dt_));
            }
            // TFSF will determine the correct polarization
            tfsfArr_.push_back(std::make_shared<parallelTFSFReal>(gridComm_, IP.tfsfLoc_[tt], IP.tfsfSize_[tt], IP.tfsfTheta_[tt], IP.tfsfPhi_[tt], IP.tfsfPsi_[tt], IP.tfsfCircPol_[tt], IP.tfsfEllipticalKratio_[tt], d_[0], dt_, pul, Ex_, Ey_, Ez_, Hx_, Hy_, Hz_) );
        }
    }
    // Polarization matters here since the z field always forms the continous box for the spatial offset (TE uses H, TM uses E)
    if(IP.fluxName_.size() > 0)
    {
        DIRECTION propDir;
        if(tfsfArr_.size() > 0 && (tfsfArr_.back()->theta() == 0.0 || tfsfArr_.back()->theta() == M_PI ) )
            propDir = DIRECTION::Z;
        else if(tfsfArr_.size() > 0 && (tfsfArr_.back()->quadrant() == 2 || tfsfArr_.back()->quadrant() == 4 ) )
            propDir = DIRECTION::X;
        else if(tfsfArr_.size() > 0 && (tfsfArr_.back()->quadrant() == 1 || tfsfArr_.back()->quadrant() == 3 ) )
            propDir = DIRECTION::Y;
        else if(srcArr_.size() > 0 && (srcArr_.back()->sz()[0] >= srcArr_.back()->sz()[1] ) && (srcArr_.back()->sz()[2] >= srcArr_.back()->sz()[1] ) )
            propDir = DIRECTION::Y;
        else if(srcArr_.size() > 0 && (srcArr_.back()->sz()[0] <= srcArr_.back()->sz()[1] ) && (srcArr_.back()->sz()[0] <= srcArr_.back()->sz()[2] ) )
            propDir = DIRECTION::X;
        else if(srcArr_.size() > 0)
            propDir = DIRECTION::Z;
        else
            throw std::logic_error("Constructing a flux with no source, it will be 0.");

        double theta = 0; double phi = 0; double psi = 0; double alpha = 0;
        if(tfsfArr_.size() > 0)
        {
            theta = std::atan( std::abs( std::tan( tfsfArr_.back()->theta() ) ) );
            phi   = std::atan( std::abs( std::tan( tfsfArr_.back()->phiPreFact() ) ) );
            psi   = std::atan( std::abs( std::tan( tfsfArr_.back()->psiPreFact() ) ) );
            alpha = tfsfArr_.back()->alpha();
        }
        for(int ff = 0; ff < IP.fluxLoc_.size(); ff ++)
            fluxArr_.push_back(std::make_shared<parallelFluxDTCReal>(gridComm_, IP.fluxName_[ff], IP.fluxWeight_[ff], Ex_, Ey_, Ez_, Hx_, Hy_, Hz_, IP.fluxLoc_[ff], IP.fluxSz_[ff], IP.fluxCrossSec_[ff], IP.fluxSave_[ff], IP.fluxLoad_[ff], IP.fluxTimeInt_[ff], IP.fluxFreqList_[ff], propDir, d_, dt_, theta, phi, psi, alpha, IP.fluxIncdFieldsFilename_[ff], IP.fluxSI_[ff], IP.I0_, IP.a_) );
    }
    // Construct all DTC based on types (all it changes is the list of fields it passes)
    for(int dd = 0; dd < IP.dtcType_.size(); dd++)
    {
        std::vector<pgrid_ptr> fields;
        // Fill the fields vector with the appropriate field values
        if(IP.dtcType_[dd] == DTCTYPE::EX)
            fields = {{ Ex_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::EY)
            fields = {{ Ey_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::EZ)
            fields = {{ Ez_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::HX)
            fields = {{ Hx_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::HY)
            fields = {{ Hy_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::HZ)
            fields = {{ Hz_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::PX)
            fields = lorPx_;
        else if(IP.dtcType_[dd] == DTCTYPE::PY)
            fields = lorPy_;
        else if(IP.dtcType_[dd] == DTCTYPE::PZ)
            fields = lorPz_;
        else if(IP.dtcType_[dd] == DTCTYPE::MX)
            fields = lorMx_;
        else if(IP.dtcType_[dd] == DTCTYPE::MY)
            fields = lorMy_;
        else if(IP.dtcType_[dd] == DTCTYPE::MZ)
            fields = lorMz_;
        else if(IP.dtcType_[dd] == DTCTYPE::EPOW && Hz_ && Ez_)
            fields = {{ Ex_, Ey_, Ez_}};
        else if(IP.dtcType_[dd] == DTCTYPE::HPOW && Hz_ && Ez_)
            fields = {{ Hx_, Hy_, Hz_}};
        else if(IP.dtcType_[dd] == DTCTYPE::EPOW && Hz_)
            fields = {{ Ex_,Ey_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::HPOW && Hz_)
            fields = {{ Hz_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::EPOW && Ez_)
            fields = {{ Ez_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::HPOW && Ez_)
            fields = {{ Hx_,Hy_ }};
        else
            throw std::logic_error("DTC TYPE IS NOT DEFINED");
        coustructDTC(IP.dtcClass_[dd], fields, IP.dtcSI_[dd], IP.dtcLoc_[dd], IP.dtcSz_[dd], IP.dtcName_[dd], IP.dtcOutBMPFxnType_[dd], IP.dtcOutBMPOutType_[dd], IP.dtcType_[dd], IP.dtcFreqList_[dd], IP.dtcTimeInt_[dd], IP.a_, IP.I0_, IP.tMax_);
    }
    // Initialze all detectors to time 0
    for(auto& dtc : dtcArr_)
        dtc->output(tcur_);
    for(auto& dtc : dtcFreqArr_)
        dtc->output(tcur_);
    for(auto& flux : fluxArr_)
        flux->fieldIn(tcur_);

    E_incd_.push_back(0.0);
    E_pl_incd_.push_back(0.0);
    H_incd_.push_back(0.0);
    H_mn_incd_.push_back(0.0);
}

// Same as the real version but uses the complex versions of everything
parallelFDTDFieldCplx::parallelFDTDFieldCplx(parallelProgramInputs &IP, std::shared_ptr<mpiInterface> gridComm) :
    parallelFDTDFieldBase<cplx>(IP, gridComm)
{
    // If there is an Hz field set up all TE mode functions otherwise set them to do nothing
    if(Hz_)
    {
        // initialize the PMLs
        if(magMatInPML_)
        {
            HzPML_   = std::make_shared<parallelCPMLCplx>(gridComm_, weights_, Bz_, Ex_, Ey_, POLARIZATION::HZ, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Hz_, objArr_);
        }
        else
        {
            HzPML_   = std::make_shared<parallelCPMLCplx>(gridComm_, weights_, Hz_, Ex_, Ey_, POLARIZATION::HZ, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Hz_, objArr_);
        }
        if(dielectricMatInPML_)
        {
            ExPML_   = std::make_shared<parallelCPMLCplx>(gridComm_, weights_, Dx_, Hy_, Hz_, POLARIZATION::EX, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Ex_, objArr_);
            EyPML_   = std::make_shared<parallelCPMLCplx>(gridComm_, weights_, Dy_, Hz_, Hx_, POLARIZATION::EY, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Ex_, objArr_);
        }
        else
        {
            ExPML_   = std::make_shared<parallelCPMLCplx>(gridComm_, weights_, Ex_, Hy_, Hz_, POLARIZATION::EX, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Ey_, objArr_);
            EyPML_   = std::make_shared<parallelCPMLCplx>(gridComm_, weights_, Ey_, Hz_, Hx_, POLARIZATION::EY, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Ey_, objArr_);
        }
        initializeList(phys_Ex_, ExPML_,  true, std::array<int,3>( {{ 0, -1,  0 }} ), std::array<int,3>( {{n_vec_[0]  , n_vec_[1]+1, n_vec_[2]+1 }} ), d_[1], upEx_, upDx_);
        initializeList(phys_Ey_, EyPML_,  true, std::array<int,3>( {{ 0,  0, -1 }} ), std::array<int,3>( {{n_vec_[0]+1, n_vec_[1]  , n_vec_[2]+1 }} ), d_[0], upEy_, upDy_);
        initializeList(phys_Hz_, HzPML_, false, std::array<int,3>( {{ 1,  0,  0 }} ), std::array<int,3>( {{n_vec_[0]  , n_vec_[1]  , n_vec_[2]   }} ), d_[0], upHz_, upBz_);

        if(!Ez_)
        {
            upExFxn_ = &FDTDCompUpdateFxnCplx::OneCompCurlK;
            upEyFxn_ = &FDTDCompUpdateFxnCplx::OneCompCurlJ;
        }
        else
        {
            upExFxn_ = &FDTDCompUpdateFxnCplx::TwoCompCurl;
            upEyFxn_ = &FDTDCompUpdateFxnCplx::TwoCompCurl;
        }

        upHzFxn_ = &FDTDCompUpdateFxnCplx::TwoCompCurl;

        updateExPML_ = [](pml_ptr pml){pml->updateGrid();};
        updateEyPML_ = [](pml_ptr pml){pml->updateGrid();};
        updateHzPML_ = [](pml_ptr pml){pml->updateGrid();};

        upLorPxFxn_ = &FDTDCompUpdateFxnCplx::UpdateLorPol;
        upLorPyFxn_ = &FDTDCompUpdateFxnCplx::UpdateLorPol;
        upLorMzFxn_ = &FDTDCompUpdateFxnCplx::UpdateLorMag;

        D2ExFxn_ = &FDTDCompUpdateFxnCplx::DtoE;
        D2EyFxn_ = &FDTDCompUpdateFxnCplx::DtoE;

        B2HzFxn_ = &FDTDCompUpdateFxnCplx::BtoH;

        if(IP.periodic_ && gridComm_->size() > 1)
        {
            pbcEx_ = &FDTDCompUpdateFxnCplx::applyPBC;
            pbcEy_ = &FDTDCompUpdateFxnCplx::applyPBC;
            pbcHz_ = &FDTDCompUpdateFxnCplx::applyPBC;

            yExPBC_ = ln_vec_[1]+1;
            yEyPBC_ = ln_vec_[1]+1;
            yHzPBC_ = ln_vec_[1]+1;

            if(gridComm_->size()-1 == gridComm_->rank())
            {
                yHzPBC_ -= 1;
                yEyPBC_ -= 1;
            }
        }
        else if(IP.periodic_)
        {
            pbcEx_ = &FDTDCompUpdateFxnCplx::applyPBC1Proc;
            pbcEy_ = &FDTDCompUpdateFxnCplx::applyPBC1Proc;
            pbcHz_ = &FDTDCompUpdateFxnCplx::applyPBC1Proc;

            yHzPBC_ = ln_vec_[1];
            yExPBC_ = ln_vec_[1]+1;
            yEyPBC_ = ln_vec_[1];
        }
        else
        {
            pbcEx_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
            pbcEy_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
            pbcHz_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
        }

        // Transfer fields should do nothing if only on one processor
        if(gridComm_->size() > 1)
        {
            transferEx_ = [=](){Ex_->transferDat();};
            transferEy_ = [=](){Ey_->transferDat();};
            transferHz_ = [=](){Hz_->transferDat();};
        }
        else
        {
            transferEx_ = [](){return;};
            transferEy_ = [](){return;};
            transferHz_ = [](){return;};
        }
    }
    else
    {
        upHzFxn_ = [](std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr){return;};

        upExFxn_ = [](std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr){return;};;
        upEyFxn_ = [](std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr){return;};;

        updateExPML_ = [](pml_ptr pml){return;};
        updateEyPML_ = [](pml_ptr pml){return;};
        updateHzPML_ = [](pml_ptr pml){return;};

        upLorPxFxn_ = []( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, cplx*, std::shared_ptr<Obj> ){return;};
        upLorPyFxn_ = []( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, cplx*, std::shared_ptr<Obj> ){return;};
        upLorMzFxn_ = []( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, cplx*, std::shared_ptr<Obj> ){return;};

        D2ExFxn_ = []( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj> ){return;};
        D2EyFxn_ = []( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj> ){return;};
        B2HzFxn_ = []( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj> ){return;};

        pbcEx_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
        pbcEy_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
        pbcHz_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};

        transferEx_ = [](){return;};
        transferEy_ = [](){return;};
        transferHz_ = [](){return;};
    }
    // If there is an Ez field set up all TM functions otherwise set them to do nothing
    if(Ez_)
    {
        // Initialize PMLs

        if(magMatInPML_)
        {
            HxPML_   = std::make_shared<parallelCPMLCplx>(gridComm_, weights_, Bx_, Ey_, Ez_, POLARIZATION::HX, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Hx_, objArr_);
            HyPML_   = std::make_shared<parallelCPMLCplx>(gridComm_, weights_, By_, Ez_, Ex_, POLARIZATION::HY, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Hx_, objArr_);
        }
        else
        {
            HxPML_   = std::make_shared<parallelCPMLCplx>(gridComm_, weights_, Hx_, Ey_, Ez_, POLARIZATION::HX, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Hy_, objArr_);
            HyPML_   = std::make_shared<parallelCPMLCplx>(gridComm_, weights_, Hy_, Ez_, Ex_, POLARIZATION::HY, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Hy_, objArr_);
        }
        if(dielectricMatInPML_)
        {
            EzPML_   = std::make_shared<parallelCPMLCplx>(gridComm_, weights_, Dz_, Hx_, Hy_, POLARIZATION::EZ, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Ez_, objArr_);
        }
        else
        {
            EzPML_   = std::make_shared<parallelCPMLCplx>(gridComm_, weights_, Ez_, Hx_, Hy_, POLARIZATION::EZ, pmlThickness_, IP.pmlM_, IP.pmlMa_, IP.pmlAMax_, d_, dt_, phys_Ez_, objArr_);
        }

        initializeList(phys_Ez_, EzPML_,  true, std::array<int,3>( {{-1,  0,  0 }} ), std::array<int,3>( {{n_vec_[0]+1, n_vec_[1]+1, n_vec_[2]   }} ), d_[0], upEz_, upDz_);
        initializeList(phys_Hx_, HxPML_, false, std::array<int,3>( {{ 0,  1,  0 }} ), std::array<int,3>( {{n_vec_[0]+1, n_vec_[1]  , n_vec_[2]   }} ), d_[1], upHx_, upBx_);
        initializeList(phys_Hy_, HyPML_, false, std::array<int,3>( {{ 0,  0,  1 }} ), std::array<int,3>( {{n_vec_[0]  , n_vec_[1]+1, n_vec_[2]+1 }} ), d_[0], upHy_, upBy_);

        if(!Hz_)
        {
            upHxFxn_ = &FDTDCompUpdateFxnCplx::OneCompCurlK;
            upHyFxn_ = &FDTDCompUpdateFxnCplx::OneCompCurlJ;
        }
        else
        {
            upHxFxn_ = &FDTDCompUpdateFxnCplx::TwoCompCurl;
            upHyFxn_ = &FDTDCompUpdateFxnCplx::TwoCompCurl;
        }

        upEzFxn_ = &FDTDCompUpdateFxnCplx::TwoCompCurl;

        upLorMxFxn_ = &FDTDCompUpdateFxnCplx::UpdateLorMag;
        upLorMyFxn_ = &FDTDCompUpdateFxnCplx::UpdateLorMag;

        upLorPzFxn_ = &FDTDCompUpdateFxnCplx::UpdateLorPol;

        B2HxFxn_ = &FDTDCompUpdateFxnCplx::BtoH;
        B2HyFxn_ = &FDTDCompUpdateFxnCplx::BtoH;

        D2EzFxn_ = &FDTDCompUpdateFxnCplx::DtoE;

        updateHxPML_ = [](pml_ptr pml){pml->updateGrid();};
        updateHyPML_ = [](pml_ptr pml){pml->updateGrid();};
        updateEzPML_ = [](pml_ptr pml){pml->updateGrid();};

        if(IP.periodic_ && gridComm_->size() > 1)
        {
            pbcHx_ = &FDTDCompUpdateFxnCplx::applyPBC;
            pbcHy_ = &FDTDCompUpdateFxnCplx::applyPBC;
            pbcEz_ = &FDTDCompUpdateFxnCplx::applyPBC;

            yHxPBC_ = ln_vec_[1]+1;
            yHyPBC_ = ln_vec_[1]+1;
            yEzPBC_ = ln_vec_[1]+1;

            if(gridComm_->size()-1 == gridComm_->rank())
                yHxPBC_ = ln_vec_[1];
        }
        else if(IP.periodic_)
        {
            pbcHx_ = &FDTDCompUpdateFxnCplx::applyPBC1Proc;
            pbcHy_ = &FDTDCompUpdateFxnCplx::applyPBC1Proc;
            pbcEz_ = &FDTDCompUpdateFxnCplx::applyPBC1Proc;
            yHxPBC_ = ln_vec_[1];
            yHyPBC_ = ln_vec_[1]+1;
            yEzPBC_ = ln_vec_[1]+1;
        }
        else
        {
            pbcHx_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
            pbcHy_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
            pbcEz_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
        }

        // Transfer fields should do nothing if only on one processor
        if(gridComm_->size() > 1)
        {
            transferHx_ = [=](){Hx_->transferDat();};
            transferHy_ = [=](){Hy_->transferDat();};
            transferEz_ = [=](){Ez_->transferDat();};

        }
        else
        {
            transferEz_ = [](){return;};
            transferHx_ = [](){return;};
            transferHy_ = [](){return;};
        }
    }
    else
    {
        upHxFxn_ = [](std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr){return;};;
        upHyFxn_ = [](std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr){return;};

        upEzFxn_ = [](std::array<int,8>&, std::array<double,2>&, pgrid_ptr, pgrid_ptr, pgrid_ptr){return;};

        upLorMxFxn_ = []( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, cplx*, std::shared_ptr<Obj> ){return;};
        upLorMyFxn_ = []( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, cplx*, std::shared_ptr<Obj> ){return;};
        upLorPzFxn_ = []( std::array<int,8>&, pgrid_ptr, std::vector<pgrid_ptr>&, std::vector<pgrid_ptr>&, cplx*, std::shared_ptr<Obj> ){return;};

        B2HxFxn_ = []( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj> ){return;};
        B2HyFxn_ = []( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj> ){return;};
        D2EzFxn_ = []( std::array<int,8>&, pgrid_ptr, pgrid_ptr, std::vector<pgrid_ptr>&, std::shared_ptr<Obj> ){return;};

        updateHxPML_ = [](pml_ptr pml){return;};
        updateHyPML_ = [](pml_ptr pml){return;};
        updateEzPML_ = [](pml_ptr pml){return;};

        pbcHx_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
        pbcHy_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};
        pbcEz_ = [](pgrid_ptr, std::array<double,3>&, int, int, int, int, int, int, int, double&, double&, double&){return;};

        transferEz_ = [](){return;};
        transferHx_ = [](){return;};
        transferHy_ = [](){return;};
    }
    // Construct all soft sources
    for(int ss = 0; ss < IP.srcPol_.size(); ss++)
    {
        // Make the pulse (including all pulses to be used)
        std::vector<std::shared_ptr<PulseBase>> pul;
        for(int pp = 0; pp < IP.srcPulShape_[ss].size(); pp ++)
        {
            if(IP.srcPulShape_[ss][pp] == PLSSHAPE::CONTINUOUS)
                pul.push_back(std::make_shared<PulseCont>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp], dt_));
            else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::BH)
                pul.push_back(std::make_shared<PulseBH>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp], dt_));
            else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::RECT)
                pul.push_back(std::make_shared<PulseRect>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp], dt_));
            else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::GAUSSIAN)
                pul.push_back(std::make_shared<PulseGauss>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp], dt_));
            else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::RICKER)
                pul.push_back(std::make_shared<PulseRicker>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp], dt_));
            else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::RAMP_CONT)
                pul.push_back(std::make_shared<PulseRampCont>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp], dt_));
        }
        //  Make the source act on any of the fields
        if(IP.srcPol_[ss] == POLARIZATION::EX)
        {
            if(int( round(IP.srcPhi_[ss]) ) % 90 == 0)
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalCplx>(gridComm_, pul, Ex_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            else
            {
                srcArr_.push_back( std::make_shared<parallelSourceObliqueCplx >(gridComm_, pul, Ex_, POLARIZATION::EX, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceObliqueCplx >(gridComm_, pul, Ey_, POLARIZATION::EY, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
            }
        }
        else if(IP.srcPol_[ss] == POLARIZATION::EY)
        {
            if(int( round(IP.srcPhi_[ss]) ) % 90 == 0)
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalCplx>(gridComm_, pul, Ey_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            else
            {
                srcArr_.push_back( std::make_shared<parallelSourceObliqueCplx >(gridComm_, pul, Ex_, POLARIZATION::EX, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceObliqueCplx >(gridComm_, pul, Ey_, POLARIZATION::EY, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
            }
        }
        else if(IP.srcPol_[ss] == POLARIZATION::EZ)
        {
            if(int( round(IP.srcPhi_[ss]) ) % 90 == 0)
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalCplx>(gridComm_, pul, Ez_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            else
            {
                srcArr_.push_back( std::make_shared<parallelSourceObliqueCplx >(gridComm_, pul, Ez_, POLARIZATION::EZ, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
            }
        }
        else if(IP.srcPol_[ss] == POLARIZATION::HX)
        {
            if(int( round(IP.srcPhi_[ss]) ) % 90 == 0)
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalCplx>(gridComm_, pul, Hx_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            else
            {
                srcArr_.push_back( std::make_shared<parallelSourceObliqueCplx >(gridComm_, pul, Hy_, POLARIZATION::HY, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceObliqueCplx >(gridComm_, pul, Hx_, POLARIZATION::HX, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
            }
        }
        else if(IP.srcPol_[ss] == POLARIZATION::HY)
        {
            if(int( round(IP.srcPhi_[ss]) ) % 90 == 0)
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalCplx>(gridComm_, pul, Hy_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            else
            {
                srcArr_.push_back( std::make_shared<parallelSourceObliqueCplx >(gridComm_, pul, Hx_, POLARIZATION::HX, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceObliqueCplx >(gridComm_, pul, Hy_, POLARIZATION::HY, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
            }
        }
        else if(IP.srcPol_[ss] == POLARIZATION::HZ)
        {
            if(int( round(IP.srcPhi_[ss]) ) % 90 == 0)
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalCplx>(gridComm_, pul, Hz_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            else
            {
                srcArr_.push_back( std::make_shared<parallelSourceObliqueCplx >(gridComm_, pul, Hz_, POLARIZATION::HZ, dt_, IP.srcLoc_[ss], IP.srcSz_[ss], IP.srcPhi_[ss], IP.srcTheta_[ss] ) );
            }
        }
        else
        {
            double axRat = IP.srcEllipticalKratio_[ss];
            double psi = IP.srcPsi_[ss];
            double psiPrefactCalc = psi;
            double alphaOff = 0.0;
            double prefactor_k_ = 1.0;
            double prefactor_j_ = 1.0;
            double c = pow(axRat, 2.0);

            // phi/psi control the light polarization angle
            psiPrefactCalc = 0.5 * asin( sqrt( ( pow(cos(2.0*psi),2.0)*4.0*c + pow( (1.0+c)*sin(2.0*psi), 2.0) ) / pow(1.0+c, 2.0) ) );
            alphaOff = acos( ( (c - 1.0)*sin(2.0*psi) ) / sqrt( pow(cos(2.0*psi),2.0)*4.0*c + pow( (1.0+c)*sin(2.0*psi), 2.0) ) );
            if(std::abs( std::tan(psi) ) > 1)
                psiPrefactCalc = M_PI/2.0 - psiPrefactCalc;
            if(IP.srcPol_[ss] == POLARIZATION::R)
                alphaOff *= -1.0;

            if( isamin_(IP.srcSz_[ss].size(), IP.srcSz_[ss].data(), 1)-1 == 0 )
            {
                prefactor_j_ *= -1.0 * cos(psiPrefactCalc);
                prefactor_k_ *= -1.0 * sin(psiPrefactCalc);
            }
            else if( isamin_(IP.srcSz_[ss].size(), IP.srcSz_[ss].data(), 1)-1 == 1 )
            {
                prefactor_j_ *=  -1.0 * cos(psiPrefactCalc);
                prefactor_k_ *=         sin(psiPrefactCalc);
            }
            else if( isamin_(IP.srcSz_[ss].size(), IP.srcSz_[ss].data(), 1)-1 == 2 )
            {
                prefactor_j_ *= -1.0*sin(psiPrefactCalc);
                prefactor_k_ *=      cos(psiPrefactCalc);
            }

            for(auto pulse : pul)
                pulse->modE0(prefactor_j_);
            std::vector<std::shared_ptr<PulseBase>> phaseOffPul;
            for(int pp = 0; pp < IP.srcPulShape_[ss].size(); pp ++)
            {
                if(IP.srcPulShape_[ss][pp] == PLSSHAPE::CONTINUOUS)
                    phaseOffPul.push_back(std::make_shared<PulseCont>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp]*prefactor_k_, dt_, alphaOff));
                else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::BH)
                    phaseOffPul.push_back(std::make_shared<PulseBH>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp]*prefactor_k_, dt_, alphaOff));
                else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::RECT)
                    phaseOffPul.push_back(std::make_shared<PulseRect>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp]*prefactor_k_, dt_, alphaOff));
                else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::GAUSSIAN)
                    phaseOffPul.push_back(std::make_shared<PulseGauss>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp]*prefactor_k_, dt_, alphaOff));
                else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::RICKER)
                    phaseOffPul.push_back(std::make_shared<PulseRicker>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp]*prefactor_k_, dt_, alphaOff));
                else if(IP.srcPulShape_[ss][pp] == PLSSHAPE::RAMP_CONT)
                    phaseOffPul.push_back(std::make_shared<PulseRampCont>(IP.srcFxn_[ss][pp], IP.srcEmax_[ss][pp]*prefactor_k_, dt_, alphaOff));
            }
            if(isamin_(IP.srcSz_[ss].size(), IP.srcSz_[ss].data(), 1)-1 == 0 )
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalCplx>(gridComm_, pul, Ey_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceNormalCplx>(gridComm_, phaseOffPul, Ez_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            if(isamin_(IP.srcSz_[ss].size(), IP.srcSz_[ss].data(), 1)-1 == 1 )
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalCplx>(gridComm_, pul, Ex_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceNormalCplx>(gridComm_, phaseOffPul, Ez_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
            if(isamin_(IP.srcSz_[ss].size(), IP.srcSz_[ss].data(), 1)-1 == 2 )
            {
                srcArr_.push_back( std::make_shared<parallelSourceNormalCplx>(gridComm_, pul, Ex_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
                srcArr_.push_back( std::make_shared<parallelSourceNormalCplx>(gridComm_, phaseOffPul, Ey_, dt_, IP.srcLoc_[ss], IP.srcSz_[ss] ) );
            }
        }
    }

    // Construct all TFSF surfaces
    for(int tt = 0; tt < IP.tfsfSize_.size(); tt++)
    {
        if(IP.tfsfSize_[tt][0] != 0.0 || IP.tfsfSize_[tt][1] != 0.0)
        {
        // Make the pulse (including all pulses to be used)
            std::vector<std::shared_ptr<PulseBase>> pul;
            for(int pp = 0; pp < IP.tfsfPulShape_[tt].size(); pp ++)
            {
                if(IP.tfsfPulShape_[tt][pp] == PLSSHAPE::CONTINUOUS)
                    pul.push_back(std::make_shared<PulseCont>(IP.tfsfPulFxn_[tt][pp], IP.tfsfEmax_[tt][pp], dt_));
                else if(IP.tfsfPulShape_[tt][pp] == PLSSHAPE::BH)
                    pul.push_back(std::make_shared<PulseBH>(IP.tfsfPulFxn_[tt][pp], IP.tfsfEmax_[tt][pp], dt_));
                else if(IP.tfsfPulShape_[tt][pp] == PLSSHAPE::RECT)
                    pul.push_back(std::make_shared<PulseRect>(IP.tfsfPulFxn_[tt][pp], IP.tfsfEmax_[tt][pp], dt_));
                else if(IP.tfsfPulShape_[tt][pp] == PLSSHAPE::GAUSSIAN)
                    pul.push_back(std::make_shared<PulseGauss>(IP.tfsfPulFxn_[tt][pp], IP.tfsfEmax_[tt][pp], dt_));
                else if(IP.tfsfPulShape_[tt][pp] == PLSSHAPE::RICKER)
                    pul.push_back(std::make_shared<PulseRicker>(IP.tfsfPulFxn_[tt][pp], IP.tfsfEmax_[tt][pp], dt_));
                else if(IP.tfsfPulShape_[tt][pp] == PLSSHAPE::RAMP_CONT)
                    pul.push_back(std::make_shared<PulseRampCont>(IP.tfsfPulFxn_[tt][pp], IP.tfsfEmax_[tt][pp], dt_));
            }
            // TFSF will determine the correct polarization
            tfsfArr_.push_back(std::make_shared<parallelTFSFCplx>(gridComm_, IP.tfsfLoc_[tt], IP.tfsfSize_[tt], IP.tfsfTheta_[tt], IP.tfsfPhi_[tt], IP.tfsfPsi_[tt], IP.tfsfCircPol_[tt], IP.tfsfEllipticalKratio_[tt], d_[0], dt_, pul, Ex_, Ey_, Ez_, Hx_, Hy_, Hz_) );
        }
    }
    // Polarization matters here since the z field always forms the continous box for the spatial offset (TE uses H, TM uses E)
    if(Hz_ && Ez_ && IP.fluxName_.size() > 0)
    {
        // Flux needs to know about how to calculate incident light.
        DIRECTION propDir;
        if(tfsfArr_.size() > 0 && (tfsfArr_.back()->theta() == 0.0 || tfsfArr_.back()->theta() == M_PI ) )
            propDir = DIRECTION::Z;
        else if(tfsfArr_.size() > 0 && (tfsfArr_.back()->quadrant() == 2 || tfsfArr_.back()->quadrant() == 4 ) )
            propDir = DIRECTION::X;
        else if(tfsfArr_.size() > 0 && (tfsfArr_.back()->quadrant() == 1 || tfsfArr_.back()->quadrant() == 3 ) )
            propDir = DIRECTION::Y;
        else if(srcArr_.size() > 0 && (srcArr_.back()->sz()[0] >= srcArr_.back()->sz()[1] ) && (srcArr_.back()->sz()[2] >= srcArr_.back()->sz()[1] ) )
            propDir = DIRECTION::Y;
        else if(srcArr_.size() > 0 && (srcArr_.back()->sz()[0] <= srcArr_.back()->sz()[1] ) && (srcArr_.back()->sz()[0] <= srcArr_.back()->sz()[2] ) )
            propDir = DIRECTION::X;
        else if(srcArr_.size() > 0)
            propDir = DIRECTION::Z;
        else
            throw std::logic_error("Constructing a flux with no source, it will be 0.");

        double theta = 0; double phi = 0; double psi = 0; double alpha = 0;
        if(tfsfArr_.size() > 0)
        {
            theta = std::atan( std::abs( std::tan( tfsfArr_.back()->theta() ) ) );
            phi   = std::atan( std::abs( std::tan( tfsfArr_.back()->phiPreFact() ) ) );
            psi   = std::atan( std::abs( std::tan( tfsfArr_.back()->psiPreFact() ) ) );
            alpha = tfsfArr_.back()->alpha();
        }
        for(int ff = 0; ff < IP.fluxLoc_.size(); ff ++)
        {
            // Flux set by wavelength or frequency?
            fluxArr_.push_back(std::make_shared<parallelFluxDTCCplx>(gridComm_, IP.fluxName_[ff], IP.fluxWeight_[ff], Ex_, Ey_, Ez_, Hx_, Hy_, Hz_, IP.fluxLoc_[ff], IP.fluxSz_[ff], IP.fluxCrossSec_[ff], IP.fluxSave_[ff], IP.fluxLoad_[ff], IP.fluxTimeInt_[ff], IP.fluxFreqList_[ff], propDir, d_, dt_, theta, phi, psi, alpha, IP.fluxIncdFieldsFilename_[ff], IP.fluxSI_[ff], IP.I0_, IP.a_) );
        }
    }
    // Construct all DTC based on types (all it changes is the list of fields it passes)
    for(int dd = 0; dd < IP.dtcType_.size(); dd++)
    {
        std::vector<pgrid_ptr> fields;
        // Fill the fields vector with the appropriate field values
        if(IP.dtcType_[dd] == DTCTYPE::EX)
            fields = {{ Ex_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::EY)
            fields = {{ Ey_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::EZ)
            fields = {{ Ez_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::HX)
            fields = {{ Hx_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::HY)
            fields = {{ Hy_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::HZ)
            fields = {{ Hz_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::PX)
            fields = lorPx_;
        else if(IP.dtcType_[dd] == DTCTYPE::PY)
            fields = lorPy_;
        else if(IP.dtcType_[dd] == DTCTYPE::PZ)
            fields = lorPz_;
        else if(IP.dtcType_[dd] == DTCTYPE::MX)
            fields = lorMx_;
        else if(IP.dtcType_[dd] == DTCTYPE::MY)
            fields = lorMy_;
        else if(IP.dtcType_[dd] == DTCTYPE::MZ)
            fields = lorMz_;
        else if(IP.dtcType_[dd] == DTCTYPE::EPOW && Hz_ && Ez_)
            fields = {{ Ex_, Ey_, Ez_}};
        else if(IP.dtcType_[dd] == DTCTYPE::HPOW && Hz_ && Ez_)
            fields = {{ Hx_, Hy_, Hz_}};
        else if(IP.dtcType_[dd] == DTCTYPE::EPOW && Hz_)
            fields = {{ Ex_,Ey_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::HPOW && Hz_)
            fields = {{ Hz_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::EPOW && Ez_)
            fields = {{ Ez_ }};
        else if(IP.dtcType_[dd] == DTCTYPE::HPOW && Ez_)
            fields = {{ Hx_,Hy_ }};
        else
            throw std::logic_error("DTC TYPE IS NOT DEFINED");
        coustructDTC(IP.dtcClass_[dd], fields, IP.dtcSI_[dd], IP.dtcLoc_[dd], IP.dtcSz_[dd], IP.dtcName_[dd], IP.dtcOutBMPFxnType_[dd], IP.dtcOutBMPOutType_[dd], IP.dtcType_[dd], IP.dtcFreqList_[dd], IP.dtcTimeInt_[dd], IP.a_, IP.I0_, IP.tMax_);
    }
    // Initialze all detectors to time 0
    for(auto& dtc : dtcArr_)
        dtc->output(tcur_);
    for(auto& dtc : dtcFreqArr_)
        dtc->output(tcur_);
    for(auto& flux : fluxArr_)
        flux->fieldIn(tcur_);

    E_incd_.push_back(0.0);
    E_pl_incd_.push_back(0.0);
    H_incd_.push_back(0.0);
    H_mn_incd_.push_back(0.0);
}

void parallelFDTDFieldReal::coustructDTC(DTCCLASS c, std::vector<pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, GRIDOUTFXN fxn, GRIDOUTTYPE txtType, DTCTYPE type, std::vector<double> freqList, double timeInterval, double a, double I0, double t_max)
{
    if(c == DTCCLASS::BIN)
        dtcArr_.push_back( std::make_shared<parallelDetectorBINReal>( grid, SI, loc, sz, out_name, type, timeInterval, a, I0, dt_) );
    else if(c == DTCCLASS::TXT)
        dtcArr_.push_back( std::make_shared<parallelDetectorTXTReal>(  grid, SI, loc, sz, out_name, type, timeInterval, a, I0, dt_ ) );
    else if(c == DTCCLASS::COUT)
        dtcArr_.push_back( std::make_shared<parallelDetectorCOUTReal>(grid, SI, loc, sz, out_name, type, timeInterval, a, I0, dt_) );
    else if(c == DTCCLASS::FREQ)
        dtcFreqArr_.push_back(std::make_shared<parallelDetectorFREQReal>(out_name, grid, loc, sz, type, static_cast<int>(std::floor(timeInterval/dt_+0.5) ), freqList, d_, dt_, SI, I0, a) );
    else
        throw std::logic_error("The detector class is undefined.");
}

void parallelFDTDFieldCplx::coustructDTC(DTCCLASS c, std::vector<pgrid_ptr> grid, bool SI, std::array<int,3> loc, std::array<int,3> sz, std::string out_name, GRIDOUTFXN fxn, GRIDOUTTYPE txtType, DTCTYPE type, std::vector<double> freqList, double timeInterval, double a, double I0, double t_max)
{
    if(c == DTCCLASS::BIN)
        dtcArr_.push_back( std::make_shared<parallelDetectorBINCplx>( grid, SI, loc, sz, out_name, type, timeInterval, a, I0, dt_) );
    else if(c == DTCCLASS::TXT)
        dtcArr_.push_back( std::make_shared<parallelDetectorTXTCplx>( grid, SI, loc, sz, out_name, type, timeInterval, a, I0, dt_ ) );
    else if(c == DTCCLASS::COUT)
        dtcArr_.push_back( std::make_shared<parallelDetectorCOUTCplx>(grid, SI, loc, sz, out_name, type, timeInterval, a, I0, dt_) );
    else if(c == DTCCLASS::FREQ)
        dtcFreqArr_.push_back(std::make_shared<parallelDetectorFREQCplx>(out_name, grid, loc, sz, type, static_cast<int>(std::floor(timeInterval/dt_+0.5) ),freqList, d_, dt_, SI, I0, a) );
    else
        throw std::logic_error("The detector class is undefined.");
}