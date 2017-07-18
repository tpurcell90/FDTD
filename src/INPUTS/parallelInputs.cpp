#include <INPUTS/parallelInputs.hpp>
parallelProgramInputs::parallelProgramInputs(boost::property_tree::ptree IP,std::string fn) :
    periodic_(IP.get<bool>("CompCell.PBC", false) ),
    pol_(string2pol(IP.get<std::string>("CompCell.pol") ) ),
    filename_(fn),
    res_(IP.get<int>("CompCell.res") ),
    size_( as_ptArr<double>( IP, "CompCell.size") ),
    courant_(IP.get<double>("CompCell.courant", 0.5) ),
    a_(IP.get<double>("CompCell.a",1e-7) ),
    tMax_(IP.get<double>("CompCell.tLim") ),
    I0_(IP.get<double>("CompCell.I0", a_ * EPS0 * SPEED_OF_LIGHT) ),
    cplxFields_(false),
    saveFreqField_(false),
    // pmlThickness_( as_ptArr<double>(IP, "PML.thickness") ),
    k_point_( as_ptArr<double>( IP, "k-point") ),

    pmlAMax_( IP.get<double>("PML.aMax",0.25) ),
    pmlMa_( IP.get<double>("PML.ma", 1.0) ),
    pmlM_( IP.get<double>("PML.m", 3.0) ),

    inputMapSlicesX_(as_vector<double>(IP, "CompCell.InputMaps_x") ),
    inputMapSlicesY_(as_vector<double>(IP, "CompCell.InputMaps_y") ),
    inputMapSlicesZ_(as_vector<double>(IP, "CompCell.InputMaps_z") ),

    srcPol_( std::vector<POLARIZATION>(IP.get_child("SourceList").size(), POLARIZATION::EX) ),
    srcFxn_( std::vector<std::vector<std::vector<double>>>(IP.get_child("SourceList").size(), std::vector<std::vector<double>>() ) ),
    srcLoc_( std::vector<std::array<int,3>>(IP.get_child("SourceList").size() ) ),
    srcSz_ ( std::vector<std::array<int,3>>(IP.get_child("SourceList").size() ) ),
    srcPhi_( std::vector<double>(IP.get_child("SourceList").size(), 0.0) ),
    srcTheta_( std::vector<double>(IP.get_child("SourceList").size(), 0.0) ),
    srcPulShape_( std::vector<std::vector<PLSSHAPE>>(IP.get_child("SourceList").size(), std::vector<PLSSHAPE>() ) ),
    srcEmax_( std::vector<std::vector<double>>(IP.get_child("SourceList").size(), std::vector<double>() ) ),
    srcEllipticalKratio_( std::vector<double>(IP.get_child("SourceList").size(), 0.0) ),
    srcPsi_( std::vector<double>(IP.get_child("SourceList").size(), 0.0) ),

    tfsfLoc_( std::vector<std::array<int,3>>(IP.get_child("TFSF").size() ) ),
    tfsfSize_( std::vector<std::array<int,3>>(IP.get_child("TFSF").size() ) ),
    tfsfTheta_( std::vector<double>(IP.get_child("TFSF").size(), 0.0) ),
    tfsfPhi_( std::vector<double>(IP.get_child("TFSF").size(), 0.0) ),
    tfsfPsi_( std::vector<double>(IP.get_child("TFSF").size(), 0.0) ),
    tfsfPulFxn_( std::vector<std::vector<std::vector<double>>>(IP.get_child("TFSF").size(),  std::vector<std::vector<double>>() ) ),
    tfsfPulShape_( std::vector<std::vector<PLSSHAPE>>(IP.get_child("TFSF").size(),  std::vector<PLSSHAPE>() ) ),
    tfsfEmax_( std::vector<std::vector<double>>(IP.get_child("TFSF").size(),  std::vector<double>() ) ),
    tfsfCircPol_( std::vector<POLARIZATION>(IP.get_child("TFSF").size(), POLARIZATION::EX) ),
    tfsfEllipticalKratio_( std::vector<double>(IP.get_child("TFSF").size(), 0.0) ),

    dtcSI_( std::vector<bool>(IP.get_child("DetectorList").size(), false) ),
    dtcClass_( std::vector<DTCCLASS>(IP.get_child("DetectorList").size(), DTCCLASS::COUT) ),
    dtcLoc_( std::vector<std::array<int,3>>(IP.get_child("DetectorList").size() ) ),
    dtcSz_( std::vector<std::array<int,3>>(IP.get_child("DetectorList").size() ) ),
    dtcName_( std::vector<std::string>(IP.get_child("DetectorList").size(), std::string() ) ),
    dtcType_( std::vector<DTCTYPE>(IP.get_child("DetectorList").size(), DTCTYPE::EX ) ),
    dtcTimeInt_( std::vector<double>(IP.get_child("DetectorList").size(), 0.0) ),
    dtcNFreq_( std::vector<int>(IP.get_child("DetectorList").size(), -1) ),
    dtcfCen_( std::vector<double>(IP.get_child("DetectorList").size(), -1.0) ),
    dtcfWidth_( std::vector<double>(IP.get_child("DetectorList").size(), -1.0) ),
    dtcLamL_( std::vector<double>(IP.get_child("DetectorList").size(), -1.0) ),
    dtcLamR_( std::vector<double>(IP.get_child("DetectorList").size(), -1.0) ),

    fluxXOff_( std::vector<int>(IP.get_child("FluxList").size(), 0) ),
    fluxYOff_( std::vector<int>(IP.get_child("FluxList").size(), 0) ),
    fluxTimeInt_( std::vector<int>(IP.get_child("FluxList").size(), 0) ),
    fluxLoc_( std::vector<std::array<int,3>>(IP.get_child("FluxList").size() ) ),
    fluxSz_( std::vector<std::array<int,3>>(IP.get_child("FluxList").size() ) ),
    fluxWeight_( std::vector<double>(IP.get_child("FluxList").size(), 0.0) ),
    fluxName_( std::vector<std::string>(IP.get_child("FluxList").size(), std::string()) ),
    fluxFCen_( std::vector<double>(IP.get_child("FluxList").size(), 0.0) ),
    fluxFWidth_( std::vector<double>(IP.get_child("FluxList").size(), 0.0) ),
    fluxLamL_( std::vector<double>(IP.get_child("FluxList").size(), 0.0) ),
    fluxLamR_( std::vector<double>(IP.get_child("FluxList").size(), 0.0) ),
    fluxSI_( std::vector<bool>(IP.get_child("FluxList").size(), false) ),
    fluxCrossSec_( std::vector<bool>(IP.get_child("FluxList").size(), false) ),
    fluxSave_( std::vector<bool>(IP.get_child("FluxList").size(), false) ),
    fluxLoad_( std::vector<bool>(IP.get_child("FluxList").size(), false) ),
    fluxIncdFieldsFilename_( std::vector<std::string>(IP.get_child("FluxList").size(), std::string() ) ),
    fluxNFreq_( std::vector<int>( IP.get_child("FluxList").size(), 0 ) )
{
    std::array<double,3> pmlThickness = as_ptArr<double>( IP, "PML.thickness");
    for(int ii = 0; ii < 3; ++ii)
        pmlThickness_[ii] = find_pt(pmlThickness[ii]);
    if(periodic_)
        for(int kk = 0; kk < k_point_.size(); kk++)
            if(k_point_[kk] != 0)
                cplxFields_= true;


    int ii = 0;
    for (auto& iter : IP.get_child("SourceList") )
    {
        // PLSSHAPE prof = string2prof(iter.second.get<std::string>("profile"));
        srcPol_[ii]      = string2pol(iter.second.get<std::string>("pol"));

        boost::property_tree::ptree& PulseList = iter.second.get_child("PulseList");
        std::vector<std::vector<double>> pulFxn_(PulseList.size(), std::vector<double>());
        std::vector<PLSSHAPE> pulShape_(PulseList.size(), PLSSHAPE::GAUSSIAN);
        std::vector<double> pulEmax_(PulseList.size(), 0.0);
        int pp = 0;
        for(auto& pul : PulseList )
        {
            pulShape_[pp] = string2prof(pul.second.get<std::string>("profile"));
            pulEmax_[pp]  = pul.second.get<double>("Field_Intensity",1.0) * a_ * EPS0 * SPEED_OF_LIGHT / I0_;
            std::vector<double> fxn;
            switch (pulShape_[pp])
            {
                case PLSSHAPE::GAUSSIAN:
                    fxn.push_back(pul.second.get<double>("fcen"));
                    fxn.push_back(1.0/pul.second.get<double>("fwidth"));
                    fxn.push_back(pul.second.get<double>("cutoff"));
                    fxn.push_back(pul.second.get<double>("t_0", fxn[1]*fxn[2]));
                break;
                case PLSSHAPE::BH:
                    fxn.push_back(pul.second.get<double>("fcen"));
                    fxn.push_back(1.0/pul.second.get<double>("fwidth", 1.0/pul.second.get<double>("tau") ) );
                    fxn.push_back(pul.second.get<double>("t_0", fxn[1]*fxn[2]));
                    fxn.push_back(pul.second.get<double>("BH1"));
                    fxn.push_back(pul.second.get<double>("BH2"));
                    fxn.push_back(pul.second.get<double>("BH3"));
                    fxn.push_back(pul.second.get<double>("BH4"));
                break;
                case PLSSHAPE::RECT :
                    fxn.push_back(pul.second.get<double>("tau"));
                    fxn.push_back(pul.second.get<double>("t_0"));
                    fxn.push_back(pul.second.get<double>("n", 30) * 2.0);
                    fxn.push_back(pul.second.get<double>("fcen"));
                break;
                case PLSSHAPE::CONTINUOUS:
                    fxn.push_back(pul.second.get<double>("fcen"));
                break;
                case PLSSHAPE::RAMP_CONT:
                    fxn.push_back(pul.second.get<double>("fcen"));
                    fxn.push_back(pul.second.get<double>("ramp_val"));
                break;
                case PLSSHAPE::RICKER:
                    fxn.push_back(pul.second.get<double>("fcen"));
                    fxn.push_back(pul.second.get<double>("fwidth"));
                    fxn.push_back(pul.second.get<double>("cutoff"));
                break;
                default:
                    throw std::logic_error("This pulse shape is undefined in source ");
                break;
            }
            pulFxn_[pp] = fxn;
            ++pp;
        }
        srcPulShape_[ii] = pulShape_;
        srcEmax_[ii] = pulEmax_;
        srcFxn_[ii] = pulFxn_;

        srcEllipticalKratio_[ii] = iter.second.get<double>("ellpiticalKRat", 1.0);
        std::array<double,3> tempSz = as_ptArr<double>(iter.second, "size");
        for(int cc = 0; cc < 3; ++cc )
            srcSz_[ii][cc] = find_pt(tempSz[cc]) + 1;

        srcPhi_[ii] = iter.second.get<double>("phi", 90);
        if( (srcPol_[ii] == POLARIZATION::R || srcPol_[ii] == POLARIZATION::L) && isamin_(srcSz_[ii].size(), srcSz_[ii].data(), 1)-1 == 0 )
            srcPsi_[ii] = M_PI * iter.second.get<double>("psi", 0) / 180.0;
        else if( (srcPol_[ii] == POLARIZATION::R || srcPol_[ii] == POLARIZATION::L) && isamin_(srcSz_[ii].size(), srcSz_[ii].data(), 1)-1 == 1 )
            srcPsi_[ii] = M_PI * iter.second.get<double>("psi", 0) / 180.0;
        else if( (srcPol_[ii] == POLARIZATION::R || srcPol_[ii] == POLARIZATION::L) && isamin_(srcSz_[ii].size(), srcSz_[ii].data(), 1)-1 == 2 )
            srcPsi_[ii] = M_PI * iter.second.get<double>("psi", 0) / 180.0;
        int i = 0;
        for(auto& loc : as_ptArr<double>(iter.second, "loc") )
        {
            if(srcPhi_[ii] == 90 || srcPhi_[ii] == 180 || srcPhi_[ii] == 270 || srcPhi_[ii] == 0 )
            {
                srcLoc_[ii][i] = find_pt(loc + size_[i]/2.0 - tempSz[i]/2.0);
            }
            else
            {
                srcLoc_[ii][i] = find_pt(loc + size_[i]/2.0);
            }
            ++i;
        }
        ++ii;
    }

    ii = 0;
    for (auto& iter : IP.get_child("TFSF"))
    {
        int i = 0;
        for(auto& sz : as_ptArr<double>(iter.second, "size") )
        {
            tfsfSize_[ii][i] = find_pt(sz) + 1;
            ++i;
        }
        i = 0;
        for(auto& loc : as_ptArr<double>(iter.second, "loc") )
        {
            tfsfLoc_[ii][i] = find_pt(loc + size_[i]/2.0) - (tfsfSize_[ii][i] - (tfsfSize_[ii][i] % 2) ) / 2 ;
            ++i;
        }
        tfsfTheta_[ii] = M_PI*iter.second.get<double>("theta" , 90.0)/180.0;
        if(tfsfTheta_[ii] < 0 || tfsfTheta_[ii] > M_PI)
            throw std::logic_error("The theta value for the TFSF surfaces must be between 0 and 180 degrees.");
        tfsfPhi_[ii] = M_PI*iter.second.get<double>("phi" , 90.0)/180.0;
        tfsfPsi_[ii] = M_PI*iter.second.get<double>("psi" , 90.0)/180.0;
        tfsfCircPol_[ii] = string2pol( iter.second.get<std::string>("circPol", "Ex") );

        boost::property_tree::ptree& PulseList = iter.second.get_child("PulseList");
        std::vector<std::vector<double>> pulFxn_(PulseList.size(), std::vector<double>());
        std::vector<PLSSHAPE> pulShape_(PulseList.size(), PLSSHAPE::GAUSSIAN);
        std::vector<double> pulEmax_(PulseList.size(), 0.0);
        int pp = 0;
        for(auto & pul : PulseList )
        {
            pulShape_[pp] = string2prof(pul.second.get<std::string>("profile"));
            pulEmax_[pp]  = pul.second.get<double>("Field_Intensity",1.0) * a_ * EPS0 * SPEED_OF_LIGHT / I0_;
            std::vector<double> fxn;
            switch (pulShape_[pp])
            {
                case PLSSHAPE::GAUSSIAN:
                    fxn.push_back(pul.second.get<double>("fcen"));
                    fxn.push_back(1.0/pul.second.get<double>("fwidth"));
                    fxn.push_back(pul.second.get<double>("cutoff"));
                    fxn.push_back(pul.second.get<double>("t_0", fxn[1]*fxn[2]));
                break;
                case PLSSHAPE::BH:
                    fxn.push_back(pul.second.get<double>("fcen"));
                    fxn.push_back(1.0/pul.second.get<double>("fwidth", 1.0/pul.second.get<double>("tau") ) );
                    fxn.push_back(pul.second.get<double>("t_0"));
                    fxn.push_back(pul.second.get<double>("BH1"));
                    fxn.push_back(pul.second.get<double>("BH2"));
                    fxn.push_back(pul.second.get<double>("BH3"));
                    fxn.push_back(pul.second.get<double>("BH4"));
                break;
                case PLSSHAPE::RECT :
                    fxn.push_back(pul.second.get<double>("tau"));
                    fxn.push_back(pul.second.get<double>("t_0"));
                    fxn.push_back(pul.second.get<double>("n", 30) * 2.0);
                    fxn.push_back(pul.second.get<double>("fcen"));
                break;
                case PLSSHAPE::CONTINUOUS:
                    fxn.push_back(pul.second.get<double>("fcen"));
                break;
                case PLSSHAPE::RAMP_CONT:
                    fxn.push_back(pul.second.get<double>("fcen"));
                    fxn.push_back(pul.second.get<double>("ramp_val"));
                break;
                case PLSSHAPE::RICKER:
                    fxn.push_back(pul.second.get<double>("fcen"));
                    fxn.push_back(pul.second.get<double>("fwidth"));
                    fxn.push_back(pul.second.get<double>("cutoff"));
                break;
                default:
                    throw std::logic_error("This pulse shape is undefined in source ");
                break;
            }
            pulFxn_[pp] = fxn;
            ++pp;
        }
        tfsfEmax_[ii] = pulEmax_;
        tfsfPulShape_[ii] = pulShape_;
        tfsfPulFxn_[ii] = pulFxn_;
        tfsfEllipticalKratio_[ii] = iter.second.get<double>("ellpiticalKRat", 1.0);
    }

    int qq = 0;
    std::array<std::array<double,3>,3> uVecs;
    for(ii = 0; ii < uVecs.size(); ++ii)
    {
        uVecs[ii] = {{ 0.0, 0.0, 0.0}};
        uVecs[ii][ii] = 1.0;
    }
    std::vector<double> szCell = {{ size_[0], size_[1], size_[2] }};
    objArr_.push_back(std::make_shared<block>(std::vector<double>(1,1.0), std::vector<double>(1,1.0), szCell, std::array<double,3>({{0.0,0.0,0.0}}) , uVecs ) );
    for (auto& iter : IP.get_child("ObjectList"))
    {
        std::shared_ptr<Obj> obj = ptreeToObject(iter);
        objArr_.push_back(obj);
    }
    ii = 0;
    int jj = 0;
    int dd = 0;
    for (auto& iter : IP.get_child("DetectorList"))
    {
        // dtcOutType t = string2out(iter.second.get<std::string>("type"));
        // dtcType_.push_back(t);
        dtcType_[dd] = string2out(iter.second.get<std::string>("type") );
        std::string out_name = iter.second.get<std::string>("fname") + "_field_" + std::to_string(ii) + ".dat";
        boost::filesystem::path p(out_name.c_str());
        boost::filesystem::create_directories(p.remove_filename());
        ++ii;
        dtcName_[dd] = out_name;
        dtcClass_[dd] = string2dtcclass(iter.second.get<std::string>("dtc_class" ,"cout"));
        dtcSI_[dd] = iter.second.get<bool>("SI", true);
        dtcTimeInt_[dd] = iter.second.get<double>("Time_Interval", courant_/static_cast<double>(res_));

        std::array<double,3> tempSz = as_ptArr<double>(iter.second, "size");
        for(int cc = 0; cc < 3; ++cc )
            dtcSz_[dd][cc] = find_pt(tempSz[cc]) + 1;

        int i = 0;
        for(auto& loc : as_ptArr<double>(iter.second, "loc") )
        {
            dtcLoc_[dd][i] = find_pt(loc + size_[i]/2.0 - tempSz[i]/2.0);
            ++i;
        }

        dtcNFreq_[dd] = iter.second.get<int>("nfreq", -1),
        dtcfCen_[dd] = iter.second.get<double>("fcen"  , -1.0),
        dtcfWidth_[dd] = iter.second.get<double>("fwidth", -1.0),
        dtcLamL_[dd] = iter.second.get<double>("lamL", -1.0),
        dtcLamR_[dd] = iter.second.get<double>("lamR", -1.0),
        ++dd;
    }
    dd = 0;
    for (auto& iter : IP.get_child("FluxList"))
    {
        fluxName_[dd] = iter.second.get<std::string>("name");

        boost::filesystem::path p(iter.second.get<std::string>("name").c_str());
        boost::filesystem::create_directories(p.remove_filename());

        std::array<double,3> tempSz = as_ptArr<double>(iter.second, "size");
        for(int cc = 0; cc < 3; ++cc )
            fluxSz_[dd][cc] = find_pt(tempSz[cc]) + 1;

        int i = 0;
        for(auto& loc : as_ptArr<double>(iter.second, "loc") )
        {
            fluxLoc_[dd][i] = find_pt(loc + size_[i]/2.0 - tempSz[i]/2);
            ++i;
        }
        fluxWeight_[dd] = iter.second.get<double>("weight", 1.0);

        fluxTimeInt_[dd] = static_cast<int>( std::floor(iter.second.get<double>("Time_Interval", courant_/static_cast<double>(res_)) / (courant_/static_cast<double>(res_)) + 0.50) );
        if(fluxTimeInt_[dd] <= 0)
        {
            fluxTimeInt_[dd]=1;
        }

        fluxFCen_[dd] = iter.second.get<double>("fcen",-1.0);
        fluxFWidth_[dd] = iter.second.get<double>("fwidth",-1.0);
        fluxLamL_[dd] = iter.second.get<double>("lamL",-1.0);
        fluxLamR_[dd] = iter.second.get<double>("lamR",-1.0);
        fluxNFreq_[dd] = iter.second.get<int>("nfreq",-1);


        fluxSI_[dd] = iter.second.get<bool>("SI", false);
        fluxCrossSec_[dd] = iter.second.get<bool>("cross_sec", false);

        fluxIncdFieldsFilename_[dd] = iter.second.get<std::string>("incd_fileds", "");
        fluxSave_[dd] = iter.second.get<bool>("save", false);
        fluxLoad_[dd] = iter.second.get<bool>("load", false);

        if(fluxLoad_[dd] && fluxIncdFieldsFilename_[dd] == "")
            throw std::logic_error("Trying to load in file without a valid path, in the " + std::to_string(dd) +" flux detector");

        dd++;
    }

}

std::tuple<std::vector<double>,std::vector<double>>  parallelProgramInputs::getMater(std::string mat)
{
    if(mat.compare("Au") == 0 || mat.compare("au") == 0 || mat.compare("AU") == 0)
    {
        return std::make_tuple( getMetal(AU_MAT_), std::vector<double>(1,1.0) );
    }
    else if(mat.compare("Ag") == 0 || mat.compare("ag") == 0 || mat.compare("AG") == 0)
    {
        return std::make_tuple( getMetal(AG_MAT_), std::vector<double>(1,1.0) );
    }
    else if(mat.compare("Al") == 0 || mat.compare("al") == 0 || mat.compare("AL") == 0)
    {
        return std::make_tuple( getMetal(AL_MAT_), std::vector<double>(1,1.0) );
    }
    else if(mat.compare("Cu") == 0 || mat.compare("cu") == 0 || mat.compare("CU") == 0)
    {
        return std::make_tuple( getMetal(CU_MAT_), std::vector<double>(1,1.0) );
    }
    else if(mat.compare("Be") == 0 || mat.compare("be") == 0 || mat.compare("BE") == 0)
    {
        return std::make_tuple( getMetal(BE_MAT_), std::vector<double>(1,1.0) );
    }
    else if(mat.compare("Cr") == 0 || mat.compare("cr") == 0 || mat.compare("CR") == 0)
    {
        return std::make_tuple( getMetal(CR_MAT_), std::vector<double>(1,1.0) );
    }
    else if(mat.compare("Ni") == 0 || mat.compare("ni") == 0 || mat.compare("NI") == 0)
    {
        return std::make_tuple( getMetal(NI_MAT_), std::vector<double>(1,1.0) );
    }
    else if(mat.compare("Pt") == 0 || mat.compare("pt") == 0 || mat.compare("PT") == 0)
    {
        return std::make_tuple( getMetal(PT_MAT_), std::vector<double>(1,1.0) );
    }
    else if(mat.compare("Pd") == 0 || mat.compare("pd") == 0 || mat.compare("PD") == 0)
    {
        return std::make_tuple( getMetal(PD_MAT_), std::vector<double>(1,1.0) );
    }
    else if(mat.compare("Ti") == 0 || mat.compare("ti") == 0 || mat.compare("TI") == 0)
    {
        return std::make_tuple( getMetal(TI_MAT_), std::vector<double>(1,1.0) );
    }
    else if(mat.compare("W") == 0 || mat.compare("w") == 0)
    {
        return std::make_tuple( getMetal(W_MAT_), std::vector<double>(1,1.0) );
    }
    else if(mat.compare("TiO2") == 0)
    {
        std::vector<double> mater(1,6.2001);
        return std::make_tuple( mater, std::vector<double>(1,1.0) );
    }
    else if(mat.compare("SiO2") == 0)
    {
        std::vector<double> mater(1,2.1025);
        return std::make_tuple( mater, std::vector<double>(1,1.0) );
    }
    else if(mat.compare("CdSe") == 0)
    {
        std::vector<double> mater(1,6.20);
        return std::make_tuple( mater, std::vector<double>(1,1.0) );
    }
    else if(mat.compare("CdSeQD") == 0)
    {
        std::vector<double> mater = getMetal(CDSE_QD_MAT_);
        mater[0] = 6.20;
        return std::make_tuple( mater, std::vector<double>(1,1.0) );
    }
    else if(mat.compare("vac") == 0 || mat.compare("Vac") == 0 || mat.compare("VAC") == 0)
    {
        std::vector<double> mater(1,1.0);
        return std::make_tuple( mater, std::vector<double>(1,1.0) );
    }
    else
    {
        throw std::logic_error("The material name " + mat + " is not a predefined material. Please give me all the parameters or select from one of our predefined materials");
    }
}

inline double parallelProgramInputs::ev2FDTD(double eV)
{
    return eV / 4.135666e-15 * a_ / SPEED_OF_LIGHT; // Energy -> freq (SI) -> freq(FDTD)
}

std::vector<double> parallelProgramInputs::getMetal(std::vector<double> params)
{
    std::vector<double> mater = {1.0};
    double wp = ev2FDTD(params[0]);
    for(int ii = 0; ii < (params.size()-1)/3; ii++)
    {
        double f   = params[3*ii+1];
        double GAM = ev2FDTD(params[3*ii+2]);
        double OMG = ev2FDTD(params[3*ii+3]);
        if(OMG == 0.0)
            OMG = ev2FDTD(1.0e-20);
        mater.push_back(f*pow(wp/OMG,2));
        mater.push_back(GAM*M_PI);
        mater.push_back(OMG*2.0*M_PI);
    }
    return mater;
}

POLARIZATION parallelProgramInputs::string2pol(std::string p)
{
    if(p.compare("Ex") == 0)
        return POLARIZATION::EX;
    else if(p.compare("Ey") == 0)
        return POLARIZATION::EY;
    else if(p.compare("Ez") == 0)
        return POLARIZATION::EZ;
    else if(p.compare("Hx") == 0)
        return POLARIZATION::HX;
    else if(p.compare("Hy") == 0)
        return POLARIZATION::HY;
    else if(p.compare("Hz") == 0)
        return POLARIZATION::HZ;
    else if(p.compare("L") == 0)
        return POLARIZATION::L;
    else if(p.compare("R") == 0)
        return POLARIZATION::R;
    else
        throw std::logic_error("POLARIZATION undefined");

}

DIRECTION parallelProgramInputs::string2dir(std::string dir)
{
    if((dir.compare("x") == 0) || (dir.compare("X") == 0))
        return DIRECTION::X;
    else if((dir.compare("y") == 0) || (dir.compare("Y") == 0))
        return DIRECTION::Y;
    else if ((dir.compare("z") == 0) || (dir.compare("Z") == 0))
        return DIRECTION::Z;
    else
        throw std::logic_error("I would appricate it if you stick to the standard X,Y,Z DIRECTIONs. While it's fun to invent new ones, it is very hard to do calculations if I don't even understand what dimension I am in. Please try again!");

}

DISTRIBUTION parallelProgramInputs::string2dist(std::string dist)
{
    if( (dist.compare("Gaussian") == 0) || (dist.compare("gaussian") == 0))
        return DISTRIBUTION::GAUSSIAN;
    else if( (dist.compare("delta_fxn") == 0) || (dist.compare("Delta_fxn") == 0) || (dist.compare("Delta_Fxn") == 0) || (dist.compare("delta_Fxn") == 0) || (dist.compare("DeltaFxn") == 0))
        return DISTRIBUTION::DELTAFXN;
    else if( (dist.compare("skew_normal") == 0) || (dist.compare("Skew_Normal") == 0) || (dist.compare("Skew_normal") == 0) || (dist.compare("skew_Normal") == 0) )
        return DISTRIBUTION::SKEW_NORMAL;
    else if( (dist.compare("chi_squared") == 0) || (dist.compare("Chi_Squared") == 0) || (dist.compare("Chi_squared") == 0) || (dist.compare("chi_Squared") == 0))
        return DISTRIBUTION::CHI_SQUARED;
    else
        throw std::logic_error("The distribution type " + dist + " is not defined, please define it for me.");

}
SHAPE parallelProgramInputs::string2shape(std::string s)
{
    if(s.compare("sphere") == 0)
        return SHAPE::SPHERE;
    else if (s.compare("block") == 0)
        return SHAPE::BLOCK;
    else if (s.compare("triangle_prism") == 0)
        return SHAPE::TRIANGLE_PRISM;
    else if (s.compare("trapezoid_prism") == 0)
        return SHAPE::TRAPEZOIDAL_PRISM;
    else if (s.compare("ters_tip") == 0)
        return SHAPE::TERS_TIP;
    else if (s.compare("ellipsoid") == 0)
        return SHAPE::ELLIPSOID;
    else if (s.compare("cylinder") == 0)
        return SHAPE::CYLINDER;
    else if (s.compare("cone") == 0)
        return SHAPE::CONE;
    else if (s.compare("parabolic_ters_tip") == 0)
        return SHAPE::PARABOLIC_TERS_TIP;
    else
        throw std::logic_error("Shape undefined");
}

DTCTYPE parallelProgramInputs::string2out(std::string t)
{
    if(t.compare("Ex") == 0)
        return DTCTYPE::EX;
    else if(t.compare("Ey") == 0)
        return DTCTYPE::EY;
    else if(t.compare("Ez") == 0)
        return DTCTYPE::EZ;
    else if(t.compare("Hx") == 0)
        return DTCTYPE::HX;
    else if(t.compare("Hy") == 0)
        return DTCTYPE::HY;
    else if(t.compare("Hz") == 0)
        return DTCTYPE::HZ;
    else if(t.compare("Px") == 0)
        return DTCTYPE::PX;
    else if(t.compare("Py") == 0)
        return DTCTYPE::PY;
    else if(t.compare("Pz") == 0)
        return DTCTYPE::PZ;
    else if(t.compare("E_pow") == 0)
        return DTCTYPE::EPOW;
    else if(t.compare("H_pow") == 0)
        return DTCTYPE::HPOW;
    else
        throw std::logic_error("I'm afraid I don't understand this detector type please, either give me a new definition in the string2out function of Inputs.cpp or try again. I hope you go for the former I always like new definitions");
}

DTCCLASS parallelProgramInputs::string2dtcclass(std::string c)
{
    if(c.compare("bin") == 0)
        return DTCCLASS::BIN;
    else if(c.compare("txt") == 0)
        return DTCCLASS::TXT;
    else if(c.compare("cout") == 0)
        return DTCCLASS::COUT;
    else if(c.compare("freq") == 0)
        return DTCCLASS::FREQ;
    else
        throw std::logic_error("I'm afraid I don't understand this detector class please, either give me a new definition in the string2dtcclass function of Inputs.cpp or try again. I hope you go for the former I always like new definitions");
}

PLSSHAPE parallelProgramInputs::string2prof(std::string p)
{
    if(p.compare("gaussian") == 0)
        return PLSSHAPE::GAUSSIAN;
    else if(p.compare("BH") == 0)
        return PLSSHAPE::BH;
    else if(p.compare("rectangle") == 0)
        return PLSSHAPE::RECT;
    else if(p.compare("continuous") == 0)
        return PLSSHAPE::CONTINUOUS;
    else if(p.compare("ricker") == 0)
        return PLSSHAPE::RICKER;
    else if(p.compare("ramped_cont") == 0)
        return PLSSHAPE::RAMP_CONT;
    else
        throw std::logic_error("Pulse sahpe undefined");
}


std::shared_ptr<Obj> parallelProgramInputs::ptreeToObject(boost::property_tree::ptree::value_type &iter)
{
    //  Get information common amongst all objects
    std::shared_ptr<Obj> out;
    double ang;

    std::array<double,3> loc = as_ptArr<double>(iter.second, "loc");

    // Parse the material information
    std::string material = iter.second.get<std::string>("material");
    std::vector<double> mater;
    std::vector<double> magMater;
    if(material.compare("custom") == 0)
    {
        mater.push_back(iter.second.get<double>("eps",1.00));
        boost::property_tree::ptree& pols = iter.second.get_child("pols");
        for(auto& iter2 : pols)
        {
            mater.push_back(iter2.second.get<double>("sigma"));
            mater.push_back(iter2.second.get<double>("gamma")*M_PI);
            mater.push_back(iter2.second.get<double>("omega")*2*M_PI);
        }
        magMater.push_back( iter.second.get<double>("mu", 0.0) );
        boost::property_tree::ptree& magPols = iter.second.get_child("magPols");
        for(auto& iter2 : magPols)
        {
            magMater.push_back(iter2.second.get<double>("sigma"));
            magMater.push_back(iter2.second.get<double>("gamma")*M_PI);
            magMater.push_back(iter2.second.get<double>("omega")*2*M_PI);
        }
    }
    else
    {
        std::tie(mater, magMater) = getMater(material);
    }
    std::vector<double> geo_param;
    //Make unit vectors from user inputs
    std::array<std::array<double,3>,3> unitVecs;
    boost::property_tree::ptree& uvecs = iter.second.get_child("unit_vectors");
    int cc = 0;
    for(auto& iter2 : uvecs)
    {
        unitVecs[cc] = as_ptArr<double>(iter2.second, "uvec");
        ++cc;
    }
    //If user inputs are not there use angles
    if(uvecs.size() == 0)
    {
        double orTheta = iter.second.get<double>("orTheta", 0.0) * M_PI / 180.0;
        double orPhi   = iter.second.get<double>("orPhi", 0.0) * M_PI / 180.0;
        if(loc.size() == 2)
        {
            unitVecs[0] = {{ cos(orPhi), -1.0*sin(orPhi), 0}};
            unitVecs[1] = {{ sin(orPhi),      cos(orPhi), 0}};
            unitVecs[2] = {{        0.0,             0.0, 0}};
        }
        else
        {
            // First rotate by theta about x axis then phi around z axis
            unitVecs[0] = {{ cos(orPhi)             , -1.0*        sin(orPhi),  0.0 }};
            unitVecs[1] = {{ cos(orTheta)*sin(orPhi), cos(orTheta)*cos(orPhi), -1.0*sin(orTheta) }};
            unitVecs[2] = {{ sin(orTheta)*sin(orPhi), sin(orTheta)*cos(orPhi),      cos(orTheta) }};
        }
    }
    SHAPE s = string2shape(iter.second.get<std::string>("shape"));
    switch (s)
    {
        case(SHAPE::BLOCK):
            geo_param = as_vector<double>(iter.second, "size");
            geo_param.push_back(iter.second.get<double>("rad_curve", 0.0) );
            if(geo_param.back() == 0.0)
                out = std::make_shared<block>(mater, magMater, geo_param, loc, unitVecs);
            else
                out = std::make_shared<rounded_block>(mater, magMater, geo_param, loc, unitVecs);
        break;
        case (SHAPE::ELLIPSOID):
            geo_param = as_vector<double>(iter.second, "size");
            out = std::make_shared<ellipsoid>(mater, magMater, geo_param, loc, unitVecs);
        break;
        case (SHAPE::SPHERE):
            geo_param = { iter.second.get<double>("radius",0.0) };
            out = std::make_shared<sphere>(mater, magMater, geo_param, loc, unitVecs);
        break;
        case (SHAPE::CYLINDER):
            geo_param = { iter.second.get<double>("radius",0.0), iter.second.get<double>("length",0.0) };
            out = std::make_shared<cylinder>(mater, magMater, geo_param, loc, unitVecs);
        break;
        case (SHAPE::CONE):
        {
            double rad1 =  iter.second.get<double>("radius1",0.0);
            double rad2 =  iter.second.get<double>("radius2",0.0);
            double len  =  iter.second.get<double>("height",0.0);
            double raiseAng = atan( (rad2 - rad1) / len);
            geo_param = { rad2, rad1, len  };
            out = std::make_shared<cone>(mater, magMater, geo_param, loc, unitVecs);
            break;
        }
        case (SHAPE::TRAPEZOIDAL_PRISM):
            geo_param = {{ iter.second.get<double>("base1", 0.0), iter.second.get<double>("base2", 0.0), iter.second.get<double>("height", 0.0), iter.second.get<double>("length", 0.0) }};
            out = std::make_shared<trapezoid_prism>(mater, magMater, geo_param, loc, unitVecs);
        break;
        case (SHAPE::TRIANGLE_PRISM):
            geo_param = {{ iter.second.get<double>("base", 0.0), iter.second.get<double>("height", 0.0), iter.second.get<double>("length", 0.0) }};
            out = std::make_shared<isosceles_tri_prism>(mater, magMater, geo_param, loc, unitVecs);
        break;
        case (SHAPE::TERS_TIP):
            geo_param = {{ 2.0*iter.second.get<double>("rad_curve", 0.0), iter.second.get<double>("base", 0.0),  iter.second.get<double>("length", 0.0) }};
            out = std::make_shared<ters_tip>(mater, magMater, geo_param, loc, unitVecs);
        break;
        case (SHAPE::PARABOLIC_TERS_TIP):
            geo_param = {{ 1.0/(2.0*iter.second.get<double>("rad_curve", 0.0) ) , iter.second.get<double>("length", 0.0) }};
            out = std::make_shared<parabolic_ters_tip>(mater, magMater, geo_param, loc, unitVecs);
        break;
        default:
            throw std::logic_error("You probably should specify a supported shape. You'll like your results better that way.");
        break;
    }
    return out;
}
void stripComments(std::string& filename)
{
    //Open input and output file
    std::string newfn = "stripped_" + filename;
    std::fstream inputfile;
    inputfile.open(filename);
    std::ofstream inputcopy;
    inputcopy.open(newfn);

    //search for '//', delete everything following, print remainder to new file
    std::string line;
    int found, found2;
    while (getline(inputfile,line))
    {
        found  = line.find('/');
        found2 = line.find('/', found+1);
        if (found != line.npos && found2 == found+1)
            inputcopy << line.erase(found, line.length()) << std::endl;
        else
            inputcopy << line << std::endl;
    }
    inputcopy.close();
    //update filename;
    filename = newfn;
}

