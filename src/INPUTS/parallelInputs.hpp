#ifndef PRALLEL_FDTD_INPUTS
#define PRALLEL_FDTD_INPUTS

#include <src/OBJECTS/Obj.hpp>
#include <src/UTIL/FDTD_consts.hpp>
#include <src/UTIL/dielectric_params.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iterator>

struct EnergyLevelDiscriptor
{
    DISTRIBUTION EDist_; //!< Distribution type for the energy level
    int nstates_; //!< number of states in the level
    int levDescribed_; //!< how many levels are described by this distribution
    std::vector<double> energyStates_; //!< energy of all the states in the level
    std::vector<double> weights_; //!< weights for the levels
};

/**
 * @brief input parameters class for the parallel FDTD class
 * @details sores all the values necessary to convert the input files into a fdtd grid
 */
class parallelProgramInputs
{
public:
    bool periodic_; //!< if true use PBC

    POLARIZATION pol_; //!< polarization of the grid (EX, EY, HZ = TE; HX, HY, EZ = TM)

    std::string filename_; //!< filename of the input file
    int res_; //!< number of grid points per unit length

    double courant_; //!< Courant factor of the cell
    double a_; //!< the unit length of the calculations
    double tMax_; //!< the max time of the calculation
    double I0_; //!< the unit current of the cell

    bool cplxFields_; //!< if true use complex fields
    bool saveFreqField_; //!< if true save the flux fields

    std::array<int,3> pmlThickness_; //!< thickness of the PMLs in all directions
    std::array<double,3> k_point_; //!< k_point vector of the light
    std::array<double,3> size_; //!< size of the cell in units of the unit length
    double pmlAMax_; //!< max a value for CPML
    double pmlMa_; //!< scaling factor of a for the CPML
    double pmlM_; //!< scaling factor for m for the CPML

    std::vector<double> inputMapSlicesX_; //!< list of slices in the YZ plane
    std::vector<double> inputMapSlicesY_; //!< list of slices in the XZ plane
    std::vector<double> inputMapSlicesZ_; //!< list of slices in the XY plane

    std::vector<POLARIZATION> srcPol_; //!<polarization of all sources
    std::vector<std::vector<std::vector<double>>> srcFxn_; //!< pulse function parameters for the sources
    std::vector<std::array<int,3>> srcLoc_; //!<location of all sources
    std::vector<std::array<int,3>> srcSz_; //!<sizes of all all sources
    std::vector<double> srcPhi_; //!< angle of incidence for the source detector
    std::vector<double> srcTheta_; //!< angle of incidence for the source detector
    std::vector<std::vector<PLSSHAPE>> srcPulShape_; //!< pulse shape of all sources
    std::vector<std::vector<double>> srcEmax_; //!< max E incd field of all sources
    std::vector<double> srcEllipticalKratio_; //!< ratio between the long and short axis for current source surfaces
    std::vector<double> srcPsi_; //!< angle of polarization for circular and elliptical light

    std::vector<std::array<int,3>> tfsfLoc_;//!< location of the tfsf lower left point
    std::vector<std::array<int,3>> tfsfSize_; //!< size of the TFSF region
    std::vector<double> tfsfTheta_; //!< Phi angle of the TFSF surface
    std::vector<double> tfsfPhi_; //!< Phi angle of the TFSF surface
    std::vector<double> tfsfPsi_; //!< Phi angle of the TFSF surface
    std::vector<std::vector<std::vector<double>>> tfsfPulFxn_; //!< pulse function parameters for the sources
    std::vector<std::vector<PLSSHAPE>> tfsfPulShape_; //!< pulse shape of all sources
    std::vector<std::vector<double>> tfsfEmax_; //!< max E incd field of all sources
    std::vector<POLARIZATION> tfsfCircPol_; //!< sets the TFSF surface to be L or R polarized
    std::vector<double> tfsfEllipticalKratio_; //!< ratio between the long and short axis for TFSF surfaces

    std::vector<std::shared_ptr<Obj>> objArr_; //!< array of all objects in the grid

    std::vector<DTCCLASS> dtcClass_; //!< describer of output file type for all detectors
    std::vector<bool> dtcSI_; //!< if true use SI units
    std::vector<std::array<int,3>> dtcLoc_; //!< location of all detectors
    std::vector<std::array<int,3>> dtcSz_; //!< sizes of all detectors
    std::vector<std::string> dtcName_; //!< file name for all detectors
    std::vector<DTCTYPE> dtcType_; //!< whether the dtc stores the Ex, Ey, Ez, Hx, Hy, Hz fields or the H or E power
    std::vector<double> dtcTimeInt_; //!< time interval for all detectors
    std::vector<GRIDOUTFXN> dtcOutBMPFxnType_; //!< what function should bmp converter use
    std::vector<GRIDOUTTYPE> dtcOutBMPOutType_; //!< how to output the values for the detector ina text file
    std::vector<std::vector<double>> dtcFreqList_; //!< center frequency

    std::vector<int> fluxXOff_; //!< the x location offset of the fields
    std::vector<int> fluxYOff_; //!< the y location offset of the fields
    std::vector<int> fluxTimeInt_; //!< time interval used to intake the fields
    std::vector<std::array<int,3>> fluxLoc_; //!< Location of the lower left corner of the flux surface
    std::vector<std::array<int,3>> fluxSz_; //!< size of the flux surface
    std::vector<double> fluxWeight_; //!< weight of the flux surface
    std::vector<std::string> fluxName_; //!< file name for the output file
    std::vector<std::vector<double>> fluxFreqList_; //!<  center frequency
    std::vector<bool> fluxSI_; //!<  use SI units
    std::vector<bool> fluxCrossSec_; //!< calculate the cross-section?
    std::vector<bool> fluxSave_; //!< save the fields?
    std::vector<bool> fluxLoad_; //!< load the fields?
    std::vector<std::string> fluxIncdFieldsFilename_; //!< incident file names

    /**
     * @brief      Constructs the input parameter object
     *
     * @param[in]  IP    boost property tree generated from the input json file
     * @param[in]  fn    The filename of the input file
     */
    parallelProgramInputs(boost::property_tree::ptree IP,std::string fn);

    /**
     * @brief      Gets the dielectric parameters for a material.
     *
     * @param[in]  mat   String identifier of the material
     *
     * @return     The dielectric parameters.
     */
    std::vector<double> getDielectricParams(std::string mat);

    /**
     * @brief      converts a string to GRIDOUTFXN
     *
     * @param[in]  f     String identifier to a GRIDOUTFXN
     *
     * @return     GRIDOUTFXN from that input string
     */
    GRIDOUTFXN string2GRIDOUTFXN (std::string f);

    /**
     * @brief      converts a string to GRIDOUTTYPE
     *
     * @param[in]  t     String identifier to a GRIDOUTTYPE
     *
     * @return     GRIDOUTTYPE from that input string
     */
    GRIDOUTTYPE string2GRIDOUTTYPE (std::string t);

    /**
     * @brief      converts a string to POLARIZATION
     *
     * @param[in]  p     String identifier to a POLARIZATION
     *
     * @return     POLARIZATION from that input string
     */
    POLARIZATION string2pol(std::string p);

    /**
     * @brief      converts a string to SHAPE
     *
     * @param[in]  s     String identifier to a SHAPE
     *
     * @return     SHAPE from that input string
     */
    SHAPE string2shape(std::string s);

    /**
     * @brief      converts a string to DTCTYPE
     *
     * @param[in]  t     String identifier to a DTCTYPE
     *
     * @return     DTCTYPE from that input string
     */
    DTCTYPE string2out(std::string t);

    /**
     * @brief      converts a string to DTCCLASS
     *
     * @param[in]  c     String identifier to a DTCCLASS
     *
     * @return     DTCCLASS from that input string
     */
    DTCCLASS string2dtcclass(std::string c);

    /**
     * @brief      converts a string to PLSSHAPE
     *
     * @param[in]  p     String identifier to a PLSSHAPE
     *
     * @return     PLSSHAPE from that input string
     */
    PLSSHAPE string2prof(std::string p);

    /**
     * @brief      converts a string to DIRECTION
     *
     * @param[in]  dir   String identifier to a DIRECTION
     *
     * @return     DIRECTION from that input string
     */
    DIRECTION string2dir(std::string dir);

    /**
     * @brief      Gets the material parameters for a given material
     *
     * @param[in]  mat   String identifying the material
     *
     * @return     The material parameters
     */
    std::tuple<std::vector<double>, std::vector<double>> getMater(std::string mat);

    /**
     * @brief      Converts the metal parameters from eV based units to FDTD based units
     *
     * @param[in]  params  The parameters for the metallic material in eV
     *
     * @return     The parameters for the metallic material in FDTD units
     */
    std::vector<double> getMetal(std::vector<double> params);

    /**
     * @brief      Converts eV units to FDTD units
     *
     * @param[in]  eV    The value in eV
     *
     * @return     The value in FDTD frequency units
     */
    double ev2FDTD(double eV);

    /**
     * @brief      Converts the point from real space to grid points
     *
     * @param[in]  pt    Real space value
     *
     * @return     Corresponding grid point value
     */
    inline int find_pt(double pt) {return floor(pt*static_cast<double>(res_) + 0.5);}

    /**
     * @return     Maximum time of the simulation
     */
    inline double tMax() {return tMax_;}
    /**
     * @brief converts a ptree object into an Object
     * @details take values in a ptree to create a new object
     *
     * @param iter ptree object describing the object
     * @return the Object described by iter
     */
    /**
     * @brief      Constructs an object from the object list child tree
     *
     * @param      iter  boost::ptree child corresponding to the object
     *
     * @return     shared_ptr to the object
     */
    std::shared_ptr<Obj> ptreeToObject(boost::property_tree::ptree::value_type &iter);
};
/**
 * @brief      strips comments from the input file
 *
 * @param      filename  The filename of the file to strip
 */
void stripComments(std::string& filename);


/**
 * @brief      boost json to std::vector<T>
 *
 * @param[in]  pt          property tree
 * @param[in]  key         property tree key
 *
 * @tparam     T           double, int
 *
 * @return     json input as a std::vector<T>
 */
template <typename T>
std::vector<T> as_vector(boost::property_tree::ptree const &pt, boost::property_tree::ptree::key_type const &key)
{
    std::vector<T> r;
    for (auto& item : pt.get_child(key))
        r.push_back(item.second.get_value<T>());
    return r;
}

/**
 * @brief      boost json to std::vector<T>
 *
 * @param[in]  pt          property tree
 * @param[in]  key         property tree key
 * @param[in]  defaultVal  The default value
 * @param[in]  szVec       The size of the vector
 *
 * @tparam     T           double, int
 *
 * @return     json input as a std::vector<T>
 */
template <typename T>
std::vector<T> as_vector(boost::property_tree::ptree const &pt, boost::property_tree::ptree::key_type const &key, T defaultVal, int szVec)
{
    std::vector<T> r;
    try
    {
        for (auto& item : pt.get_child(key))
            r.push_back(item.second.get_value<T>());
    }
    catch(std::exception& e)
    {
        r = std::vector<T>(szVec, defaultVal);
    }
    return r;
}


/**
 * @brief      boost json to std::array<T,3>
 *
 * @param[in]  pt          property tree
 * @param[in]  key         property tree key
 * @param[in]  defaultVal  The default value
 *
 * @tparam     T           double, int
 *
 * @return     json input as a std::array<T,3>
 */
template <typename T>
std::array<T,3> as_ptArr(boost::property_tree::ptree const &pt, boost::property_tree::ptree::key_type const &key, T defaultVal=0)
{
    std::array<T, 3> r = {0,0,0};
    try
    {
        int ii = 0;
        for (auto& item : pt.get_child(key))
        {
            r[ii] = item.second.get_value<T>();
            ++ii;
        }
    }
    catch(std::exception& e)
    {
        r = {{ defaultVal, defaultVal, defaultVal}};
    }
    return r;
}

#endif