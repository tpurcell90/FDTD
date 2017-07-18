#ifndef PRALLEL_FDTD_INPUTS
#define PRALLEL_FDTD_INPUTS

#include <src/OBJECTS/Obj.hpp>
#include <src/UTIL/ml_consts.hpp>
#include <src/UTIL/dielectric_params.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iterator>
#ifdef MKL
#include <UTIL/utilities_MKL.hpp>
#elif defined ACML
#include <UTIL/utilities_acml.hpp>
#else
#include <UTIL/utilities_MKL.hpp>
#endif

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
    std::vector<int> dtcNFreq_; //!< number of frequencies to take in the range
    std::vector<double> dtcfCen_; //!< center frequency
    std::vector<double> dtcfWidth_; //!< frequency width
    std::vector<double> dtcLamL_; //!< left wavelength end point
    std::vector<double> dtcLamR_; //!< right wavelength end point

    std::vector<int> fluxXOff_; //!< the x location offset of the fields
    std::vector<int> fluxYOff_; //!< the y location offset of the fields
    std::vector<int> fluxTimeInt_; //!< time interval used to intake the fields
    std::vector<std::array<int,3>> fluxLoc_; //!< Location of the lower left corner of the flux surface
    std::vector<std::array<int,3>> fluxSz_; //!< size of the flux surface
    std::vector<double> fluxWeight_; //!< weight of the flux surface
    std::vector<std::string> fluxName_; //!< file name for the output file
    std::vector<double> fluxFCen_; //!<  center frequency
    std::vector<double> fluxFWidth_; //!< frequency width of the flux
    std::vector<double> fluxLamL_; //!< left wavelength limit
    std::vector<double> fluxLamR_; //!< right wavelength limit
    std::vector<bool> fluxSI_; //!<  use SI units
    std::vector<bool> fluxCrossSec_; //!< calculate the cross-section?
    std::vector<bool> fluxSave_; //!< save the fields?
    std::vector<bool> fluxLoad_; //!< load the fields?
    std::vector<std::string> fluxIncdFieldsFilename_; //!< incident file names
    std::vector<int> fluxNFreq_; //!< number of frequencies

    /**
     * @brief Constructs the input parameter file
     *
     * @param IP boost property tree object
     * @param fn input filename
     */
    parallelProgramInputs(boost::property_tree::ptree IP,std::string fn);
    /**
     * @brief get the dielectric parameters for a material
     * @details uses a string identifier to get the dielectric parameters for a material
     *
     * @param mat string identifier for the material
     * @return dielectric parameters
     */
    std::vector<double> getDielectricParams(std::string mat);

    /**
     * @brief Converts the string from an input file into a POLARIZATION
     *
     * @param p String type of POLARIZATION
     * @return The POLARIZATION
     */
    POLARIZATION string2pol(std::string p);
    /**
     * @brief Converts the string from an input file into a shape
     *
     * @param p String type of shape
     * @return The shape
     */
    SHAPE string2shape(std::string s);
    /**
     * @brief Converts the string from an input file into a Detector Type
     *
     * @param p String type of Detector Type
     * @return The Detector Type
     */
    DTCTYPE string2out(std::string t);
    /**
     * @brief Converts the string from an input file into a DTCCLASS
     *
     * @param c String type of DTCCLASS
     * @return The DTCCLASS
     */
    DTCCLASS string2dtcclass(std::string c);

    /**
     * @brief Converts the string from an input file into a Pulse Shape
     *
     * @param p String type of Pulse Shape
     * @return The Pulse Shape
     */
    PLSSHAPE string2prof(std::string p);
    /**
     * @brief Converts the string from an input file into a DIRECTION
     *
     * @param dir String type of DIRECTION
     * @return The DIRECTION
     */
    DIRECTION string2dir(std::string dir);

    /**
     * @brief Converts the string from an input file into a DIRECTION
     *
     * @param dist String key for each distribution
     * @return The DISTRIBUTION
     */
    DISTRIBUTION string2dist(std::string dist);

    /**
     * @brief get the material parameters from a string identifier
     * @details returns the material parameters
     *
     * @param mat string material identifier
     * @return material parameters
     */
    std::tuple<std::vector<double>, std::vector<double> > getMater(std::string mat);
    /**
     * @brief Converts the eV dielectric function parameters to FDTD units
     * @details Preforms a unit conversion on all dielectric parameters to get them in the correct units
     *
     * @param eV parameter
     * @return dielectric parameters in FDTD units
     */
    std::vector<double> getMetal(std::vector<double> params);
    /**
     * @brief converts eV to FDTDsunits
     *
     * @param eV a value in eV
     * @return a value in FDTD units
     */
    double ev2FDTD(double eV);
    /**
     * @brief Converts a distance in unit lengths to grid points
     *
     * @param pt distance in unit lengths
     * @return distance in grid points
     */
    inline int find_pt(double pt) {return floor(pt*static_cast<double>(res_) + 0.5);}
    /**
     * @return tMax
     */
    inline double tMax() {return tMax_;}
    /**
     * @brief converts a ptree object into an Object
     * @details take values in a ptree to create a new object
     *
     * @param iter ptree object describing the object
     * @return the Object described by iter
     */
    std::shared_ptr<Obj> ptreeToObject(boost::property_tree::ptree::value_type &iter);
    /**
     * @brief removes comments from input file
     */
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