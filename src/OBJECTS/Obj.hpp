#ifndef FDTD_OBJECT
#define FDTD_OBJECT

#include <math.h>
#include <UTIL/enum.hpp>
#include <UTIL/typedefs.hpp>
#ifdef MKL
#include <UTIL/utilities_MKL.hpp>
#elif defined ACML
#include <UTIL/utilities_acml.hpp>
#else
#include <UTIL/utilities_MKL.hpp>
#endif

// #include <iostream>

class Obj
{
protected:
    std::vector<double> geoParam_; //!< parameters describing the geometry of the objects
    std::vector<double> material_; //!< parameters for storing the material properties of the system
    std::vector<double> magMaterial_; //!< parameters for storing the material properties of the system
    std::vector<double> chiMaterial_; //!< parameters describing the chiral response of the material
    std::vector<double> alpha_; //!<Lorentz model pole parameter (see Taflove 2005 book ch 9)
    std::vector<double> xi_; //!<Lorentz model pole parameter (see Taflove 2005 book ch 9)
    std::vector<double> gamma_; //!<Lorentz model pole parameter (see Taflove 2005 book ch 9)

    std::vector<double> magAlpha_; //!<Lorentz model pole parameter for magnetic materials (equivlant model for normal dispersive material)
    std::vector<double> magXi_; //!<Lorentz model pole parameter for magnetic materials (equivlant model for normal dispersive material)
    std::vector<double> magGamma_; //!<Lorentz model pole parameter for magnetic materials (equivlant model for normal dispersive material)

    std::array<double,3> location_; //!< location of the center point of the object
    std::array<double,9> coordTransform_; //!< Coordinate Transform Matrix
public:

    /**
     * @brief Constructor
     *
     * @param s Shape of the object
     * @param geo size (or other geometric parameters of the object)
     * @param loc location of the object
     */
    Obj(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);
    /**
     * @brief Copy Constructor
     *
     * @param o object to be copied
     */
    Obj(const Obj &o);
    /**
     * @brief Calculates the dispersion parameters from physical ones
     * @details Show equations here
     *
     * @param dt time step
     */
    void setUpConsts (double dt);

    /**
     * @return the shape of the object
     */
    // Shape s() {return part_;}
    /**
     * @return the shape of the object
     */
    inline std::vector<double>& geo() {return geoParam_;}
    /**
     * @return the location of the object
     */
    inline std::array<double,3>& loc() {return location_;}
    /**
     * @return the material parameters of the object
     */
    inline std::vector<double>& mat() {return material_;}
    /**
     * @return the alpha value for each Lorenz-Drude Pole in the material of the object
     */
    inline std::vector<double>& alpha() {return alpha_;}
    /**
     * @return the xi value for each Lorenz-Drude Pole in the material of the object
     */
    inline std::vector<double>& xi() {return xi_;}
    /**
     * @return the gamma value for each Lorenz-Drude Pole in the material of the object
     */
    inline std::vector<double>& gamma() {return gamma_;}

    /**
     * @return the material parameters of the object
     */
    inline std::vector<double>& magMat() {return magMaterial_;}
    /**
     * @return the alpha value for each Lorenz-Drude Pole in the material of the object
     */
    inline std::vector<double>& magAlpha() {return magAlpha_;}
    /**
     * @return the xi value for each Lorenz-Drude Pole in the material of the object
     */
    inline std::vector<double>& magXi() {return magXi_;}
    /**
     * @return the gamma value for each Lorenz-Drude Pole in the material of the object
     */
    inline std::vector<double>& magGamma() {return magGamma_;}

    /**
     * @return     return the $\eps_{\infty}$
     */
    inline double epsInfty(){return material_[0];}

    inline double muInfty(){return magMaterial_[0];}
    /**
     * @return     return the $\eps_{\infty}$
     */
    inline double tellegen(){return chiMaterial_[0];}

    /**
     * @brief Determine if a point is in an object
     * @details Deterimes a point in vector form is inside the object or not
     *
     * @param v point to test if it's in an object
     * @return True if inside the object false if not
     */
    virtual bool isObj(std::array<double,3> &v, double dx) = 0;
    /**
     * @brief Determines teh distance between two points
     *
     * @param pt1 The Second point
     * @param pt2 The Second point
     * @return The distance between points
     */
    double dist(std::array<double,3> pt1, std::array<double,3> pt2);
    // void mvLoc(std::array<double,3> locOff);
    virtual SHAPE shape() = 0;
};

class sphere : public Obj
{
public:
    /**
     * @brief Constructs a sphere
     * @details Constructs a sphere using the Obj constructor
     *
     * @param mater A vector with all the material descriptors in it
     * @param geo A vector containing all the geometric properties of the sphere
     * @param loc A vector for the center point location of the sphere
     * @param ang orientation angle of the sphere
     */
    sphere(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);
    /**
     * @brief Copy Constructor
     *
     * @param o sphere to be copied
     */
    sphere(const sphere &o);
    /**
     * @brief Determines if a point is inside the object
     * @details Determines if a real space point is within the sphere
     *
     * @param v Point to check
     * @param dx spatial step of the Cell
     *
     * @return Returns true if v is in the object
     */
    bool isObj(std::array<double,3> &v, double dx);
    /**
     * @brief returns the Shape parameter of the sphere
     * @details Returns the type of object it is by mapping the object name to a enum Shape
     * @return SHAPE::SPHERE
     */
    SHAPE shape() {return SHAPE::SPHERE;}
};

class block : public Obj
{
public:
    /**
     * @brief Constructs a block
     * @details Constructs a block using the Obj constructor
     *
     * @param mater A vector with all the material descriptors in it
     * @param geo A vector containing all the geometric properties of the block
     * @param loc A vector for the center point location of the block
     * @param ang orientation angle of the block
     */
    block(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);
    /**
     * @brief Copy Constructor
     *
     * @param o block to be copied
     */
    block(const block &o);
    /**
     * @brief Determines if a point is inside the object
     * @details Determines if a real space point is within the block
     *
     * @param v Point to check
     * @param dx spatial step of the Cell
     *
     * @return Returns true if v is in the object
     */
    bool isObj(std::array<double,3> &v, double dx);
    /**
     * @brief returns the Shape parameter of the block
     * @details Returns the type of object it is by mapping the object name to a enum Shape
     * @return SHAPE::BLOCK
     */
    SHAPE shape() {return SHAPE::BLOCK;}
};
class rounded_block : public Obj
{
protected:
    std::array<std::array<double,3>,8> curveCens_;
public:
    /**
     * @brief Constructs a block
     * @details Constructs a block using the Obj constructor
     *
     * @param mater A vector with all the material descriptors in it
     * @param geo A vector containing all the geometric properties of the block
     * @param loc A vector for the center point location of the block
     * @param ang orientation angle of the block
     */
    rounded_block(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);
    /**
     * @brief Copy Constructor
     *
     * @param o block to be copied
     */
    rounded_block(const rounded_block &o);
    /**
     * @brief Determines if a point is inside the object
     * @details Determines if a real space point is within the block
     *
     * @param v Point to check
     * @param dx spatial step of the Cell
     *
     * @return Returns true if v is in the object
     */
    bool isObj(std::array<double,3> &v, double dx);
    /**
     * @brief returns the Shape parameter of the block
     * @details Returns the type of object it is by mapping the object name to a enum Shape
     * @return SHAPE::BLOCK
     */
    SHAPE shape() {return SHAPE::ROUNDED_BLOCK;}
};
class ellipsoid : public Obj
{
public:
    /**
     * @brief Constructs a ellipsoid
     * @details Constructs a ellipsoid using the Obj constructor
     *
     * @param mater A vector with all the material descriptors in it
     * @param geo A vector containing all the geometric properties of the ellipsoid
     * @param loc A vector for the center point location of the ellipsoid
     * @param ang orientation angle of the ellipsoid
     */
    ellipsoid(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);
    /**
     * @brief Copy Constructor
     *
     * @param o ellipsoid to be copied
     */
    ellipsoid(const ellipsoid &o);
    /**
     * @brief Determines if a point is inside the object
     * @details Determines if a real space point is within the ellipsoid
     *
     * @param v Point to check
     * @param dx spatial step of the Cell
     *
     * @return Returns true if v is in the object
     */
    bool isObj(std::array<double,3> &v, double dx);
    /**
     * @brief returns the Shape parameter of the ellipsoid
     * @details Returns the type of object it is by mapping the object name to a enum Shape
     * @return SHAPE::ELLIPSE
     */
    SHAPE shape() {return SHAPE::ELLIPSOID;}
};

class cone : public Obj
{
public:

    cone(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    cone(const cone &o);

    bool isObj(std::array<double,3> &v, double dx);

    SHAPE shape() {return SHAPE::CONE;}
};

class cylinder : public Obj
{
public:

    cylinder(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    cylinder(const cylinder &o);

    bool isObj(std::array<double,3> &v, double dx);

    SHAPE shape() {return SHAPE::CYLINDER;}
};

class isosceles_tri_prism : public Obj
{
public:
    /**
     * @brief Constructs a isosceles_tri_prism
     * @details Constructs a isosceles_tri_prism using the Obj constructor
     *
     * @param mater A vector with all the material descriptors in it
     * @param geo A vector containing all the geometric properties of the isosceles_tri_prism
     * @param loc A vector for the center point location of the isosceles_tri_prism
     * @param ang orientation angle of the isosceles_tri_prism
     */
    isosceles_tri_prism(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);
    /**
     * @brief Copy Constructor
     *
     * @param o isosceles_tri_prism to be copied
     */
    isosceles_tri_prism(const isosceles_tri_prism &o);
    /**
     * @brief Determines if a point is inside the object
     * @details Determines if a real space point is within the isosceles_tri_prism
     *
     * @param v Point to check
     * @param dx spatial step of the Cell
     *
     * @return Returns true if v is in the object
     */
    bool isObj(std::array<double,3> &v, double dx);
    /**
     * @brief returns the Shape parameter of the isosceles_tri_prism
     * @details Returns the type of object it is by mapping the object name to a enum Shape
     * @return SHAPE::CONE
     */
    SHAPE shape() {return SHAPE::TRIANGLE_PRISM;}
};


class trapezoid_prism : public Obj
{
public:
    /**
     * @brief Constructs a isosceles_tri_prism
     * @details Constructs a isosceles_tri_prism using the Obj constructor
     *
     * @param mater A vector with all the material descriptors in it
     * @param geo A vector containing all the geometric properties of the isosceles_tri_prism
     * @param loc A vector for the center point location of the isosceles_tri_prism
     * @param ang orientation angle of the isosceles_tri_prism
     */
    trapezoid_prism(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);
    /**
     * @brief Copy Constructor
     *
     * @param o trapezoid_prism to be copied
     */
    trapezoid_prism(const trapezoid_prism &o);
    /**
     * @brief Determines if a point is inside the object
     * @details Determines if a real space point is within the trapezoid_prism
     *
     * @param v Point to check
     * @param dx spatial step of the Cell
     *
     * @return Returns true if v is in the object
     */
    bool isObj(std::array<double,3> &v, double dx);
    /**
     * @brief returns the Shape parameter of the trapezoid_prism
     * @details Returns the type of object it is by mapping the object name to a enum Shape
     * @return SHAPE::CONE
     */
    inline SHAPE shape() {return SHAPE::TRAPEZOIDAL_PRISM;}
};


class ters_tip : public Obj
{
protected:
    std::array<double,3> radCen_; //!< location of the center point of the circle
public:
    /**
     * @brief Constructs a isosceles_tri_prism
     * @details Constructs a isosceles_tri_prism using the Obj constructor
     *
     * @param mater A vector with all the material descriptors in it
     * @param geo A vector containing all the geometric properties of the isosceles_tri_prism
     * @param loc A vector for the center point location of the isosceles_tri_prism
     * @param ang orientation angle of the isosceles_tri_prism
     */
    ters_tip(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);
    /**
     * @brief Copy Constructor
     *
     * @param o ters_tip to be copied
     */
    ters_tip(const ters_tip &o);
    /**
     * @brief Determines if a point is inside the object
     * @details Determines if a real space point is within the ters_tip
     *
     * @param v Point to check
     * @param dx spatial step of the Cell
     *
     * @return Returns true if v is in the object
     */
    bool isObj(std::array<double,3> &v, double dx);
    /**
     * @brief returns the Shape parameter of the ters_tip
     * @details Returns the type of object it is by mapping the object name to a enum Shape
     * @return SHAPE::CONE
     */
    inline SHAPE shape() {return SHAPE::TERS_TIP;}
};

class parabolic_ters_tip : public Obj
{
public:
    /**
     * @brief Constructs a isosceles_tri_prism
     * @details Constructs a isosceles_tri_prism using the Obj constructor
     *
     * @param mater A vector with all the material descriptors in it
     * @param geo A vector containing all the geometric properties of the isosceles_tri_prism
     * @param loc A vector for the center point location of the isosceles_tri_prism
     * @param ang orientation angle of the isosceles_tri_prism
     */
    parabolic_ters_tip(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);
    /**
     * @brief Copy Constructor
     *
     * @param o parabolic_ters_tip to be copied
     */
    parabolic_ters_tip(const parabolic_ters_tip &o);
    /**
     * @brief Determines if a point is inside the object
     * @details Determines if a real space point is within the parabolic_ters_tip
     *
     * @param v Point to check
     * @param dx spatial step of the Cell
     *
     * @return Returns true if v is in the object
     */
    bool isObj(std::array<double,3> &v, double dx);
    /**
     * @brief returns the Shape parameter of the parabolic_ters_tip
     * @details Returns the type of object it is by mapping the object name to a enum Shape
     * @return SHAPE::CONE
     */
    inline SHAPE shape() {return SHAPE::PARABOLIC_TERS_TIP;}
};

#endif

