#ifndef FDTD_OBJECT
#define FDTD_OBJECT

#include <math.h>
#include <UTIL/typedefs.hpp>

// #include <iostream>

class Obj
{
protected:
    std::array<std::array<double,3>,3> unitVec_; //!< vectors describing object's coordinate sys relative to main grids
    std::vector<double> geoParam_; //!< parameters describing the geometry of the objects
    std::vector<double> material_; //!< parameters for storing the material properties of the system
    std::vector<double> magMaterial_; //!< parameters for storing the material properties of the system
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
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       List of parameters describing the object geometry
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises
     */
    Obj(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc, std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     Object to be copied
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

    /**
     * @return     returns $\mu_{\infty}$
     */
    inline double muInfty(){return magMaterial_[0];}

    /**
     * @return     the unit vectors of the object's coordinate system
     */
    inline std::array<std::array<double,3>,3> unitVec() {return unitVec_;}

    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    virtual bool isObj(std::array<double,3> v, double dx) = 0;

    /**
     * @brief      Returns the distance between two points
     *
     * @param[in]  pt1   The point 1
     * @param[in]  pt2   The point 2
     *
     * @return     The distance between the points
     */
    double dist(std::array<double,3> pt1, std::array<double,3> pt2);

    /**
     * @return     the shape of the object
     */
    virtual SHAPE shape() = 0;
};

class sphere : public Obj
{
public:
    /**
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       Vector containing {radius}
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises
     */
    sphere(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     sphere to be copied
     */
    sphere(const sphere &o);

    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    bool isObj(std::array<double,3> v, double dx);

    /**
     * @return     SHAPE::SPHERE
     */
    SHAPE shape() {return SHAPE::SPHERE;}
};

class hemisphere : public Obj
{
public:
    /**
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       Vector containing {radius}
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises; The first vector determines the direction of the normal vector of the symmetry plane
     */
    hemisphere(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     hemisphere to be copied
     */
    hemisphere(const hemisphere &o);

    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    bool isObj(std::array<double,3> v, double dx);

    /**
     * @return     SHAPE::HEMISPHERE
     */
    SHAPE shape() {return SHAPE::HEMISPHERE;}
};

class block : public Obj
{
public:
    /**
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       Vector Containing {edge size in direction of unit vec 0, edge size in direction of unit vec 1, edge size in direction of unit vec 2}
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises
     */
    block(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     Block to be copied
     */
    block(const block &o);

    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    bool isObj(std::array<double,3> v, double dx);

    /**
     * @return     SHAPE::BLOCK
     */
    SHAPE shape() {return SHAPE::BLOCK;}

};
class rounded_block : public Obj
{
protected:
    std::array<std::array<double,3>,8> curveCens_;
public:
    /**
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       Vector Containing {edge size in direction of unit vec 0, edge size in direction of unit vec 1, edge size in direction of unit vec 2, radius of curvature of the corners}
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises
     */
    rounded_block(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     rounded_block to be copied
     */
    rounded_block(const rounded_block &o);

    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    bool isObj(std::array<double,3> v, double dx);

    /**
     * @return     SHAPE::ROUNDED_BLOCK
     */
    SHAPE shape() {return SHAPE::ROUNDED_BLOCK;}

};
class ellipsoid : public Obj
{
public:
    /**
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       Vector Containing {axis length in direction of unit vec 0, axis length in direction of unit vec 1, axis length in direction of unit vec 2}
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises
     */
    ellipsoid(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     ellipsoid to be copied
     */
    ellipsoid(const ellipsoid &o);
    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    bool isObj(std::array<double,3> v, double dx);

    /**
     * @return     SHAPE::ELLIPSOID
     */
    SHAPE shape() {return SHAPE::ELLIPSOID;}

};

class hemiellipsoid : public Obj
{
public:
    /**
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       Vector Containing {axis length in direction of unit vec 0, axis length in direction of unit vec 1, axis length in direction of unit vec 2}
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises; The first axis is the direction of the normal vector to the symmetry plane of the hemiellipsoid
     */
    hemiellipsoid(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     hemiellipsoid to be copied
     */
    hemiellipsoid(const hemiellipsoid &o);
    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    bool isObj(std::array<double,3> v, double dx);

    /**
     * @return     SHAPE::HEMIELLIPSOID
     */
    SHAPE shape() {return SHAPE::HEMIELLIPSOID;}

};
class cone : public Obj
{
public:

    /**
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       Vector Containing {radius of base 1, radius of base 2, length of the cone}
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises
     */
    cone(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     Cone to be copied
     */
    cone(const cone &o);

    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    bool isObj(std::array<double,3> v, double dx);

    /**
     * @return     SHAPE::CONE
     */
    SHAPE shape() {return SHAPE::CONE;}
};

class cylinder : public Obj
{
public:
    /**
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       Vector Containing {radius of the cylinder, length of the cylinder}
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises
     */
    cylinder(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     Cylinder to be copied
     */
    cylinder(const cylinder &o);

    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    bool isObj(std::array<double,3> v, double dx);

    /**
     * @return     SHAPE::CYLINDER
     */
    SHAPE shape() {return SHAPE::CYLINDER;}
};

class isosceles_tri_prism : public Obj
{
public:
    /**
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       Vector Containing { base of the triangle base, height of the triangle base, length of the prism}
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises
     */
    isosceles_tri_prism(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     isosceles triangular prism to be copied
     */
    isosceles_tri_prism(const isosceles_tri_prism &o);

    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    bool isObj(std::array<double,3> v, double dx);

    /**
     * @return     SHAPE::TRIANGLE_PRISM
     */
    SHAPE shape() {return SHAPE::TRIANGLE_PRISM;}
};


class trapezoid_prism : public Obj
{
public:
    /**
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       Vector Containing { base1 of the trapezoid base, base2 of the trapezoid base, height of the trapezoid base, length of the prism}
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises
     */
    trapezoid_prism(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     trapezoidal prism to be copied
     */
    trapezoid_prism(const trapezoid_prism &o);

    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    bool isObj(std::array<double,3> v, double dx);

    /**
     * @return     SHAPE::TRAPEZOIDAL_PRISM
     */
    inline SHAPE shape() {return SHAPE::TRAPEZOIDAL_PRISM;}
};


class ters_tip : public Obj
{
protected:
    std::array<double,3> radCen_; //!< location of the center point of the circle
public:
    /**
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       Vector Containing { radius of curvature of the tip point, radius of the base of the tip, tip length}
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises
     */
    ters_tip(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     TERS Tip to be copied
     */
    ters_tip(const ters_tip &o);

    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    bool isObj(std::array<double,3> v, double dx);

    /**
     * @return     SHAPE::TERS_TIP
     */
    inline SHAPE shape() {return SHAPE::TERS_TIP;}
};

class parabolic_ters_tip : public Obj
{
public:
    /**
     * @brief      Constructor
     *
     * @param[in]  mater     Vector containing all electric field dispersive material parameters
     * @param[in]  magMater  Vector containing all magnetic field dispersive material parameters
     * @param[in]  geo       Vector Containing { radius of curvature of the tip, the tip length }
     * @param[in]  loc       The location of the center of the object
     * @param[in]  unitVec   An array of vectors describing the coordinate transform to make the objects orientation in the grids along the x,y,z axises
     */
    parabolic_ters_tip(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc,std::array<std::array<double,3>,3> unitVec);

    /**
     * @brief      Copy Constructor
     *
     * @param[in]  o     Parabolic TERS Tip to be copied
     */
    parabolic_ters_tip(const parabolic_ters_tip &o);

    /**
     * @brief      Determines if a point is inside the object
     *
     * @param[in]  v     real space point
     * @param[in]  dx    grid spacing
     *
     * @return     True if point is in the object, False otherwise.
     */
    bool isObj(std::array<double,3> v, double dx);

    /**
     * @return     SHAPE::PARABOLIC_TERS_TIP
     */
    inline SHAPE shape() {return SHAPE::PARABOLIC_TERS_TIP;}
};

#endif

