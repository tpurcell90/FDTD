#include "Obj.hpp"

Obj::Obj(const Obj &o) :
    coordTransform_(o.coordTransform_),
    material_(o.material_),
    magMaterial_(o.magMaterial_),
    geoParam_(o.geoParam_),
    location_(o.location_),
    alpha_(o.alpha_),
    xi_(o.xi_),
    gamma_(o.gamma_),
    magAlpha_(o.magAlpha_),
    magXi_(o.magXi_),
    magGamma_(o.magGamma_)
{}


Obj::Obj(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc, std::array<std::array<double,3>,3> unitVec) :
    geoParam_(geo),
    material_(mater),
    magMaterial_(magMater),
    location_(loc),
    alpha_(std::vector<double>((material_.size()-1)/3, 0.0)),
    xi_(std::vector<double>((material_.size()-1)/3, 0.0)),
    gamma_(std::vector<double>((material_.size()-1)/3, 0.0)),
    magAlpha_(std::vector<double>((magMaterial_.size()-1)/3, 0.0)),
    magXi_(std::vector<double>((magMaterial_.size()-1)/3, 0.0)),
    magGamma_(std::vector<double>((magMaterial_.size()-1)/3, 0.0))
{
    if(unitVec.size() != location_.size())
        throw std::logic_error("The size of the unit vector descriptor vector must be the same as the dimensionality of the system");
    for(int ii = 0; ii < 3; ++ii)
    {
        if(unitVec[ii].size() != unitVec.size())
            throw std::logic_error("The size of each unit vector must be the same as the dimensionality of the system");
        for(int jj = 0; jj < 3; ++jj)
        {
            std::array<double,3>cartCoord = {{ 0, 0, 0 }};
            cartCoord[jj] = 1.0;
            coordTransform_[(ii)*3+jj] = ddot_(unitVec[ii].size(), unitVec[ii].data(), 1, cartCoord.data(), 1) / ( std::sqrt( std::accumulate(unitVec[ii].begin(),unitVec[ii].end(),0.0, [](double a, double b){return a + pow(b,2.0); } ) ) * std::sqrt( std::accumulate(unitVec[ii].begin(),unitVec[ii].end(),0.0, [](double a, double b){return a + pow(b,2.0); } ) ) );
        }
    }
}

sphere::sphere(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc, std::array<std::array<double,3>,3> unitVec) :
    Obj(mater, magMater, geo, loc, unitVec)
{}
cylinder::cylinder(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc, std::array<std::array<double,3>,3> unitVec) :
    Obj(mater, magMater, geo, loc, unitVec)
{}
cone::cone(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc, std::array<std::array<double,3>,3> unitVec) :
    Obj(mater, magMater, geo, loc, unitVec)
{}
ellipsoid::ellipsoid(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc, std::array<std::array<double,3>,3> unitVec) :
    Obj(mater, magMater, geo, loc, unitVec)
{}
block::block(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc, std::array<std::array<double,3>,3> unitVec) :
    Obj(mater, magMater, geo, loc, unitVec)
{}
rounded_block::rounded_block(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc, std::array<std::array<double,3>,3> unitVec) :
    Obj(mater, magMater, geo, loc, unitVec)
{
    double radCurv = geoParam_[geoParam_.size()-1];
    curveCens_[0] = {     geoParam_[0]/2.0 - radCurv,      geoParam_[1]/2.0 - radCurv,      geoParam_[2]/2.0 - radCurv};
    curveCens_[1] = {-1.0*geoParam_[0]/2.0 + radCurv,      geoParam_[1]/2.0 - radCurv,      geoParam_[2]/2.0 - radCurv};
    curveCens_[2] = {     geoParam_[0]/2.0 - radCurv, -1.0*geoParam_[1]/2.0 + radCurv,      geoParam_[2]/2.0 - radCurv};
    curveCens_[3] = {-1.0*geoParam_[0]/2.0 + radCurv, -1.0*geoParam_[1]/2.0 + radCurv,      geoParam_[2]/2.0 - radCurv};
    curveCens_[4] = {     geoParam_[0]/2.0 - radCurv,      geoParam_[1]/2.0 - radCurv, -1.0*geoParam_[2]/2.0 + radCurv};
    curveCens_[5] = {-1.0*geoParam_[0]/2.0 + radCurv,      geoParam_[1]/2.0 - radCurv, -1.0*geoParam_[2]/2.0 + radCurv};
    curveCens_[6] = {     geoParam_[0]/2.0 - radCurv, -1.0*geoParam_[1]/2.0 + radCurv, -1.0*geoParam_[2]/2.0 + radCurv};
    curveCens_[7] = {-1.0*geoParam_[0]/2.0 + radCurv, -1.0*geoParam_[1]/2.0 + radCurv, -1.0*geoParam_[2]/2.0 + radCurv};
}
isosceles_tri_prism::isosceles_tri_prism(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc, std::array<std::array<double,3>,3> unitVec) :
    Obj(mater, magMater, geo, loc, unitVec)
{}
trapezoid_prism::trapezoid_prism(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc, std::array<std::array<double,3>,3> unitVec) :
    Obj(mater, magMater, geo, loc, unitVec)
{}
ters_tip::ters_tip(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc, std::array<std::array<double,3>,3> unitVec) :
    Obj(mater, magMater, geo, loc, unitVec)
{
    radCen_ = {{ 0.0, 0.0, -1.0*geo[2]/2.0 + geo[0]/2.0 }};
}
parabolic_ters_tip::parabolic_ters_tip(std::vector<double> mater, std::vector<double> magMater, std::vector<double> geo, std::array<double,3> loc, std::array<std::array<double,3>,3> unitVec) :
    Obj(mater, magMater, geo, loc, unitVec)
{}

sphere::sphere(const sphere &o) : Obj(o)
{}
cone::cone(const cone &o) : Obj(o)
{}
cylinder::cylinder(const cylinder &o) : Obj(o)
{}
block::block(const block &o) : Obj(o)
{}
rounded_block::rounded_block(const rounded_block &o) : Obj(o), curveCens_(o.curveCens_)
{}
ellipsoid::ellipsoid(const ellipsoid &o) : Obj(o)
{}
isosceles_tri_prism::isosceles_tri_prism(const isosceles_tri_prism &o) : Obj(o)
{}
trapezoid_prism::trapezoid_prism(const trapezoid_prism &o) : Obj(o)
{}
ters_tip::ters_tip(const ters_tip &o) : Obj(o), radCen_(o.radCen_)
{}
parabolic_ters_tip::parabolic_ters_tip(const parabolic_ters_tip &o) : Obj(o)
{}

void Obj::setUpConsts (double dt)
{
    for(int ii = 0; ii < (material_.size()-1)/3; ++ii)
    {
        double sig = material_[3*ii+1];
        double gam = material_[3*ii+2];
        double omg = material_[3*ii+3];
        alpha_[ii] = ((2-pow(omg*dt,2.0))   / (1+gam*dt));
           xi_[ii] = ((gam*dt -1)           / (1+gam*dt));
        gamma_[ii] = ((sig*pow(omg*dt,2.0)) / (1+gam*dt));
    }
    for(int ii = 0; ii < (magMaterial_.size()-1)/3; ++ii)
    {
        double sig = magMaterial_[3*ii+1];
        double gam = magMaterial_[3*ii+2];
        double omg = magMaterial_[3*ii+3];
        magAlpha_[ii] = ((2-pow(omg*dt,2.0))   / (1+gam*dt));
           magXi_[ii] = ((gam*dt -1)           / (1+gam*dt));
        magGamma_[ii] = ((sig*pow(omg*dt,2.0)) / (1+gam*dt));
    }
}

bool sphere::isObj(std::array<double,3> &v, double dx)
{
    if(dist(v,location_) > geoParam_[0] + dx/1.0e6)
    {
        return false;
    }
    return true;
}
bool block::isObj(std::array<double,3> &v, double dx)
{
    if(v.size() != location_.size())
        throw std::logic_error("When checking if a point is inside an object the point and object must have same dimensionality.");
    std::array<double,3> v_trans = {{0.0, 0.0, 0.0}};
    std::array<double,3> v_cen = {{0.0, 0.0, 0.0}};
    for(int ii = 0; ii < v_cen.size(); ++ii)
        v_cen[ii] = v[ii]-location_[ii];
    dgemv_('T', v_cen.size(), v_cen.size(), 1.0, coordTransform_.data(), v_cen.size(), v_cen.data(), 1, 0.0, v_trans.data(), 1 );
    for(int ii = 0; ii < v.size(); ii++)
    {
        if((v_trans[ii] > geoParam_[ii]/2.0 + dx/1.0e6) || (v_trans[ii] < -1.0*geoParam_[ii]/2.0 - dx/1.0e6))
            return false;
    }
    return true;
}
bool rounded_block::isObj(std::array<double,3> &v, double dx)
{
    if(v.size() != location_.size())
        throw std::logic_error("When checking if a point is inside an object the point and object must have same dimensionality.");
    std::array<double,3> v_trans = {{0.0, 0.0, 0.0}};
    std::array<double,3> v_cen = {{0.0, 0.0, 0.0}};
    for(int ii = 0; ii < v_cen.size(); ++ii)
        v_cen[ii] = v[ii]-location_[ii];
    dgemv_('T', v_cen.size(), v_cen.size(), 1.0, coordTransform_.data(), v_cen.size(), v_cen.data(), 1, 0.0, v_trans.data(), 1 );
    double radCurv = geoParam_.back();
    for(int ii = 0; ii < v.size(); ii++)
    {
        if((v_trans[ii] > geoParam_[ii]/2.0 + dx/1.0e6) || (v_trans[ii] < -1.0*geoParam_[ii]/2.0 - dx/1.0e6))
            return false;
    }
    if( ( (v_trans[0] > geoParam_[0]/2.0 - radCurv + dx/1.0e6) || (v_trans[0] < -1.0*geoParam_[0]/2.0 + radCurv - dx/1.0e6) ) && ( (v_trans[1] > geoParam_[1]/2.0 - radCurv + dx/1.0e6) || (v_trans[1] < -1.0*geoParam_[1]/2.0 + radCurv - dx/1.0e6) ) && ( (v_trans[2] > geoParam_[2]/2.0 - radCurv + dx/1.0e6) || (v_trans[2] < -1.0*geoParam_[2]/2.0 + radCurv - dx/1.0e6) ) )
    {
        for(int cc = 0; cc < curveCens_.size(); cc++)
            if(dist(v_trans,curveCens_[cc]) < radCurv + dx/1.0e6)
                return true;
        return false;
    }
    return true;
}
bool ellipsoid::isObj(std::array<double,3> &v, double dx)
{
    if(v.size() != location_.size())
        throw std::logic_error("When checking if a point is inside an object the point and object must have same dimensionality.");
    std::array<double,3> v_trans = {{0.0, 0.0, 0.0}};
    std::array<double,3> v_cen = {{0.0, 0.0, 0.0}};
    for(int ii = 0; ii < v_cen.size(); ++ii)
        v_cen[ii] = v[ii]-location_[ii];
    dgemv_('T', v_cen.size(), v_cen.size(), 1.0, coordTransform_.data(), v_cen.size(), v_cen.data(), 1, 0.0, v_trans.data(), 1 );
    double ptSum = 0.0;
    for( int ii = 0; ii < v_trans.size(); ++ii)
        ptSum += pow( (v_trans[ii]-dx/1.0e6) / (geoParam_[ii]/2.0), 2.0 );

    if(1.0 < ptSum )
        return false;
    return true;
}
bool isosceles_tri_prism::isObj(std::array<double,3> &v, double dx)
{
    if(v.size() != location_.size())
        throw std::logic_error("When checking if a point is inside an object the point and object must have same dimensionality.");
    std::array<double,3> v_trans = {{0.0, 0.0, 0.0}};
    std::array<double,3> v_cen = {{0.0, 0.0, 0.0}};
    for(int ii = 0; ii < v_cen.size(); ++ii)
        v_cen[ii] = v[ii]-location_[ii];
    dgemv_('T', v_cen.size(), v_cen.size(), 1.0, coordTransform_.data(), v_cen.size(), v_cen.data(), 1, 0.0, v_trans.data(), 1 );
    // The last compenent is the length of the prism
    if( (v_trans[1] < -1.0*geoParam_[1]/2.0 + dx*1e-6) || (v_trans[1] > geoParam_[1]/2.0 + dx*1e-6) )
        return false;
    else if( (v_trans[0] > -1.0*(v_trans[1]-geoParam_[1]/2.0) * geoParam_[0]/(2.0*geoParam_[1]) - dx*1e-6 ) || (v_trans[0] < (v_trans[1]-geoParam_[1]/2.0)*geoParam_[0]/(2.0*geoParam_[1]) + dx*1e-6 ) )
        return false;
    else if( (v_trans[2] < -1.0*geoParam_[2]/2.0 + dx*1e-6) || (v_trans[2] > geoParam_[2]/2.0 + dx*1e-6) )
        return false;
    return true;
}

bool trapezoid_prism::isObj(std::array<double,3> &v, double dx)
{
    if(v.size() != location_.size())
        throw std::logic_error("When checking if a point is inside an object the point and object must have same dimensionality.");
    std::array<double,3> v_trans = {{0.0, 0.0, 0.0}};
    std::array<double,3> v_cen = {{0.0, 0.0, 0.0}};
    for(int ii = 0; ii < v_cen.size(); ++ii)
        v_cen[ii] = v[ii]-location_[ii];
    dgemv_('T', v_cen.size(), v_cen.size(), 1.0, coordTransform_.data(), v_cen.size(), v_cen.data(), 1, 0.0, v_trans.data(), 1 );
    // Last component the length of the prism
    if( (v_trans[1] < -1.0*geoParam_[2]/2.0 - dx*1e-6) || (v_trans[1] > geoParam_[2]/2.0 + dx*1e-6) )
        return false;
    else if( ( v_trans[0]-dx*1.e-15 < (geoParam_[0]-geoParam_[1])/(2.0*geoParam_[2]) * v_trans[1] - (geoParam_[0] +geoParam_[1])/ 4.0 ) || ( v_trans[0]-dx*1.e-15 > -1.0 * (geoParam_[0]-geoParam_[1])/(2.0*geoParam_[2]) * v_trans[1] + (geoParam_[0] +geoParam_[1])/ 4.0 ) )
        return false;
    else if( (v_trans[2] < -1.0*geoParam_[3]/2.0 - dx*1e-6) || (v_trans[2] > geoParam_[3]/2.0 + dx*1e-6) )
        return false;
    return true;
}

bool cone::isObj(std::array<double,3> &v, double dx)
{
    if(v.size() != location_.size())
        throw std::logic_error("When checking if a point is inside an object the point and object must have same dimensionality.");
    std::array<double,3> v_trans = {{0.0, 0.0, 0.0}};
    std::array<double,3> v_cen = {{0.0, 0.0, 0.0}};
    for(int ii = 0; ii < v_cen.size(); ++ii)
        v_cen[ii] = v[ii]-location_[ii];
    dgemv_('T', v_cen.size(), v_cen.size(), 1.0, coordTransform_.data(), v_cen.size(), v_cen.data(), 1, 0.0, v_trans.data(), 1 );
    if( (v_trans[2] < -1.0*geoParam_[1]/2.0 - dx*1e-6) || (v_trans[2] > geoParam_[1]/2.0 + dx*1e-6) )
        return false;
    else if( ( pow( v_trans[0], 2.0 ) + pow( v_trans[1], 2.0 ) ) * pow( cos(geoParam_[0]), 2.0 ) - pow( v_trans[2] * sin(geoParam_[0]), 2.0) > dx*1e-6 )
        return false;
    return true;
}

bool cylinder::isObj(std::array<double,3> &v, double dx)
{
    if(v.size() != location_.size())
        throw std::logic_error("When checking if a point is inside an object the point and object must have same dimensionality.");
    std::array<double,3> v_trans = {{0.0, 0.0, 0.0}};
    std::array<double,3> v_cen = {{0.0, 0.0, 0.0}};
    for(int ii = 0; ii < v_cen.size(); ++ii)
        v_cen[ii] = v[ii]-location_[ii];
    dgemv_('T', v_cen.size(), v_cen.size(), 1.0, coordTransform_.data(), v_cen.size(), v_cen.data(), 1, 0.0, v_trans.data(), 1 );

    if( (v_trans[2] < -1.0*geoParam_[1]/2.0 - dx*1e-6) || (v_trans[2] > geoParam_[1]/2.0 + dx*1e-6) )
        return false;
    else if( dist(std::array<double,3>( {{ v_trans[0], v_trans[1], 0.0 }} ), std::array<double,3>({{0,0,0}}) ) > geoParam_[0] + 1.0e-6*dx )
        return false;
    return true;
}

bool parabolic_ters_tip::isObj(std::array<double,3> &v, double dx)
{
    if(v.size() != location_.size())
        throw std::logic_error("When checking if a point is inside an object the point and object must have same dimensionality.");
    std::array<double,3> v_trans = {{0.0, 0.0, 0.0}};
    std::array<double,3> v_cen = {{0.0, 0.0, 0.0}};
    for(int ii = 0; ii < v_cen.size(); ++ii)
        v_cen[ii] = v[ii]-location_[ii];
    dgemv_('T', v_cen.size(), v_cen.size(), 1.0, coordTransform_.data(), v_cen.size(), v_cen.data(), 1, 0.0, v_trans.data(), 1 );
    v_trans[1] += geoParam_[1]/2.0;
    if( (v_trans[2] > geoParam_[1] - dx*1e-6 ) || (v_trans[2] < dx*1e-6 + geoParam_[0] * pow( dist( std::array<double,3>( {{ v_trans[0], v_trans[1], 0 }} ), std::array<double,3>({{ 0, 0, 0 }}) ), 2.0 ) ) )
        return false;
    return true;
}


bool ters_tip::isObj(std::array<double,3> &v, double dx)
{
    if(v.size() != location_.size())
        throw std::logic_error("When checking if a point is inside an object the point and object must have same dimensionality.");
    std::array<double,3> v_trans = {{0.0, 0.0, 0.0}};
    std::array<double,3> v_cen = {{0.0, 0.0, 0.0}};
    for(int ii = 0; ii < v_cen.size(); ++ii)
        v_cen[ii] = v[ii]-location_[ii];
    dgemv_('T', v_cen.size(), v_cen.size(), 1.0, coordTransform_.data(), v_cen.size(), v_cen.data(), 1, 0.0, v_trans.data(), 1 );
    if( dist(v_trans, radCen_) < geoParam_[0]/2.0 )
        return true;

    if( (v_trans[2] < geoParam_[0] - geoParam_[1]/2.0 - dx*1e-6) || (v_trans[2] > geoParam_[1]/2.0 + dx*1e-6) )
        return false;
    else if( ( pow( v_trans[0], 2.0 ) + pow( v_trans[1], 2.0 ) ) * pow( cos(geoParam_[0]), 2.0 ) - pow( v_trans[2] * sin(geoParam_[0]), 2.0)  > dx*1e-6 )
        return false;

    return true;
}

double Obj::dist(std::array<double,3> pt1, std::array<double,3> pt2)
{
    double sum = 0;
    for(int cc = 0; cc < pt1.size(); cc ++)
        sum += pow((pt1[cc]-pt2[cc]),2);
    return sqrt(sum);
}