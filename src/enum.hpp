#ifndef FDTD_ENUMS
#define FDTD_ENUMS


namespace my_enums
{
    enum Polarization {EX,EY,EZ,HX,HY,HZ};
    enum OupuptsData {field, flux};
    enum ProfType {gaussian, continuous};
    enum Shape {sphere,block,ellipse,cone,cylinder};
}

#endif
