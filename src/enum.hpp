#ifndef FDTD_ENUMS
#define FDTD_ENUMS


namespace my_enums
{
    enum Polarization {Ex,Ey,Ez,Hx,Hy,Hz};
    enum OupuptsData {field, flux};
    enum ProfType {gaussian, continuous};
    enum Shape {sphere,block,ellipse,cone,cylinder};
}

#endif
