#ifndef FDTD_ENUMS
#define FDTD_ENUMS

namespace my_enums
{
    enum Polarization {EX,EY,EZ,HX,HY,HZ};
    enum dtcOutType {field, flux};
    enum plsShape {gaussian, continuous};
    enum Shape {sphere,block,ellipse,cone,cylinder};
    enum PMLTopBot {TOP, BOT};
}
#endif
