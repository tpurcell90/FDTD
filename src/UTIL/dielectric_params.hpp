#ifndef FDTD_DIELEC_CONST
#define FDTD_DIELEC_CONST

#include <vector>
#include <cmath>
    // Source: Aleksandar D. Rakic´ , Aleksandra B. Djurisˇ ic´ , Jovan M. Elazar, and Marian L. Majewski 1998 Vol. 37, No. 22 y APPLIED OPTICS 5271
    /**
     * List of all the Dielectric parameters in Rakic´ Applied Optics 1998, just a library
     */
    const std::vector<double>      AG_MAT_ = {9.01 ,   0.845,0.048,0.0,   0.065,3.886,0.816,   0.124,0.452,4.481,    0.011,0.065,8.185,   0.840,0.916,9.083,   5.646,2.419,20.29};
    const std::vector<double>      AU_MAT_ = {9.03 ,   0.760,0.053,0.0,   0.024,0.241,0.415,   0.010,0.345,0.830,    0.071,0.870,2.969,   0.601,2.494,4.304,   4.384,2.214,13.32};
    const std::vector<double>      AL_MAT_ = {14.98,   0.523,0.047,0.0,   0.227,0.333,0.162,   0.050,0.312,1.544,    0.166,1.351,1.808,   0.030,3.382,3.473};
    const std::vector<double>      CU_MAT_ = {10.83,   0.575,0.030,0.0,   0.061,0.378,0.291,   0.104,1.056,2.957,    0.723,2.213,5.300,   0.638,4.305,11.18};
    const std::vector<double>      BE_MAT_ = {18.51,   0.084,0.035,0.0,   0.031,1.664,0.100,   0.140,3.395,1.032,    0.530,4.454,3.183,   0.130,1.802,4.604};
    const std::vector<double>      CR_MAT_ = {10.75,   0.168,0.047,0.0,   0.151,3.175,0.121,   0.150,1.305,0.543,    1.149,2.676,1.970,   0.825,1.335,8.775};
    const std::vector<double>      NI_MAT_ = {15.92,   0.096,0.048,0.0,   0.100,4.511,0.174,   0.135,1.334,0.582,    0.106,2.178,1.597,   0.729,6.292,6.089};
    const std::vector<double>      PT_MAT_ = {9.59 ,   0.333,0.080,0.0,   0.191,0.517,0.780,   0.659,1.838,1.314,    0.547,3.668,3.141,   3.576,8.517,9.249};
    const std::vector<double>      PD_MAT_ = {9.72 ,   0.330,0.008,0.0,   0.649,2.950,0.336,   0.121,0.555,0.501,    0.638,4.621,1.659,   0.453,3.236,5.715};
    const std::vector<double>      TI_MAT_ = {7.29 ,   0.148,0.082,0.0,   0.899,2.276,0.777,   0.393,2.518,1.545,    0.187,1.663,2.509,   0.001,1.762,19.43};
    const std::vector<double>       W_MAT_ = {13.22,   0.206,0.064,0.0,   0.054,0.530,1.004,   0.166,1.281,1.917,    0.706,3.332,3.580,   2.590,5.836,7.498};
    const std::vector<double> CDSE_QD_MAT_ = {2.03253, 0.65, 0.140,2.03253};
    const std::vector<double>  PBS_QD_MAT_ = {1.45864, 5.00, 0.26,1.45864};

#endif
