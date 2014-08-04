#ifndef FDTD_PML
#define FDTD_PML

#include "Source.hpp"
#include "Grid.hpp"
// #include <assert.h>
// #include <iomanip>
// #include <iostream>
#include <memory>
// #include <random>
// #include <stdexcept>
// #include <string>
#include <vector>
#include <map>

// #include <complex>
// typedef std::complex<double> cplx;
enum Direction {X,Y,Z};

template <typename T> class UPML
{
protected:
	int thickness_;
	Direction d_;
	//PMLTopBot tb_;
	double m_;
	double R0_;
	std::vector<double> kappa_;
	std::vector<double> sigma_;
	double sigmaMax_;
	double kappaMax_;

public:
	std::shared_ptr<Grid2D<T>> Dx_,Dy_,Dz_,Bx_,By_,Bz_,Dx_end_,Dy_end_,Dz_end_,Bx_end_,By_end_,Bz_end_;

	UPML(int thickness, Direction d, double m, double R0, int nx, double dx, double dy, Polarization pol) : thickness_(thickness), d_(d), m_(m), R0_(R0)
	{
		//double sigma0 = -log(R0_)*log(g_)/(2.0*dx*pow(g,static_cast<double>(thickness_)/dx) - 1.00); // eta should be in the denominator eta = sqrt(mu_0*Material(1,2)/epsilon_0/Material(1,1));
		//double sigmaMax_ = -(m_+1)*log(R0_)/(2*thickness_*dx); // eta should be included
		//double sigmaMax_ = (m+1.0)/(150.0*M_PI*dx);
		sigmaMax_ = 0.8*(m_+1)/dx;
		kappaMax_ = 1.0;
		if(d == X)
		{
			if(pol == EX || pol == EY || pol == HZ)
			{
		        Dx_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
		        Dy_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
		        Bz_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
		        Dx_end_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
		        Dy_end_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
		        Bz_end_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
		        Bx_ = nullptr;
		        By_ = nullptr;
		        Dz_ = nullptr;
		        Bx_end_ = nullptr;
		        By_end_ = nullptr;
		        Dz_end_ = nullptr;
			}
			else
			{
		        Bx_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
		        By_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
		        Dz_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
		        Dx_ = nullptr;
		        Dy_ = nullptr;
		        Bz_ = nullptr;
		        Bx_end_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
		        By_end_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
		        Dz_end_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
		        Dx_end_ = nullptr;
		        Dy_end_ = nullptr;
		        Bz_end_ = nullptr;
			}
		}
		else if (d == Y)
		{
			if(pol == EX || pol == EY || pol == HZ)
			{
		        Dx_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
		        Dy_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
		        Bz_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
		        Bx_ = nullptr;
		        By_ = nullptr;
		        Dz_ = nullptr;
		        Dx_end_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
		        Dy_end_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
		        Bz_end_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
		        Bx_end_ = nullptr;
		        By_end_ = nullptr;
		        Dz_end_ = nullptr;
			}
			else
			{
		        Bx_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
		        By_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
		       	Dz_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
		        Dx_ = nullptr;
		        Dy_ = nullptr;
		        Bz_ = nullptr;
		        Bx_end_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
		        By_end_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
		       	Dz_end_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
		        Dx_end_ = nullptr;
		        Dy_end_ = nullptr;
		        Bz_end_ = nullptr;
			}

		}
		else
			 throw std::logic_error("While yes we could have a thrid dimension to run, I have yet to be implimented to do such a thing. So please accept this error as my sincerest appology.");
	}
	// Accessor Functions
	int thickness(){return thickness_;}
	Direction d(){return d_;}
	double kappa(int x){return kappaMax_;}
	double sigma(double x)
	{
		if(x <= thickness_ -1)
			return sigmaMax_ * pow((static_cast<double>(thickness_-1) - x) / static_cast<double>(thickness_-1) , m_);
		else
			return 0.0;
	}
};

#endif