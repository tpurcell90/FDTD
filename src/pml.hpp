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
enum PMLTopBot {TOP, BOT};

template <typename T> class UPML
{
protected:
	int thickness_;
	Direction d_;
	//PMLTopBot tb_;
	int m_;
	double R0_;
	std::vector<double> kappa_;
	std::vector<double> sigma_;

public:
	std::shared_ptr<Grid2D<T>> Dx_,Dy_,Dz_,Bx_,By_,Bz_;

	UPML(int thickness, Direction d, int m, double R0, int nx, double dx, double dy, Polarization pol) : thickness_(thickness), d_(d), m_(m), R0_(R0)
	{
		double sigmaMax = -(m_+1)*log(R0_)/(2*thickness_*dx); // eta should be in the denominator eta = sqrt(mu_0*Material(1,2)/epsilon_0/Material(1,1));
		double kappaMax = 1.0;
		if(d == X)
		{
			if(pol == EX || pol == EY || pol == HZ)
			{
		        Dx_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
		        Dy_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
		        Bz_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
		        Bx_ = nullptr;
		        By_ = nullptr;
		        Dz_ = nullptr;
			}
			else
			{
		        Bx_ = std::make_shared<Grid2D<double>>(thickness,nx,dx,dy);
		        By_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
		        Dz_ = std::make_shared<Grid2D<double>>(thickness,nx-1,dx,dy);
		        Dx_ = nullptr;
		        Dy_ = nullptr;
		        Bz_ = nullptr;
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
			}
			else
			{
		        Bx_ = std::make_shared<Grid2D<double>>(nx-1,thickness,dx,dy);
		        By_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
		       	Dz_ = std::make_shared<Grid2D<double>>(nx,thickness,dx,dy);
		        Dx_ = nullptr;
		        Dy_ = nullptr;
		        Bz_ = nullptr;
			}

		}
		else
			 throw std::logic_error("While yes we could have a thrid dimension to run, I have yet to be implimented to do such a thing. So please accept this error as my sincerest appology.");
		for(int ii = 0; ii < thickness_; ii++)
		{
			//seperate out the points for the different fields
			kappa_.push_back(1 +(kappaMax - 1) * pow((static_cast<double>(thickness_) - static_cast<double>(ii))/static_cast<double>(thickness_),m_));
			sigma_.push_back(sigmaMax * pow((static_cast<double>(thickness_) - static_cast<double>(ii))/static_cast<double>(thickness_),m_));
		}
	}
	// Accessor Functions
	int thickness(){return thickness_;}
	Direction d(){return d_;}
	double kappa(int x){return kappa_[x];}
	double sigma(int x){return sigma_[x];}
	//PMLTopBot tb(){return tb_;}


/*k_Fz_1 = zeros(ny,1);   k_Fz_2 = zeros(ny,1);   
k_Ez_1 = zeros(nx,1);   k_Ez_2 = zeros(nx,1);   
k_Gx_1 = zeros(ny-1,1); k_Gx_2 = zeros(ny-1,1);
k_Hx_1 = zeros(nx,1);   k_Hx_2 = zeros(nx,1);   
k_Gy_1 = zeros(nx-1,1); k_Gy_2 = zeros(nx-1,1); 
k_Hy_1 = zeros(ny,1);   k_Hy_2 = zeros(ny,1);
// TE coefficients
k_Hz_1 = zeros(nx,1);   k_Hz_2 = zeros(nx,1);   
k_Ex_1 = zeros(nx,1);   k_Ex_2 = zeros(nx,1);   
k_Ey_1 = zeros(ny,1);   k_Ey_2 = zeros(ny,1);


//// Field transformation coefficients in PML areas for TM ans TE modes
// Along x-axis
sigma_max = -(m+1)*log(R_err)/(2*eta*PML*dx);
sigma_x = sigma_max*((PML-(1:PML)+1)/PML).^m;
ka_x = 1+(ka_max-1)*((PML-(1:PML)+1)/PML).^m;
// TM Modes
k_Ez_1(1:PML) = (2*epsilon_0*ka_x-sigma_x*dt)./(2*epsilon_0*ka_x+sigma_x*dt);
k_Ez_1(end-(1:PML)+1) = k_Ez_1(1:PML);
k_Ez_2(1:PML) = 2./(2*epsilon_0*ka_x + sigma_x*dt);
k_Ez_2(end-(1:PML)+1) = k_Ez_2(1:PML);
k_Hx_1(1:PML) = (2*epsilon_0*ka_x+sigma_x*dt)/(2*epsilon_0*mu_0);
k_Hx_1(end-(1:PML)+1) = k_Hx_1(1:PML);
k_Hx_2(1:PML) = (2*epsilon_0*ka_x-sigma_x*dt)/(2*epsilon_0*mu_0);
k_Hx_2(end-(1:PML)+1) = k_Hx_2(1:PML);

// TE Modes
k_Hz_1(1:PML) = (2*epsilon_0*ka_x-sigma_x*dt)./(2*epsilon_0*ka_x+sigma_x*dt);
k_Hz_1(end-(1:PML)+1) = k_Hz_1(1:PML);
k_Hz_2(1:PML) = 2*epsilon_0./(2*epsilon_0*ka_x+sigma_x*dt)/mu_0;
k_Hz_2(end-(1:PML)+1) = k_Hz_2(1:PML);
k_Ex_1(1:PML) = (2*epsilon_0*ka_x+sigma_x*dt)/epsilon_0;
k_Ex_1(end-(1:PML)+1) = k_Ex_1(1:PML);
k_Ex_2(1:PML) = (2*epsilon_0*ka_x-sigma_x*dt)/epsilon_0;
k_Ex_2(end-(1:PML)+1) = k_Ex_2(1:PML);

// TM  Modes
sigma_x = sigma_max*((PML-(1:PML)+0.5)/PML).^m;
ka_x = 1+(ka_max-1)*((PML-(1:PML)+0.5)/PML).^m;
k_Gy_1(1:PML) = (2*epsilon_0*ka_x-sigma_x*dt)./(2*epsilon_0*ka_x+sigma_x*dt);
k_Gy_1(end-(1:PML)+1) = k_Gy_1(1:PML);
k_Gy_2(1:PML) = 2*epsilon_0*dt./(2*epsilon_0*ka_x+sigma_x*dt)/dx;
k_Gy_2(end-(1:PML)+1) = k_Gy_2(1:PML);

// Along y-axis
sigma_max = -(m+1)*log(R_err)/(2*eta*PML*dy);
sigma_y = sigma_max*((PML-(1:PML)+1)/PML).^m;
ka_y = 1+(ka_max-1)*((PML-(1:PML)+1)/PML).^m;
// TM Modes
k_Fz_1(1:PML) = (2*epsilon_0*ka_y-sigma_y*dt)./(2*epsilon_0*ka_y+sigma_y*dt);
k_Fz_1(end-(1:PML)+1) = k_Fz_1(1:PML);
k_Fz_2(1:PML) = 2*epsilon_0*dt./(2*epsilon_0*ka_y+sigma_y*dt);
k_Fz_2(end-(1:PML)+1) = k_Fz_2(1:PML);
k_Hy_1(1:PML) = (2*epsilon_0*ka_y+sigma_y*dt)/(2*epsilon_0*mu_0);
k_Hy_1(end-(1:PML)+1) = k_Hy_1(1:PML);
k_Hy_2(1:PML) = (2*epsilon_0*ka_y-sigma_y*dt)/(2*epsilon_0*mu_0);
k_Hy_2(end-(1:PML)+1) = k_Hy_2(1:PML);

// TE Modes
k_Ey_1(1:PML) = (2*epsilon_0*ka_y+sigma_y*dt)/epsilon_0;
k_Ey_1(end-(1:PML)+1) = k_Ey_1(1:PML);
k_Ey_2(1:PML) = (2*epsilon_0*ka_y-sigma_y*dt)/epsilon_0;
k_Ey_2(end-(1:PML)+1) = k_Ey_2(1:PML);

//TM Modes
sigma_y = sigma_max*((PML-(1:PML)+0.5)/PML).^m;
ka_y = 1+(ka_max-1)*((PML-(1:PML)+0.5)/PML).^m;
k_Gx_1(1:PML) = (2*epsilon_0*ka_y-sigma_y*dt)./(2*epsilon_0*ka_y+sigma_y*dt);
k_Gx_1(end-(1:PML)+1) = k_Gx_1(1:PML);
k_Gx_2(1:PML) = 2*epsilon_0*dt./(2*epsilon_0*ka_y+sigma_y*dt)/dy;
k_Gx_2(end-(1:PML)+1) = k_Gx_2(1:PML);


*/

};

#endif