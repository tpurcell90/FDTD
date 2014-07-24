#include "pml.hpp"
#include "Source.hpp"
#include "enum.hpp"
//#include <math>

using namespace std;

UPML::UPML(int thickness, Direction d, int m, double R0, int nx, double dx, double dy)
{
	thickness_ = thickness;
	d_ = d;
	m_ = m;
	R0_ = R0;
	double sigmaMax = -(m_+1)*log(R0_)/(2*thickness_*dx); // eta should be in the denominator eta = sqrt(mu_0*Material(1,2)/epsilon_0/Material(1,1));
	double kappaMax = 1.0;
	for(int ii = 0; ii < thickness_; ii++)
	{
		kappa_.push_back(1 +(kappaMax - 1) * pow(static_cast<double>(ii)/static_cast<double>(thickness_),m_));
		sigma_.push_back(sigmaMax * pow(static_cast<double>(ii)/static_cast<double>(thickness_),m_));
	}
	pml_ = make_shared<Grid2D<double>>(nx, thickness_,dx,dy);
}
double UPML::F_out(int x, int y, std::shared_ptr<Grid2D<double>> field, Polarization pol)
{
	return 0.0;
}