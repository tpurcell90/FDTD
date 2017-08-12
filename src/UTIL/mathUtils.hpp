#ifndef SRC_MATHUTILS_HPP
#define SRC_MATHUTILS_HPP

#include <type_traits>
#include <stdexcept>
#include <cassert>
#ifdef MKL
#include <UTIL/utilities_MKL.hpp>
#elif defined ACML
#include <UTIL/utilities_acml.hpp>
#else
#include <UTIL/utilities_MKL.hpp>
#endif

template <typename T> T conj(const T& a) { throw std::logic_error("conj");}

template<> inline double conj<double> (const double& a) { return a; }

template<> inline cplx conj(const cplx& a) { return std::conj(a); }


#endif