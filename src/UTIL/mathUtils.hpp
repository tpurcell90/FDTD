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

// AXPY
// template<class T, class U, typename Type,
// class = typename std::enable_if <(std::is_pointer<T>::value) && (std::is_pointer<U>::value) >::type >
// void ax_plus_y_n(const Type& a, const T p, const size_t n, U q)
// {
//   std::transform(p, p+n, q, q, [&a](decltype(*p) i, decltype(*q) j) { return j+a*i; });
// }

// template<> inline void ax_plus_y_n(const double& a, const double* p, const size_t n, double* q)
// {
//   daxpy_(n, a, p, 1, q, 1);
// }

// template<> inline void ax_plus_y_n(const cplx& a, const cplx* p, const size_t n, cplx* q)
// {
//   zaxpy_(n, a, p, 1, q, 1);
// }

// template<> inline void ax_plus_y_n(const double& a, const cplx* p, const size_t n, cplx* q)
// {
//   zaxpy_(n, static_cast<cplx>(a), p, 1, q, 1);
// }

// DOT
// template<class T, class U,
//          class = typename std::enable_if< (std::is_pointer<T>::value) &&
//                                           (std::is_pointer<U>::value) >::type >
// auto dot_product_(const T p, const size_t n, const U q) -> decltype(*p * *q)
// {
//   using ResultType = decltype(*p * *q);
//   return std::inner_product(p, p+n, q, static_cast<ResultType>(0.0), std::plus<ResultType>(), [](decltype(*p) i, decltype(*q) j) { return conj(i)*j; });
// }

// template<> inline double dot_product_(const double* p, const size_t n, const double* q)
// {
//   return ddot_(n, p, 1, q, 1);
// }

// template<> inline cplx dot_product_(const cplx* p, const size_t n, const cplx* q)
// {
//   return zdotc_(n, p, 1, q, 1);
// }

#endif