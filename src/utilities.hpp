/**
 * \file "utilities.hpp"
 * \author  J. Szekely
 */

#ifndef FDTD_UTILITIES
#define FDTD_UTILITIES

#include <complex>
#include <memory>

//BLAS
extern "C"
{
 void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
    const double* alpha, const double* a, const int* lda, const double* b, const int* ldb,
    const double* beta, double* c, const int* ldc);

  void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*);

  void zheev_(const char*, const char*, const int*, std::complex<double>*, const int*, double*, std::complex<double>*, const int*, double*, int*);

  double ddot_(const int*, const double*, const int*, const double*, const int*);

  void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);

  void zaxpy_(const int*, const std::complex<double>*, const std::complex<double>*, const int*, const std::complex<double>*, const int*);

  void zcopy_(const int*, const std::complex<double>*, const int*, std::complex<double>*, const int*);

  void zgemm3m_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
               const std::complex<double>* alpha, const std::complex<double>* a, const int* lda, const std::complex<double>* b, const int* ldb,
               const std::complex<double>* beta, std::complex<double>* c, const int* ldc);

  void zgemv_(const char*, const int*, const int*, const std::complex<double>*, const std::complex<double>*, const int*, const std::complex<double>*, const int*,
             const std::complex<double>*, std::complex<double>*, const int*);
  int izamax_(const int*, const std::complex<double>*, const int*);

  int izamin_(const int*, const std::complex<double>*, const int*);

}

//LAPACK
extern "C"
{
  void dgesvd_(const char*, const char*, const int*, const int*, double*, const int*, double*, double*, const int*, double*, const int*,  double*, const int*, int*);

  void dsyevr_(const char*, const char*, const char*, const int*, double*, const int*, const double*, const double*, const int*, const int*, const double*,  int*, double*, double*, const int*, int*, double*, int*, int*, int*, int*);
}

//AlignmentTool interface
namespace
{
  void dgemm_(const char* transa, const char* transb, const int m, const int n, const int k,
              const double alpha, const double* a, const int lda, const double* b, const int ldb,
              const double beta, double* c, const int ldc)
              { ::dgemm_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }

  void dgemm_(const char* transa, const char* transb, const int m, const int n, const int k,
              const double alpha, const std::unique_ptr<double []>& a, const int lda, const std::unique_ptr<double []>& b, const int ldb,
              const double beta, std::unique_ptr<double []>& c, const int ldc)
             { ::dgemm_(transa,transb,&m,&n,&k,&alpha,a.get(),&lda,b.get(),&ldb,&beta,c.get(),&ldc); }


  void dsyev_(const char* a, const char* b, const int c, double* d, const int e, double* f, double* g, const int h, int& i)
             { ::dsyev_(a,b,&c,d,&e,f,g,&h,&i);}

  void dsyev_(const char* a, const char* b, const int c, std::unique_ptr<double []>& d, const int e,
             std::unique_ptr<double []>& f, std::unique_ptr<double []>& g, const int h, int& i)
             { ::dsyev_(a,b,&c,d.get(),&e,f.get(),g.get(),&h,&i);}

  void dgesvd_(const char* a, const char* b, const int c, const int d, double* e, const int f, double* g, double* h, const int i, double* j, const int k,
              double* l, const int m, int& n) { ::dgesvd_(a,b,&c,&d,e,&f,g,h,&i,j,&k,l,&m,&n); }


  void mkl_domatcopy_(const char* ordering, const char* trans, const int r, const int c, const double alpha,
                    const double* A, const int nr, double* B, const int nc)
                    {mkl_domatcopy_(ordering,trans,r,c,alpha,A,nr,B,nc);}

  double ddot_(const int a, const double* b, const int c, const double* d, const int e) { return ::ddot_(&a,b,&c,d,&e);}

  double ddot_(const int a, const std::unique_ptr<double []>& b, const int c, const std::unique_ptr<double []>& d, const int e)
             { return ::ddot_(&a,b.get(),&c,d.get(),&e); }

  void dsyevr_(const char* jobz, const char* range, const char* uplo, const int n, double* a, const int lda, double vl,
                 double vu, const int il, const int iu, const double abstol, int m, double* w, double* z, const int ldz,
                int* isuppz, double* work, int lwork, int* iwork, int liwork, int& info)
                { ::dsyevr_(jobz, range, uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);}

  void daxpy_(const int a, const double b, const double* c, const int d, double* e, const int f) { ::daxpy_(&a,&b,c,&d,e,&f); }

  void zaxpy_(const int a, const std::complex<double> b, const std::complex<double>* c, const int d, std::complex<double>* e, const int f) { ::zaxpy_(&a,&b,c,&d,e,&f); }

  void zcopy_(const int a, const std::complex<double>* b, const int c, std::complex<double>* d, const int e) { ::zcopy_(&a,b,&c,d,&e);}

  void zgemm3m_(const char* transa, const char* transb, const int m, const int n, const int k,
                 const std::complex<double> alpha, const std::complex<double>* a, const int lda, const std::complex<double>* b, const int ldb,
                 const std::complex<double> beta, std::complex<double>* c, const int ldc) { ::zgemm3m_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }

  void zgemm3m_(const char* transa, const char* transb, const int m, const int n, const int k,
                 const std::complex<double> alpha, const std::unique_ptr<std::complex<double>[]>& a, const int lda, const std::unique_ptr<std::complex<double>[]>& b, const int ldb,
                 const std::complex<double> beta, std::unique_ptr<std::complex<double>[]>& c, const int ldc)
                 { ::zgemm3m_(transa,transb,&m,&n,&k,&alpha,a.get(),&lda,b.get(),&ldb,&beta,c.get(),&ldc); }

  void zgemv_(const char* a, const int b, const int c, const std::complex<double> d, const std::complex<double>* e, const int f, const std::complex<double>* g, const int h,
             const std::complex<double> i, std::complex<double>* j, const int k) { ::zgemv_(a,&b,&c,&d,e,&f,g,&h,&i,j,&k); }

  void zgemv_(const char* a, const int b, const int c, const std::complex<double> d, const std::unique_ptr<std::complex<double> []>& e, const int f,
             const std::unique_ptr<std::complex<double> []>& g, const int h, const std::complex<double> i, std::unique_ptr<std::complex<double> []>& j, const int k)
             { ::zgemv_(a,&b,&c,&d,e.get(),&f,g.get(),&h,&i,j.get(),&k); }

  int izamax_(const int n, const std::complex<double> *x, const int inc)
              { return ::izamax_(&n,x,&inc);}

  int izamin_(const int n, const std::complex<double> *x, const int inc)
              { return ::izamin_(&n,x,&inc);}

  void zheev_(const char* a, const char* b, const int c, std::complex<double>* d, const int e, double* f, std::complex<double>* g, const int h, double* i, int& j)
             { ::zheev_(a,b,&c,d,&e,f,g,&h,i,&j); }

  void zheev_(const char* a, const char* b, const int c, std::unique_ptr<std::complex<double> []>& d, const int e,
             std::unique_ptr<double []>& f, std::unique_ptr<std::complex<double> []>& g, const int h, std::unique_ptr<double[]>& i, int& j)
             { ::zheev_(a,b,&c,d.get(),&e,f.get(),g.get(),&h,i.get(),&j); }

// Sparse routines
 void mkl_dcsrsymv_(const char* uplo, const int m, const double* a, const int* ia, const int* ja, double* x, double* y)
                {mkl_dcsrsymv_(uplo, m, a, ia, ja, x, y); };

 void mkl_dscrgemv_(const char* transa, const int m, const double* a, const int* ia, const int* ja, double* x, double* y)
                {mkl_dscrgemv_(transa, m, a, ia, ja, x, y); };

}

#endif