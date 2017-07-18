/**
 * \file "utilities.hpp"
 * \author  J. Szekely
 */

// To see documentation look up functions here: https://software.intel.com/en-us/mkl-reference-manual-for-c
#ifndef FDTD_UTILITIES
#define FDTD_UTILITIES

#include <UTIL/typedefs.hpp>

//BLAS
extern "C"
{
    void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha, const double* a, const int* lda, const double* b, const int* ldb, const double* beta, double* c, const int* ldc);

    void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*);

    void zheev_(const char*, const char*, const int*, cplx*, const int*, double*, cplx*, const int*, double*, int*);

    double ddot_(const int*, const double*, const int*, const double*, const int*);

    void saxpy_(const int*, const int*, const int*, const int*, int*, const int*);

    void scopy_(const int*, const int*, const int*, int*, const int*);

    void sscal_(const int*, const int*, int*, const int*);

    void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);

    void dcopy_(const int*, const double*, const int*, double*, const int*);

    void dscal_(const int*, const double*, double*, const int*);

    void mkl_domatcopy_(const char*, const char*, const int *, const int *, const double *, const double* , const int *, double* , const int *);

    void zaxpy_(const int*, const cplx*, const cplx*, const int*, const cplx*, const int*);

    void zcopy_(const int*, const cplx*, const int*, cplx*, const int*);

    void zscal_(const int*, const cplx*, cplx*, const int*);

    void zgemm3m_(const char* transa, const char* transb, const int* m, const int* n, const int* k, const cplx* alpha, const cplx* a, const int* lda, const cplx* b, const int* ldb, const cplx* beta, cplx* c, const int* ldc);

    void zhemm_(const char* side, const char* uplo, const int* m, const int* n, const cplx *alpha, const cplx *a, const int* lda, const cplx *b, const int* ldb, const cplx *beta, const cplx *c, const int* ldc);

    void zgemv_(const char*, const int*, const int*, const cplx*, const cplx*, const int*, const cplx*, const int*, const cplx*, cplx*, const int*);

    void mkl_zomatcopy_(const char*, const char*, const int *, const int *, const cplx *, const cplx* , const int *, cplx* , const int *);

    int izamax_(const int*, const cplx*, const int*);

    int izamin_(const int*, const cplx*, const int*);

    int idamax_(const int*, const double*, const int*);

    int idamin_(const int*, const double*, const int*);

    int isamax_(const int*, const int*, const int*);

    int isamin_(const int*, const int*, const int*);

    #ifndef ZDOT_RETURN
        void zdotc_(cplx*, const int*, const cplx*, const int*, const cplx*, const int*);
    #else
        cplx zdotc_(const int*, const cplx*, const int*, const cplx*, const int*);
    #endif
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
    /**
     * @brief      Wrapper for dgemm_: c<-$\alpha a*b + /beta*c
     *
     * @param[in]  transa  'T' or 't' for Transpose, 'C' or c' for Hermitian Conjugate, 'N' or 'n' for neither on Matrix A
     * @param[in]  transb  'T' or 't' for Transpose, 'C' or c' for Hermitian Conjugate, 'N' or 'n' for neither on Matrix B
     * @param[in]  m       Specifies the number of rows of the matrix op(A) and of the matrix C. The value of m must be at least zero.
     * @param[in]  n       Specifies the number of columns of the matrix op(A) and the number of rows of the matrix op(B).
     * @param[in]  k       Specifies the number of columns of the matrix op(A) and the number of rows of the matrix op(B).
     * @param[in]  alpha   Specifies the scalar alpha.
     * @param[in]  a       Array size lda by ka, where ka is k when transa = 'N' or 'n', and is m otherwise. Before entry with transa = 'N' or 'n', the leading m-by-k part of the array a must contain the matrix A, otherwise the leading k-by-m part of the array a must contain the matrix A.
     * @param[in]  lda     specifies the leading dimension of a as declared in the calling (sub)program.
     * @param[in]  b       Array, size ldb by kb, where kb is n when transa = 'N' or 'n', and is k otherwise. Before entry with transa = 'N' or 'n', the leading k-by-n part of the array b must contain the matrix B, otherwise the leading n-by-k part of the array b must contain the matrix B.
     * @param[in]  ldb     specifies the leading dimension of a as declared in the calling (sub)program.
     * @param[in]  beta    Specifies the scalar beta.
     * @param[in]  c       Array, size ldc by n. Before entry, the leading m-by-n part of the array c must contain the matrix C, except when beta is equal to zero, in which case c need not be set on entry.
     * @param[in]  ldc     Specifies the leading dimension of c as declared in the calling (sub)program.
     */
    void dgemm_(const char* transa, const char* transb, const int m, const int n, const int k, const double alpha, const double* a, const int lda, const double* b, const int ldb, const double beta, double* c, const int ldc)
        { ::dgemm_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }
    /**
     * @brief      Wrapper for dgemm_: c<-$\alpha a*b + /beta*c
     *
     * @param[in]  transa  'T' or 't' for Transpose, 'C' or c' for Hermitian Conjugate, 'N' or 'n' for neither on Matrix A
     * @param[in]  transb  'T' or 't' for Transpose, 'C' or c' for Hermitian Conjugate, 'N' or 'n' for neither on Matrix B
     * @param[in]  m       Specifies the number of rows of the matrix op(A) and of the matrix C. The value of m must be at least zero.
     * @param[in]  n       Specifies the number of columns of the matrix op(A) and the number of rows of the matrix op(B).
     * @param[in]  k       Specifies the number of columns of the matrix op(A) and the number of rows of the matrix op(B).
     * @param[in]  alpha   Specifies the scalar alpha.
     * @param[in]  a       Array size lda by ka, where ka is k when transa = 'N' or 'n', and is m otherwise. Before entry with transa = 'N' or 'n', the leading m-by-k part of the array a must contain the matrix A, otherwise the leading k-by-m part of the array a must contain the matrix A.
     * @param[in]  lda     specifies the leading dimension of a as declared in the calling (sub)program.
     * @param[in]  b       Array, size ldb by kb, where kb is n when transa = 'N' or 'n', and is k otherwise. Before entry with transa = 'N' or 'n', the leading k-by-n part of the array b must contain the matrix B, otherwise the leading n-by-k part of the array b must contain the matrix B.
     * @param[in]  ldb     specifies the leading dimension of a as declared in the calling (sub)program.
     * @param[in]  beta    Specifies the scalar beta.
     * @param[in]  c       Array, size ldc by n. Before entry, the leading m-by-n part of the array c must contain the matrix C, except when beta is equal to zero, in which case c need not be set on entry.
     * @param[in]  ldc     Specifies the leading dimension of c as declared in the calling (sub)program.
     */
    void dgemm_(const char* transa, const char* transb, const int m, const int n, const int k, const double alpha, const std::unique_ptr<double []>& a, const int lda, const std::unique_ptr<double []>& b, const int ldb, const double beta, std::unique_ptr<double []>& c, const int ldc)
        { ::dgemm_(transa,transb,&m,&n,&k,&alpha,a.get(),&lda,b.get(),&ldb,&beta,c.get(),&ldc); }


    /**
     * @brief       The routine computes all the eigenvalues and, optionally, the eigenvectors of a square real symmetric matrix A.
     *
     * @param[in]  jobz   Must be 'N' or 'V'. If jobz = 'N', then only eigenvalues are computed. If jobz = 'V', then eigenvalues and eigenvectors are computed.
     * @param[in]  uplo   Must be 'U' or 'L'. If uplo = 'U', a stores the upper triangular part of A. If uplo = 'L', a stores the lower triangular part of A.
     * @param[in]  n      The order of the matrix A (n ≥ 0).
     * @param      a      a (size max(1, lda*n)) is an array containing either upper or lower triangular part of the symmetric matrix A, as specified by uplo.
     * @param[in]  lda    The leading dimension of the array a.
     * @param      w      Array, size at least max(1, n). If INFO = 0, contains the eigenvalues of the matrix A in ascending order.
     * @param      WORK   array, dimension (MAX(1,LWORK)) On exit, if i = 0, g(1) returns the optimal LWORK.
     * @param[in]  LWORK  The length of the array WORK.  LWORK >= max(1,3*N-1). For optimal efficiency, WORK >= (NB+2)*N, where NB is the blocksize for DSYTRD returned by ILAENV. If WORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the LWORK array, returns this value as the first entry of the LWORK array, and no error message related to WORK is issued by XERBLA.
     * @param      INFO   If INFO=0, the execution is successful. If INFO = -i, the i-th parameter had an illegal value. If INFO = i, then the algorithm failed to converge; i indicates the number of elements of an intermediate tridiagonal form which did not converge to zero.
     */
    void dsyev_(const char* jobz, const char* uplo, const int n, double* a, const int lda, double* w, double* WORK, const int LWORK, int& INFO)
        { ::dsyev_(jobz,uplo,&n,a,&lda,w, WORK, &LWORK,&INFO);}

    /**
     * @brief       The routine computes all the eigenvalues and, optionally, the eigenvectors of a square real symmetric matrix A.
     *
     * @param[in]  jobz   Must be 'N' or 'V'. If jobz = 'N', then only eigenvalues are computed. If jobz = 'V', then eigenvalues and eigenvectors are computed.
     * @param[in]  uplo   Must be 'U' or 'L'. If uplo = 'U', a stores the upper triangular part of A. If uplo = 'L', a stores the lower triangular part of A.
     * @param[in]  n      The order of the matrix A (n ≥ 0).
     * @param      a      a (size max(1, lda*n)) is an array containing either upper or lower triangular part of the symmetric matrix A, as specified by uplo.
     * @param[in]  lda    The leading dimension of the array a.
     * @param      w      Array, size at least max(1, n). If INFO = 0, contains the eigenvalues of the matrix A in ascending order.
     * @param      WORK   array, dimension (MAX(1,LWORK)) On exit, if i = 0, g(1) returns the optimal LWORK.
     * @param[in]  LWORK  The length of the array WORK.  LWORK >= max(1,3*N-1). For optimal efficiency, WORK >= (NB+2)*N, where NB is the blocksize for DSYTRD returned by ILAENV. If WORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the LWORK array, returns this value as the first entry of the LWORK array, and no error message related to WORK is issued by XERBLA.
     * @param      INFO   If INFO=0, the execution is successful. If INFO = -i, the i-th parameter had an illegal value. If INFO = i, then the algorithm failed to converge; i indicates the number of elements of an intermediate tridiagonal form which did not converge to zero.
     */
    void dsyev_(const char* jobz, const char* uplo, const int n, std::unique_ptr<double []>& a, const int lda, std::unique_ptr<double []>& w, std::unique_ptr<double []>& WORK, const int LWORK, int& INFO)
        { ::dsyev_(jobz,uplo,&n,a.get(),&lda,w.get(), WORK.get(), &LWORK,&INFO);}

    /**
     * @brief      wrapper for dgesvd_
     *
     * @param[in]  jobu    Must be 'A', 'S', 'O', or 'N'. Specifies options for computing all or part of the matrix U. If jobu = 'A', all m columns of U are returned in the array u; if jobu = 'S', the first min(m, n) columns of U (the left singular vectors) are returned in the array u; if jobu = 'O', the first min(m, n) columns of U (the left singular vectors) are overwritten on the array a; if jobu = 'N', no columns of U (no left singular vectors) are computed.
     * @param[in]  jobvt   Must be 'A', 'S', 'O', or 'N'. Specifies options for computing all or part of the matrix VT/VH. If jobvt = 'A', all n rows of VT/VH are returned in the array vt; if jobvt = 'S', the first min(m,n) rows of VT/VH (the right singular vectors) are returned in the array vt; if jobvt = 'O', the first min(m,n) rows of VT/VH) (the right singular vectors) are overwritten on the array a; if jobvt = 'N', no rows of VT/VH (no right singular vectors) are computed. jobvt and jobu cannot both be 'O'.
     * @param[in]  m       The number of rows of the matrix A (m ≥ 0).
     * @param[in]  n       The number of columns in A (n ≥ 0).
     * @param      a       Arrays: a(size at least max(1, lda*n) for column major layout and max(1, lda*m) for row major layout) is an array containing the m-by-n matrix A.
     * @param[in]  lda     The leading dimension of the array a. Must be at least max(1, m) for column major layout and at least max(1, n) for row major layout .
     * @param      s       array, dimension (min(M,N)) The singular values of A, sorted so that S(i) >= S(i+1).
     * @param      u       U is DOUBLE PRECISION array, dimension (LDU,UCOL) (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'. If JOBU = 'A', U contains the M-by-M orthogonal matrix U; if JOBU = 'S', U contains the first min(m,n) columns of U (the left singular vectors, stored columnwise); if JOBU = 'N' or 'O', U is not referenced.
     * @param[in]  ldu     LDU is INTEGER The leading dimension of the array U.  LDU >= 1; if JOBU = 'S' or 'A', LDU >= M.
     * @param      vt      VT is DOUBLE PRECISION array, dimension (LDVT,N) If JOBVT = 'A', VT contains the N-by-N orthogonal matrix V**T; if JOBVT = 'S', VT contains the first min(m,n) rows of V**T (the right singular vectors, stored rowwise); if JOBVT = 'N' or 'O', VT is not referenced.
     * @param[in]  ldvt    LDVT is INTEGER The leading dimension of the array VT.  LDVT >= 1; if JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
     * @param      superb  If ?bdsqr does not converge (indicated by the return value info > 0), on exit superb(0:min(m,n)-2) contains the unconverged superdiagonal elements of an upper bidiagonal matrix B whose diagonal is in s (not necessarily sorted). B satisfies A = u*B*VT (real flavors) or A = u*B*VH (complex flavors), so it has the same singular values as A, and singular vectors related by u and vt.
     * @param[in]  lwork    LWORK is INTEGER The dimension of the array WORK. LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code): - PATH 1  (M much larger than N, JOBU='N')  - PATH 1t (N much larger than M, JOBVT='N') LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)) for the other paths For good performance, LWORK should generally be larger. If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.
     * @param      info    INFO is INTEGER = 0:  successful exit. < 0:  if INFO = -i, the i-th argument had an illegal value. > 0:  if DBDSQR did not converge, INFO specifies how many superdiagonals of an intermediate bidiagonal form B did not converge to zero. See the description of WORK above for details.
     */
    void dgesvd_(const char* jobu, const char* jobvt, const int m, const int n, double* a, const int lda, double* s, double* u, const int ldu, double* vt, const int ldvt, double* superb, const int lwork , int& info) {return ::dgesvd_(jobu, jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, superb, &lwork , &info); }

    /**
     * @brief      mkl_domatcopy_ wrapper
     *
     * @param[in]  ordering  Ordering of the matrix storage. If ordering = 'R' or 'r', the ordering is row-major. If ordering = 'C' or 'c', the ordering is column-major.
     * @param[in]  trans     Parameter that specifies the operation type. If trans = 'N' or 'n', op(A)=A and the matrix A is assumed unchanged on input. If trans = 'T' or 't', it is assumed that A should be transposed. If trans = 'C' or 'c', it is assumed that A should be conjugate transposed. If trans = 'R' or 'r', it is assumed that A should be only conjugated. If the data is real, then trans = 'R' is the same as trans = 'N', and trans = 'C' is the same as trans = 'T'.
     * @param[in]  r         he number of rows in the source matrix.
     * @param[in]  c         The number of columns in the source matrix.
     * @param[in]  alpha     This parameter scales the input matrix by alpha.
     * @param[in]  A         Array
     * @param[in]  nr        Distance between the first elements in adjacent columns (in the case of the column-major order) or rows (in the case of the row-major order) in the source matrix; measured in the number of elements. This parameter must be at least max(1,rows) if ordering = 'C' or 'c', and max(1,cols) otherwise.
     * @param[in]  B         Array
     * @param[in]  nc        Distance between the first elements in adjacent columns (in the case of the column-major order) or rows (in the case of the row-major order) in the destination matrix; measured in the number of elements. To determine the minimum value of ldb on output, consider the following guideline: If ordering = 'C' or 'c', then If trans = 'T' or 't' or 'C' or 'c', this parameter must be at least max(1,cols) If trans = 'N' or 'n' or 'R' or 'r', this parameter must be at least max(1,rows) If ordering = 'R' or 'r', then If trans = 'T' or 't' or 'C' or 'c', this parameter must be at least max(1,rows) If trans = 'N' or 'n' or 'R' or 'r', this parameter must be at least max(1,cols)
     */
    void mkl_domatcopy_(const char* ordering, const char* trans, const int r, const int c, const double alpha, const double* A, const int nr, double* B, const int nc) {::mkl_domatcopy_(ordering,trans,&r,&c,&alpha,A,&nr,B,&nc);}

    /**
     * @brief      wrapper for ddot
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  x     Array, size at least (1+(n-1)*abs(incx)).
     * @param[in]  incx  INTEGER. Specifies the increment for the elements of x.
     * @param[in]  y     Array, size at least (1+(n-1)*abs(incy)).
     * @param[in]  incy  INTEGER. Specifies the increment for the elements of y.
     *
     * @return     Contains the result of the dot product of x and y, if n is positive. Otherwise,  0.
     */
    double ddot_(const int n, const double* x, const int incx, const double* y, const int incy) { return ::ddot_(&n,x,&incx,y,&incy);}

    /**
     * @brief      wrapper for ddot
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  x     Array, size at least (1+(n-1)*abs(incx)).
     * @param[in]  incx  INTEGER. Specifies the increment for the elements of x.
     * @param[in]  y     Array, size at least (1+(n-1)*abs(incy)).
     * @param[in]  incy  INTEGER. Specifies the increment for the elements of y.
     *
     * @return     Contains the result of the dot product of x and y, if n is positive. Otherwise,  0.
     */
    double ddot_(const int n, const std::unique_ptr<double []>& x, const int incx, const std::unique_ptr<double []>& y, const int incy)
        { return ::ddot_(&n,x.get(),&incx,y.get(),&incy); }

    /**
     * @brief      { function_description }
     *
     * @param[in]  jobz    Must be 'N' or 'V'. If jobz = 'N', then only eigenvalues are computed. If jobz = 'V', then eigenvalues and eigenvectors are computed.
     * @param[in]  range   Must be 'A' or 'V' or 'I'. If range = 'A', the routine computes all eigenvalues. If range = 'V', the routine computes eigenvalues w[i] in the half-open interval: vl < w[i] ≤ vu. If range = 'I', the routine computes eigenvalues with indices il to iu.
     * @param[in]  uplo    Must be 'U' or 'L'. If uplo = 'U', a stores the upper triangular part of A. If uplo = 'L', a stores the lower triangular part of A.
     * @param[in]  n       The order of the matrix A (n ≥ 0).
     * @param[in]  a       a (size max(1, lda*n)) is an array containing either upper or lower triangular part of the symmetric matrix A, as specified by uplo.
     * @param[in]  lda     The leading dimension of the array a.
     * @param[in]  vl      If range = 'V', the lower and upper bounds of the interval to be searched for eigenvalues.Constraint: vl< vu.If range = 'A' or 'I', vl and vu are not referenced.
     * @param[in]  vu      The upper bound for vl
     * @param[in]  il      If range = 'I', the indices in ascending order of the smallest and largest eigenvalues to be returned. Constraint: 1 ≤ il ≤ iu ≤ n, if n > 0; il=1 and iu=0, if n = 0. If range = 'A' or 'V', il and iu are not referenced.
     * @param[in]  iu      upper bound for il
     * @param[in]  abstol  If jobz = 'V', the eigenvalues and eigenvectors output have residual norms bounded by abstol, and the dot products between different eigenvectors are bounded by abstol. If abstol < n *eps*||T||, then n *eps*||T|| is used instead, where eps is the machine precision, and ||T|| is the 1-norm of the matrix T. The eigenvalues are computed to an accuracy of eps*||T|| irrespective of abstol. If high relative accuracy is important, set abstol to ?lamch('S').
     * @param[in]  m       The total number of eigenvalues found, 0 ≤ m ≤ n. If range = 'A', m = n, if range = 'I', m = iu-il+1, and if range = 'V' the exact value of m is not known in advance.
     * @param[in]  w       Arrays: w, size at least max(1, n), contains the selected eigenvalues in ascending order, stored in w[0] to w[m - 1];
     * @param[in]  z       z(size max(1, ldz*m) for column major layout and max(1, ldz*n) for row major layout) . If jobz = 'V', then if info = 0, the first m columns of z contain the orthonormal eigenvectors of the matrix A corresponding to the selected eigenvalues, with the i-th column of z holding the eigenvector associated with w[i - 1]. If jobz = 'N', then z is not referenced.
     * @param[in]  ldz     The leading dimension of the output array z. Constraints: ldz ≥ 1 if jobz = 'N' and ldz ≥ max(1, n) for column major layout and ldz ≥ max(1, m) for row major layout if jobz = 'V'.
     * @param[in]  isuppz  Array, size at least 2 *max(1, m). The support of the eigenvectors in z, i.e., the indices indicating the nonzero elements in z. The i-th eigenvector is nonzero only in elements isuppz[2i - 2] through isuppz[2i - 1]. Referenced only if eigenvectors are needed (jobz = 'V') and all eigenvalues are needed, that is, range = 'A' or range = 'I' and il = 1 and iu = n.
     * @param[in]  work    outputs I can't find definitions for on the website
     * @param[in]  lwork   outputs I can't find definitions for on the website
     * @param[in]  iwork   outputs I can't find definitions for on the website
     * @param[in]  liwork  outputs I can't find definitions for on the website
     * @param[in]  info    If info=0, the execution is successful. If info = -i, the i-th parameter had an illegal value. If info = i, then the algorithm failed to converge; i indicates the number of elements of an intermediate tridiagonal form which did not converge to zero.
     */
    void dsyevr_(const char* jobz, const char* range, const char* uplo, const int n, double* a, const int lda, double vl, double vu, const int il, const int iu, const double abstol, int m, double* w, double* z, const int ldz, int* isuppz, double* work, int lwork, int* iwork, int liwork, int& info)
        { ::dsyevr_(jobz, range, uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);}

    /**
     * @brief      wrapper for saxpy_
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  a     Specifies the scalar a.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  incx  INTEGER. Specifies the increment for the elements of x
     * @param      y     Array, size at least (1 + (n-1)*abs(incy)).
     * @param[in]  incy  INTEGER. Specifies the increment for the elements of y.
     */
    void saxpy_(const int n, const int a, const int* x, const int incx, int* y, const int incy) { ::saxpy_(&n,&a,x,&incx,y,&incy); }

    /**
     * @brief      wrapper for scopy_
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  incx  INTEGER. Specifies the increment for the elements of x
     * @param      y     Array, size at least (1 + (n-1)*abs(incy)).
     * @param[in]  incy  INTEGER. Specifies the increment for the elements of y.
     */
    void scopy_(const int n, const int* x, const int incx, int* y, const int incy) { ::scopy_(&n,x,&incx,y,&incy);}

   /**
     * @brief      wrapper for saxpy_
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  a     Specifies the scalar a.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  incx  INTEGER. Specifies the increment for the elements of x
     */
    void sscal_(const int n, const int a, int* x, const int incx) {::sscal_(&n,&a,x,&incx);}

   /**
     * @brief      wrapper for saxpy_
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  a     Specifies the scalar a.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  incx  INTEGER. Specifies the increment for the elements of x
     * @param      y     Array, size at least (1 + (n-1)*abs(incy)).
     * @param[in]  incy  INTEGER. Specifies the increment for the elements of y.
     */
    void daxpy_(const int n, const double a, const double* x, const int incx, double* y, const int incy) { ::daxpy_(&n,&a,x,&incx,y,&incy); }

    /**
     * @brief      wrapper for daxpy_
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  incx  INTEGER. Specifies the increment for the elements of x
     * @param      y     Array, size at least (1 + (n-1)*abs(incy)).
     * @param[in]  incy  INTEGER. Specifies the increment for the elements of y.
     */
    void dcopy_(const int n, const double* x, const int incx, double* y, const int incy) { ::dcopy_(&n,x,&incx,y,&incy);}

   /**
     * @brief      wrapper for saxpy_
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  a     Specifies the scalar a.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  incx  INTEGER. Specifies the increment for the elements of x
     */
    void dscal_(const int n, const double a, double* x, const int incx) {::dscal_(&n,&a,x,&incx);}

    /**
     * @brief      wrapper for zaxpy_
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  a     Specifies the scalar a.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  incx  INTEGER. Specifies the increment for the elements of x
     * @param      y     Array, size at least (1 + (n-1)*abs(incy)).
     * @param[in]  incy  INTEGER. Specifies the increment for the elements of y.
     */
    void zaxpy_(const int n, const cplx a, const cplx* x, const int incx, cplx* y, const int incy) { ::zaxpy_(&n,&a,x,&incx,y,&incy); }

    /**
     * @brief      wrapper for daxpy_
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  incx  INTEGER. Specifies the increment for the elements of x
     * @param      y     Array, size at least (1 + (n-1)*abs(incy)).
     * @param[in]  incy  INTEGER. Specifies the increment for the elements of y.
     */
    void zcopy_(const int n, const cplx* x, const int incx, cplx* y, const int incy) { ::zcopy_(&n,x,&incx,y,&incy);}

   /**
     * @brief      wrapper for saxpy_
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  a     Specifies the scalar a.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  incx  INTEGER. Specifies the increment for the elements of x
     */
    void zscal_(const int n, const cplx a, cplx* x, const int incx) {::zscal_(&n,&a,x,&incx);}

    /**
     * @brief      Wrapper for zgemm3m_: c<-$\alpha a*b + /beta*c (Similar to zgemm_ but different algorithm to do it)
     *
     * @param[in]  transa  'T' or 't' for Transpose, 'C' or c' for Hermitian Conjugate, 'N' or 'n' for neither on Matrix A
     * @param[in]  transb  'T' or 't' for Transpose, 'C' or c' for Hermitian Conjugate, 'N' or 'n' for neither on Matrix B
     * @param[in]  m       Specifies the number of rows of the matrix op(A) and of the matrix C. The value of m must be at least zero.
     * @param[in]  n       Specifies the number of columns of the matrix op(A) and the number of rows of the matrix op(B).
     * @param[in]  k       Specifies the number of columns of the matrix op(A) and the number of rows of the matrix op(B).
     * @param[in]  alpha   Specifies the scalar alpha.
     * @param[in]  a       Array size lda by ka, where ka is k when transa = 'N' or 'n', and is m otherwise. Before entry with transa = 'N' or 'n', the leading m-by-k part of the array a must contain the matrix A, otherwise the leading k-by-m part of the array a must contain the matrix A.
     * @param[in]  lda     specifies the leading dimension of a as declared in the calling (sub)program.
     * @param[in]  b       Array, size ldb by kb, where kb is n when transa = 'N' or 'n', and is k otherwise. Before entry with transa = 'N' or 'n', the leading k-by-n part of the array b must contain the matrix B, otherwise the leading n-by-k part of the array b must contain the matrix B.
     * @param[in]  ldb     specifies the leading dimension of a as declared in the calling (sub)program.
     * @param[in]  beta    Specifies the scalar beta.
     * @param[in]  c       Array, size ldc by n. Before entry, the leading m-by-n part of the array c must contain the matrix C, except when beta is equal to zero, in which case c need not be set on entry.
     * @param[in]  ldc     Specifies the leading dimension of c as declared in the calling (sub)program.
     */
    void zgemm3m_(const char* transa, const char* transb, const int m, const int n, const int k, const cplx alpha, const cplx* a, const int lda, const cplx* b, const int ldb, const cplx beta, cplx* c, const int ldc)
        { ::zgemm3m_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }

    /**
     * @brief      Wrapper for zgemm3m_: c<-$\alpha a*b + /beta*c (Similar to zgemm_ but different algorithm to do it)
     *
     * @param[in]  transa  'T' or 't' for Transpose, 'C' or c' for Hermitian Conjugate, 'N' or 'n' for neither on Matrix A
     * @param[in]  transb  'T' or 't' for Transpose, 'C' or c' for Hermitian Conjugate, 'N' or 'n' for neither on Matrix B
     * @param[in]  m       Specifies the number of rows of the matrix op(A) and of the matrix C. The value of m must be at least zero.
     * @param[in]  n       Specifies the number of columns of the matrix op(A) and the number of rows of the matrix op(B).
     * @param[in]  k       Specifies the number of columns of the matrix op(A) and the number of rows of the matrix op(B).
     * @param[in]  alpha   Specifies the scalar alpha.
     * @param[in]  a       Array size lda by ka, where ka is k when transa = 'N' or 'n', and is m otherwise. Before entry with transa = 'N' or 'n', the leading m-by-k part of the array a must contain the matrix A, otherwise the leading k-by-m part of the array a must contain the matrix A.
     * @param[in]  lda     specifies the leading dimension of a as declared in the calling (sub)program.
     * @param[in]  b       Array, size ldb by kb, where kb is n when transa = 'N' or 'n', and is k otherwise. Before entry with transa = 'N' or 'n', the leading k-by-n part of the array b must contain the matrix B, otherwise the leading n-by-k part of the array b must contain the matrix B.
     * @param[in]  ldb     specifies the leading dimension of a as declared in the calling (sub)program.
     * @param[in]  beta    Specifies the scalar beta.
     * @param[in]  c       Array, size ldc by n. Before entry, the leading m-by-n part of the array c must contain the matrix C, except when beta is equal to zero, in which case c need not be set on entry.
     * @param[in]  ldc     Specifies the leading dimension of c as declared in the calling (sub)program.
     */
    void zgemm3m_(const char* transa, const char* transb, const int m, const int n, const int k, const cplx alpha, const std::unique_ptr<cplx[]>& a, const int lda, const std::unique_ptr<cplx[]>& b, const int ldb, const cplx beta, std::unique_ptr<cplx[]>& c, const int ldc)
        { ::zgemm3m_(transa,transb,&m,&n,&k,&alpha,a.get(),&lda,b.get(),&ldb,&beta,c.get(),&ldc); }

    void zhemm_(const char* side, const char* uplo, const int m, const int n, const cplx alpha, const cplx *a, const int lda, const cplx *b, const int ldb, const cplx beta, const cplx *c, const int ldc)
        { ::zhemm_(side,uplo,&m,&n,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }


    void zhemm_(const char* side, const char* uplo, const int m, const int n, const cplx alpha, const std::unique_ptr<cplx[]>& a, const int lda, const std::unique_ptr<cplx[]>& b, const int ldb, const cplx beta, const std::unique_ptr<cplx[]>& c, const int ldc)
        { ::zhemm_(side,uplo,&m,&n,&alpha,a.get(),&lda,b.get(),&ldb,&beta,c.get(),&ldc); }

    /**
     * @brief      zgemv_ wrapper
     *
     * @param[in]  trans  CHARACTER*1. Specifies the operation: if trans= 'N' or 'n', then y := alpha*A*x + beta*y; if trans= 'T' or 't', then y := alpha*A'*x + beta*y; if trans= 'C' or 'c', then y := alpha *conjg(A')*x + beta*y.
     * @param[in]  m      INTEGER. Specifies the number of rows of the matrix A. The value of m must be at least zero.
     * @param[in]  n      INTEGER. Specifies the number of columns of the matrix A. The value of n must be at least zero.
     * @param[in]  alpha  Specifies the scalar alpha.
     * @param[in]  a      Array, size (lda, n). Before entry, the leading m-by-n part of the array a must contain the matrix of coefficients.
     * @param[in]  lda    Specifies the leading dimension of a as declared in the calling (sub)program.
     * @param[in]  x      Array, size at least (1+(n-1)*abs(incx)) when trans= 'N' or 'n' and at least (1+(m - 1)*abs(incx)) otherwise. Before entry, the incremented array x must contain the vector x.
     * @param[in]  incx   Specifies the increment for the elements of x. The value of incx must not be zero.
     * @param[in]  beta   Specifies the scalar beta. When beta is set to zero, then y need not be set on input.
     * @param      y      Array, size at least (1 +(m - 1)*abs(incy)) when trans= 'N' or 'n' and at least (1 +(n - 1)*abs(incy)) otherwise. Before entry with non-zero beta, the incremented array y must contain the vector y.
     * @param[in]  incy   Specifies the increment for the elements of y. The value of incy must not be zero.
     */
    void zgemv_(const char* trans, const int m, const int n, const cplx alpha, cplx* a, const int lda, cplx* x, const int incx, const cplx beta, cplx* y, const int incy)
        { ::zgemv_(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }

    /**
     * @brief      zgemv_ wrapper
     *
     * @param[in]  trans  CHARACTER*1. Specifies the operation: if trans= 'N' or 'n', then y := alpha*A*x + beta*y; if trans= 'T' or 't', then y := alpha*A'*x + beta*y; if trans= 'C' or 'c', then y := alpha *conjg(A')*x + beta*y.
     * @param[in]  m      INTEGER. Specifies the number of rows of the matrix A. The value of m must be at least zero.
     * @param[in]  n      INTEGER. Specifies the number of columns of the matrix A. The value of n must be at least zero.
     * @param[in]  alpha  Specifies the scalar alpha.
     * @param[in]  a      Array, size (lda, n). Before entry, the leading m-by-n part of the array a must contain the matrix of coefficients.
     * @param[in]  lda    Specifies the leading dimension of a as declared in the calling (sub)program.
     * @param[in]  x      Array, size at least (1+(n-1)*abs(incx)) when trans= 'N' or 'n' and at least (1+(m - 1)*abs(incx)) otherwise. Before entry, the incremented array x must contain the vector x.
     * @param[in]  incx   Specifies the increment for the elements of x. The value of incx must not be zero.
     * @param[in]  beta   Specifies the scalar beta. When beta is set to zero, then y need not be set on input.
     * @param      y      Array, size at least (1 +(m - 1)*abs(incy)) when trans= 'N' or 'n' and at least (1 +(n - 1)*abs(incy)) otherwise. Before entry with non-zero beta, the incremented array y must contain the vector y.
     * @param[in]  incy   Specifies the increment for the elements of y. The value of incy must not be zero.
     */
    void zgemv_(const char* trans, const int m, const int n, const cplx alpha, const std::unique_ptr<cplx []>& a, const int lda, const std::unique_ptr<cplx []>& x, const int incx, const cplx beta, std::unique_ptr<cplx []>& y, const int incy)
        { ::zgemv_(trans, &m, &n, &alpha, a.get(), &lda, x.get(), &incx, &beta, y.get(), &incy); }

    /**
     * @brief      mkl_zomatcopy_ wrapper
     *
     * @param[in]  ordering  Ordering of the matrix storage. If ordering = 'R' or 'r', the ordering is row-major. If ordering = 'C' or 'c', the ordering is column-major.
     * @param[in]  trans     Parameter that specifies the operation type. If trans = 'N' or 'n', op(A)=A and the matrix A is assumed unchanged on input. If trans = 'T' or 't', it is assumed that A should be transposed. If trans = 'C' or 'c', it is assumed that A should be conjugate transposed. If trans = 'R' or 'r', it is assumed that A should be only conjugated. If the data is real, then trans = 'R' is the same as trans = 'N', and trans = 'C' is the same as trans = 'T'.
     * @param[in]  r         he number of rows in the source matrix.
     * @param[in]  c         The number of columns in the source matrix.
     * @param[in]  alpha     This parameter scales the input matrix by alpha.
     * @param[in]  A         Array
     * @param[in]  nr        Distance between the first elements in adjacent columns (in the case of the column-major order) or rows (in the case of the row-major order) in the source matrix; measured in the number of elements. This parameter must be at least max(1,rows) if ordering = 'C' or 'c', and max(1,cols) otherwise.
     * @param[in]  B         Array
     * @param[in]  nc        Distance between the first elements in adjacent columns (in the case of the column-major order) or rows (in the case of the row-major order) in the destination matrix; measured in the number of elements. To determine the minimum value of ldb on output, consider the following guideline: If ordering = 'C' or 'c', then If trans = 'T' or 't' or 'C' or 'c', this parameter must be at least max(1,cols) If trans = 'N' or 'n' or 'R' or 'r', this parameter must be at least max(1,rows) If ordering = 'R' or 'r', then If trans = 'T' or 't' or 'C' or 'c', this parameter must be at least max(1,rows) If trans = 'N' or 'n' or 'R' or 'r', this parameter must be at least max(1,cols)
     */
    void mkl_zomatcopy_(const char* ordering, const char* trans, const int r, const int c, const cplx alpha, const cplx* A, const int nr, cplx* B, const int nc)
        {::mkl_zomatcopy_(ordering,trans,&r,&c,&alpha,A,&nr,B,&nc);}

    /**
     * @brief      izamax_ wrapper
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  inc   INTEGER. Specifies the increment for the elements of x
     *
     * @return     INTEGER Contains the position of vector element that has the largest absolute value such that x(index) has the largest absolute value.
     */
    int izamax_(const int n, const cplx *x, const int inc)
        { return ::izamax_(&n,x,&inc);}

    /**
     * @brief      izamin_ wrapper
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  inc   INTEGER. Specifies the increment for the elements of x
     *
     * @return     INTEGER. Indicates the position of vector element with the smallest absolute value such that x(index) has the smallest absolute value.
     */
    int izamin_(const int n, const cplx *x, const int inc)
        { return ::izamin_(&n,x,&inc);}

    /**
     * @brief      idamax_ wrapper
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  inc   INTEGER. Specifies the increment for the elements of x
     *
     * @return     INTEGER Contains the position of vector element that has the largest absolute value such that x(index) has the largest absolute value.
     */
    int idamax_(const int n, const double *x, const int inc)
        { return ::idamax_(&n,x,&inc);}

    /**
     * @brief      idamin_ wrapper
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  inc   INTEGER. Specifies the increment for the elements of x
     *
     * @return     INTEGER. Indicates the position of vector element with the smallest absolute value such that x(index) has the smallest absolute value.
     */
    int idamin_(const int n, const double *x, const int inc)
        { return ::idamin_(&n,x,&inc);}

    /**
     * @brief      isamax_ wrapper
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  inc   INTEGER. Specifies the increment for the elements of x
     *
     * @return     INTEGER Contains the position of vector element that has the largest absolute value such that x(index) has the largest absolute value.
     */
    int isamax_(const int n, const int *x, const int inc)
        { return ::isamax_(&n,x,&inc);}
    /**
     * @brief      isamin_ wrapper
     *
     * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
     * @param[in]  x     Array, size at least (1 + (n-1)*abs(incx)).
     * @param[in]  inc   INTEGER. Specifies the increment for the elements of x
     *
     * @return     INTEGER. Indicates the position of vector element with the smallest absolute value such that x(index) has the smallest absolute value.
     */
    int isamin_(const int n, const int *x, const int inc)
        { return ::isamin_(&n,x,&inc);}

    /**
     * @brief      wrapper for zheev_
     *
     * @param[in]  jobz   Must be 'N' or 'V'. If jobz = 'N', then only eigenvalues are computed. If jobz = 'V', then eigenvalues and eigenvectors are computed.
     * @param[in]  uplo   Must be 'U' or 'L'. If uplo = 'U', a stores the upper triangular part of A. If uplo = 'L', a stores the lower triangular part of A.
     * @param[in]  n      The order of the matrix A (n ≥ 0).
     * @param      a      a (size max(1, lda*n)) is an array containing either upper or lower triangular part of the Hermitian matrix A, as specified by uplo.
     * @param[in]  lda    The leading dimension of the array a. Must be at least max(1, n).
     * @param      w      Array, size at least max(1, n).
     * @param      WORK   (workspace/output) COMPLEX*16 array, dimension (LWORK) On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
     * @param[in]  LWORK  (input) INTEGER The length of the array WORK.  LWORK >= max(1,2*N-1). For optimal efficiency, LWORK >= (NB+1)*N, where NB is the blocksize for ZHETRD returned by ILAENV.  If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.
     * @param      RWORK  (workspace) DOUBLE PRECISION array, dimension (max(1, 3*N-2))
     * @param      info   INTEGER = 0:  successful exit < 0:  if INFO = -i, the i-th argument had an illegal value > 0:  if INFO = i, the algorithm failed to converge; i   off-diagonal elements of an intermediate tridiagonal form did not converge to zero.
     */
    void zheev_(const char* jobz, const char* uplo, const int n, cplx* a, const int lda,double* w, cplx* WORK, const int LWORK,double* RWORK, int& info )
        { ::zheev_( jobz, uplo, &n, a, &lda, w, WORK, &LWORK, RWORK, &info); }


    /**
     * @brief      wrapper for zheev_
     *
     * @param[in]  jobz   Must be 'N' or 'V'. If jobz = 'N', then only eigenvalues are computed. If jobz = 'V', then eigenvalues and eigenvectors are computed.
     * @param[in]  uplo   Must be 'U' or 'L'. If uplo = 'U', a stores the upper triangular part of A. If uplo = 'L', a stores the lower triangular part of A.
     * @param[in]  n      The order of the matrix A (n ≥ 0).
     * @param      a      a (size max(1, lda*n)) is an array containing either upper or lower triangular part of the Hermitian matrix A, as specified by uplo.
     * @param[in]  lda    The leading dimension of the array a. Must be at least max(1, n).
     * @param      w      Array, size at least max(1, n).
     * @param      WORK   (workspace/output) COMPLEX*16 array, dimension (LWORK) On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
     * @param[in]  LWORK  (input) INTEGER The length of the array WORK.  LWORK >= max(1,2*N-1). For optimal efficiency, LWORK >= (NB+1)*N, where NB is the blocksize for ZHETRD returned by ILAENV.  If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.
     * @param      RWORK  (workspace) DOUBLE PRECISION array, dimension (max(1, 3*N-2))
     * @param      info   INTEGER = 0:  successful exit < 0:  if INFO = -i, the i-th argument had an illegal value > 0:  if INFO = i, the algorithm failed to converge; i   off-diagonal elements of an intermediate tridiagonal form did not converge to zero.
     */
    void zheev_(const char* jobz, const char* uplo, const int n, std::unique_ptr<cplx []>& a, const int lda, std::unique_ptr<double[]>& w, std::unique_ptr<cplx []>& WORK, const int LWORK, std::unique_ptr<double[]>& RWORK, int& info )
        { ::zheev_( jobz, uplo, &n, a.get(), &lda, w.get(), WORK.get(), &LWORK, RWORK.get(), &info ); }
    #ifndef ZDOT_RETURN
        /**
         * @brief      wrapper for zdotc_
         *
         * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
         * @param[in]  x     Array, size at least (1+(n-1)*abs(incx)).
         * @param[in]  incx  INTEGER. Specifies the increment for the elements of x.
         * @param[in]  y     Array, size at least (1+(n-1)*abs(incy)).
         * @param[in]  incy  INTEGER. Specifies the increment for the elements of y.
         *
         * @return     Contains the result of the dot product of x and y, if n is positive. Otherwise,  0.
         */
        cplx zdotc_(const int n, const cplx* x, const int incx, const cplx* y, const int incy)
        {
            cplx a;
            ::zdotc_(&a,&n,x,&incx,y,&incy);
            return a;
        }
        /**
         * @brief      wrapper for zdotc_
         *
         * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
         * @param[in]  x     Array, size at least (1+(n-1)*abs(incx)).
         * @param[in]  incx  INTEGER. Specifies the increment for the elements of x.
         * @param[in]  y     Array, size at least (1+(n-1)*abs(incy)).
         * @param[in]  incy  INTEGER. Specifies the increment for the elements of y.
         *
         * @return     Contains the result of the dot product of x and y, if n is positive. Otherwise,  0.
         */
        cplx zdotc_(const int n, const std::unique_ptr<cplx []>& x, const int incx, const std::unique_ptr<cplx []>& y, const int incy)
        {
            cplx a;
            ::zdotc_(&a,&n,x.get(),&incx,y.get(),&incy);
            return a;
         }
    #else
        /**
         * @brief      wrapper for zdotc_
         *
         * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
         * @param[in]  x     Array, size at least (1+(n-1)*abs(incx)).
         * @param[in]  incx  INTEGER. Specifies the increment for the elements of x.
         * @param[in]  y     Array, size at least (1+(n-1)*abs(incy)).
         * @param[in]  incy  INTEGER. Specifies the increment for the elements of y.
         *
         * @return     Contains the result of the dot product of x and y, if n is positive. Otherwise,  0.
         */
        cplx zdotc_(const int n, const cplx* x, const int incx, const cplx* y, const int incy) { return ::zdotc_(&n,x,&incx,y,&incy); }
        /**
         * @brief      wrapper for zdotc_
         *
         * @param[in]  n     INTEGER. Specifies the number of elements in vectors x and y.
         * @param[in]  x     Array, size at least (1+(n-1)*abs(incx)).
         * @param[in]  incx  INTEGER. Specifies the increment for the elements of x.
         * @param[in]  y     Array, size at least (1+(n-1)*abs(incy)).
         * @param[in]  incy  INTEGER. Specifies the increment for the elements of y.
         *
         * @return     Contains the result of the dot product of x and y, if n is positive. Otherwise,  0.
         */
        cplx zdotc_(const int n, const std::unique_ptr<cplx []>& x, const int incx, const std::unique_ptr<cplx []>& y, const int incy)
            { return ::zdotc_(&n,x.get(),&incx,y.get(),&incy); }
    #endif
    // Sparse routines
    // void mkl_dcsrsymv_(const char* uplo, const int m, const double* a, const int* ia, const int* ja, double* x, double* y)
    //             {mkl_dcsrsymv_(uplo, m, a, ia, ja, x, y); };

    // void mkl_dscrgemv_(const char* transa, const int m, const double* a, const int* ia, const int* ja, double* x, double* y)
    //             {mkl_dscrgemv_(transa, m, a, ia, ja, x, y); };

}

#endif