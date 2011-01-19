/*
 * Copyright (C) 2010 Philipp Wähnert
 *               2010-2011 Stefan Handschuh
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __BLASINTERFACE_HPP
#define __BLASINTERFACE_HPP

// TODO: Replace Complex by standard std::complex<T>
#include <complex>

// TODO: Uncomment namespace if all parts of the TensorCalculus library
//       are "namespace-ready", i.e. itself in this namespace
namespace TensorCalculus {

enum BlasTranspose
{
  NoTrans = 'N',
  Trans = 'T',
  ConjTrans = 'C'
};

enum BlasUpLo
{
  Upper = 'U',
  Lower = 'L'
};

enum BlasDiag
{
  NonUnit = 'N',
  Unit = 'U'
};

enum BlasSide
{
  Left = 'L',
  Right = 'R'
};

// Empty template for the general case
// The Blas interface declares the functions
// for T = double and complex
template<typename T>
struct Blas {
  static void copy(long int n, const T* x, long int incx, T* y, long int incy) {
    while (n-- > 0) {
      *y = *x;
      y += incy;
      x += incx;
    }
  }

  static void scal(long int n, T alpha, T* x, long int incx) {
    while (n-- > 0) {
      *x *= alpha;
      x += incx;
    }
  }

  static const T dot(long int n, const T* x, long int incx, const T* y, long int incy) {
    T result = T();
    while (n-- > 0) {
      result += (*x) * (*y);
      x += incx;
      y += incy;
    }
    return result;
  }

  static void axpy(long int n, T alpha, const T* x, long int incx, T* y, long int incy) {
    while (n-- > 0) {
      *y += alpha * (*x);
      x += incx;
      y += incy;
    }
  }

  static const T nrm2(long int n, const T* x, long int incx) {
    T result = 0.0;
    while (n-- > 0) {
      result += (*x) * (*x);
      x += incx;
    }
    return std::sqrt(result);
  }

  static const long int iamax(long int n, const T* x, const long int incx) {
	T max = std::abs(x[0]);

	long int index = 0;

	while (n-- > 1) {
		if (std::abs(x[n]) > max) {
			max = std::abs(x[n]);
			index = n;
		}
	}
	return index;
  }

};

// Forward declarations of the used BLAS routines
extern "C" void dcopy_(const int*, const double*, const int*,
                       double*, const int*);

extern "C" void daxpy_(const int*, const double*, const double*,
                       const int*, double*, const int*);

extern "C" void dscal_(const int*, const double*, double*, const int*);

extern "C" double ddot_(const int*, const double*, const int*,
                        const double*, const int*);

extern "C" double dnrm2_(const int*, const double*, const int*);

extern "C" void dgemm_(const char* , const char*, const int*,
                       const int*, const int*, const double*,
                       const double*, const int*, const double*,
                       const int*, const double*, double*, const int*);

extern "C" void dgemv_(const char* , const int*, const int*, const double*,
                       const double*, const int*,
                       const double*, const int*, const double*,
                       double*, const int*);

extern "C" void dger_(const int*, const int*, const double*, const double*,
                      const int*, const double*, const int*,
                      double*, const int*);

extern "C" double dasum_(const int*, const double*, const int*);

extern "C" void dsyr_(const char* uplo, const int* n, const double* alpha,
                      const double* x, const int* incx,
                      double* a, const int* lda);

extern "C" double dtrmm_(const char*, const char*, const char*, const char*,
                         const int*, const int*, const double*, const double*, const int*,
                         double*, const int*);

extern "C" void dsyev_(const char* jobz, const char* uplo, const int* n, double* a,
		               const int* lda, double* w, double* workspace, const int* lwork,
		               int* info);

extern "C" long idamax_(const long*, const double*, const int*);


template<>
struct Blas<double>
{

  static void copy(int n, const double* x, int incx, double* y, int incy)
  {
    dcopy_(&n, x, &incx, y, &incy);
  }

  static void axpy(int n, double alpha, const double* x, int incx,
                   double* y, int incy)
  {
    daxpy_(&n, &alpha, x, &incx, y, &incy);
  }

  static void scal(int n, double alpha, double* x, int incx)
  {
    dscal_(&n, &alpha, x, &incx);
  }

  static double dot(int n, const double* x, int incx,
                    const double* y, int incy)
  {
    return ddot_(&n, x, &incx, y, &incy);
  }

  static double nrm2(int n, const double* x, int incx)
  {
    return dnrm2_(&n, x, &incx);
  }

  static void gemm(char transa, char transb, int m, int n, int k, double alpha,
                   const double* A, int lda, const double* B, int ldb,
                   double beta, double* C, int ldc)
  {
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb,
           &beta, C, &ldc);
  }

  static void gemm(char transa, char transb, int m, int n, int k, double alpha,
                   const double* A, const double* B, double beta, double* C)
  {
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, transa == 'N' ? &m : &k, B,
    	   transb == 'N' ? &k : &n, &beta, C, &m);
  }

  static void gemv(char trans, int m, int n, double alpha, const double* A,
                   int lda, const double* X, int incx, double beta,
                   double* Y, int incy)
  {
    dgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
  }

  static void ger(int m, int n, double alpha, const double* x, int incx,
                  const double* y, int incy, double* A, int lda)
  {
    dger_(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
  }

  static double asum(int n, const double* x, int incx)
  {
    return dasum_(&n, x, &incx);
  }

  static void syr(char uplo, int n, double alpha, const double* x, int incx,
                  double* a, int lda)
  {
    dsyr_(&uplo, &n, &alpha, x, &incx, a, &lda);
  }
  
  static void trmm(char side, char uplo, char transA, char diag, int m,
                     int n, double alpha, const double* A, int lda, double* B, int ldb)
    {
      dtrmm_(&side, &uplo, &transA, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
    }


  static int syev(const char jobz, const char uplo, const int n, double* a,
  		               const int lda, double* w, double* workspace, const int lwork)
  {
	int info;

	dsyev_(&jobz, &uplo, &n, a, &lda, w, workspace, &lwork, &info);
	return info;
  }

  static long iamax(const long n, const double* x, const int incx) {
	return idamax_(&n, x, &incx);
  }

};

// TODO: Replace Complex by std::complex<double>!

extern "C" void zcopy_(const int*, const std::complex<double>*, const int*,
                       std::complex<double>*, const int*);

extern "C" void zaxpy_(const int*, const std::complex<double>*, const std::complex<double>*, const int*,
                       std::complex<double>*, const int*);

extern "C" void zscal_(const int*, const std::complex<double>*, std::complex<double>*, const int*);

extern "C" double dznrm2_(const int*, const std::complex<double>*, const int*);

extern "C" void zgemm_(const char* , const char*, const int*, const int*,
                       const int*, const std::complex<double>*, const std::complex<double>*, const int*,
                       const std::complex<double>*, const int*, const std::complex<double>*, std::complex<double>*,
                       const int*);

extern "C" void zgemv_(const char* , const int*, const int*, const std::complex<double>*,
                       const std::complex<double>*, const int*, const std::complex<double>*, const int*,
                       const std::complex<double>*, std::complex<double>*, const int*);

extern "C" double dzasum_(const int*, const std::complex<double>*, const int*);

extern "C" void zdotc_(std::complex<double>*, const int*, const std::complex<double>*, const int*,
                       const std::complex<double>*, const int*);

extern "C" void zdotu_(std::complex<double>*, const int*, const std::complex<double>*, const int*,
                       const std::complex<double>*, const int*);

extern "C" void zgeru_(const int*, const int*, const std::complex<double>*, const std::complex<double>*,
                       const int*, const std::complex<double>*, const int*, std::complex<double>*,
                       const int*);

extern "C" long izamax_(const long*, const std::complex<double>*, const int*);

// Specialization of the Blas template for Complex
template<>
struct Blas<std::complex<double> >
{

  static void copy(int n, const std::complex<double>* x, int incx, std::complex<double>* y, int incy)
  {
    zcopy_(&n, x, &incx, y, &incy);
  }

  static void axpy(int n, std::complex<double> alpha, const std::complex<double>* x, int incx,
                   std::complex<double>* y, int incy)
  {
    zaxpy_(&n, &alpha, x, &incx, y, &incy);
  }

  static void scal(int n, std::complex<double> alpha, std::complex<double>* x, int incx)
  {
    zscal_(&n, &alpha, x, &incx);
  }

  static double nrm2(int n, const std::complex<double>* x, int incx)
  {
    return dznrm2_(&n, x, &incx);
  }

  static void gemm(char transa, char transb, int m, int n, int k, std::complex<double> alpha,
                   const std::complex<double>* A, int lda, const std::complex<double>* B, int ldb,
                   std::complex<double> beta, std::complex<double>* C, int ldc)
  {
    zgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb,
           &beta, C, &ldc);
  }

  static void gemv(char trans, int m, int n, std::complex<double> alpha,
                   const std::complex<double>* A, int lda, const std::complex<double>* x, int incx,
                   std::complex<double> beta, std::complex<double>* y, int incy)
  {
    zgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }

  static double asum(int n, const std::complex<double>* x, int incx)
  {
    return dzasum_(&n, x, &incx);
  }

  static std::complex<double> dotc(int n, const std::complex<double>* x, int incx,
                      const std::complex<double>* y, int incy)
  {
    std::complex<double> result;
    zdotc_(&result, &n, x, &incx, y, &incy);
    return result;
  }

  static std::complex<double> dot(int n, const std::complex<double>* x, int incx,
                     const std::complex<double>* y, int incy)
  {
    std::complex<double> result;
    zdotu_(&result, &n, x, &incx, y, &incy);
    return result;
  }

  static void ger(int m, int n, const std::complex<double> alpha, const std::complex<double>* x,
                  int incx, const std::complex<double>* y, int incy,
                  std::complex<double>* A, int lda)
  {
    zgeru_(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
  }

  static long iamax(const long n, const std::complex<double> *x, const int incx) {
	return izamax_(&n, x, &incx);
  }

};

} // namespace TensorCalculus

#endif // __BLASINTERFACE_HPP
