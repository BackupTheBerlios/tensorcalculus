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

#ifndef __LAPACKINTERFACE_HPP
#define __LAPACKINTERFACE_HPP

#include <complex>
#include <cstring>

extern "C" int ilaenv_ (const int* ispec, const char* name, const char* opts,
                        const int* n1, const int* n2, const int* n3, const int* n4,
                        const int name_length, const int opts_length);

extern "C" void dgetri_(const int*, double*, const int*, const int*, double*,
                        const int*, int*);

extern "C" void dgetrf_(const int*, const int*, double*, const int*,
                        int*, int*);

extern "C" void dgeqp3_(const int*, const int*, double*, const int*, int*,
                        double*, double*, const int*, int*);

extern "C" void dorgqr_(const int*, const int*, const int*, double*, const int*,
                        const double*, double*, const int*, int*);

extern "C" void dsyev_(const char*, const char*, const int*, double*,
                       const int*, double*, double*, const int*, int*);

extern "C" void dgeqrf_(const int*, const int*, double*, const int*, double*,
                        double*, const int*, int*);

extern "C" void dgesvd_(const char*, const char*, const int*, const int*,
                        double*, const int*, double*, double*, const int*,
                        double*, const int*, double*, const int*, int*);

extern "C" void dpotrf_(const char*, const int*, double*, const int*, int*);

extern "C" void dpotri_(const char*, const int*, double*, const int*, int*);

extern "C" void sorgqr_(const int*, const int*, const int*, double*, const int*,
                        double*, double*, const int*, int*);

extern "C" void dgesv_(const int*, const int*, const double*, const int*,
		               const int*, const double*, const int*, const int*);

extern "C" void dstevx_ (const char* jobz, const char* range, const int* n,
                         double* d, double* e,
                         const double* vl, const double* vu,
                         const int* il, const int* iu,
                         const double* abstol,
                         int* m, double* w, double* z, const int* ldz,
                         double* work, int* iwork,
                         int* ifail, int* info);

extern "C" void dsyevx_ (const char* jobz, const char* range, const char* uplo,
                         const int* n, double* a, const int* lda,
                         const double* vl, const double* vu,
                         const int* il, const int* iu,
                         const double* abstol,
                         int* m, double* w, double* z, const int* ldz,
                         double* work, const int* lwork, int* iwork,
                         int* ifail, int* info);

extern "C" void dsygvx_ (const int* itype, const char* jobz, const char* range, const char* uplo,
                         const int* n, double* a, const int* lda, double* b, const int* ldb,
                         const double* vl, const double* vu,
                         const int* il, const int* iu,
                         const double* abstol,
                         int* m, double* w, double* z, const int* ldz,
                         double* work, const int* lwork, int* iwork,
                         int* ifail, int* info);

extern "C" double dlansy_ (const char * norm, const char * uplo, const int * n, const double * a,
		                   const int * lda, double * work);

extern "C" void dgecon_ (const char * norm, const int * n, const double * a, const int * lda,
		                 const double * anorm, double * rcond, double * work, int * iwork,
		                 int * info);

extern "C" double dlange_ (const char * norm, const int * m, const int * n, const double * a,
		                   const int * lda, double * work);

namespace TensorCalculus {

template<typename T>
struct Lapack;

template<>
struct Lapack<double> {

  static int ilaenv(const int ispec, char* name, const char* opts,
                    const int n1, const int n2, const int n3, const int n4) {
    name[0] = 'D';
    return ilaenv_(&ispec, name, opts, &n1, &n2, &n3, &n4,
                   std::strlen(name), std::strlen(opts));
  }

  static int getri(int n, double* A, int lda, const int* ipiv,
                   double* work, int lwork)
  {
    int info;
    dgetri_(&n, A, &lda, ipiv, work, &lwork, &info);
    return info;
  }

  static int getrf(int m, int n, double* A, int lda, int* ipiv)
  {
    int info;
    dgetrf_(&m, &n, A, &lda, ipiv, &info);
    return info;
  }

  static int geqp3(int m, int n, double* A, int lda, int* jpvt,
                   double* tau, double* work, int lwork)
  {
    int info;
    dgeqp3_(&m, &n, A, &lda, jpvt, tau, work, &lwork, &info);
    return info;
  }

  static int orgqr(int m, int n, int k, double* A, int lda,
                   const double* tau, double* work, int lwork)
  {
    int info;
    dorgqr_(&m, &n, &k, A, &lda, tau, work, &lwork, &info);
    return info;
  }

  static int syev(char job, char uplo, int n, double* A, int lda,
                  double* w, double* work, int lwork)
  {
    int info;
    dsyev_(&job, &uplo, &n, A, &lda, w, work, &lwork, &info);
    return info;
  }

  static int geqrf(int m, int n, double* A, int lda, double* tau,
                   double* work, int lwork)
  {
    int info;
    dgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);
    return info;
  }

  static int gesvd(char jobu, char jobvt, int m, int n, double* A,
                   int lda, double* s, double* u, int ldu,
                   double* vt, int ldvt, double* work, int lwork)
  {
    int info;
    dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, &info);
    return info;
  }

  static int potrf(char uplo, int n, double* A, int lda)
  {
    int info;
    dpotrf_(&uplo, &n, A, &lda, &info);
    return info;
  }

  static int potri(char uplo, int n, double* A, int lda)
  {
    int info;
    dpotri_(&uplo, &n,A, &lda, &info);
    return info;
  }

  static int sorgqr(int m, int n, int k, double* A, int lda,
                    double* tau, double* work, int lwork)
  {
    int info;
    sorgqr_(&m, &n, &k, A, &lda, tau, work, &lwork, &info);
    return info;
  }

  static int gesv(const int n, const int m, const double* A,
		           const int lda, int* ipiv, const double* B,
		           const int ldb)
  {
	int info;

	dgesv_(&n, &m, A, &lda, ipiv, B, &ldb, &info);
	return info;
  }

  static int syevx(char jobz, char range, char uplo,
                   int n, double* a, int lda,
                   double vl, double vu, int il, int iu, double abstol,
                   int& m, double* w, double* z, int ldz,
                   double* work, int lwork, int* iwork, int* ifail)
  {
    int info;
    dsyevx_(&jobz, &range, &uplo,
            &n, a, &lda,
            &vl, &vu, &il, &iu, &abstol,
            &m, w, z, &ldz,
            work, &lwork, iwork, ifail, &info);
    return info;
  }

  static int stevx(char jobz, char range,
                   int n, double* d, double* s,
                   double vl, double vu, int il, int iu, double abstol,
                   int& m, double* w, double* z, int ldz,
                   double* work, int* iwork, int* ifail)
  {
    int info;
    dstevx_(&jobz, &range,
            &n, d, s,
            &vl, &vu, &il, &iu, &abstol,
            &m, w, z, &ldz,
            work, iwork, ifail, &info);
    return info;
  }

  static int sygvx(int itype, char jobz, char range, char uplo,
                   int n, double* a, int lda, double* b, int ldb,
                   double vl, double vu,
                   int il, int iu,
                   double abstol,
                   int& m, double* w, double* z, int ldz,
                   double* work, int lwork, int* iwork, int* ifail) {
    int info;
    dsygvx_(&itype, &jobz, &range, &uplo,
            &n, a, &lda, b, &ldb,
            &vl, &vu,
            &il, &iu,
            &abstol,
            &m, w, z, &ldz,
            work, &lwork, iwork, ifail, &info);
    return info;
  }

  static double lansy(const char norm, const char uplo, const int n,
		            const double * a, const int lda, double * work)
  {
	return dlansy_(&norm, &uplo, &n, a, &lda, work);
  }

  static double gecon(const char norm, const int n, const double * a, const int lda,
                       const double anorm, double * work, int * iwork)
  {
	int info;

	double rcond;

	dgecon_(&norm, &n, a, &lda, &anorm, &rcond, work, iwork, &info);
	if (info < 0) {
	  // do something
	}
	return rcond;
  }

  static double lange(const char norm, const int m, const int n, const double * a,
          const int lda, double * work) {
	return dlange_(&norm, &m, &m, a, &lda, work);
  }

};

} // namespace TensorCalculus

extern "C" void zgetri_(const int*, std::complex<double>*, const int*, const int*,
                        std::complex<double>*, const int*, int*);

extern "C" void zgetrf_(const int*, const int*, std::complex<double>*, const int*,
                        int*, int*);

extern "C" void zgeqrf_(const int*, const int*, std::complex<double>*, const int*,
                        std::complex<double>*, std::complex<double>*, const int*, int*);

extern "C" void zgesvd_(const char*, const char*, const int*, const int*,
                        std::complex<double>*, const int*, double*, std::complex<double>*, const int*,
                        std::complex<double>*, const int*, std::complex<double>*, const int*,
                        double*, int*);

extern "C" void zgeqp3_(const int*, const int*, std::complex<double>*, const int*, int*,
                        std::complex<double>*, std::complex<double>*, const int*, double*, int*);

extern "C" void zungqr_(const int*, const int*, const int*, std::complex<double>*,
                        const int*, const std::complex<double>*, std::complex<double>*,
                        const int*, int*);

namespace TensorCalculus {

template<>
struct Lapack< std::complex<double> > {

  static int ilaenv(const int ispec, char* name, const char* opts,
                    const int n1, const int n2, const int n3, const int n4) {
    name[0] = 'Z';
    return ilaenv_(&ispec, name, opts, &n1, &n2, &n3, &n4,
                   std::strlen(name), std::strlen(opts));
  }

  static int getri(int n, std::complex<double>* A, int lda, const int* ipiv,
                   std::complex<double>* work, int lwork)
  {
    int info;
    zgetri_(&n, A, &lda, ipiv, work, &lwork, &info);
    return info;
  }

  static int getrf(int m, int n, std::complex<double>* A, int lda, int* ipiv)
  {
    int info;
    zgetrf_(&m, &n, A, &lda, ipiv, &info);
    return info;
  }

  static int geqrf(int m, int n, std::complex<double>* A, int lda, std::complex<double>* tau,
                   std::complex<double>* work, int lwork)
  {
    int info;
    zgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);
    return info;
  }

  static int ungqr(int m, int n, int k, std::complex<double>* A, int lda,
                   const std::complex<double>* tau, std::complex<double>* work, int lwork)
  {
    int info;
    zungqr_(&m, &n, &k, A, &lda, tau, work, &lwork, &info);
    return info;
  }

  static int gesvd(char jobu, char jobvt, int m, int n, std::complex<double>* A,
                   int lda, double* s, std::complex<double>* u, int ldu, std::complex<double>* vt,
                   int ldvt, std::complex<double>* work, int lwork, double* rwork)
  {
    int info;
    zgesvd_(&jobu, &jobvt, &m, &n, A, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, rwork, &info);
    return info;
  }

  static int geqp3(int m, int n, std::complex<double>* A, int lda, int* jpvt,
                   std::complex<double>* tau, std::complex<double>* work, int lwork, double* rwork)
  {
    int info;
    zgeqp3_(&m, &n, A, &lda, jpvt, tau, work, &lwork, rwork, &info);
    return info;
  }

};

} // namespace TensorCalculus

#endif // __LAPACKINTERFACE_HPP
