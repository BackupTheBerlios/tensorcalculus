/*
 * Copyright (C) 2011 Henry Auer, Stefan Handschuh
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

#ifndef __MATRIX_OPERATORS_HPP
#define __MATRIX_OPERATORS_HPP

#include <vector>
#include "Vector/VectorOperators.hpp"
#include "BlasInterface.hpp"
#include "LapackInterface2.hpp"


namespace TensorCalculus {
	template<typename T>
	struct minvert_result_t {
		T norm_sqr;

		T norm_inv_sqr;

		bool approx;
	};

  
  template<typename T>
  void transpose (const int m, const int n, T* Ain, T* Aout){
   // getestet
   #pragma omp parallel for
   for(int i=0; i<m; i++){
    Blas<T>::copy(n, &Ain[i], m, &Aout[i*n], 1);
   }
  }

  template<typename T>
    void transpose (const int m, const int n, T* Ain, const int ldain, T* Aout, const int ldaout){
     // getestet
     #pragma omp parallel for
     for(int i=0; i<m; i++){
      Blas<T>::copy(n, &Ain[i], ldain, &Aout[i*ldaout], 1);
     }
    }


  
  template<typename T>
    void transpose (const int m, const int n, const T alpha, T* Ain, T* Aout){
     // todo: optimieren
     // getestet
     Blas<T>::scal(m*n,alpha, Ain, 1);
     transpose(m,n,Ain,Aout);
    }
  
  template<typename T>
  T checkSVD(int m, int n, const std::vector<T> &A,
                      int lda, T* s, T* u, int ldu,
                      T* vt, int ldvt){
			T info;
			int min = std::min(m,n);
			std::vector<T> B(m*n);
			for(int i=0; i<(min); i++){
				Blas<T>::ger(m,n,s[i],&u[i*ldu],1,&vt[i],ldvt,&B[0],m);
			}
			//std::cout<<A<<std::endl<<B<<std::endl;
			using namespace VectorOperators;
			info = l2_norm(A-B);
				return info;
  }

  template<typename T>

  void setSizesForSVD (const char jobu, const char jobvt, int m, int n, std::vector<T> &A,
        int lda, std::vector<T> &S, std::vector<T> &U, int ldu, std::vector<T> &VT, int ldvt,
        std::vector<T> &work, int &lwork){
     int min = std::min(m,n);
     int uSize;
     int vtSize;
     if(jobu == 'S' || jobu == 's'){
        uSize  = m * min;
        if(U.size()<uSize){
           U.resize(uSize);
        }
     }
     if(jobu == 'A' || jobu == 'a'){
       uSize  = m * m;
       if(U.size()<uSize){
          U.resize(uSize);
       }
    }
    if(jobvt == 'S' || jobvt == 's'){
       vtSize = min * n;
       if(VT.size()<vtSize){
          VT.resize(vtSize);
       }
    }
    if(jobvt == 'A' || jobvt == 'a'){
       vtSize = n * n;
       if(VT.size()<vtSize){
          VT.resize(vtSize);
       }
    }
    if(S.size()<min){
       S.resize(min);
    }
    Lapack<T>::gesvd(jobu, jobvt, m, n, &A[0], m, &S[0], &U[0], m, &VT[0], min, &work[0], -1);
    if(work[0] > lwork){
       lwork = work[0];
       work.resize(lwork);
    }
  }

  template<typename T>
  std::vector<T> giveQ (int m, int n, T* a, int lda, T* tau){
     //returns matrix q form elementary reflectors produced by dgeqrf
     std::vector<T> Q(lda*n);
     Blas<T>::copy(lda*n, a, 1, &Q[0], 1);
     int k = std::min(m,n);
     std::vector<T> work(1);
     int lwork = -1;
     Lapack<T>::orgqr(m, n, k, &Q[0], lda, tau, &work[0], lwork);
     lwork = work[0];
     work.resize(lwork);
     Lapack<T>::orgqr(m, n, k, &Q[0], lda, tau, &work[0], lwork);

    return(Q);
  }

  template<typename T>
  void printMatrix(const int m, const int n, T* A, int lda){
     // print vector as columnwise matrix
     //std::cout.setf(std::ios::scientific, std::ios::floatfield);
     //std::cout.precision(4);

     for (int i=0; i< m; i++){
        for(int j=0; j< n; j++){
           std::cout<< A[i+j*lda]<<"   ";
        }
        std::cout<<std::endl;
     }
     //std::cout.unsetf(std::ios::scientific, std::ios::floatfield);
  }

  template<typename T>
      void printMatrix(const int m, const int n, T* A){
       printMatrix(m, n, A, m);
    }

  template<typename T>
   void printMatrix(const int m, const int n, T* A, int lda, T eps){
    // print vector as columnwise matrix
    // if A[index] < eps --> print 0
    T value;
    for (int i=0; i< m; i++){
       for(int j=0; j< n; j++){
          value = A[i+j*lda];
          if(value < eps){
             std::cout<< 0 <<"   ";
          }else{
             std::cout<< value <<"   ";
          }
       }
       std::cout<<std::endl;
    }
   }

  template<typename T>
    void printMatrix(const int m, const int n, T* A, T eps){
     printMatrix(m, n, A, m, eps);
  }
/* use invert_matrix
  template <typename T>
  const int invertA (int m, int n, std::vector<T> &A, int lda){
     //tested for n x n matrices
     int info;
     int min = std::min(m,n);
     std::vector<int> ipiv(min);

     info = Lapack<T>::getrf(m, n, &A[0], lda, &ipiv[0]);

     if(info){
        std::cout<<" A is singular"<<std::endl;
        return(info);
     }

     std::vector<T> work(1);
     int lwork = -1;

     Lapack<T>::getri(min, &A[0], lda, &ipiv[0], &work[0],lwork);

     lwork = work[0];
     work.resize(lwork);

     info = Lapack<T>::getri(min, &A[0], lda, &ipiv[0], &work[0],lwork);

     if(info){
       std::cout<<" A is singular"<<std::endl;
       return(info);
     }

     return(info);
  }
*/

  /**
   * This method performs an adaptive cross approximation of the given matrix
   * until a certain rank is reached.
   *
   * m row count
   * n col count
   * rows.size() >= rank*n
   * cols.size() >= rank*m
   *
   * cols and rows will be completely overwritten
   */
  template<typename T>
  int aca(std::vector<T> &matrix, int m, int n, int rank, std::vector<T> &cols, std::vector<T> &rows) {
	  long max_index;

	  int col, row;

	  T max;

	  for (int c = 0; c < rank; c++) {
		  max_index = Blas<T>::iamax(m*n, &matrix[0], 1)-1;

		  max = matrix[max_index];
		  if (std::abs(max) < 1e-20) {
			  return c;
		  }
		  col = max_index / m;
		  row = max_index - col* m;
		  Blas<T>::copy(m, &matrix[col*m], 1, &cols[c*m], 1);
		  Blas<T>::scal(m, 1.0/max, &cols[c*m], 1);
		  Blas<T>::copy(n, &matrix[row], m, &rows[c*n], 1);
		  Blas<T>::ger(m, n, -1.0, &cols[c*m], 1, &rows[c*n], 1, &matrix[0], m);
	  }
	  return rank;
  }

  /**
   * This method performs an adaptive cross approximation of the given matrix
   * until a certain epsilon is reached.
   *
   * m row count
   * n col count
   * rows.size() >= min(n,m)*n
   * cols.size() >= min(n,m)*m
   *
   * cols and rows will be completely overwritten
   */
  template<typename T>
  int aca(std::vector<T> &matrix, int m, int n, T eps, std::vector<T> &cols, std::vector<T> &rows) {
	  long max_index;

	  int col, row, rank = std::min(m,n), mn = m*n;

	  T max;

	  for (int c = 0; c < rank; c++) {
		  max_index = Blas<T>::iamax(mn, &matrix[0], 1)-1;

		  max = matrix[max_index];
		  if (std::abs(max) < 1e-80) {
			  return c;
		  }
		  col = max_index / m;
		  row = max_index - col* m;
		  Blas<T>::copy(m, &matrix[col*m], 1, &cols[c*m], 1);
		  Blas<T>::scal(m, 1.0/max, &cols[c*m], 1);
		  Blas<T>::copy(n, &matrix[row], m, &rows[c*n], 1);
		  Blas<T>::ger(m, n, -1.0, &cols[c*m], 1, &rows[c*n], 1, &matrix[0], m);

		  if (Blas<T>::nrm2(mn, &matrix[0], 1) < eps) {
			  return c+1;
		  }
	  }
	  return rank;
  }


  /*
   * m row count
   * n col count
   */
  template<typename T>
  int invert_matrix(std::vector<T> & matrix, const int m, const int n, minvert_result_t<T> * minvert_result = 0) {
	  T COND_BORDER = 1e-19;

	  int min = std::min(n, m);

	  std::vector<T> workspace(4*min);

	  T h_norm = Lapack<T>::lange('1', m, n, &matrix[0], m, 0); // 1-norm

	  std::vector<int> ipiv(min);

	  std::vector<T> matrix_copy;

	  if (minvert_result != 0) { // we might need the original matrix
		  matrix_copy = std::vector<T>(matrix);
	  }
	  Lapack<T>::getrf(m, n, &matrix[0], m, &ipiv[0]);// LU factorization

	  std::vector<int> iwork(min);

	  T cond = Lapack<T>::gecon('1', min, &matrix[0], min, h_norm, &workspace[0], &iwork[0]);

	  T workspaceLength = 0;

	  if (cond > COND_BORDER) { // otherwise the matrix is almost singular
	  	Lapack<T>::getri(min, &matrix[0], min, &ipiv[0], &workspaceLength, -1);

	  	if (4*min < workspaceLength) {
	  	  workspace.resize((int) workspaceLength);
	  	}
	  	if (minvert_result != 0) {
	  		minvert_result->approx = false;
	  	}
		return Lapack<T>::getri(min, &matrix[0], min, &ipiv[0], &workspace[0], (int) workspaceLength); // calculate the inverse
	  } else if (minvert_result != 0) { // compute pseudo-inverse with SVD
		std::vector<T> S(min);
		std::vector<T> U(n*n);
		std::vector<T> VT(m*m);

		Lapack<T>::gesvd('S', 'S', m, n, &matrix_copy[0], m, &S[0], &U[0], m,
		                 &VT[0], min, &workspaceLength, -1);

		if (workspace.size() < workspaceLength) {
			workspace.resize((int) workspaceLength);
		}

		Lapack<T>::gesvd('S', 'S', m, n, &matrix_copy[0], m, &S[0], &U[0], m,
		                 &VT[0], min, &workspace[0], (int)workspaceLength);

		Blas<T>::scal(n*m, 0, &matrix[0], 1);// reset the matrix

		T norm_svd_sqr = 0.0, norm_inv_sqr = 0.0;

		int newRank;

		for (newRank = 0; newRank < min; newRank++) {
			norm_svd_sqr += S[newRank]*S[newRank];
			norm_inv_sqr += 1.0 / (S[newRank]*S[newRank]);
			if (1.0 / std::sqrt(norm_svd_sqr*norm_inv_sqr) < COND_BORDER) {
				newRank--;
				break;
			}
		}
		/* now we have figured out the maximal rank that fulfills our condition number */

		minvert_result->norm_sqr = norm_svd_sqr;
		minvert_result->norm_inv_sqr = norm_inv_sqr;
		minvert_result->approx = true;

		if (newRank < 1 || newRank == min) {
		  return -n*m;
		}

		for (int l = 0; l < newRank; l++) {
		  Blas<T>::ger(m, n, 1.0/S[l], &U[l*m], 1, &VT[l], m, &matrix[0], m);
		}
		return 0;
	  } else {
		  std::cout << "condition " << cond << std::endl;
		return -n*m;
	  }
  }

  template<typename T>
  int invert_matrix(std::vector<T> &matrix, const int n, minvert_result_t<T> * minvert_result = 0) {
	return invert_matrix(matrix, n, n, minvert_result);
  }

  /*
   * n = rowcount
   * m = colcount
   */
  template<typename T>
  void delete_rows(std::vector<T> &matrix, int n, int m, const std::vector<int>& rows) {
	  int i = rows.size()-1;

	  std::vector<int> min(rows);
	  std::sort(min.begin(), min.end());

	  int TSIZE = sizeof(T);

	  for (int k = i; k > -1; k--) {
		  /* remove column j and row j */
		  const int j = min[k];

		  for (int l = m-1; l > -1; l--) {
			  std::memmove(&matrix[l*n+j], &matrix[l*n+j+1] , TSIZE*(n-j-1+(n-1)*(m-l-1)));
		  }
		  n--;
	  }
	  matrix.resize(n*m);
  }

  /*
   * n = rowcount
   * m = colcount
   */
  template<typename T>
  void delete_cols(std::vector<T> &matrix, int n, int m, const std::vector<int>& cols) {
	  int i = cols.size()-1;

	  std::vector<int> min(cols);
	  std::sort(min.begin(), min.end());

	  int TSIZE = sizeof(T);

	  for (int k = i; k > -1; k--) {
		  const int j = min[k];

		  if (j < m-1) {
			  std::memmove(&matrix[n*j], &matrix[n*(j+1)] , n*(m-j-1)*TSIZE);
		  }
		  m--;
	  }
	  matrix.resize(n*m);
  }

  /*
   * Complexity (n-1)*indices.size()+indices.size() = n*indices.size()
   */
  template<typename T>
  void delete_matrix_indices(std::vector<T> &matrix, int n, int m, std::vector<int>& indices) {
	  int i = indices.size()-1;

	  std::sort(indices.begin(), indices.end());

	  int TSIZE = sizeof(T);

	  for (int k = i; k > -1; k--) {
	      /* remove column j  */
		  const int j = indices[k];

		  if (j < m-1) {
			  std::memmove(&matrix[n*j], &matrix[n*(j+1)] , n*(m-j-1)*TSIZE);
		  }
		  m--;

	  }
	  for (int k = i; k > -1; k--) {
		  /* remove row j */
		  const int j = indices[k];

		  for (int l = m-1; l > -1; l--) {
			  std::memmove(&matrix[l*n+j], &matrix[l*n+j+1] , TSIZE*(n-j-1+(n-1)*(m-l-1)));
		  }
		  n--;
	  }
	  matrix.resize(n*m);
  }

  template<typename T>
  int invert_gramian_matrix(std::vector<T> &matrix, int n, std::vector<int> & dependentIndices) {
	  /*
	   * Checking the linearly independent rows using LU-factorization
	   */

	  std::vector<int> ipiv(n);

	  std::vector<T> matrix_copy(matrix); // this can be handled more cleverly; maybe if there is some free time, FIXME

	  Lapack<T>::getrf(n, n, &matrix[0], n, &ipiv[0]); // LU factorization

	  /* ROW-interchange, not column-interchange
	  std::vector<int> interchange(n);

	  for (int k = 1; k < n; k++) {
		  interchange[k] = k;
	  }
	  for (int k = 0; k < n-1; k++) {
		  std::swap(interchange[ipiv[k]-1], interchange[k]);
	  }
	  */

	  for (int k = 0; k < n; k++) {
		  if (std::abs(matrix[0])*1e-15 > std::abs(matrix[k*n+k]) && indexOf(dependentIndices, k) == -1) {
			  dependentIndices.push_back(k);
		  }
	  }

	  T workspaceLength;

	  /*
	   * If the matrix is regular, we can invert it in its current form, otherwise,
	   * we have to delete the linearly dependent rows/columns.
	   */
	  int dependencySize = dependentIndices.size();

	  if (dependencySize > 0) {
		  delete_matrix_indices(matrix_copy, n, n, dependentIndices); // since this matrix is already permuted

		  n -= dependencySize; // since n is only a copy, we can change it
		  std::memcpy(&matrix[0], &matrix_copy[0], n*n*sizeof(T));
		  matrix.resize(n*n);
		  Lapack<T>::getrf(n, n, &matrix[0], n, &ipiv[0]); // LU factorization
	  }

	  Lapack<T>::getri(n, &matrix[0], n, &ipiv[0], &workspaceLength, -1);

	  std::vector<T> workspace((int) workspaceLength);
	  return Lapack<T>::getri(n, &matrix[0], n, &ipiv[0], &workspace[0], (int) workspaceLength); // calculate the inverse
  }



  /*
   * Uses SVD to determine the matrix rank.
   *
   * This method should be only used for testing purposes.
   */
  template<typename T>
  int getMatrixRank(const std::vector<T> &matrix, const int m, const int n) {
	int min = std::min(n, m); // the rank is at least

	std::vector<T> matrix_copy(matrix);

	std::vector<T> S(min);

	T wl;

	Lapack<T>::gesvd('N', 'N', m, n, &matrix_copy[0], m, &S[0], 0, m, 0, min, &wl, -1);

	std::vector<T> work((int) wl);

	Lapack<T>::gesvd('N', 'N', m, n, &matrix_copy[0], m, &S[0], 0, m, 0, min, &work[0], (int) wl);

	std::vector<T> partialSums(min);

	partialSums[m-1] = S[min-1]*S[min-1];

	for (int n = min-2; n > -1; n--) {
		partialSums[n] = S[n]*S[n] + partialSums[n-1];
	}
	using namespace VectorOperators;
	std::cout << S << std::endl;
	T sum = 0;

	for (int n = 0; n < min-1; n++) {
		if (sum + partialSums[n+1] == sum) {
			return n;
		} else {
			sum += S[n]*S[n];
	  	}
	}
	return min;

  }
} // namespace TensorCalculus

#endif // __MATRIX_OPERATORS_HPP
