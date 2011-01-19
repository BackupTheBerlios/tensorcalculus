/*
 * Copyright (C) 2010 Philipp Wähnert
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

#ifndef __SIMPLEX_HPP
#define __SIMPLEX_HPP

#include <vector>
#include <algorithm>
#include <iterator>
#include <stdexcept>
// #include <functional>
#include <ostream>
#include <cmath>

#include "BlasInterface.hpp"
#include "LapackInterface2.hpp"
#include "Vector/VectorOperators.hpp"
#include "Utilities/Factorials.hpp"

namespace TensorCalculus {

  class SimplexError : std::invalid_argument {
  public:
    SimplexError(const std::string& what) : std::invalid_argument(what) { }
  };

  template<typename T>
  class Simplex {
  private:
    int dimension;
    std::vector<T> A;
    std::vector<T> x0;
    std::vector< std::vector<T> > base;
    T determinant;
    // T area;

    void calc_determinant/*_and_area*/();
    void calc_base_vectors(/*const std::vector<int>& border*/);

    void fill_points(const std::vector< std::vector<T> >& thePoints) {
      x0 = thePoints[0];
      A.reserve(dimension * dimension);
      for (int i = 1; i <= dimension; ++i) {
        const std::vector<T>& xi = thePoints[i];
        for (int j = 0; j < dimension; ++j) {
          A.push_back(xi[j] - x0[j]);
        }
      }
    }

    Simplex(int dimension, const std::vector<T>& A, const std::vector<T>& x0, T determinant/*, T area*/)
      : dimension(dimension),
        A(A), x0(x0),
        determinant(determinant)/*,
        area(area)*/
    {
      calc_base_vectors();
    }

  public:
    Simplex(const std::vector< std::vector<T> >& thePoints)
      : dimension(thePoints.size() - 1)
    {
      fill_points(thePoints);
      calc_determinant/*_and_area*/();
      calc_base_vectors(/*std::vector<int>(dimension + 1, 0)*/);
    }

    // Simplex(const std::vector< std::vector<T> >& thePoints, const std::vector<int>& border)
    //   : dimension(thePoints.size() - 1)
    // {
    //   fill_points(thePoints);
    //   calc_determinant_and_area();
    //   calc_base_vectors(border);
    // }

    /**
     * @brief   Transforms the coordinates of the standard \f$d\f$-Simplex into actual coordinates
     */
    std::vector<T> transform(const std::vector<T>& standard) const {
#ifdef RANGE_CHECKS_ON
      if (standard.size() != dimension) {
        throw SimplexError("Dimensions mismatch");
      }
#endif // RANGE_CHECKS_ON
      std::vector<T> result(dimension);
      Blas<T>::gemv('N', dimension, dimension,
                    1.0, &A[0], dimension,
                    &standard[0], 1,
                    0.0, &result[0], 1);
      Blas<T>::axpy(dimension, 1.0, &x0[0], 1, &result[0], 1);
      return result;
    }

    // std::vector<Simplex> dissect() const;

    T get_determinant() const { return determinant; }
    // T get_area() const { return area; }
    T get_dimension() const { return dimension; }
    
    const std::vector<T>& get_matrix() const { return A; }
    const std::vector<T>& get_origin() const { return x0; }
    const std::vector<T>& get_base(int index) const { return base[index]; }
    
    void print(std::ostream& stream) const;
    
    std::vector< Simplex<T> > cross_product(const Simplex<T>& b) const;
  };

  template<typename T>
  void Simplex<T>::calc_determinant/*_and_area*/() {
    if (dimension == 2) { // fast calculations for d = 2
      determinant = (A[0] * A[3]) - (A[2] * A[1]);
      // using std::abs;
      // area = 0.5 * abs(determinant);
    } else { // else use the LU factorization to obtain the determinant
      std::vector<T> tempA(A);
      std::vector<int> pivots(dimension);
      if (Lapack<T>::getrf(dimension, dimension, &tempA[0], dimension, &pivots[0]) != 0) {
        throw SimplexError("Error in LU decomposition");
      }

      determinant = 1;
      int diagonal = 0;
      for (int i = 0; i < dimension; ++i) {
        determinant *= tempA[diagonal];
        diagonal += dimension + 1;
      }

      // using std::abs;
      // area = abs(determinant) / factorial<T>(dimension);
    }
  }

  // template<typename T>
  // std::vector< Simplex<T> > Simplex<T>::dissect() const  {
  //   if (dimension == 2) {
  //     std::vector<Simplex> result;
  //     result.reserve(4);

  //     using VectorOperators::operator /=;
  //     std::vector<T> new_A = A;
  //     new_A /= 2.0;
  //     result.push_back(Simplex(dimension, new_A, x0, determinant/4.0, area/4.0));

  //     std::vector<T> mid_x0_x1(2);
  //     mid_x0_x1[0] = new_A[0] + x0[0];
  //     mid_x0_x1[1] = new_A[1] + x0[1];
  //     result.push_back(Simplex(dimension, new_A, mid_x0_x1, determinant/4.0, area/4.0));

  //     std::vector<T> mid_x0_x2(2);
  //     mid_x0_x2[0] = new_A[2] + x0[0];
  //     mid_x0_x2[1] = new_A[3] + x0[1];
  //     result.push_back(Simplex(dimension, new_A, mid_x0_x2, determinant/4.0, area/4.0));

  //     new_A[0] = new_A[2] - new_A[0];
  //     new_A[1] = new_A[3] - new_A[1];
  //     result.push_back(Simplex(dimension, new_A, mid_x0_x1, determinant/4.0, area/4.0));

  //     return result;
  //   } else {
  //     throw SimplexError("Not yet implemented - how do you dissect a higher dimensional simplex?");
  //   }
  // }

  template<typename T>
  inline void Simplex<T>::print(std::ostream& stream) const {
    int x0_index = 0;
    int A_index = 0;
    stream << "( ( ";
    using VectorOperators::operator <<;
    stream << x0 << ") ";

    for (int i = 0; i < dimension; ++i) {
      stream << "( ";
      for (int j = 0; j < dimension; ++j) {
        stream << A[A_index] + x0[x0_index] << " ";
        ++A_index;
        ++x0_index;
      }
      stream << ") ";
      x0_index = 0;
    }
    stream << determinant << " )";
  }

  template<typename T>
  std::ostream& operator << (std::ostream& stream, const Simplex<T>& simplex) {
    simplex.print(stream);
    return stream;
  }
  
  template<typename T>
  void Simplex<T>::calc_base_vectors(/*const std::vector<int>& border*/)
  {
    std::vector<T> Ainv(A);
    std::vector<int> pivots(dimension);
    if (Lapack<T>::getrf(dimension, dimension, &Ainv[0], dimension, &pivots[0]) != 0) {
      throw std::runtime_error("Error in \"?getrf\"");
    }
    double dbl_lwork;
    if (Lapack<T>::getri(dimension, &Ainv[0], dimension, &pivots[0], &dbl_lwork, -1) != 0) {
      throw std::runtime_error("Error in workspace query in \"?getri\"");
    }
    int lwork = static_cast<int>(dbl_lwork);
    std::vector<double> work(lwork);
    if (Lapack<T>::getri(dimension, &Ainv[0], dimension, &pivots[0], &work[0], lwork) != 0) {
      throw std::runtime_error("Error in \"?getri\"");
    }
        
    base.assign(dimension + 1, std::vector<T>(dimension, 0));
    for (int i = 1; i <= dimension; ++i) {
      // if (border[i] == 0) {
        Blas<T>::copy(dimension, &Ainv[i - 1], dimension, &base[i][0], 1);
      // }
      // if (border[0] == 0) {
        Blas<T>::axpy(dimension, -1.0, &Ainv[i - 1], dimension, &base[0][0], 1);
      // }
    }
  }
  
  template<typename T>
  inline std::vector< Simplex<T> > cross_product(const Simplex<T>& a, const Simplex<T>& b)
  {
    return a.cross_product(b);
  }
  
  template<typename T>
  std::vector< Simplex<T> > Simplex<T>::cross_product(const Simplex<T>& b) const
  {
    if (get_dimension() != 2 || b.get_dimension() != 2) {
      throw SimplexError("Not yet implemented");
    }
    // Triangulation from:
    // Gelfand, Kapranov, Zelevinsky: Discriminants, cross_resultants, and multidimensional determinants, p. 250
    std::vector< Simplex<T> > cross_result;
    cross_result.reserve(6);
    
    std::vector<T> A_new(4*dimension*dimension);
    std::vector<T> x0_new(2*dimension);
    { // 1)
      T trafo[4] = { 1, -1, 0, -1 };
      std::vector<T> A_trafo(dimension*dimension);
      Blas<T>::gemm('N', 'N', dimension, dimension, dimension,
                    1.0, &A[0], dimension, &trafo[0], dimension,
                    0.0, &A_trafo[0], dimension);
      for (int i = 0; i < dimension; ++i) {
        Blas<T>::copy(dimension, &b.A[i*dimension], 1, &A_new[i*2*dimension + dimension], 1);
        Blas<T>::copy(dimension, &b.A[i*dimension], 1, &A_new[i*2*dimension + 2*dimension*dimension + dimension], 1);
        Blas<T>::copy(dimension, &A_trafo[i*dimension], 1, &A_new[i*2*dimension + 2*dimension*dimension], 1);
      }
      Blas<T>::copy(dimension, &x0[0], 1, &x0_new[0], 1);
      Blas<T>::axpy(dimension, 1, &A[dimension], 1, &x0_new[0], 1);
      Blas<T>::copy(dimension, &b.x0[0], 1, &x0_new[dimension], 1);
      cross_result.push_back(Simplex<T>(2*dimension, A_new, x0_new, (dimension % 2 == 0 ? 1 : -1) * (-1) * determinant * b.determinant));
    }
    { // 2)
      A_new.assign(4*dimension*dimension, 0);
#ifndef NDEBUG
      x0_new.assign(2*dimension, 0);
#endif
      T trafo[4] = { 1, -1, 0, -1 };
      std::vector<T> A_trafo(dimension*dimension);
      Blas<T>::gemm('N', 'N', dimension, dimension, dimension,
                    1.0, &b.A[0], dimension, &trafo[0], dimension,
                    0.0, &A_trafo[0], dimension);
      for (int i = 0; i < dimension; ++i) {
        Blas<T>::copy(dimension, &A[i*dimension], 1, &A_new[i*2*dimension], 1);
        Blas<T>::copy(dimension, &A[i*dimension], 1, &A_new[i*2*dimension + 2*dimension*dimension], 1);
        Blas<T>::copy(dimension, &A_trafo[i*dimension], 1, &A_new[i*2*dimension + 2*dimension*dimension + dimension], 1);
      }
      Blas<T>::copy(dimension, &x0[0], 1, &x0_new[0], 1);
      Blas<T>::copy(dimension, &b.x0[0], 1, &x0_new[dimension], 1);
      Blas<T>::axpy(dimension, 1, &b.A[dimension], 1, &x0_new[dimension], 1);
      cross_result.push_back(Simplex<T>(2*dimension, A_new, x0_new, determinant * (-1) * b.determinant));
    }
    { // 3)
      A_new.assign(4*dimension*dimension, 0);
#ifndef NDEBUG
      x0_new.assign(2*dimension, 0);
#endif
      T trafo1[4] = { -1, 0, -1, 1 };
      std::vector<T> A1_trafo(dimension*dimension);
      Blas<T>::gemm('N', 'N', dimension, dimension, dimension,
                    1.0, &A[0], dimension, &trafo1[0], dimension,
                    0.0, &A1_trafo[0], dimension);
      T trafo2[4] = { -1, 1, -1, 0 };
      std::vector<T> A2_trafo(dimension*dimension);
      Blas<T>::gemm('N', 'N', dimension, dimension, dimension,
                    1.0, &b.A[0], dimension, &trafo2[0], dimension,
                    0.0, &A2_trafo[0], dimension);
      for (int i = 0; i < dimension; ++i) {
        Blas<T>::copy(dimension, &A1_trafo[i*dimension], 1, &A_new[i*2*dimension], 1);
        Blas<T>::copy(dimension, &A1_trafo[i*dimension], 1, &A_new[i*2*dimension + 2*dimension*dimension], 1);
        Blas<T>::copy(dimension, &A2_trafo[i*dimension], 1, &A_new[i*2*dimension + 2*dimension*dimension + dimension], 1);
      }
      Blas<T>::copy(dimension, &x0[0], 1, &x0_new[0], 1);
      Blas<T>::axpy(dimension, 1, &A[0], 1, &x0_new[0], 1);
      Blas<T>::copy(dimension, &b.x0[0], 1, &x0_new[dimension], 1);
      Blas<T>::axpy(dimension, 1, &b.A[0], 1, &x0_new[dimension], 1);
      cross_result.push_back(Simplex<T>(2*dimension, A_new, x0_new, (-1) * determinant * b.determinant));
    }
    { // 4)
      A_new.assign(4*dimension*dimension, 0);
      T trafo1[4] = { -1, 1, -1, 0 };
      std::vector<T> A1_trafo(dimension*dimension);
      Blas<T>::gemm('N', 'N', dimension, dimension, dimension,
                    1.0, &A[0], dimension, &trafo1[0], dimension,
                    0.0, &A1_trafo[0], dimension);
      T trafo2[4] = { -1, 0, -1, 1 };
      std::vector<T> A2_trafo(dimension*dimension);
      Blas<T>::gemm('N', 'N', dimension, dimension, dimension,
                    1.0, &b.A[0], dimension, &trafo2[0], dimension,
                    0.0, &A2_trafo[0], dimension);
      for (int i = 0; i < dimension; ++i) {
        Blas<T>::copy(dimension, &A2_trafo[i*dimension], 1, &A_new[i*2*dimension + dimension], 1);
        Blas<T>::copy(dimension, &A2_trafo[i*dimension], 1, &A_new[i*2*dimension + 2*dimension*dimension + dimension], 1);
        Blas<T>::copy(dimension, &A1_trafo[i*dimension], 1, &A_new[i*2*dimension + 2*dimension*dimension], 1);
      }
      cross_result.push_back(Simplex<T>(2*dimension, A_new, x0_new, (dimension % 2 == 0 ? 1 : -1) * determinant * (-1) * b.determinant));
    }
    { // 5)
      A_new.assign(4*dimension*dimension, 0);
#ifndef NDEBUG
      x0_new.assign(2*dimension, 0);
#endif
      T trafo1[4] = { 0, 1, 1, 0 };
      std::vector<T> A1_trafo(dimension*dimension);
      Blas<T>::gemm('N', 'N', dimension, dimension, dimension,
                    1.0, &A[0], dimension, &trafo1[0], dimension,
                    0.0, &A1_trafo[0], dimension);
      T trafo2[4] = { 0, -1, 1, -1 };
      std::vector<T> A2_trafo(dimension*dimension);
      Blas<T>::gemm('N', 'N', dimension, dimension, dimension,
                    1.0, &b.A[0], dimension, &trafo2[0], dimension,
                    0.0, &A2_trafo[0], dimension);
      for (int i = 0; i < dimension; ++i) {
        Blas<T>::copy(dimension, &A2_trafo[i*dimension], 1, &A_new[i*2*dimension + dimension], 1);
        Blas<T>::copy(dimension, &A2_trafo[i*dimension], 1, &A_new[i*2*dimension + 2*dimension*dimension + dimension], 1);
        Blas<T>::copy(dimension, &A1_trafo[i*dimension], 1, &A_new[i*2*dimension + 2*dimension*dimension], 1);
      }
      Blas<T>::copy(dimension, &x0[0], 1, &x0_new[0], 1);
      Blas<T>::copy(dimension, &b.x0[0], 1, &x0_new[dimension], 1);
      Blas<T>::axpy(dimension, 1, &b.A[dimension], 1, &x0_new[dimension], 1);
      cross_result.push_back(Simplex<T>(2*dimension, A_new, x0_new, (dimension % 2 == 0 ? 1 : -1) * (-1) * determinant * b.determinant));
    }
    { // 6)
      A_new.assign(4*dimension*dimension, 0);
#ifndef NDEBUG
      x0_new.assign(2*dimension, 0);
#endif
      T trafo1[4] = { 0, -1, 1, -1 };
      std::vector<T> A1_trafo(dimension*dimension);
      Blas<T>::gemm('N', 'N', dimension, dimension, dimension,
                    1.0, &A[0], dimension, &trafo1[0], dimension,
                    0.0, &A1_trafo[0], dimension);
      T trafo2[4] = { 0, 1, 1, 0 };
      std::vector<T> A2_trafo(dimension*dimension);
      Blas<T>::gemm('N', 'N', dimension, dimension, dimension,
                    1.0, &b.A[0], dimension, &trafo2[0], dimension,
                    0.0, &A2_trafo[0], dimension);
      for (int i = 0; i < dimension; ++i) {
        Blas<T>::copy(dimension, &A1_trafo[i*dimension], 1, &A_new[i*2*dimension], 1);
        Blas<T>::copy(dimension, &A1_trafo[i*dimension], 1, &A_new[i*2*dimension + 2*dimension*dimension], 1);
        Blas<T>::copy(dimension, &A2_trafo[i*dimension], 1, &A_new[i*2*dimension + 2*dimension*dimension + dimension], 1);
      }
      Blas<T>::copy(dimension, &x0[0], 1, &x0_new[0], 1);
      Blas<T>::axpy(dimension, 1, &A[dimension], 1, &x0_new[0], 1);
      Blas<T>::copy(dimension, &b.x0[0], 1, &x0_new[dimension], 1);
      cross_result.push_back(Simplex<T>(2*dimension, A_new, x0_new, determinant * (-1) * b.determinant));
    }
    return std::vector< Simplex<T> >(cross_result);
  }
  
} // namespace TensorCalculus

#endif // __SIMPLEX_HPP
