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

#ifndef __MESHSIMPLEX_HPP
#define __MESHSIMPLEX_HPP

#include <istream> // std::istream
#include <ostream> // std::ostream
#include <vector>  // std::vector
#include <map>     // std::map
#include <utility> // std::pair, std::make_pair

#include "Mesh/MeshFwd.hpp"
#include "Mesh/MeshPoint.hpp"
#include "BlasInterface.hpp"
#include "LapackInterface2.hpp"

namespace TensorCalculus {

  template<typename T>
  class Mesh<T>::Simplex {
  public:
    Simplex(std::istream& stream, const int dimension, const std::vector<Mesh<T>::Point>& points);
    
    const int getPoint(int i) const;
    const int countPoints() const;
    const std::vector<int>& getPoints() const;
    std::vector<int> getCommonPoints(const Simplex& simplex) const;

    const std::vector<T>& getBaseVector(int i) const;
    std::vector<T> getRank1ApproxBaseVectors(const int k = 1) const;
    
    const T getDeterminant() const;
    const T getDiameter() const;
    
    void print(std::ostream& stream, const std::vector<Point>& points) const;
    
  private:
    std::vector<int> points_indices;
    std::map<int, std::vector<T> > base_vectors;
    std::pair<int, int> longest_edge;
    T determinant;
    T diameter;
    
    void load_from_stream(std::istream& stream, const int dimension, const std::vector<Mesh<T>::Point>& points);
    void calculate_determinant_and_base_vectors(const std::vector<Mesh<T>::Point>& points);
    void find_longest_edge(const std::vector<Mesh<T>::Point>& points);
  };
  
  template<typename T>
  Mesh<T>::Simplex::Simplex(std::istream& stream, const int dimension, const std::vector<Mesh<T>::Point>& points) {
    load_from_stream(stream, dimension, points);
    find_longest_edge(points);
  }
  
  template<typename T>
  const int Mesh<T>::Simplex::getPoint(int i) const {
    return points_indices[i];
  }
  
  template<typename T>
  const int Mesh<T>::Simplex::countPoints() const {
    return points_indices.size();
  }
  
  template<typename T>
  const std::vector<int>& Mesh<T>::Simplex::getPoints() const {
    return points_indices;
  }
  
  template<typename T>
  std::vector<int> Mesh<T>::Simplex::getCommonPoints(const Simplex& simplex) const {
    std::vector<int> result;
    std::set_intersection(points_indices.begin(), points_indices.end(),
                          simplex.points_indices.begin(), simplex.points_indices.end(),
                          std::back_inserter(result));
    return result;
  }
  
  template<typename T>
  const std::vector<T>& Mesh<T>::Simplex::getBaseVector(int i) const {
    typename std::map<int, std::vector<T> >::const_iterator ibase = base_vectors.find(i);
    if (ibase != base_vectors.end()) {
      return ibase->second;
    } else {
      throw std::runtime_error("Point doesn't belong to the vertex set of this simplex");
      // return std::vector<T>(points_indices.size() - 1, 0);
    }
  }
  
  template<typename T>
  std::vector<T> Mesh<T>::Simplex::getRank1ApproxBaseVectors(const int k) const {
    const int d = points_indices.size() - 1;
    std::vector<T> base_vectors_matrix((d + 1)*d);
    for (int i = 0; i <= d; ++i) {
      Blas<T>::copy(d, &getBaseVector(points_indices[i])[0], 1, &base_vectors_matrix[i], d + 1);
    }
    std::vector<T> A((d + 1)*(d + 1));
    Blas<T>::gemm('N', 'T', d + 1, d + 1, d, 
                  1, &base_vectors_matrix[0], d + 1, &base_vectors_matrix[0], d + 1,
                  0, &A[0], d + 1);
    T vu, vl;
    int m;
    std::vector<T> evl(d + 1);
    std::vector<T> evc((d + 1)*(d + 1));
    T work_temp;
    std::vector<int> iwork(5*(d + 1));
    std::vector<int> fail(d + 1);
    int info = Lapack<T>::syevx('V', 'I', 'U', d + 1, &A[0], d + 1, vl, vu, d + 2 - k, d + 1, -1,
                                m, &evl[0], &evc[0], d + 1, &work_temp, -1, &iwork[0], &fail[0]);
    if (info != 0) {
      throw std::runtime_error("Error in workspace size calculation");
    }
    int lwork = static_cast<int>(work_temp);
    std::vector<T> work(lwork);
    info = Lapack<T>::syevx('V', 'I', 'U', d + 1, &A[0], d + 1, vu, vl, d + 2 - k, d + 1, -1,
                            m, &evl[0], &evc[0], d + 1, &work[0], lwork, &iwork[0], &fail[0]);
    if (info != 0) {
      throw std::runtime_error("Error in eigenvalue/-vector calculation");
    }
    using std::sqrt;
    std::vector<T> result((d + 1)*k, 0);
    for (int l = 0; l < k; ++l) {
      Blas<T>::axpy((d + 1), sqrt(evl[l]), &evc[(d+1)*l], 1, &result[(d+1)*l], 1);
    }
    return result;
  }
  
  template<typename  T>
  const T Mesh<T>::Simplex::getDeterminant() const {
    return determinant;
  }
  
  template<typename T>
  const T Mesh<T>::Simplex::getDiameter() const {
    return diameter;
  }
  
  template<typename T>
  void Mesh<T>::Simplex::print(std::ostream& stream, const std::vector<Mesh<T>::Point>& points) const {
    stream << "( ";
    for (unsigned int i = 0; i < points_indices.size(); ++i) {
      points[points_indices[i]].print(stream);
      stream << " ";
    }
    stream << determinant << ' ' << longest_edge.first << ' ' << longest_edge.second << ' ' << diameter << " )";
  }

  template<typename T>
  void Mesh<T>::Simplex::load_from_stream(std::istream& stream, const int dimension, const std::vector<Mesh<T>::Point>& points) {
    points_indices.reserve(dimension + 1);
    for (int i = 0; i <= dimension; ++i) {
      int point_index;
      if (!(stream >> point_index)) throw StreamFormatError();
      points_indices.push_back(point_index);
    }
    std::sort(points_indices.begin(), points_indices.end());
    calculate_determinant_and_base_vectors(points);
  }

  template<typename T>
  void Mesh<T>::Simplex::calculate_determinant_and_base_vectors(const std::vector<Mesh<T>::Point>& points) {
    const int dimension = points_indices.size() - 1;
    if (dimension == 1) {
      const T A = points[points_indices[1]].getCoordinate(0) - points[points_indices[0]].getCoordinate(0);
      // const T x0 = points[points_indices[0]].getCoordinate(0);
      
      determinant = A;
      const T base_vector_value = static_cast<T>(1) / A;
      base_vectors.insert(std::make_pair(points_indices[0], std::vector<T>(1, -base_vector_value)));
      base_vectors.insert(std::make_pair(points_indices[0], std::vector<T>(1, base_vector_value)));
    } else if (dimension == 2) {
      std::vector<T> A(2 * 2);
      A[0] = points[points_indices[1]].getCoordinate(0) - points[points_indices[0]].getCoordinate(0);
      A[1] = points[points_indices[1]].getCoordinate(1) - points[points_indices[0]].getCoordinate(1);
      A[2] = points[points_indices[2]].getCoordinate(0) - points[points_indices[0]].getCoordinate(0);
      A[3] = points[points_indices[2]].getCoordinate(1) - points[points_indices[0]].getCoordinate(1);
      determinant = A[0] * A[3] - A[2] * A[1];
      
      std::vector<T> base_vector(2);
      base_vector[0] = (A[1] - A[3]) / determinant;
      base_vector[1] = (A[2] - A[0]) / determinant;
      base_vectors.insert(std::make_pair(points_indices[0], base_vector)); 
    
      base_vector[0] = A[3] / determinant;
      base_vector[1] = -A[2] / determinant;
      base_vectors.insert(std::make_pair(points_indices[1], base_vector));

      base_vector[0] = -A[1] / determinant;
      base_vector[1] = A[0] / determinant;
      base_vectors.insert(std::make_pair(points_indices[2], base_vector));

    } else {
      std::vector<T> A(dimension * dimension);
      const T* px0 = points[points_indices[0]].getData();
      for (int j = 0; j < dimension; ++j) {
        const T* pxj = points[points_indices[j + 1]].getData();
        Blas<T>::copy(dimension, pxj, 1, &A[j*dimension], 1);
        Blas<T>::axpy(dimension, -1, px0, 1, &A[j*dimension], 1);
      }
      {
        std::vector<int> pivots(dimension);
        if (Lapack<T>::getrf(dimension, dimension, &A[0], dimension, &pivots[0]) != 0) {
          throw std::runtime_error("Error in \"?getrf\"");
        }

        determinant = 1;
        for (int j = 0, offset = 0; j < dimension; ++j, offset += dimension + 1) {
          determinant *= A[offset];
        }

        {
          double dbl_lwork;
          if (Lapack<T>::getri(dimension, &A[0], dimension, &pivots[0], &dbl_lwork, -1) != 0) {
            throw std::runtime_error("Error in workspace query in \"?getri\"");
          }
          int lwork = static_cast<int>(dbl_lwork);
          std::vector<double> work(lwork);
          if (Lapack<T>::getri(dimension, &A[0], dimension, &pivots[0], &work[0], lwork) != 0) {
            throw std::runtime_error("Error in \"?getri\"");
          }
          
          std::vector<T> first_base_vector(dimension);
          std::vector<T> base_vector(dimension);
          for (int j = 0; j < dimension; ++j) {
            if (!points[points_indices[j + 1]].isBorder()) {
              Blas<T>::copy(dimension, &A[j], dimension, &base_vector[0], 1);
              base_vectors.insert(make_pair(j + 1, base_vector));
            }
            Blas<T>::axpy(dimension, -1, &A[j * dimension], 1, &base_vector[0], 1);
          }
          if (!points[points_indices[0]].isBorder()) {
            base_vectors.insert(make_pair(0, first_base_vector));
          }
        }
      }
    }
  }

  template<typename T>
  void Mesh<T>::Simplex::find_longest_edge(const std::vector<Mesh<T>::Point>& points)
  {
    const int dimension = points_indices.size() - 1;
    T length = 0;
    std::vector<T> difference(dimension);
    for (int i = 0; i < points_indices.size(); ++i) {
      const std::vector<T>& first = points[points_indices[i]].getCoordinates();
      for (int j = i+1; j < points_indices.size(); ++j) {
        const std::vector<T>& second = points[points_indices[j]].getCoordinates();
        Blas<T>::copy(dimension, &first[0], 1,
                                 &difference[0], 1);
        Blas<T>::axpy(dimension, -1, &second[0], 1, &difference[0], 1);
        const T new_length = Blas<T>::nrm2(dimension, &difference[0], 1);
        if (new_length > length) {
          longest_edge.first = i;
          longest_edge.second = j;
          length = new_length;
        }
      }
    }
    diameter = length;
  }

} // namespace TensorCalculus


#endif // __MESHSIMPLEX_HPP
