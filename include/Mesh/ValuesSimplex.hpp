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

#ifndef __VALUESSIMPLEX_HPP
#define __VALUESSIMPLEX_HPP

#include <vector>
#include <stdexcept>
#include <ostream>
#include <limits>

#include "BlasInterface.hpp"
#include "Mesh/Simplex.hpp"
// #include "VectorOperators.hpp"

namespace TensorCalculus {

  // template<typename T>
  // class Mesh;
  
  template<typename T>
  class ValuesSimplex {
  public:
    ValuesSimplex(const std::vector< std::vector<T> >& values,
                  // const typename Mesh<T>::Simplex& original_simplex,
                  T factor)
      : domain_dim(values.size() - 1),
        value_dim(values.at(0).size()),
        // original_simplices(1, &original_simplex),
        factor(factor)
    {
      fill_matrix_and_origin(values);
    }
    
    // explicit ValuesSimplex(const Simplex<T>& simplex)
    //   : num_vertices(simplex.get_dimension() + 1),
    //     dimension(simplex.get_dimension()),
    //     matrix(simplex.get_matrix()),
    //     origin(simplex.get_origin()),
    //     factor(simplex.get_determinant())
    // { }

    std::vector<T> transform(const std::vector<T>& standard) const
    {
      std::vector<T> transformed(value_dim);
      transform(standard, transformed);
      return transformed;
    }

    const int getDomainDim() const { return domain_dim; }
    const int getValueDim() const { return value_dim; }
    
    void setFector(T new_factor) { factor = new_factor; } 
    const T getFactor() const { return factor; }
    // const std::vector<const typename Mesh<T>::Simplex*>& getOriginalSimplices() const { return original_simplices; }
    
    void transform(const std::vector<T>& standard, std::vector<T>& transformed) const;
    std::vector<ValuesSimplex> decompose() const;
    std::vector<ValuesSimplex> cross_product(const ValuesSimplex& s) const;

    void print(std::ostream& stream) const;

  private:
    ValuesSimplex(int domain_dim, int value_dim,
                  const std::vector<T>& matrix, const std::vector<T>& origin,
                  // const std::vector<const typename Mesh<T>::Simplex*>& original_simplices,
                  T factor)
      : domain_dim(domain_dim),
        value_dim(value_dim),
        matrix(matrix),
        origin(origin),
        // original_simplices(original_simplices),
        factor(factor)
    { }
  
    int domain_dim; // number of vertices in the simplex
    int value_dim; // dimension of the values attached to the vertices
    std::vector<T> matrix;
    std::vector<T> origin;
    // std::vector<const typename Mesh<T>::Simplex*> original_simplices;
    T factor;
    
    void fill_matrix_and_origin(const std::vector< std::vector<T> >& values);
    std::pair<int, int> find_longest_edge() const;
  };
  
  template<typename T>
  void ValuesSimplex<T>::fill_matrix_and_origin(const std::vector< std::vector<T> >& values) {
#if defined(RANGE_CHECKS_ON) || defined(_DEBUG) || !defined(NDEBUG)
    for (unsigned int i = 1; i < values.size(); ++i) {
      if (values[i].size() != static_cast<unsigned int>(value_dim)) {
        throw std::invalid_argument("Values must have same size");
      }
    }
#endif

    matrix.resize(value_dim * domain_dim);
    origin = values[0];
    for (int i = 0; i < domain_dim; ++i) {
      Blas<T>::copy(value_dim, &(values[i+1][0]), 1, &matrix[i*value_dim], 1);
      Blas<T>::axpy(value_dim, -1, &(origin[0]), 1, &matrix[i*value_dim], 1);
    }
  }

  template<typename T>
  void ValuesSimplex<T>::transform(const std::vector<T>& standard, std::vector<T>& transformed) const {
#ifdef RANGE_CHECKS_ON
    if (standard.size() != domain_dim) {
      throw std::invalid_argument("Coordinate and transformation matrix size doesn't match");
    }
    if (transformed.size() != value_dim) {
      throw std::invalid_argument("Result vector and transformation matrix size doesn't match");
    }
#endif // RANGE_CHECKS_ON
    Blas<T>::gemv('N', value_dim, domain_dim,
                  1.0, &matrix[0], value_dim,
                  &standard[0], 1,
                  0.0, &transformed[0], 1);
    Blas<T>::axpy(value_dim, 1.0, &origin[0], 1, &transformed[0], 1);
    // return result;
  }
  
  template<typename T>
  std::vector< ValuesSimplex<T> > ValuesSimplex<T>::decompose() const {
    if (domain_dim == 1) {
      std::vector< ValuesSimplex<T> > result;
      result.reserve(2);
      
      std::vector<T> new_matrix = matrix;
      Blas<T>::scal(value_dim, 0.5, &new_matrix[0], 1);
      std::vector<T> new_origin = origin;
      const T new_factor = 0.5;
      
      result.push_back(ValuesSimplex(domain_dim, value_dim, new_matrix, new_origin, new_factor));
      
      Blas<T>::axpy(value_dim, 1, &new_matrix[0], 1, &new_origin[0], 1);
      result.push_back(ValuesSimplex(domain_dim, value_dim, new_matrix, new_origin, new_factor));
      
      return result;
    } else if (domain_dim == 2) {
      std::vector< ValuesSimplex<T> > result;
      result.reserve(4);
      
      std::vector<T> new_matrix = matrix;
      scale<T>(new_matrix, 0.5);
      std::vector<T> new_origin = origin;
      const T new_factor = 0.25 * factor;
      
      result.push_back(ValuesSimplex(domain_dim, value_dim, new_matrix, new_origin, /*original_simplices,*/ new_factor));

      Blas<T>::axpy(value_dim, 1, &new_matrix[0], 1, &new_origin[0], 1);
      result.push_back(ValuesSimplex(domain_dim, value_dim, new_matrix, new_origin, /*original_simplices,*/ new_factor));
      
      new_origin = origin;
      Blas<T>::axpy(value_dim, 1, &new_matrix[value_dim], 1, &new_origin[0], 1);
      result.push_back(ValuesSimplex(domain_dim, value_dim, new_matrix, new_origin, /*original_simplices,*/ new_factor));
      
      Blas<T>::axpy(value_dim, 1, &new_matrix[0], 1, &new_origin[0], 1);
      scale<T>(new_matrix, -1);
      result.push_back(ValuesSimplex(domain_dim, value_dim, new_matrix, new_origin, /*original_simplices,*/ -new_factor));

      return result;
    } else {
      std::pair<int, int> longest_edge = find_longest_edge(); // note: first < second
      const T new_factor = 0.5 * factor;
      
      std::vector< ValuesSimplex<T> > result;
      result.reserve(2);

      if (longest_edge.first == 0) {
        std::vector<T> new_matrix = matrix;
        std::vector<T> new_origin = origin;
        
        Blas<T>::scal(value_dim, 0.5, &new_matrix[(longest_edge.second - 1) * value_dim], 1);
        result.push_back(ValuesSimplex(domain_dim, value_dim, new_matrix, new_origin, new_factor));
        
        for (int i = 0; i < longest_edge.second - 1; ++i) {
          Blas<T>::axpy(value_dim, -1, &new_matrix[(longest_edge.second - 1) * value_dim], 1,
                                       &new_matrix[i * value_dim], 1);
        }
        for (int i = longest_edge.second; i < domain_dim; ++i) {
          Blas<T>::axpy(value_dim, -1, &new_matrix[(longest_edge.second - 1) * value_dim], 1,
                                       &new_matrix[i * value_dim], 1);
        }
        Blas<T>::axpy(value_dim, 1, &new_matrix[(longest_edge.second - 1) * value_dim], 1,
                                    &new_origin[0], 1);
        result.push_back(ValuesSimplex(domain_dim, value_dim, new_matrix, new_origin, new_factor));
        
        return result;
      } else {
        std::vector<T> new_matrix = matrix;
        const T new_factor = 0.5 * factor;

        Blas<T>::scal(value_dim, 0.5, &new_matrix[(longest_edge.first - 1) * value_dim], 1);
        Blas<T>::axpy(value_dim, 0.5, &new_matrix[(longest_edge.second - 1) * value_dim], 1,
                                      &new_matrix[(longest_edge.first - 1) * value_dim], 1);
        result.push_back(ValuesSimplex(domain_dim, value_dim, new_matrix, origin, new_factor));
        
        Blas<T>::copy(value_dim, &new_matrix[(longest_edge.first - 1) * value_dim], 1,
                                 &new_matrix[(longest_edge.second - 1) * value_dim], 1);
        Blas<T>::copy(value_dim, &matrix[(longest_edge.first - 1) * value_dim], 1,
                                 &new_matrix[(longest_edge.first - 1) * value_dim], 1);
        result.push_back(ValuesSimplex(domain_dim, value_dim, new_matrix, origin, new_factor));
        
        return result;
      }           
    }
  }

  template<typename T>
  struct transformation {
    T factors[2][2];
    T matrix[2][4];
    T origin[2][2];
  };

  template<typename T>
  std::vector< ValuesSimplex<T> > ValuesSimplex<T>::cross_product(const ValuesSimplex<T>& s) const
  {
    if (s.domain_dim == 2) {
      if (domain_dim == 2) {
        transformation<T> trafos[6] =
          { { { { 0, 1 }, { 1, 1 } }, { {  1, -1,  0, -1 }, {  1,  0,  0,  1 } }, { { 0, 1 }, { 0, 0 } } },
            { { { 1, 1 }, { 0, 1 } }, { {  1,  0,  0,  1 }, {  1, -1,  0, -1 } }, { { 0, 0 }, { 0, 1 } } },
            { { { 1, 1 }, { 0, 1 } }, { { -1,  0, -1,  1 }, { -1,  1, -1,  0 } }, { { 1, 0 }, { 1, 0 } } },
            { { { 0, 1 }, { 1, 1 } }, { { -1,  1, -1,  0 }, { -1,  0, -1,  1 } }, { { 1, 0 }, { 1, 0 } } },
            { { { 0, 1 }, { 1, 1 } }, { {  0,  1,  1,  0 }, {  0, -1,  1, -1 } }, { { 0, 0 }, { 0, 1 } } },
            { { { 1, 1 }, { 0, 1 } }, { {  0, -1,  1, -1 }, {  0,  1,  1,  0 } }, { { 0, 1 }, { 0, 0 } } } };
        
        std::vector< ValuesSimplex<T> > result;
        result.reserve(6);
        std::vector<T> new_matrix((value_dim + s.value_dim) * (domain_dim + s.domain_dim));
        std::vector<T> new_origin( value_dim + s.value_dim );
        
        for (int i = 0; i < 6; ++i) {
          Blas<T>::gemm('N', 'N', value_dim, domain_dim, domain_dim,
                        trafos[i].factors[0][0], &matrix[0], value_dim, &trafos[i].matrix[0][0], domain_dim,
                        0, &new_matrix[0], value_dim + s.value_dim);
          Blas<T>::gemm('N', 'N', value_dim, domain_dim, domain_dim,
                        trafos[i].factors[0][1], &matrix[0], value_dim, &trafos[i].matrix[0][0], domain_dim,
                        0, &new_matrix[(value_dim + s.value_dim)*domain_dim], value_dim + s.value_dim);
    
          Blas<T>::gemm('N', 'N', s.value_dim, domain_dim, domain_dim,
                        trafos[i].factors[1][0], &s.matrix[0], s.value_dim, &trafos[i].matrix[1][0], domain_dim,
                        0, &new_matrix[value_dim], value_dim + s.value_dim);
          Blas<T>::gemm('N', 'N', s.value_dim, domain_dim, domain_dim,
                        trafos[i].factors[1][1], &s.matrix[0], s.value_dim, &trafos[i].matrix[1][0], domain_dim,
                        0, &new_matrix[(value_dim + s.value_dim)*domain_dim + value_dim], value_dim + s.value_dim);
    
          Blas<T>::copy(value_dim, &origin[0], 1, &new_origin[0], 1);
          Blas<T>::copy(s.value_dim, &s.origin[0], 1, &new_origin[value_dim], 1);
          
          Blas<T>::gemv('N', value_dim, domain_dim,
                        1, &matrix[0], value_dim, &trafos[i].origin[0][0], 1,
                        1, &new_origin[0], 1);
          Blas<T>::gemv('N', s.value_dim, domain_dim,
                        1, &s.matrix[0], s.value_dim, &trafos[i].origin[1][0], 1,
                        1, &new_origin[value_dim], 1);
                        
          // std::vector<const typename Mesh<T>::Simplex*> new_original_simplices;
          // new_original_simplices.reserve(original_simplices.size() + s.original_simplices.size());
          // new_original_simplices.assign(original_simplices.begin(), original_simplices.end());
          // new_original_simplices.insert(new_original_simplices.end(), s.original_simplices.begin(), s.original_simplices.end());
          result.push_back(ValuesSimplex<T>(domain_dim + s.domain_dim, value_dim + s.value_dim,
                                            new_matrix, new_origin, /*new_original_simplices,*/
                                            - factor * s.factor)); 
        }
        return result;
      } else if (domain_dim == 1) {
        return s.cross_product(*this);
      }
    } else if (s.domain_dim == 1) {
      if (domain_dim == 2) {
        std::vector< ValuesSimplex<T> > result;
        result.reserve(3);

        std::vector<T> new_origin(value_dim + s.value_dim);
        Blas<T>::copy(value_dim, &origin[0], 1, &new_origin[0], 1);
        Blas<T>::copy(s.value_dim, &s.origin[0], 1, &new_origin[value_dim], 1);
        
        std::vector<T> new_matrix((value_dim + s.value_dim) * 3);
        Blas<T>::copy(value_dim, &matrix[0], 1, &new_matrix[0], 1);
        Blas<T>::copy(value_dim, &matrix[value_dim], 1, &new_matrix[value_dim + s.value_dim], 1);
        Blas<T>::copy(s.value_dim, &s.matrix[0], 1, &new_matrix[2 * (value_dim + s.value_dim) + value_dim], 1);
        
        result.push_back(ValuesSimplex(domain_dim + s.domain_dim, value_dim + s.value_dim,
                                       new_matrix, new_origin, factor * s.factor));
        
        Blas<T>::axpy(s.value_dim, 1, &s.matrix[0], 1, &new_origin[value_dim], 1);
        
        Blas<T>::copy(value_dim, &matrix[value_dim], 1, &new_matrix[2 * (value_dim + s.value_dim)], 1);
        Blas<T>::scal(s.value_dim, -1, &new_matrix[2 * (value_dim + s.value_dim) + value_dim], 1);
        result.push_back(ValuesSimplex(domain_dim + s.domain_dim, value_dim + s.value_dim,
                                       new_matrix, new_origin, - factor * s.factor));
        
        Blas<T>::copy(value_dim, &matrix[0], 1, &new_matrix[2 * (value_dim + s.value_dim)], 1);
        Blas<T>::axpy(s.value_dim, -1, &s.matrix[0], 1, &new_matrix[(value_dim + s.value_dim) + value_dim], 1);
        result.push_back(ValuesSimplex(domain_dim + s.domain_dim, value_dim + s.value_dim,
                                       new_matrix, new_origin, - factor * s.factor));
        
        return result;
      } else if (domain_dim == 1) {
        std::vector<T> new_origin(value_dim + s.value_dim);
        Blas<T>::copy(value_dim, &origin[0], 1, &new_origin[0], 1);
        Blas<T>::copy(s.value_dim, &s.origin[0], 1, &new_origin[value_dim], 1);
        
        std::vector<T> new_matrix((value_dim + s.value_dim) * 2);
        Blas<T>::copy(value_dim, &matrix[0], 1, &new_matrix[0], 1);
        Blas<T>::copy(value_dim, &matrix[0], 1, &new_matrix[value_dim + s.value_dim], 1);
        Blas<T>::copy(s.value_dim, &s.matrix[0], 1, &new_matrix[value_dim + s.value_dim + value_dim], 1);
        
        std::vector< ValuesSimplex<T> > result;
        result.reserve(2);
        
        result.push_back(ValuesSimplex(domain_dim + s.domain_dim, value_dim + s.value_dim,
                                       new_matrix, new_origin, factor * s.factor));
                                       
        Blas<T>::scal(value_dim, 0.0, &new_matrix[0], 1);
        Blas<T>::copy(s.value_dim, &s.matrix[0], 1, &new_matrix[value_dim], 1);
        result.push_back(ValuesSimplex(domain_dim + s.domain_dim, value_dim + s.value_dim,
                                       new_matrix, new_origin, - factor * s.factor));
        
        return result;
      }
    }
    throw std::runtime_error("Cross product only implemented for dimensions up to two");
  }
  
  template<typename T>
  std::vector< ValuesSimplex<T> > cross_product(const ValuesSimplex<T>& a, const ValuesSimplex<T>& b)
  {
    return a.cross_product(b);
  }
  
  template<typename T>
  void ValuesSimplex<T>::print(std::ostream& stream) const
  {
    stream << "( ";
    for (int j = 0; j < value_dim; ++j) {
      stream << origin[j] << " ";
    }
    stream << ") ";
    std::vector<T> values(value_dim);
    for (int i = 0; i < domain_dim; ++i) {
      Blas<T>::copy(value_dim, &matrix[i * value_dim], 1, &values[0], 1);
      Blas<T>::axpy(value_dim, 1, &origin[0], 1, &values[0], 1);
      stream << "( ";
      for (int j = 0; j < value_dim; ++j) {
        stream << values[j] << " ";
      }
      stream << ")";
    }
  }
  
  template<typename T>
  std::pair<int, int> ValuesSimplex<T>::find_longest_edge() const
  {
    int i0 = 0;
    int i1 = 1;
    T max = Blas<T>::nrm2(value_dim, &matrix[0], 1);
    for (int i = 1; i < domain_dim; ++i) {
      const T new_max = Blas<T>::nrm2(value_dim, &matrix[i * value_dim], 1);
      if (new_max > max) {
        max = new_max;
        i1 = i + 1;
      }
    }
    std::vector<T> difference(value_dim);
    for (int i = 0; i < domain_dim; ++i) {
      for (int j = i + 1; j < domain_dim; ++j) {
        Blas<T>::copy(value_dim, &matrix[i * value_dim], 1, &difference[0], 1);
        Blas<T>::axpy(value_dim, -1, &matrix[j * value_dim], 1, &difference[0], 1);
        const T new_max = Blas<T>::nrm2(value_dim, &difference[0], 1);
        if (new_max > max) {
          max = new_max;
          i0 = i + 1;
          i1 = j + 1;
        }
      }
    }
    return std::make_pair(i0, i1); // fact: i0 < i1
  }

} // namespace TensorCalculus

#endif // __VALUESSIMPLEX_HPP
