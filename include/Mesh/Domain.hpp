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

#ifndef __DOMAIN_HPP
#define __DOMAIN_HPP

#include <vector>

#include "Mesh/Mesh.hpp"
#include "Mesh/MeshFunction.hpp"
#include "Mesh/ValuesSimplex.hpp"
#include "Utilities/Factorials.hpp"

namespace TensorCalculus {

  template<typename T>
  class Domain {
  public:
    class Simplex {
    public:
      Simplex(const std::vector< std::vector<T> >& values, T factor, int original_simplex);
              //const std::vector<const typename Mesh<T>::Simplex*>& original_simplices);
      
      std::vector<Simplex> cross_product(const Simplex& s) const;
      
      //const std::vector<const typename Mesh<T>::Simplex*>& getOriginalSimplices() const { return original_simplices; }
      const std::vector<int>& getOriginalSimplices() const { return original_simplices; }
      const ValuesSimplex<T>& getValuesSimplex() const { return values_simplex; }
    private:
      Simplex(const ValuesSimplex<T>& values_simplex,
              const std::vector<int>& original_simplices)
              // const std::vector<const typename Mesh<T>::Simplex*>& original_simplices)
        : values_simplex(values_simplex), original_simplices(original_simplices) { }
    
      ValuesSimplex<T> values_simplex;
      std::vector<int> original_simplices;
      // std::vector<const typename Mesh<T>::Simplex*> original_simplices;
    };
  
    Domain(const Mesh<T>& mesh);
    Domain(const Mesh<T>& mesh, const MeshFunction<T>& mesh_function);
    
    Domain cross_product(const Domain& d) const;
    
    const unsigned int countSimplices() const { return simplices.size(); }
    const Domain::Simplex& getSimplex(int index) const { return simplices[index]; }
    const int getDimension() const { return dimension; }
    
    const T getArea() const;
  private:
    Domain(int dimension, const std::vector<Domain::Simplex>& simplices)
      : dimension(dimension), simplices(simplices) { }
  
    int dimension;
    std::vector<Domain::Simplex> simplices;
  };
  
  template<typename T>
  Domain<T>::Simplex::Simplex(const std::vector< std::vector<T> >& values, T factor, int original_simplex)
                              // const std::vector<const typename Mesh<T>::Simplex*>& original_simplices)
    : values_simplex(values, factor), original_simplices(std::vector<int>(1, original_simplex)) { }
  
  template<typename T>
  std::vector<typename Domain<T>::Simplex> Domain<T>::Simplex::cross_product(const Domain::Simplex& s) const
  {
    // std::vector<const typename Mesh<T>::Simplex*> common_original_simplices;
    std::vector<int> common_original_simplices;
    common_original_simplices.reserve(original_simplices.size() + s.original_simplices.size());
    common_original_simplices.assign(original_simplices.begin(), original_simplices.end());
    common_original_simplices.insert(common_original_simplices.end(),
                                     s.original_simplices.begin(), s.original_simplices.end());
                                     
    std::vector< ValuesSimplex<T> > values_simplices = values_simplex.cross_product(s.values_simplex);
    std::vector<Domain::Simplex> result;
    result.reserve(values_simplices.size());
    for (int i = 0; i < values_simplices.size(); ++i) {
      result.push_back(Domain::Simplex(values_simplices[i], common_original_simplices));
    }
    return result;
  }
  
  template<typename T>
  Domain<T> Domain<T>::cross_product(const Domain<T>& d) const
  {
    std::vector<Domain::Simplex> new_simplices;
    new_simplices.reserve(simplices.size() * d.simplices.size() * 6);
    for (int i = 0; i < simplices.size(); ++i) {
      for (int j = 0; j < d.simplices.size(); ++j) {
        const std::vector<Domain::Simplex> parts = simplices[i].cross_product(d.simplices[j]);
        new_simplices.insert(new_simplices.end(),
                             parts.begin(), parts.end());
      }
    }
    return Domain<T>(dimension + d.dimension, new_simplices);
  }
  
  template<typename T>
  Domain<T> cross_product(const Domain<T>& u, const Domain<T>& v)
  {
    return u.cross_product(v);
  }
  
  template<typename T>
  Domain<T>::Domain(const Mesh<T>& mesh)
    : dimension(mesh.getDimension())
  {
    simplices.reserve(mesh.countSimplices());
    for (int i = 0; i < mesh.countSimplices(); ++i) {
      const typename Mesh<T>::Simplex& simplex = mesh.getSimplex(i);
      std::vector< std::vector<T> > values;
      values.reserve(mesh.getDimension() + 1);
      for (int j = 0; j <= mesh.getDimension(); ++j) {
        values.push_back(mesh.getPoint(j).getCoordinates());
      }
      simplices.push_back(Domain::Simplex(values, simplex.getDeterminant(), i));
                          //std::vector<const typename Mesh<T>::Simplex*>(1, &simplex)));
    }
  }

  template<typename T>
  Domain<T>::Domain(const Mesh<T>& mesh, const MeshFunction<T>& mesh_function)
    : dimension(mesh.getDimension())
  {
    simplices.reserve(mesh.countSimplices());
    for (int i = 0; i < mesh.countSimplices(); ++i) {
      const typename Mesh<T>::Simplex& simplex = mesh.getSimplex(i);
      std::vector< std::vector<T> > values;
      values.reserve(mesh.getDimension() + 1);
      for (int j = 0; j <= mesh.getDimension(); ++j) {
        values.push_back(mesh_function.getValues(simplex.getPoint(j)));
      }
      simplices.push_back(Domain::Simplex(values, simplex.getDeterminant(), i));
                          // std::vector<const typename Mesh<T>::Simplex*>(1, &simplex)));
    }
  }
  
  template<typename T>
  const T Domain<T>::getArea() const
  {
    using std::abs;
    const T fac = factorial<T>(getDimension());
    T value = 0;
    for (int i = 0; i < simplices.size(); ++i) {
      value += abs(simplices[i].getValuesSimplex().getFactor());
    }
    value /= fac;
    return value;
  }

} // namespace TensorCalculus

#endif // __DOMAIN_HPP
