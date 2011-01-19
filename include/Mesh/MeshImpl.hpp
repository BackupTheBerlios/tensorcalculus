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

#ifndef __MESHIMPL_HPP
#define __MESHIMPL_HPP

#include "Mesh/MeshFwd.hpp"     // Mesh
#include "Mesh/MeshSimplex.hpp" // Mesh::Simplex
#include "Mesh/MeshPoint.hpp"   // Mesh::Point
#include "Utilities/Factorials.hpp"  // factorial
#include "Utilities/Utilities.hpp"   // StreamFormatError

namespace TensorCalculus {
  
  template<typename T>
  Mesh<T>::Mesh(std::istream& stream) {
    load_from_stream(stream);
    fill_mass_matrix();
  }

  template<typename T>
  const int Mesh<T>::getDimension() const {
    return dimension;
  }
  
  template<typename T>
  const int Mesh<T>::countPoints() const {
    return points.size();
  }
  
  template<typename T>
  const typename Mesh<T>::Point& Mesh<T>::getPoint(int index) const {
    return points[index];
  }
  
  template<typename T>
  const int Mesh<T>::countSimplices() const {
    return simplices.size();
  }
  
  template<typename T>
  const typename Mesh<T>::Simplex& Mesh<T>::getSimplex(int index) const {
    return simplices[index];
  }
  
  template<typename T>
  const T Mesh<T>::get_mass_matrix(unsigned int i, unsigned int j) const {
    using std::swap;
    if (i > j) swap(i, j);
      
    const std::map<int, T>& row = mass_matrix[i];
    const typename std::map<int, T>::const_iterator it = row.find(j);
    if (it != row.end()) {
      return it->second;
    } else {
      return 0;
    }
  }
  
  template<typename T>
  void Mesh<T>::load_from_stream(std::istream& stream) {
    std::string line;
    int num_points;
    if (!getline(stream, line)) throw StreamFormatError();
    {
      std::stringstream sstream(line);
      if (!(sstream >> num_points)) throw StreamFormatError();
      if (!(sstream >> dimension)) throw StreamFormatError();
      // points.reserve(num_points);
      // for (int i = 0; i < num_points; ++i) {
      //   points.push_back(std::vector<T>(dimension));
      // }
    }
    {
      std::map<int, Point> points_map;
      for (int i = 0; i < num_points; ++i) {
        if (!getline(stream, line)) throw StreamFormatError();
        std::stringstream sstream(line);
        
        int index;
        if (!(sstream >> index)) throw StreamFormatError();
        
        Point point(sstream, dimension);
        std::pair<typename std::map<int, Point>::iterator, bool> answer
          = points_map.insert(std::make_pair(index, point));
        if (!answer.second) throw StreamFormatError("Point occures doubly in the list");
      }
      points.reserve(num_points);
      for (int i = 0; i < num_points; ++i) {
        const typename std::map<int, Point>::iterator ipoint = points_map.find(i);
        if (ipoint == points_map.end()) throw StreamFormatError("Point not in list");
        points.push_back(ipoint->second);
      }
    }

    int num_simplices;
    {
      if (!getline(stream, line)) throw StreamFormatError();
      std::stringstream sstream(line);
      if (!(sstream >> num_simplices)) throw StreamFormatError();
      // simplices.reserve(num_simplices);
      // for (int i = 0; i < num_simplices; ++i) {
      //   simplices.push_back(std::vector<int>(dimension + 1));
      // }
    }
    {
      std::map<int, Simplex> simplices_map;
      for (int i = 0; i < num_simplices; ++i) {
        if (!getline(stream, line)) throw StreamFormatError();
        std::stringstream sstream(line);
  
        int index;
        if (!(sstream >> index)) throw StreamFormatError();
        
        Simplex simplex(sstream, dimension, points);
        std::pair<typename std::map<int, Simplex>::iterator, bool> answer
          = simplices_map.insert(std::make_pair(index, simplex));
        if (!answer.second) throw StreamFormatError("Simplex occures doubly in the list");
      }
      simplices.reserve(num_simplices);
      for (int i = 0; i < num_simplices; ++i) {
        const typename std::map<int, Simplex>::iterator isimplex = simplices_map.find(i);
        if (isimplex == simplices_map.end()) throw StreamFormatError("Simplex not in list");
        simplices.push_back(isimplex->second);
      }
    }
  }
  
  template<typename T>
  void Mesh<T>::fill_mass_matrix()
  {
    const double fac = factorial<T>(dimension + 2);
    
    mass_matrix.resize(points.size());

    for (unsigned int simplex_index = 0; simplex_index < simplices.size(); ++simplex_index) {
      const Simplex& simplex = simplices[simplex_index];
      using std::abs;
      const T abs_det = abs(simplex.getDeterminant());
      
      for (int point_index = 0; point_index <= dimension; ++point_index) {
        const unsigned int point = simplex.getPoint(point_index);
        mass_matrix[point][point] += 2.0*abs_det / fac;
        for (int point_index2 = point_index + 1; point_index2 <= dimension; ++point_index2) {
          const int point2 = simplex.getPoint(point_index2);

          int i = point;
          int j = point2;
          using std::swap;
          if (i > j) swap(i, j);
          
          mass_matrix[i][j] += abs_det / fac;
        }
      }
    }
  }

  template<typename T>
  void Mesh<T>::print(std::ostream& stream) const {
    for (unsigned int i = 0; i < simplices.size() - 1; ++i) {
      simplices[i].print(stream, points);
      stream << '\n';
    }
    simplices.back().print(stream, points);
  }

  template<typename T>
  std::ostream& operator << (std::ostream& stream, const Mesh<T>& mesh)
  {
    mesh.print(stream);
    return stream;
  }  
  
} // namespace TensorCalculus


#endif // __MESHIMPL_HPP
