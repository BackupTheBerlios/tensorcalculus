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

#ifndef __MESHPOINT_HPP
#define __MESHPOINT_HPP

#include <istream>
#include <vector>

#include "Mesh/MeshFwd.hpp"
#include "Utilities/Utilities.hpp"

namespace TensorCalculus {

  template<typename T>
  class Mesh<T>::Point {
  public:
    Point(std::istream& stream, const int dimension);
    
    const int getDimension() const;
    const T getCoordinate(int i) const;
    const std::vector<T>& getCoordinates() const;
    
    T* getData();
    const T* getData() const;
    
    const bool isBorder() const;
    
    void print(std::ostream& stream) const;
  private:
    std::vector<T> coords;
    bool border;
    
    void load_from_stream(std::istream& stream, const int dimension);
  };

  template<typename T>
  inline Mesh<T>::Point::Point(std::istream& stream, const int dimension) {
    load_from_stream(stream, dimension);
  }

  template<typename T>
  inline const T Mesh<T>::Point::getCoordinate(int i) const {
    return coords[i];
  }

  template<typename T>
  inline const std::vector<T>& Mesh<T>::Point::getCoordinates() const {
    return coords;
  }
  
  template<typename T>
  inline T* Mesh<T>::Point::getData() {
    return &coords[0];
  }
  
  template<typename T>
  inline const T* Mesh<T>::Point::getData() const {
    return &coords[0];
  }
  
  template<typename T>
  inline const bool Mesh<T>::Point::isBorder() const {
    return border;
  }
  
  template<typename T>
  inline const int Mesh<T>::Point::getDimension() const {
    return coords.size();
  }
  
  template<typename T>
  void Mesh<T>::Point::load_from_stream(std::istream& stream, const int dimension) {
    coords.reserve(dimension);
    for (int i = 0; i < dimension; ++i) {
      T coord;
      if (!(stream >> coord)) throw StreamFormatError();
      coords.push_back(coord);
    }
    if (!(stream >> border)) throw StreamFormatError();
  }
  
  template<typename T>
  void Mesh<T>::Point::print(std::ostream& stream) const
  {
    stream << "( ";
    for (unsigned int i = 0; i < coords.size(); ++i) {
      stream << coords[i] << " ";
    }
    stream << ")";
  }

} // namespace TensorCalculus

#endif // __MESHPOINT_HPP
