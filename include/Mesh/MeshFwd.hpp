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

#ifndef __MESHFWD_HPP
#define __MESHFWD_HPP

#include <vector>
#include <istream>
#include <ostream>
#include <map>

namespace TensorCalculus {

  template<typename T>
  class Mesh {
  public:
    class Point;
    class Simplex;
  
    explicit Mesh(std::istream& stream);
    
    const int getDimension() const;
    
    const int countPoints() const;
    const Point& getPoint(int index) const;
    
    const int countSimplices() const;
    const Simplex& getSimplex(int index) const;
    
    const T get_mass_matrix(unsigned int i, unsigned int j) const;
    
    void print(std::ostream& stream) const;
  private:
    int dimension;
    std::vector<Point> points;
    std::vector<Simplex> simplices;
    std::vector< std::map<int, T> > mass_matrix;

    void load_from_stream(std::istream& stream);
    void fill_mass_matrix();
    void calculate_determinant_and_base_vectors();
  };
  
  template<typename T>
  std::ostream& operator << (std::ostream& stream, const Mesh<T>& mesh);
  
} // namespace TensorCalculus

#endif // __MESHFWD_HPP
