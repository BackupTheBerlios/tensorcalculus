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

#include <iostream>
#include <vector>
#include "Mesh/Simplex.hpp"
#include "Vector/VectorOperators.hpp"

using namespace TensorCalculus;

int main() {
  std::vector< std::vector<double> > points;
  points.reserve(3);

  std::vector<double> point(2);
  point[0] = 0.0;
  point[1] = 0.0;
  // point[2] = 0.0;
  points.push_back(point);

  point[0] = 1.0;
  point[1] = 0.0;
  // point[2] = 0.0;
  points.push_back(point);

  point[0] = 0.0;
  point[1] = 1.0;
  // point[2] = 0.0;
  points.push_back(point);

  // point[0] = 1.0;
  // point[1] = 1.0;
  // point[2] = 1.0;
  // points.push_back(point);

  Simplex<double> s(points);

  std::cout << "Simplex: " << s << '\n';
  // std::cout << "Area: " << s.get_area() << '\n';
  std::cout << "Determinant: " << s.get_determinant() << std::endl;

  std::vector< Simplex<double> > s2 = cross_product(s, s);
  for (unsigned int i = 0; i < s2.size(); ++i) {
    std::cout << s2[i] << '\n';
  }
  
  // std::vector< Simplex<double> > d = s.dissect();
  // for (int i = 0; i < d.size(); ++i) {
  //   std::cout << d[i] << '\n';
  // }

  std::vector<double> h(2);
  h[0] = 1.0/3.0;
  h[1] = 1.0/3.0;
  // h[2] = 1 - h[0] - h[1];
  std::vector<double> c = s.transform(h);
  {
    using VectorOperators::operator <<;
    std::cout << c << std::endl;
  }
  return 0;
}
