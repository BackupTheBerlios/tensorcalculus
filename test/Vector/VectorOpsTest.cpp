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

#include <vector>
#include "Vector/VectorOperators.hpp"
#include "TestUtilities.hpp"
#include "Constants.hpp"
#include <cmath>

namespace TensorCalculus {

  template<>
  struct Constants<double> {
    static const double epsilon;
  };

}

const double TensorCalculus::Constants<double>::epsilon = 1e-8;

using namespace TensorCalculus;
using namespace VectorOperators;

int main() {
  std::vector<double> x(3);
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  
  std::vector<double> y(3);
  y[0] = 2.0;
  y[1] = 3.0;
  y[2] = 4.0;

  add(2.0, x, y);
  TEST_EQUAL(y[0], 4.0);
  TEST_EQUAL(y[1], 7.0);
  TEST_EQUAL(y[2], 10.0);

  scale(y, 0.5);
  TEST_EQUAL(y[0], 2.0);
  TEST_EQUAL(y[1], 3.5);
  TEST_EQUAL(y[2], 5.0);
  
  double alpha = innerProduct(x, y);
  TEST_EQUAL(alpha, 24.0);

//   std::vector<double> z = add(1.0, x, -1.0, y);
//   TEST_EQUAL(z[0], -1.0);
//   TEST_EQUAL(z[1], -1.5);
//   TEST_EQUAL(z[2], -2.0);

  alpha = l2_norm(x);
  TEST_EQUAL_FLOAT(alpha, std::sqrt(14.0), Constants<double>::epsilon);

  alpha = max_norm(x);
  TEST_EQUAL(alpha, 3.0);

  std::vector<double> z(3);
  pointwiseProduct(x, y, z);
  TEST_EQUAL(z[0], 2.0);
  TEST_EQUAL(z[1], 7.0);
  TEST_EQUAL(z[2], 15.0); 

  l2_normalize(z);
  TEST_EQUAL_FLOAT(z[0], 2.0/16.673332000533065, Constants<double>::epsilon);
  TEST_EQUAL_FLOAT(z[1], 7.0/16.673332000533065, Constants<double>::epsilon);
  TEST_EQUAL_FLOAT(z[2], 15.0/16.673332000533065, Constants<double>::epsilon);
  // (sqrt (+ 4 (* 7 7) (* 15 15)))

  alpha = x * y;
  TEST_EQUAL(alpha, 24.0);

  y = x;
  TEST(equals(x, y));
  y *= -1.0;
  TEST_EQUAL(y[0], -x[0]);
  TEST_EQUAL(y[1], -x[1]);
  TEST_EQUAL(y[2], -x[2]);
  TEST(!equals(x, y));

  std::cout << "Tests passed!" << std::endl;
  return 0;
}
