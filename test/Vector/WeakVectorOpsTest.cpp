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

#include "Vector/WeakVector.hpp"
#include "Vector/WeakVectorOperators.hpp"
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

int main() {
  double a[] = { 1.0, 2.0, 3.0 };
  WeakVector<double> aw(a, 3);
  
  double b[] = { 1.0, 1.0, 1.0 };
  WeakVector<const double> bw(b, 3);

  std::vector<double> x = getVectorCopy(bw);
  TEST_EQUAL(x[0], b[0]);
  TEST_EQUAL(x[1], b[1]);
  TEST_EQUAL(x[2], b[2]);

  add(2.0, bw, aw);
  TEST_EQUAL(aw[0], 3.0);
  TEST_EQUAL(aw[1], 4.0);
  TEST_EQUAL(aw[2], 5.0);

  add(-1.0, x, aw);
  TEST_EQUAL(aw[0], 2.0);
  TEST_EQUAL(aw[1], 3.0);
  TEST_EQUAL(aw[2], 4.0);

  add(2.0, bw, x);
  TEST_EQUAL(x[0], 3.0);
  TEST_EQUAL(x[1], 3.0);
  TEST_EQUAL(x[2], 3.0);

  scale(aw, 2.0);
  TEST_EQUAL(aw[0], 4.0);
  TEST_EQUAL(aw[1], 6.0);
  TEST_EQUAL(aw[2], 8.0);
  TEST_EQUAL(a[0], 4.0);
  TEST_EQUAL(a[1], 6.0);
  TEST_EQUAL(a[2], 8.0);

  double alpha = innerProduct(aw, bw);
  TEST_EQUAL(alpha, 18.0);
  
  alpha = innerProduct(bw, x);
  TEST_EQUAL(alpha, 9.0);
  alpha = innerProduct(x, bw);
  TEST_EQUAL(alpha, 9.0);

  double c[] = { 0.0, 0.0, 0.0 };
  WeakVector<double> cw(c, 3);
  pointwiseProduct(aw, bw, cw);
  TEST_EQUAL(cw[0], 4.0);
  TEST_EQUAL(cw[1], 6.0);
  TEST_EQUAL(cw[2], 8.0);

  alpha = l2_norm(aw);
  TEST_EQUAL(alpha, std::sqrt(16.0 + 36.0 + 64.0));

  l2_normalize(aw);
  TEST_EQUAL(aw[0], 4.0/std::sqrt(16.0 + 36.0 + 64.0));
  TEST_EQUAL(aw[1], 6.0/std::sqrt(16.0 + 36.0 + 64.0));
  TEST_EQUAL(aw[2], 8.0/std::sqrt(16.0 + 36.0 + 64.0));

  // aw *= std::sqrt(16.0 + 36.0 + 64.0);
  // TEST_EQUAL(aw[0], 4.0);
  // TEST_EQUAL(aw[1], 6.0);
  // TEST_EQUAL(aw[2], 8.0);
  aw[0] = 4.0;
  aw[1] = 6.0;
  aw[2] = 8.0;

  WeakVector<double>::const_iterator i_max = max_element(aw);
  TEST_EQUAL(*i_max, 8.0);
  WeakVector<double>::const_iterator i_min = min_element(aw);
  TEST_EQUAL(*i_min, 4.0);

  alpha = max_norm(aw);
  TEST_EQUAL(alpha, 8.0);

  std::cout << "Tests passed!" << std::endl;
  return 0;
}
