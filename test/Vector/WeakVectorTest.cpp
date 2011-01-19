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
#include "Vector/WeakVector.hpp"
#include "TestUtilities.hpp"

using namespace TensorCalculus;

int main() {
  double x[] = { 1.0, 2.0, 3.0 };
  WeakVector<double> wx(x, sizeof(x)/sizeof(double));
  TEST_EQUAL(wx.size(), 3);
  TEST_EQUAL(wx[0], 1.0);
  TEST_EQUAL(wx[1], 2.0);
  TEST_EQUAL(wx[2], 3.0);

  WeakVector<double> wy(wx);
  TEST_EQUAL(wy[0], 1.0);
  TEST_EQUAL(wy[1], 2.0);
  TEST_EQUAL(wy[2], 3.0);

  std::vector<double> y(3, 1.0);
  WeakVector<const double> wz(y);
  TEST_EQUAL(wz[0], 1.0);
  TEST_EQUAL(wz[1], 1.0);
  TEST_EQUAL(wz[2], 1.0);

  double a[] = { 0.0, 0.0, 0.0 };
  WeakVector<double> aw(a, 3);
  aw = wz;
  TEST_EQUAL(aw[0], 1.0);
  TEST_EQUAL(aw[1], 1.0);
  TEST_EQUAL(aw[2], 1.0);

  std::vector<double> b(3, 2.0);
  aw = b;
  TEST_EQUAL(aw[0], 2.0);
  TEST_EQUAL(aw[1], 2.0);
  TEST_EQUAL(aw[2], 2.0);

  std::cout << "Tests passed!" << std::endl;
  return 0;
}
