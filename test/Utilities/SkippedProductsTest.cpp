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
#include "Vector/VectorOperators.hpp"
#include "Utilities/SkippedProducts.hpp"

int x_data[3][3] = { { 1, 2, 3 },
                     { 2, 3, 4 },
                     { 3, 4, 5 } };

int main()
{
  std::vector< std::vector<int> > factors;
  factors.reserve(3);
  for (int i = 0; i < 3; ++i) {
    factors.push_back(std::vector<int>(&x_data[i][0], &x_data[i][3]));
  }
  TensorCalculus::SkippedProducts< std::vector<int> > sp(factors);
  
  using TensorCalculus::VectorOperators::operator <<;
  for (int i = 0; i < sp.countSkippedProducts(); ++i) {
    std::cout << sp[i] << '\n';
  }
  std::cout << sp.getFullProduct() << std::endl;

  sp.clear();
  
  using TensorCalculus::VectorOperators::operator <<;
  std::cout << "Empty:\n";
  for (int i = 0; i < sp.countSkippedProducts(); ++i) {
    std::cout << sp[i] << '\n';
  }
  std::cout << sp.getFullProduct() << std::endl;

  for (int i = 0; i < 3; ++i) {
    sp.pushFactor(std::vector<int>(&x_data[i][0], &x_data[i][3]));
    
    std::cout << "Pushed " << i+1 << "-th element\n";
    for (int i = 0; i < sp.countSkippedProducts(); ++i) {
      std::cout << sp[i] << '\n';
    }
    std::cout << sp.getFullProduct() << std::endl;
  }
  return 0;
}
