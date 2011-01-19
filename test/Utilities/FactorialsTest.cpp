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

#include "Utilities/Factorials.hpp"

using namespace TensorCalculus;

int main() {
  std::vector< std::vector<int> > p = orderedPartition(3, 3);
  for (unsigned int i = 0; i < p.size(); ++i) {
    const std::vector<int>& pi = p[i];
    for (unsigned int j = 0; j < pi.size(); ++j) {
      std::cout << pi[j];
    }
    std::cout << '\n';
  }
  return 0;
}
