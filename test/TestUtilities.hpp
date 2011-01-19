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

#ifndef __TEST_UTILITIES_HPP
#define __TEST_UTILITIES_HPP

#include <iostream>
#include <exception>
#include <cmath>

#define TEST(cond)					\
  if(!(cond)) {                                         \
    std::cerr << "Test failed: " << #cond << std::endl; \
    std::terminate();					\
  }

#define TEST_EQUAL(expr, value)                            \
  if((expr) != (value)) {                                  \
    std::cerr << "Test failed: "                           \
              << #expr << " == " << (expr) << " != "       \
              << (value) << " == " << #value << std::endl; \
    std::terminate();                                      \
  }

#define TEST_EQUAL_FLOAT(expr, value, eps)           \
  if (std::abs((expr) - (value)) >= eps) {           \
    std::cerr << "Test failed: |"                    \
              << #expr << " - " << #value << "| >= " \
              << eps << std::endl;                   \
    std::terminate();                                \
  }

#define TEST_FOR_EXCEPTION(expr, exception)                     \
  try {                                                         \
    expr;                                                       \
    std::cerr << "Test failed: Exception '"                     \
              << #exception << "' not thrown by executing '"    \
              << #expr << std::endl;                            \
    std::terminate();                                           \
  }                                                             \
  catch (exception& e) { }

#endif // __TEST_UTILITIES_HPP
