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
// #include <tr1/array>
#include <cmath>
// #include "ArrayOperators.hpp"
// #include "CurryRule.hpp"
// #include "ThreePG.hpp"
// #include "SecantMethod.hpp"
#include "Vector/VectorOperators.hpp"
#include "CG/SecantStepSize.hpp"
#include "CG/CGMethod.hpp"
#include "CG/PoljakRibiereUpdate.hpp"
#include "Utilities/Logger.hpp"

// using std::tr1::array;
// using std::size_t;

using namespace TensorCalculus;
using namespace VectorOperators;

//~ template<typename T>
//~ struct VectorTraits {
  //~ typedef std::vector<T> Vectors;
  //~ typedef T Scalars;
//~ 
  //~ void update(const Scalars alpha, const Vectors& x, Vectors& y)
  //~ {
    //~ TensorCalculus::add(alpha, x, y);
  //~ }
  //~ 
  //~ void scale(const Scalars alpha, Vectors& x)
  //~ {
    //~ TensorCalculus::scale(x, alpha);
  //~ }
//~ 
  //~ Scalars innerProduct(const Vectors& x, const Vectors& y)
  //~ {
    //~ return TensorCalculus::innerProduct(x, y);
  //~ }
//~ };
//~ 
//~ const size_t N = 2;
//~ 
//~ double f(const std::vector<double>& x) {
  //~ return (1 + std::sin(x[0])*std::sin(x[0])) * std::sin(x[1]);
//~ }

// double f(const array<double, N>& x) {
//   // return 2*(x[0] - 1.0)*(x[0] - 1.0)*(x[0] - 1.0)*(x[0] - 1.0)
//   //   + 3*(x[1] - 1.0)*(x[1] - 1.0);
//   return (1 + std::sin(x[0])*std::sin(x[0])) * std::sin(x[1]);
// }

struct df {
public:
  std::vector<double> operator () (const std::vector<double>& x) const {
    std::vector<double> result(2);
    operator () (x, result);
    return result;
  }
  
  void operator () (const std::vector<double>& x, std::vector<double>& result) const {
    result[0] = 2*std::sin(x[0])*std::cos(x[0])*std::sin(x[1]);
    result[1] = (1 + std::sin(x[0])*std::sin(x[0])) * std::cos(x[1]);
  }
};

// array<double, N> df(const array<double, N>& x) {
//   array<double, N> result;
//   // result[0] = 8*(x[0] - 1.0)*(x[0] - 1.0)*(x[0] - 1.0);
//   // result[1] = 6*(x[1] - 1.0);
//   result[0] = 2*std::sin(x[0])*std::cos(x[0])*std::sin(x[1]);
//   result[1] = (1 + std::sin(x[0])*std::sin(x[0])) * std::cos(x[1]);
//   return result;
// }

// struct Vectors {
//   typedef Vector<double> vectors;
//   typedef double scalars;
// };

double norm(const std::vector<double>& x)
{
  double result = 0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    const double xi = x[i];
    result += xi*xi;
  }
  return std::sqrt(result);
}

bool abort_criterion(int iterations, const std::vector<double>& g)
{
  return iterations > 100 || norm(g) < 1.0e-6;
}

int main() {

  CGMethod< std::vector<double> > cg_method;
  // SecantRule< std::vector<double> > stepsize(0.5);
  SecantStepSize< std::vector<double> > stepsize2(0.25, 1e-5);
  PoljakRibiereUpdate< std::vector<double> > update(true);
  std::vector<double> start(2, 0.5);
  std::cout.precision(8);
  std::cout << "Result: "
            << cg_method.minimize(start,
                                  df(),
                                  stepsize2,
                                  update,
                                  &abort_criterion,
                                  stream_logger(std::cout))
            << std::endl;
  
  // std::cout << "Result: " << cg_minimize(start,
  //                                        &df,
  //                                        secant_rule(&df, 0.5),
  //                                        poljak_ribiere_update<double>(),
  //                                        &abort_criterion) << std::endl; // ,
//                                         stream_logger(std::cout)) << std::endl;
  
//   CGMethod<
//      array<double, N>,
//      double,
//     CurryRule<array<double, N>, double, ThreePG<double> >
//   > cg(
//     CurryRule<array<double, N>, double, ThreePG<double> >(
//       ThreePG<double>(0.75, 0.75),
//       0.0000001,
//       10.0
//     )
//   );
  
//   CGMethod<
//     Vector<double>,
//     SecantMethod<Vector<double> >,
//     PoljakRibiereUpdate<Vector<double> >
//   > cg(SecantMethod<Vector<double> >(0.5));
// 
//   std::cout.precision(8);
//   Vector<double> start(2, 0.5); //  = { { 0.5, 0.5 } };
//   cg.initialize(start, 0.0000001, &df);
//   while (!cg.iterate()) {
//     std::cout << cg.value() << std::endl;
//     // std::cout << cg.get_d() << std::endl;
//   }
//   std::cout << cg.value() << std::endl;
  return 0;
}
