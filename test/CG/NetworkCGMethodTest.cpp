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
#include "Representation/DistanceFunctionGradient.hpp"
#include "Representation/CPTensorRepresentation.hpp"
#include "Representation/TensorTrainRepresentation.hpp"

#include <Utilities/Random.hpp>


using namespace TensorCalculus;
using namespace VectorOperators;



double norm(const std::vector< std::vector<double> >& x)
{
  double result = 0;
  for (int i = 0; i < x.size(); ++i) {
	for (int k = 0; k < x[i].size(); ++k) {
		const double xi = x[i][k];
		result += xi*xi;
	}

  }
  return std::sqrt(result);
}

bool abort_criterion(int iterations, const std::vector< std::vector<double> >& g)
{
  return iterations > 1 || norm(g) < 1.0e-9;
}

int main() {
  int d = 3;

  int n = 3;

  const TensorTrainRepresentation<double> tt = createRandomTensorTrain<double>(3, d, n, Random<double>());//createRandomTensorTrain(10, d, 5, Utilities<double>::rand);

  TensorTrainRepresentation<double> tt_approx = createRandomTensorTrain<double>(3, d, n, Random<double>());//createRandomTensorTrain(10, d, 5, Utilities<double>::rand);

  DistanceFunctionGradient< TensorTrainRepresentation<double> > gradient(tt, tt_approx);

  CGMethod< std::vector< std::vector<double> > > cg_method;
  SecantStepSize< std::vector< std::vector<double> > > stepsize2(1, 1);
  PoljakRibiereUpdate< std::vector< std::vector<double> > > update(true);
  std::vector< std::vector<double> > start = tt_approx.getV();

  std::cout << l2norm(tt, tt_approx) / l2norm(tt) << std::endl;

  std::cout << "Result: "
            << cg_method.minimize(start,
								  gradient,
                                  stepsize2,
                                  update,
                                  &abort_criterion,
                                  stream_logger(std::cout))
            << std::endl;
  tt_approx.setV(start);
  std::cout << l2norm(tt, tt_approx) / l2norm(tt) << std::endl;
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
