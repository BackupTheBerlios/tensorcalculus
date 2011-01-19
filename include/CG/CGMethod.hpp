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

#ifndef __CGMETHOD_HPP
#define __CGMETHOD_HPP

#include "Utilities/Logger.hpp"
#include "StandardTraits.hpp"

namespace TensorCalculus {

  template<typename VectorSpace, template <typename> class VectorSpaceTraits = StandardVectorSpaceTraits>
  class CGMethod {
  public:
    explicit CGMethod(VectorSpaceTraits<VectorSpace> traits = VectorSpaceTraits<VectorSpace>()) : traits(traits) { }
    
    template<typename Function,
             typename StepSizeRule,
             typename UpdateRule,
             typename AbortCriterion,
             typename Logger>
    typename VectorSpaceTraits<VectorSpace>::Vectors
    minimize(typename VectorSpaceTraits<VectorSpace>::Vectors& x, Function df,
             StepSizeRule& step_size_rule, UpdateRule& update_rule,
             AbortCriterion abort_criterion, Logger logger) // not const since update and scale can change traits
    {
      // typename VectorSpaceTraits<VectorSpace>::Vectors x = start;
      typename VectorSpaceTraits<VectorSpace>::Vectors g = df(x);
      typename VectorSpaceTraits<VectorSpace>::Vectors g_new = g;
      typename VectorSpaceTraits<VectorSpace>::Vectors d = g;
      traits.scale(-1, d);
      int iterations = 0;
      logger(iterations, x, g);
      while (!abort_criterion(iterations, g)) {
        typename VectorSpaceTraits<VectorSpace>::Scalars alpha = step_size_rule(df, g, x, d);
        // if (alpha == static_cast<typename VectorSpaceTraits<VectorSpace>::Scalars>(0)) {
        //   df(x, d);
        //   g = d;
        //   traits.scale(-1, d);
        //   alpha = step_size_rule(df, g, x, d);
        // }
        traits.update(alpha, d, x);
        df(x, g_new);
        traits.scale(update_rule(g_new, g), d);
        traits.update(-1, g_new, d);
        g = g_new;
        ++iterations;
        logger(iterations, x, g);
      }
      return x;
    }
    
    template<typename Function,
             typename StepSizeRule,
             typename UpdateRule,
             typename AbortCriterion>
    typename VectorSpaceTraits<VectorSpace>::Vectors
    minimize(typename VectorSpaceTraits<VectorSpace>::Vectors& x, Function df,
             StepSizeRule& step_size_rule, UpdateRule update_rule,
             AbortCriterion abort_criterion)
    {
      minimize(x, df, step_size_rule, update_rule, abort_criterion, NullLogger());
    }

  private:
    VectorSpaceTraits<VectorSpace> traits;
  }; 

} // namespace TensorCalculus

#endif // __CGMETHOD_HPP
