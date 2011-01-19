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

#ifndef __QUADRATURE_HPP
#define __QUADRATURE_HPP

#include <vector>
#include <functional>
#include <cmath>
#include <boost/progress.hpp>

#include "BlasInterface.hpp"
#include "Mesh/Domain.hpp"
#include "Quadrature/Node.hpp"
#include "StandardTraits.hpp"

namespace TensorCalculus {
  
  template<typename T,
           typename Functor,
           typename SimplexFunctor,
           typename Abscissa,
           typename Weight>
  void quadrature(const Domain<T>& domain,
                  SimplexFunctor simplex_functor,
                  Functor functor,
                  const std::vector< Node<Abscissa, Weight> >& nodes,
                  typename SimplexFunctor::result_type& result) {
    for (unsigned int i = 0; i < domain.countSimplices(); ++i) {
      const ValuesSimplex<T>& values_simplex = domain.getSimplex(i).getValuesSimplex();
      typename Functor::result_type value = functor.getZero();
      for (unsigned int j = 0; j < nodes.size(); ++j) {
        const Node<Abscissa, Weight>& node = nodes[j];
        functor(values_simplex.transform(node.getAbscissa()), node.getWeight(), value);
      }
      using std::abs;
      simplex_functor(i, value, abs(values_simplex.getFactor()), result);
    }    
  }

  template<typename T,
           typename Functor,
           typename SimplexFunctor,
           typename Abscissa,
           typename Weight>
  void productQuadrature(const Domain<T>& domainA,
                         const Domain<T>& domainB,
                         SimplexFunctor simplex_functor,
                         Functor functor,
                         const std::vector< Node<Abscissa, Weight> >& nodesA,
                         const std::vector< Node<Abscissa, Weight> >& nodesB,
                         typename SimplexFunctor::result_type& result) {
    for (unsigned int iA = 0; iA < domainA.countSimplices(); ++iA) {
      const ValuesSimplex<T>& values_simplexA = domainA.getSimplex(iA).getValuesSimplex();
      
      for (unsigned int iB = 0; iB < domainB.countSimplices(); ++iB) {
        const ValuesSimplex<T>& values_simplexB = domainB.getSimplex(iB).getValuesSimplex();
        typename Functor::result_type value = functor.getZero();
        
        for (unsigned int jA = 0; jA < nodesA.size(); ++jA) {
          const Node<Abscissa, Weight>& nodeA = nodesA[jA];
          for (unsigned int jB = 0; jB < nodesB.size(); ++jB) {
            const Node<Abscissa, Weight>& nodeB = nodesB[jB];
            functor(values_simplexA.transform(nodeA.getAbscissa()),
                    values_simplexB.transform(nodeB.getAbscissa()),
                    nodeA.getWeight() * nodeB.getWeight(), value);
          }
        }
        
        using std::abs;
        simplex_functor(iA, iB, value,
                        abs(values_simplexA.getFactor() * values_simplexB.getFactor()),
                        result);
      }
    }
  }

  template<typename T,
           typename Functor,
           typename SimplexFunctor,
           typename Abscissa,
           typename Weight,
           template<typename> class VectorSpaceTraits1,
           template<typename> class VectorSpaceTraits2>
  typename SimplexFunctor::result_type
  quadrature(const Domain<T>& domain, SimplexFunctor simplex_functor, Functor functor,
             const std::vector< Node<Abscissa, Weight> >& nodes,
             VectorSpaceTraits1<typename SimplexFunctor::result_type> simplex_functor_traits,
             VectorSpaceTraits2<typename Functor::result_type> functor_traits)
  {
    // typename Generator::Nodes nodes = generator.nodes(accuracy);
    typename SimplexFunctor::result_type result;
    for (int i = 0; i < domain.countSimplices(); ++i) {
      const ValuesSimplex<T>& values_simplex = domain.getSimplex(i).getValuesSimplex();
      // const typename Domain<T>::Simplex& simplex = domain.getSimplex(i);
      typename Functor::result_type value;
      for (int j = 0; j < nodes.size(); ++j) {
        const Node<Abscissa, Weight>& node = nodes[j];
        functor_traits.update(node.getWeight(), functor(values_simplex.transform(node.getAbscissa())), value);
      }
      using std::abs;
      simplex_functor_traits.update(abs(values_simplex.getFactor()), simplex_functor(i, value), result);
    }
    return result;

    //~ typename Functor::result_type temp = functor(0, domain.getSimplex(0).getValuesSimplex().transform(nodes[0].getAbscissa()));
    //~ typename Functor::result_type result = temp;
    //~ traits.scale(nodes[0].getWeight(), result);
    //~ {
      //~ const typename Domain<T>::Simplex& simplex = domain.getSimplex(0);
      //~ for (int j = 1; j < nodes.size(); ++j) {
        //~ const Node<Abscissa, Weight>& node = nodes[j];
        //~ traits.update(node.getWeight(),
                      //~ functor(0, simplex.getValuesSimplex().transform(node.getAbscissa())),
                      //~ result);
      //~ }
    //~ }
    //~ for (int i = 1; i < domain.countSimplices(); ++i) {
      //~ const typename Domain<T>::Simplex& simplex = domain.getSimplex(i);
//~ 
      //~ for (int j = 0; j < nodes.size(); ++j) {
        //~ const Node<Abscissa, Weight>& node = nodes[j];
        //~ traits.update(node.getWeight(),
                      //~ functor(i, simplex.getValuesSimplex().transform(node.getAbscissa())),
                      //~ result);
      //~ }
    //~ }
    //~ return result;
  }

  template<typename T,
           typename Functor,
           typename SimplexFunctor,
           typename Abscissa,
           typename Weight>
  typename SimplexFunctor::result_type
  quadrature(const Domain<T>& domain, SimplexFunctor simplex_functor, Functor functor,
             const std::vector< Node<Abscissa, Weight> >& nodes)
  {
    return quadrature(domain, simplex_functor, functor, nodes,
                      StandardVectorSpaceTraits<typename SimplexFunctor::result_type>(),
                      StandardVectorSpaceTraits<typename Functor::result_type>());
  }

  template<typename T,
           typename Functor,
           typename SimplexFunctor,
           typename AbscissaA, typename AbscissaB,
           typename WeightA, typename WeightB,
           template<typename> class VectorSpaceTraits1,
           template<typename> class VectorSpaceTraits2>
  typename SimplexFunctor::result_type
  product_quadrature(const Domain<T>& domainA, const Domain<T>& domainB,
                     SimplexFunctor simplex_functor, Functor functor,
                     const std::vector< Node<AbscissaA, WeightA> >& nodesA,
                     const std::vector< Node<AbscissaB, WeightB> >& nodesB,
             VectorSpaceTraits1<typename SimplexFunctor::result_type> simplex_functor_traits,
             VectorSpaceTraits2<typename Functor::result_type> functor_traits)
  {
    typename SimplexFunctor::result_type result = simplex_functor_traits.zero();
    boost::progress_display p(domainA.countSimplices() * domainB.countSimplices());
    for (int i = 0; i < domainA.countSimplices(); ++i) {
      const ValuesSimplex<T>& values_simplexA = domainA.getSimplex(i).getValuesSimplex();
      for (int j = 0; j < domainB.countSimplices(); ++j) {
        const ValuesSimplex<T>& values_simplexB = domainB.getSimplex(j).getValuesSimplex();
        typename Functor::result_type value = functor_traits.zero();
        for (int k = 0; k < nodesA.size(); ++k) {
          const Node<AbscissaA, WeightA>& nodeA = nodesA[k];
          for (int l = 0; l < nodesB.size(); ++l) {
            const Node<AbscissaB, WeightB>& nodeB = nodesB[l];
            functor_traits.update(nodeA.getWeight() * nodeB.getWeight(),
                                  functor(values_simplexA.transform(nodeA.getAbscissa()),
                                          values_simplexB.transform(nodeB.getAbscissa())),
                                  value);
          }
        }
        using std::abs;
        simplex_functor_traits.update(abs(values_simplexA.getFactor() * values_simplexB.getFactor()),
                                      simplex_functor(i, j, value), result);
        ++p;
      }
    }
    return result;
  }

  template<typename T,
           typename Functor,
           typename SimplexFunctor,
           typename AbscissaA, typename AbscissaB,
           typename WeightA, typename WeightB>
  typename SimplexFunctor::result_type
  product_quadrature(const Domain<T>& domainA, const Domain<T>& domainB,
                     SimplexFunctor simplex_functor, Functor functor,
                     const std::vector< Node<AbscissaA, WeightA> >& nodesA,
                     const std::vector< Node<AbscissaB, WeightB> >& nodesB)
  {
    return product_quadrature(domainA, domainB, simplex_functor, functor, nodesA, nodesB,
                              StandardVectorSpaceTraits<typename SimplexFunctor::result_type>(),
                              StandardVectorSpaceTraits<typename Functor::result_type>());
  }

  //~ template<typename T,
           //~ typename Functor,
           //~ typename Abscissa,
           //~ typename Weight>
  //~ //std::vector<T>
  //~ typename Functor::result_type
  //~ quadrature(const Domain<T>& domain, Functor function, const std::vector< Node<Abscissa, Weight> >& nodes)
  //~ {
    //~ return quadrature(domain, function, nodes,
                      //~ StandardVectorSpaceTraits<typename Functor::result_type>());
  //~ }
//~ 
  //~ template<typename T,
           //~ typename Function,
           //~ typename Abscissa,
           //~ typename Weight,
           //~ template<typename> class VectorSpaceTraits>
  //~ typename function_type<Function>::result_type
  //~ quadrature(const Domain<T>& domain, Function function, const std::vector< Node<Abscissa, Weight> >& nodes,
             //~ VectorSpaceTraits<typename function_type<Function>::result_type> traits)
  //~ {
    //~ return quadrature(domain, std::ptr_fun(function), nodes, traits);
  //~ }
//~ 
  //~ template<typename T,
           //~ typename Function,
           //~ typename Abscissa,
           //~ typename Weight>
  //~ typename function_type<Function>::result_type
  //~ quadrature(const Domain<T>& domain, Function function, const std::vector< Node<Abscissa, Weight> >& nodes)
  //~ {
    //~ return quadrature(domain, std::ptr_fun(function), nodes,
                      //~ StandardVectorSpaceTraits<typename function_type<Function>::result_type>());
  //~ }

  template<typename T, typename Function, typename SimplexFunction, typename Generator>
  T quadrature(const Domain<T>& domain, Function function, SimplexFunction simplex_function, Generator generator, typename Generator::Accuracy accuracy)
  {
    typename Generator::Nodes nodes = generator.nodes(accuracy);
    T result = 0;
    for (int i = 0; i < domain.countSimplices(); ++i) {
      const typename Domain<T>::Simplex& simplex = domain.getSimplex(i);
      std::vector<T> factors = simplex_function(simplex.getOriginalSimplices());
      std::vector<T> values(factors.size());
      for (int j = 0; j < nodes.size(); ++j) {
        const typename Generator::Node& node = nodes[j];
        std::vector<T> update = function(simplex.getValuesSimplex().transform(node.get_abscissa()));
#if defined(_DEBUG) || !defined(NDEBUG) || defined(RANGE_CHECKS_ON)
        if (update.size() != factors.size()) {
          throw std::length_error("Update and factors don't match");
        }
#endif
        Blas<T>::axpy(update.size(), node.get_weights()[0], &update[0], 1, &values[0], 1);
      }
      result += Blas<T>::dot(values.size(), &values[0], 1, &factors[0], 1);
    }
    return result;
  }

  template<typename T, typename Function, typename SimplexFunction, typename Generator>
  T product_quadrature(const Domain<T>& domain1, const Domain<T>& domain2, Function function, SimplexFunction simplex_function, Generator generator, typename Generator::Accuracy accuracy)
  {
    typename Generator::Nodes nodes = generator.nodes(accuracy);
    T result = 0;
    boost::progress_display p(domain1.countSimplices() * domain2.countSimplices());
    for (int i = 0; i < domain1.countSimplices(); ++i) {
      // std::cout << i << std::endl;
      const typename Domain<T>::Simplex& simplex1 = domain1.getSimplex(i);
      for (int j = 0; j < domain2.countSimplices(); ++j) {
        const typename Domain<T>::Simplex& simplex2 = domain2.getSimplex(j);
        const std::vector<T> factors = simplex_function(simplex1.getOriginalSimplices(), simplex2.getOriginalSimplices());
        std::vector<T> values(factors.size());
        
        for (int k = 0; k < nodes.size(); ++k) {
          const typename Generator::Node& node1 = nodes[k];
          const T weight1 = node1.get_weights()[0];
          const std::vector<T> arg1 = simplex1.getValuesSimplex().transform(node1.get_abscissa());
          
          for (int l = 0; l < nodes.size(); ++l) {
            const typename Generator::Node& node2 = nodes[l];
            const T weight2 = node2.get_weights()[0];
            const std::vector<T> arg2 = simplex1.getValuesSimplex().transform(node2.get_abscissa());
            
            std::vector<T> update = function(arg1, arg2);
#if defined(_DEBUG) || !defined(NDEBUG) || defined(RANGE_CHECKS_ON)
            if (update.size() != factors.size()) {
              throw std::length_error("Update and factors don't match");
            }
#endif
            Blas<T>::axpy(update.size(), weight1 * weight2, &update[0], 1, &values[0], 1);
          }
        }
        result += Blas<T>::dot(values.size(), &values[0], 1, &factors[0], 1);
        ++p;
      }
    }
    return result;
  }

} // namespace TensorCalculus


#endif // __QUADRATURE_HPP
