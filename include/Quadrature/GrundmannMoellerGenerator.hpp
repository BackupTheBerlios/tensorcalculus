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

#ifndef __GRUNDMANNMOELLERGENERATOR_HPP
#define __GRUNDMANNMOELLERGENERATOR_HPP

#include <vector>
#include <algorithm>
#include <stdexcept>

#include "Quadrature/Node.hpp"
#include "Utilities/Factorials.hpp"
#include "Utilities/Utilities.hpp"
#include "Vector/VectorOperators.hpp"

namespace TensorCalculus {

  class QuadratureError : public std::runtime_error {
  public:
    QuadratureError(const std::string& what) : std::runtime_error(what) { }
  };

  template<typename T>
  class GrundmannMoellerGenerator {
  public:
    // typedef unsigned int Accuracy;
    typedef Node<std::vector<T>, T> NodeType;
    typedef std::vector<NodeType> NodesType;
    
    //~ class Node {
    //~ public:
      //~ typedef std::vector<T> Abscissa;
      //~ typedef T Weight;
      //~ typedef std::vector<T> Weights;
//~ 
    //~ private:
      //~ Abscissa abscissa;
      //~ Weights weights;
//~ 
    //~ public:
      //~ Node(const Abscissa& abscissa, const Weights& weights)
        //~ : abscissa(abscissa), weights(weights) { }
//~ 
//~ #ifdef __GXX_EXPERIMENTAL_CXX0X__
      //~ Node(Node&& node)
        //~ : abscissa(std::move(node.abscissa)),
          //~ weights(std::move(node.weights))
      //~ { }
//~ 
      //~ Node& operator = (Node&& node) {
        //~ abscissa = std::move(node.abscissa);
        //~ weights = std::move(node.weights);
        //~ return *this;
      //~ }
//~ #endif
//~ 
      //~ const Abscissa& get_abscissa() const { return abscissa; }
      //~ const Weights& get_weights() const { return weights; }
    //~ };
//~ 
    //~ typedef std::vector<Node> Nodes;

  private:
    int dimension;

    std::vector< std::vector<T> > grundmann_points(const int level) const {
      const std::vector< std::vector<int> > parts = partition(level, dimension + 1);

      std::vector< std::vector<T> > result;
      result.reserve(parts.size());

      for (unsigned int i = 0; i < parts.size(); ++i) {
        std::vector<T> abscissa;
        const std::vector<int>& p = parts[i];
        abscissa.reserve(dimension);
        for (int mu = 0; mu < dimension; ++mu) {
          abscissa.push_back(static_cast<T>(2*p[mu + 1] + 1) /
                             static_cast<T>(2*level + 1 + dimension));
        }
        result.push_back(abscissa);
      }
      return result;
    }

    T grundmann_weight(const int N, const int level) const {
      T weight = ((N - level) % 2 == 0) ? 1 : -1;
      const int max = std::max(2*N + 1, N + level + dimension + 1);
      for (int j = 1; j <= max; ++j) {
        if (j <= 2*N + 1) {
          weight *= 2*level + 1 + dimension;
          if (j <= 2*N) {
            weight /= 2;
          }
        }
        if (j <= N - level) {
          weight /= j;
        }
        if (j <= N + level + dimension + 1) {
          weight /= (N + level + dimension + 1 - j + 1);
        }
      }
      return weight;
    }

  public:
    GrundmannMoellerGenerator(const int dimension) : dimension(dimension) { }

    std::vector< Node<std::vector<T>, std::vector<T> > > nodes(const int N, const int null_rules) {
      std::vector< Node<std::vector<T>, std::vector<T> > > result;
      const int reserve = choose<int>(dimension + N + 1, N);
      // std::cout << reserve << '\n';
      result.reserve(reserve);

      for (int level = 0; level <= N; ++level) {
        const std::vector< std::vector<T> > abscissae = grundmann_points(level);
        const T weight = grundmann_weight(N, level);

        std::vector<T> weights;
        weights.reserve(1 + null_rules);
        weights.push_back(weight);

        T factor = 1;
        // const max_null_level =  std::max(null_rules, N - level);
        const T alpha = 2*level + 1 + dimension;

        for (int j = 0; j < null_rules; ++j) {
          if (j < N - level) {
            factor *= -4*((N - level - j)/alpha)*((N + level + 1 + dimension - j)/alpha);
            weights.push_back((1 - factor) * weight);
          } else {
            weights.push_back(weight);
          }
        }

        for (unsigned int i = 0; i < abscissae.size(); ++i) {
          result.push_back(Node<std::vector<T>, std::vector<T> >(abscissae[i], weights));
        }
      }

      return result;
    }

    NodesType nodes(const int N) {
      NodesType result;
      const int reserve = choose<int>(dimension + N + 1, N);
      result.reserve(reserve);

      for (int level = 0; level <= N; ++level) {
        const std::vector<typename NodeType::AbscissaType> abscissae = grundmann_points(level);
        const T weight = grundmann_weight(N, level);

        for (unsigned int i = 0; i < abscissae.size(); ++i) {
          Node<std::vector<T>, T> new_node(abscissae[i], weight);
          result.push_back(new_node);
        }
      }

      return result;
    }

  };

  //~ template<typename T>
  //~ std::ostream& operator << (std::ostream& stream,
                             //~ const typename GrundmannMoellerGenerator<T>::Node& node) {
    //~ using VectorOperators::operator <<;
    //~ stream << "(" << node.get_abscissa() << ", "
           //~ << node.get_weights() << ")";
    //~ return stream;
  //~ }

} // namespace TensorCalculus

#endif // __GRUNDMANNMOELLERGENERATOR_HPP
