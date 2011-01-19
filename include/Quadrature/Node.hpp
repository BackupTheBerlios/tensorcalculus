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

#ifndef __NODE_HPP
#define __NODE_HPP

#include "Utilities/Utilities.hpp"

namespace TensorCalculus {

  template<typename Abscissa, typename Weight>
  class Node {
  public:
    typedef Abscissa AbscissaType;
    typedef Weight   WeightType;
  
    // Node(const Abscissa abscissa, const Weight weight) : abscissa(abscissa), weight(weight) { }
    Node(const AbscissaType& abscissa, const WeightType& weight) : abscissa(abscissa), weight(weight) { }
    
    const AbscissaType getAbscissa() const { return abscissa; }
    const WeightType getWeight() const { return weight; }

  private:
    AbscissaType abscissa;
    WeightType weight;
  };
  
  template<typename Abscissa, typename Weight>
  Node<Abscissa, Weight> make_node(const Abscissa abscissa, const Weight weight)
  {
    return Node<Abscissa, Weight>(abscissa, weight);
  }

} // namespace TensorCalculus

#endif // __NODE_HPP
