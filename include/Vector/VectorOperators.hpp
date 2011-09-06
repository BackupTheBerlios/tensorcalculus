/*
 * Copyright (C) 2010 Philipp Wähnert
 *               2011 Stefan Handschuh
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

/*!
 * \file VectorOperators.hpp
 * \brief Declares and implements the usual operators for the std::vector class
 * \author Philipp Wähnert
 */

#ifndef __VECTOR_OPERATORS_HPP
#define __VECTOR_OPERATORS_HPP

#include <vector>            // std::vector<T>
#include <algorithm>         // std::min, std::max, std::min_element, ...
#include <ostream>           // std::ostream
#include <iterator>          // std::ostream_iterator<T>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <cstring>           // std::memmove

#include "BlasInterface.hpp" // Blas<T>
#include "Utilities/Functional.hpp"    // abs_less<T>
#include "Constants.hpp"     // Constants<T>
#include "StandardTraits.hpp"

namespace TensorCalculus {

  /// \brief Performs \f$y:=\alpha x+y\f$
  template<typename T>
  void add(T alpha, const std::vector<T>& x, std::vector<T>& y)
  {
    if (x.size() == 0) return;
    if (y.size() == 0) {
      y.resize(x.size());
    }
#if defined(_DEBUG) || !defined(NDEBUG) || defined(RANGE_CHECKS_ON)
    if (x.size() != y.size()) {
      throw std::invalid_argument("Vectors don't match");
    }
#endif
    Blas<T>::axpy(x.size(), alpha, &x[0], 1, &y[0], 1);
  }

  template<typename T>
  void add(T alpha, const T * x, std::vector<T> &y) {
	Blas<T>::axpy(y.size(), alpha, x, 1, &y[0], 1);
  }

  template<typename T>
  void set(T value, std::vector<T> &x) {
    Blas<T>::copy(x.size(), &value[0], 0, &x[0]);
  }

  /// \brief Performs \f$x:=\alpha x\f$
  template<typename T>
  void scale(std::vector<T>& x, const T alpha)
  {
    Blas<T>::scal(x.size(), alpha, &x[0], 1);
  }

  /// \brief Returns \f$x\cdot y\f$
  template<typename T>
  const T innerProduct(const std::vector<T>& x, const std::vector<T>& y)
  {
    const typename std::vector<T>::size_type min_size =
      std::min(x.size(), y.size());

    return Blas<T>::dot(min_size, &x[0], 1, &y[0], 1);
  }

  //~ /// \brief Returns a vector containing \f$\alpha x+\beta y\f$
  //~ template<typename T>
  //~ const std::vector<T> add(T alpha, const std::vector<T>& x,
                           //~ T beta, const std::vector<T>& y)
  //~ {
    //~ const typename std::vector<T>::size_type min_size =
      //~ std::min(x.size(), y.size());
    //~ std::vector<T> result(x);
//~ 
    //~ // result *= alpha;
    //~ Blas<T>::axpy(min_size, beta, &y[0], 1, &result[0], alpha);
//~ 
    //~ return result;
  //~ }

  /// \brief Returns \f$\|x\|_2\f$
  template<typename T>
  const T l2_norm(const std::vector<T>& x)
  {
    return Blas<T>::nrm2(x.size(), &x[0], 1);
  }

  /// \brief Returns an iterator to the best element according
  ///        to the comparison function \c comp
  template<typename T, typename Compare>
  const typename std::vector<T>::iterator best_element(std::vector<T>& x,
                                                       Compare comp)
  {
    return std::max_element(x.begin(), x.end(), comp);
  }

  /// \copybrief best_element(std::vector<T>&, Compare)
  template<typename T, typename Compare>
  const typename std::vector<T>::const_iterator
  best_element(const std::vector<T>& x, Compare comp)
  {
    return std::max_element(x.begin(), x.end(), comp);
  }

  /// \brief Returns an iterator to the maximal element according
  ///        to the standard comparison function
  template<typename T>
  const typename std::vector<T>::iterator max_element(std::vector<T>& x)
  {
    return std::max_element(x.begin(), x.end());
  }

  /// \copybrief max_element(std::vector<T>&)
  template<typename T>
  const typename std::vector<T>::const_iterator
  max_element(const std::vector<T>& x)
  {
    return std::max_element(x.begin(), x.end());
  }

  /// \brief Returns an iterator to the minimal element according
  ///        to the standard comparison function
  template<typename T>
  const typename std::vector<T>::iterator min_element(std::vector<T>& x)
  {
    return std::min_element(x.begin(), x.end());
  }

  /// \copybrief min_element(std::vector<T>&)
  template<typename T>
  const typename std::vector<T>::const_iterator
  min_element(const std::vector<T>& x)
  {
    return std::min_element(x.begin(), x.end());
  }

  /// \brief Returns \f$\|x\|_\infty\f$
  template<typename T>
  const T max_norm(const std::vector<T>& x)
  {
    const typename std::vector<T>::const_iterator i =
      std::max_element(x.begin(), x.end(), abs_less<T>());

    return std::abs(*i);
  }

  template<typename T>
  const long count_nonzero(const std::vector<T>& x)
  {
    long result = 0;

    for (int n = 0, i = x.size(); n < i; n++) {
      if (x[n] != 0.0) {
        result++;
      }
    }
    return result;
  }

  /// \brief Performs \f$z:=x\odot y\f$
  template<typename T>
  void pointwiseProduct(const std::vector<T>& x, const std::vector<T>& y,
                        std::vector<T>& z)
  {
    typename std::vector<T>::size_type min_size =
      std::min(z.size(), std::min(x.size(), y.size()));

    typename std::vector<T>::const_iterator x_iter = x.begin();
    typename std::vector<T>::const_iterator y_iter = y.begin();
    typename std::vector<T>::iterator result_iter = z.begin();

    while (min_size-- > 0)
      *result_iter++ = *x_iter++ * (*y_iter++);
  }
  
  template<typename T>
  std::vector<T> concat(const std::vector<T> &vector1, const std::vector<T> &vector2)
  {
    std::vector<T> result(vector1);

    result.insert(result.end(), vector2.begin(), vector2.end());
    return result;
  }

  template<typename T>
  int indexOf(const std::vector<T> &vector, const T &value) {

    for (int n = 0, i = vector.size(); n < i; n++) {
      if (vector[n] == value) {
    	return n;
      }
    }
    return -1;
  }

  template<typename T>
  bool contains(const std::vector<T> &vector, const T &value) {
    return indexOf(vector, value) > -1;
  }

  template<typename T>
  std::vector<T> cutOut(const std::vector<T> &vector, const int start, const int end) {
    int length = end - start;

    std::vector<T> result(length);

    for (int n = 0; n < length; n++) {
      result[n] = vector[n+start];
    }
    return result;
  }

  template<typename T>
  std::vector<T> cutOut(const std::vector<T> &vector, const int start) {
    return cutOut(vector, start, vector.size()-1);
    }

  template<typename T>
  T componentProduct(const std::vector<T> &selector, const std::vector<T> &vector) {
  	T result = 1;

  	for (int n = 0, i = selector.size(); n < i; n++) {
  	  result *= vector[selector[n]];
  	}
  	return result;
  }

  template<typename T>
  T componentProduct(const std::vector<T> &vector) {
  	T result = 1;

   	for (int n = 0, i = vector.size(); n < i; n++) {
   	  result *= vector[n];
   	}
   	return result;
  }

  template<typename T>
  T componentSum(const std::vector<T> &vector) {
	T result = 0;

	for (unsigned int n = 0, i = vector.size(); n < i; n++) {
	  result += vector[n];
	}
	return result;
  }

  template<typename T>
  std::vector<T> compound(const std::vector<T> &main, const std::vector<T> &part) {
	unsigned int length = part.size();

    std::vector<T> result(length);

    for (unsigned int n = 0; n < length; n++) {
  	  result[n] = indexOf(main, part[n]);
  	}
  	return result;
  }

  template<typename T>
  struct StandardVectorSpaceTraits< std::vector<T> > {
    typedef std::vector<T> Vectors;
    typedef T Scalars;
    
    void update(const Scalars alpha, const Vectors& x, Vectors& y) const
    {
      add(alpha, x, y);
    }
    
    void scale(const Scalars alpha, Vectors& x) const
    {
      TensorCalculus::scale(x, alpha);
    }
    
    Vectors zero() const {
      return Vectors();
    }
  };

  template<typename T>
  struct StandardVectorSpaceTraits< std::vector< std::vector<T> > > {
    typedef std::vector< std::vector<T> > Vectors;
    typedef T Scalars;
    
    void update(const Scalars alpha, const Vectors& x, Vectors& y) const
    {
      StandardVectorSpaceTraits< std::vector<T> > traits;
      y.resize(x.size());
      for (unsigned int i = 0, n = x.size(); i < n; ++i) {
        traits.update(alpha, x[i], y[i]);
      }
    }
    
    void scale(const Scalars alpha, Vectors& x) const
    {
      StandardVectorSpaceTraits< std::vector<T> > traits;
      for (unsigned int i = 0, n = x.size(); i < n; ++i) {
        traits.scale(alpha, x[i]);
      }
    }
    
    Vectors zero() const {
      return Vectors();
    }
  };

  template<typename T>
  struct StandardInnerProductSpaceTraits< std::vector<T> > : public StandardVectorSpaceTraits< std::vector<T> > {
    using StandardVectorSpaceTraits< std::vector<T> >::update;
    using StandardVectorSpaceTraits< std::vector<T> >::scale;
    using StandardVectorSpaceTraits< std::vector<T> >::zero;
    
    typedef std::vector<T> Vectors;
    typedef T Scalars;
    
    Scalars innerProduct(const Vectors& x, const Vectors& y) const
    {
      return TensorCalculus::innerProduct(x, y);
    }
  };

  template<typename T>
  struct StandardInnerProductSpaceTraits< std::vector< std::vector<T> > > : public StandardVectorSpaceTraits< std::vector< std::vector<T> > > {
    using StandardVectorSpaceTraits< std::vector< std::vector<T> > >::update;
    using StandardVectorSpaceTraits< std::vector< std::vector<T> > >::scale;
    using StandardVectorSpaceTraits< std::vector< std::vector<T> > >::zero;
    
    typedef std::vector< std::vector<T> > Vectors;
    typedef T Scalars;
    
    Scalars innerProduct(const Vectors& x, const Vectors& y) const
    {
      Scalars result = 0;
      const unsigned int size = std::min(x.size(), y.size());
      StandardInnerProductSpaceTraits< std::vector<T> > traits;
      for (int i = 0; i < size; ++i) {
        result += traits.innerProduct(x[i], y[i]);
      }
      return result;
    }
  };
  
  template<typename T>
  struct StandardMonoidTraits< std::vector<T> > {
    typedef std::vector<T> MonoidType;
    
    MonoidType identityElement() const { return std::vector<T>(); }
    
    void update(const MonoidType& x, MonoidType& y) const
    {
      if (y.size() < x.size()) {
        y.resize(x.size(), 1);
      }
      for (int i = 0; i < x.size(); ++i) {
        y[i] *= x[i];
      }
    }
  };

  namespace VectorOperators {

    /// \brief Returns a vector containing \f$\frac1\alpha x\f$
    template<typename T>
    std::vector<T> operator/ (const std::vector<T>& x, const T alpha)
    {
      std::vector<T> result(x);
      scale(result, 1/alpha);
      return result;
    }

  } // namespace VectorOperators

  /// \brief Performs \f$x:=\frac1{\|x\|_2} x\f$ and returns
  ///        the old \f$\|x\|_2\f$
  template<typename T>
  const T l2_normalize(std::vector<T>& x)
  {
    const T norm = l2_norm(x);
    using namespace VectorOperators;
    x /= norm;
    return norm;
  }

  /// \brief Returns whether \c x and \c y are equal or not
  /// \remark \f$x=y\ \Leftrightarrow\ \|x-y\|_\infty<\varepsilon\f$
  template<typename T>
  bool equals(const std::vector<T>& x, const std::vector<T>& y,
              T eps = Constants<T>::epsilon)
  {
    typename std::vector<T>::size_type min_size = std::min(x.size(), y.size());
    typename std::vector<T>::const_iterator x_iter = x.begin();
    typename std::vector<T>::const_iterator y_iter = y.begin();

    while(min_size-- > 0) {
      if (std::abs(*x_iter++ - *y_iter++) >= eps) return false;
    }

    return true;
  }

  /// \brief Returns \f$\alpha\f$ if \f$x == \alpha*y\f$, 0 if there is no such \f$\alpha\f$
  /// \remark This method pays respect to the machine accuracy

  template<typename T>
  T determineFactor(const T* x, const T* y, const int n)
  {
	T factor = 0;

	int i;

	for (i = 0; i < n; i++) {
		if (y[i] != 0) {
			factor = x[i]/y[i];
			break;
		}
	}
	for (i++; i < n; i++) {
	  if (std::abs(factor - x[i] / y[i]) > 1e-14*std::abs(factor)) {
		return 0;
	  }
	}
	return factor;
  }

  template<typename T>
  void delete_entries(std::vector<T> &x, const std::vector<int> positions) {
	  std::vector<int> min(positions);
	  std::sort(min.begin(), min.end());

	  int TSIZE = sizeof(T);

	  int n = x.size();

	  for (int k = positions.size()-1; k > -1; k--) {
		  /* remove column j and row j */
		  const int j = min[k];

  		  std::memmove(&x[j], &x[j+1] , TSIZE*(n-j-1));
  		  n--;
  	  }
  	  x.resize(n);
  }

  template<typename T>
  void delete_entries_index(std::vector<T> &x, const std::vector<int> positions) {


	  int size = positions.size();
	  std::vector<int> min(size);

	  for (int k = 0; k < size; k++) {
		  min[k] = indexOf(x, positions[k]);
	  }

	  std::sort(min.begin(), min.end());

	  int TSIZE = sizeof(T);

	  int n = x.size();

	  for (int k = positions.size()-1; k > -1; k--) {
		  /* remove column j and row j */
		  const int j = min[k];

  		  std::memmove(&x[j], &x[j+1] , TSIZE*(n-j-1));
  		  n--;
  	  }
  	  x.resize(n);
  }

  // We don't want to pollute the TensorCalculus namespace with
  // these overloaded operators. Sometimes we don't want to be
  // able to sum up two vectors
  // Usage:
  //   void need_these_operators(const vector<double>& x, ...) {
  //     using namespace VectorOperators;
  //     x *= alpha;
  //     // ...
  //   }
  // Don't use "using namespace ..." in a header file on global scope
  // e.g. outside of functions.
  namespace VectorOperators {

    /// \brief Returns a vector containing \f$x+y\f$
    template<typename T>
    std::vector<T> operator+ (const std::vector<T>& x, const std::vector<T>& y)
    {
      std::vector<T> result(x);
      add(1.0, y, result);
      return result;
    }

    template<typename T>
    std::vector<T> operator+ (const std::vector<T>& x, const T * y)
    {
      std::vector<T> result(x);
      add(1.0, y, result);
      return result;
    }

    template<typename T>
    std::vector<T> operator+ (const T * y, const std::vector<T>& x)
    {
      std::vector<T> result(x);
      add(1.0, y, result);
      return result;
    }

    /// \brief Performs \f$x:=x+y\f$
    template<typename T>
    std::vector<T>& operator+= (std::vector<T>& x, const std::vector<T>& y)
    {
      add(1.0, y, x);
      return x;
    }

    /// \brief Returns a vector containing \f$x-y\f$
    template<typename T>
    std::vector<T> operator- (const std::vector<T>& x, const std::vector<T>& y)
    {
      std::vector<T> result(x);
      add(-1.0, y, result);
      return result;
    }

    template<typename T>
    std::vector<T> operator- (const std::vector<T>& x, const T * y)
    {
      std::vector<T> result(x);
      add(-1.0, y, result);
      return result;
    }

    template<typename T>
    bool operator< (const std::vector<T>& x, const std::vector<T>& y)
    {
      bool result = false;

      for (int n = 0, i = x.size(); n < i; n++)
      {
        if (x[n] > y[n]) {
          return false;
        } else if (x[n] < y[n]) {
          result = true;
        }
      }
      return result;
    }

    template<typename T>
    bool operator> (const std::vector<T>& x, const std::vector<T>& y)
    {
      return y < x;
    }

    template<typename T>
    std::vector<T> operator- (const T * y, const std::vector<T>& x)
    {
      std::vector<T> result(x);
      add(-1.0, y, result);
      return result;
    }

    /// \brief Performs \f$x:=x-y\f$
    template<typename T>
    std::vector<T>& operator-= (std::vector<T>& x, const std::vector<T>& y)
    {
      add(-1.0, y, x);
      return x;
    }

    /// \brief Pushes the contents of the vector \c x to \c out
    template<typename T>
    std::ostream& operator<< (std::ostream& out, const std::vector<T>& x)
    {
      out << "[ ";
      std::copy(x.begin(), x.end(), std::ostream_iterator<T>(out, " "));
      out << "]";
      return out;
    }

    /// \brief Pushes the contents of the vector \c x to \c out
    template<typename T>
    std::ostream& operator<< (std::ostream& out, const std::vector< std::vector<T> >& x)
    {
      out << "[ ";
      for (int n = 0, i = x.size(); n < i; n++) {
    	out << x[n];
      }
      out << "]";
      return out;
    }

    /// \brief Returns a vector containing \f$\alpha x\f$
    template<typename T>
    std::vector<T> operator* (const std::vector<T>& x, const T alpha)
    {
      std::vector<T> result(x);
      scale(result, alpha);
      return result;
    }

    /// \copybrief operator*(const std::vector<T>&, const T)
    template<typename T>
    std::vector<T> operator* (const T alpha, const std::vector<T>& x)
    {
      std::vector<T> result(x);
      scale(result, alpha);
      return result;
    }

    /// \brief Performs \f$x:=\alpha x\f$
    template<typename T>
    std::vector<T>& operator*= (std::vector<T>& x, const T alpha)
    {
      scale(x, alpha);
      return x;
    }

    /// \brief Performs \f$x:=\frac1\alpha x\f$
    template<typename T>
    std::vector<T>& operator/= (std::vector<T>& x, const T alpha)
    {
      scale(x, 1/alpha);
      return x;
    }

    /// \brief Returns \f$x\cdot y\f$
    template<typename T>
    const T operator* (const std::vector<T>& x, const std::vector<T>& y)
    {
      return innerProduct(x, y);
    }

    /// \copybrief equals(const std::vector<T>&, const std::vector<T>&, T)
    template<typename T>
    bool operator== (const std::vector<T>& x, const std::vector<T>& y)
    {
      return equals(x, y);
    }

  } // namespace VectorOperators

} // namespace TensorCalculus

#endif // __VECTOR_OPERATORS_HPP
