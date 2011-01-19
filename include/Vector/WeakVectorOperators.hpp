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

/*!
 * \file WeakVectorOperators.hpp
 * \brief Declares the usual operators for the WeakVector class
 * \author Philipp Wähnert
 */

#ifndef __WEAK_VECTOR_OPERATORS_HPP
#define __WEAK_VECTOR_OPERATORS_HPP

#include <vector>               // std::vector<T>
#include <algorithm>            // std::min_element, std::max_element, ...
#include "Vector/VectorOperators.hpp"  // operator+(const std::vector<T>&, ...
#include "Vector/WeakVector.hpp"       // WeakVector<T>
#include "BlasInterface.hpp"    // Blas<T>
#include "Utilities/Functional.hpp"       // abs_less<T>
#include "Utilities/Utilities.hpp"        // remove_const<T>

namespace TensorCalculus {

  /// \brief Returns a \c std::vector<T> containing a copy of the data
  /// \remark This method uses the BLAS interface to copy
  template<typename T>
  const std::vector<typename remove_const<T>::type>
  getVectorCopy(const WeakVector<T>& x, bool useBLAS = true)
  {
    if (useBLAS)
    {
      std::vector<typename remove_const<T>::type> result(x.size());

      Blas<typename remove_const<T>::type>::copy(x.size(), &x[0],
                                                 1, &result[0], 1);
      return result;
    } else {
      return std::vector<typename remove_const<T>::type>(x.begin(), x.end());
    }
  }

  /// \brief Performs \f$y:=\alpha x+y\f$
  template<typename S, typename T, typename U>
  void add(S alpha, const WeakVector<T>& x, const WeakVector<U>& y)
  {
    const typename WeakVector<U>::size_type min_size =
      std::min(x.size(), y.size());

    Blas<U>::axpy(min_size, alpha, &x[0], 1, &y[0], 1);
  }

  /// \copybrief TensorCalculus::add(S, const WeakVector<T>&, WeakVector<U>&)
  template<typename S, typename T>
  void add(S alpha, const std::vector<S>& x, const WeakVector<T>& y)
  {
    const typename WeakVector<T>::size_type min_size =
      std::min(x.size(), y.size());

    Blas<T>::axpy(min_size, alpha, &x[0], 1, &y[0], 1);
  }

  /// \copybrief TensorCalculus::add(S, const WeakVector<T>&, WeakVector<U>&)
  template<typename S, typename T>
  void add(S alpha, const WeakVector<T>& x, std::vector<S>& y)
  {
    const typename std::vector<S>::size_type min_size =
      std::min(x.size(), y.size());

    Blas<S>::axpy(min_size, alpha, &x[0], 1, &y[0], 1);
  }

  /// \brief Performs \f$x:=\alpha x\f$
  template<typename T>
  void scale(WeakVector<T>& x, const T alpha)
  {
    Blas<T>::scal(x.size(), alpha, &x[0], 1);
  }

  /// \brief Returns \f$x\cdot y\f$
  template<typename S, typename T>
  const T innerProduct(const WeakVector<S>& x, const WeakVector<T>& y)
  {
    const typename WeakVector<S>::size_type min_size =
      std::min(x.size(), y.size());

    return Blas<typename remove_const<S>::type>::dot(min_size, &x[0],
                                                     1, &y[0], 1);
  }

  /// \copybrief TensorCalculus::innerProduct(const WeakVector<S>&, const WeakVector<T>&)
  template<typename S, typename T>
  const T innerProduct(const WeakVector<S>& x, const std::vector<T>& y)
  {
    const typename std::vector<T>::size_type min_size =
      std::min(x.size(), y.size());

    return Blas<T>::dot(min_size, &x[0], 1, &y[0], 1);
  }

  /// \copybrief TensorCalculus::innerProduct(const WeakVector<S>&, const WeakVector<T>&)
  template<typename S, typename T>
  const S innerProduct(const std::vector<S>& x, const WeakVector<T>& y)
  {
    return innerProduct(y, x);
  }

  /// \brief Returns a vector containing \f$\alpha x+\beta y\f$
  template<typename S, typename T, typename U>
  const std::vector<S> add(S alpha, const WeakVector<T>& x,
                           S beta, const WeakVector<U>& y)
  {
    const typename std::vector<S>::size_type min_size =
      std::min(x.size(), y.size());
    std::vector<S> result(x.begin(), x.end());

    result *= alpha;
    Blas<S>::axpy(min_size, beta, &y[0], 1, &result[0], 1);

    return result;
  }

  /// \copybrief TensorCalculus::add(S, const WeakVector<T>&, S, const WeakVector<U>&)
  template<typename S, typename T>
  const std::vector<S> add(S alpha, const std::vector<S>& x,
                           S beta, const WeakVector<T>& y)
  {
    return add(alpha, WeakVector<S>(x), beta, y);
  }

  /// \copybrief TensorCalculus::add(S, const WeakVector<T>&, S, const WeakVector<U>&)
  template<typename S, typename T>
  const std::vector<S> add(S alpha, const WeakVector<T>& x,
                           S beta, const std::vector<S>& y)
  {
    return add(alpha, x, beta, WeakVector<S>(y));
  }

  /// \brief Performs \f$z:=x\odot y\f$
  template<typename S, typename T, typename U>
  void pointwiseProduct(const WeakVector<S>& x, const WeakVector<T>& y,
                        const WeakVector<U>& z)
  {
    typename WeakVector<U>::size_type min_size =
      std::min(z.size(), std::min(x.size(), y.size()));
    typename WeakVector<S>::const_iterator x_iter = x.begin();
    typename WeakVector<T>::const_iterator y_iter = y.begin();
    typename WeakVector<U>::iterator result_iter = z.begin();

    while (min_size-- > 0)
      *result_iter++ = *x_iter++ * (*y_iter++);
  }

  /// \brief Returns \f$\|x\|_2\f$
  template<typename T>
  const T l2_norm(const WeakVector<T>& x)
  {
    return Blas<typename remove_const<T>::type>::nrm2(x.size(), &x[0], 1);
  }

  /// \brief Performs \f$x:=\frac1{\|x\|_2} x\f$ and returns the old \f$\|x\|_2\f$
  template<typename T>
  const T l2_normalize(WeakVector<T>& x)
  {
    const T norm = l2_norm(x);
    x /= norm;
    return norm;
  }

  /// \brief Returns an iterator to the best element according
  ///        to the comparison function \c comp
  template<typename T, typename Compare>
  const typename WeakVector<T>::iterator best_element(const WeakVector<T>& x,
                                                      Compare comp)
  {
    return std::max_element(x.begin(), x.end(), comp);
  }

  /// \brief Returns an iterator to the maximal element according
  ///        to the standard comparison function
  template<typename T>
  const typename WeakVector<T>::iterator max_element(const WeakVector<T>& x)
  {
    return std::max_element(x.begin(), x.end());
  }

  /// \brief Returns an iterator to the minimal element according
  ///        to the standard comparison function
  template<typename T>
  const typename WeakVector<T>::iterator min_element(const WeakVector<T>& x)
  {
    return std::min_element(x.begin(), x.end());
  }

  /// \brief Returns \f$\|x\|_\infty\f$
  template<typename T>
  const T max_norm(const WeakVector<T>& x)
  {
    const typename WeakVector<T>::const_iterator i =
      std::max_element(x.begin(), x.end(), abs_less<T>());

    return abs(*i);
  }

  /// \brief Returns a vector containing \f$x+y\f$
  template<typename S, typename T>
  const std::vector<typename remove_const<S>::type>
  operator+ (const WeakVector<S>& x, const WeakVector<T>& y)
  {
    std::vector<typename remove_const<S>::type> result(x.begin(), x.end());

    add(1.0, y, result);
    return result;
  }

  /// \copybrief TensorCalculus::operator+ (const WeakVector<S>&, const WeakVector<T>&)
  template<typename S, typename T>
  const std::vector<S> operator+ (const std::vector<S>& x,
                                  const WeakVector<T>& y)
  {
    std::vector<S> result(x);
    add(1.0, y, result);
    return result;
  }

  /// \copybrief TensorCalculus::operator+ (const WeakVector<S>&, const WeakVector<T>&)
  template<typename S, typename T>
  const std::vector<T> operator+ (const WeakVector<S>& x,
                                  const std::vector<T>& y)
  {
    std::vector<T> result(x);
    add(1.0, y, result);
    return result;
  }

  /// \brief Performs \f$x:=x+y\f$
  template<typename S, typename T>
  WeakVector<S>& operator+= (WeakVector<S>& x, const WeakVector<T>& y)
  {
    add(1.0, y, x);
    return x;
  }

  /// \copybrief operator+=(WeakVector<S>&, const WeakVector<T>&)
  template<typename S, typename T>
  WeakVector<S>& operator+= (WeakVector<S>& x, const std::vector<T>& y)
  {
    add(1.0, y, x);
    return x;
  }

  /// \copybrief operator+=(WeakVector<S>&, const WeakVector<T>&)
  template<typename S, typename T>
  std::vector<S>& operator+= (std::vector<S>& x, const WeakVector<T>& y)
  {
    add(1.0, y, x);
    return x;
  }

  /// \brief Returns a vector containing \f$x-y\f$
  template<typename S, typename T>
  const std::vector<typename remove_const<S>::type>
  operator- (const WeakVector<S>& x, const WeakVector<T>& y) {
    std::vector<typename remove_const<S>::type> result(x.begin(), x.end());
    add(-1.0, y, result);
    return result;
  }

  /// \copybrief TensorCalculus::operator- (const WeakVector<S>&, const WeakVector<T>&)
  template<typename S, typename T>
  const std::vector<S> operator- (const std::vector<S>& x,
                                  const WeakVector<T>& y)
  {
    std::vector<S> result(x);
    add(-1.0, y, result);
    return result;
  }

  /// \copybrief TensorCalculus::operator- (const WeakVector<S>&, const WeakVector<T>&)
  template<typename S, typename T>
  const std::vector<T> operator- (const WeakVector<S>& x,
                                  const std::vector<T>& y)
  {
    std::vector<T> result(x);
    add(-1.0, y, result);
    return result;
  }

  /// \brief Performs \f$x:=x-y\f$
  template<typename S, typename T>
  WeakVector<S>& operator-= (WeakVector<S>& x, const WeakVector<T>& y)
  {
    add(-1.0, y, x);
    return x;
  }

  /// \copybrief operator-=(WeakVector<S>&, const WeakVector<T>&)
  template<typename S, typename T>
  WeakVector<S>& operator-= (WeakVector<S>& x, const std::vector<T>& y)
  {
    add(-1.0, y, x);
    return x;
  }

  /// \copybrief operator-=(WeakVector<S>&, const WeakVector<T>&)
  template<typename S, typename T>
  std::vector<S>& operator-= (std::vector<S>& x, const WeakVector<T>& y)
  {
    add(-1.0, y, x);
    return x;
  }

  /// \brief Pushes the contents of the WeakVector to the given std::ostream
  template<typename T>
  std::ostream& operator<< (std::ostream& out, const WeakVector<T>& x)
  {
    std::copy(x.begin(), x.end(), std::ostream_iterator<T>(out, " "));
    return out;
  }

  /// \brief Returns a vector containing \f$\alpha x\f$
  template<typename T>
  std::vector<typename remove_const<T>::type>
  operator* (const WeakVector<T>& x, const T alpha)
  {
    std::vector<typename remove_const<T>::type> result(x.begin(), x.end());

    scale(result, alpha); // calls scale(std::vector<T>&, T)
    return result;
  }

  /// \copybrief operator*(const WeakVector<T>&, const T)
  template<typename T>
  std::vector<typename remove_const<T>::type>
  operator* (const T alpha, const WeakVector<T>& x)
  {
    std::vector<typename remove_const<T>::type> result(x.begin(), x.end());

    scale(result, alpha); // calls scale(std::vector<T>&, T)
    return result;
  }

  /// \brief Performs \f$x:=\alpha x\f$
  template<typename T>
  WeakVector<T>& operator*= (WeakVector<T>& x, const T alpha)
  {
    scale(x, alpha); // calls scale(WeakVector<T>&, T)
    return x;
  }

  /// \brief Returns a vector containing \f$\frac1\alpha x\f$
  template<typename T>
  std::vector<typename remove_const<T>::type>
  operator/ (const WeakVector<T>& x, const T alpha)
  {
    std::vector<typename remove_const<T>::type> result = getVectorCopy(x);

    scale(result, 1/alpha); // calls scale(std::vector<T>&, T)
    return result;
  }

  /// \brief Performs \f$x:=\frac1\alpha x\f$
  template<typename T>
  WeakVector<T>& operator/= (WeakVector<T>& x, const T alpha)
  {
    scale(x, 1/alpha); // calls scale(WeakVector<T>&, T)
    return x;
  }

  /// \copybrief innerProduct(const WeakVector<S>&, const WeakVector<T>&)
  template<typename S, typename T>
  const S operator* (const WeakVector<S>& x, const WeakVector<T>& y)
  {
    return innerProduct(x, y);
  }

  /// \copybrief inner_product(const WeakVector<S>& x, const WeakVector<T>& y)
  template<typename S, typename T>
  const S operator* (const std::vector<S>& x, const WeakVector<T>& y)
  {
    return innerProduct(x, y);
  }

  /// \copybrief inner_product(const WeakVector<T>& x, const WeakVector<T>& y)
  template<typename S, typename T>
  const T operator* (const WeakVector<S>& x, const std::vector<T>& y)
  {
    return innerProduct(x, y);
  }

} // namespace TensorCalculus

#endif // __WEAK_VECTOR_OPERATORS_HPP
