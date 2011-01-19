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
 * \file TensorBasis.hpp
 * \brief Declares the TensorBasis class
 * \author Philipp Wähnert
 */

#ifndef __TENSOR_BASIS_HPP
#define __TENSOR_BASIS_HPP

#include <vector>          // std::vector<T>
#include <algorithm>       // std::swap
#include <stdexcept>       // std::out_of_range
#include "Vector/WeakVector.hpp"  // WeakVector<T>
//#include <BlasInterface.hpp>

namespace TensorCalculus {

  /*!
   * \class TensorBasis
   * \brief Class template holding the basis vectors
   *
   * A given \f$x=\sum_{j=1}^r\bigotimes_{\mu=1}^dx_{j\mu}\in\bigotimes_{\mu=1}^d\mathbf R^{m_\mu}\f$
   * with \f$x_{j\mu}\in\mathbf R^{m_\mu}\f$ can be mapped to a
   * \f$y=\sum_{j=1}^r\bigotimes_{\mu=1}^dy_{j\mu}\in\bigotimes_{\mu=1}^d\mathbf R^{n_\mu}\f$ by
   * \f[
   *   y_{j\mu}=\sum_{k=1}^{m_\mu}(x_{j\mu})_kb_{k\mu}
   * \f]
   * for all \f$j\in\{1,\ldots,r\}\f$ and \f$\mu\in\{1,\ldots,d\}\f$
   * and \f$b_{k\mu}\in\mathbf R^{n_\mu}\f$
   *
   * This class represents the \f$b_{i\mu}\in\mathbf R^{n_\mu}\f$ with
   * \f$\mu=1,\ldots,d\f$ and \f$k=1,\ldots,m_\mu\f$ to provide a similar
   * transformation as shown above.
   *
   * \remark The names in this class are chosen according to the
   *         above notation.
   *          - \f$d\f$ denotes the number of factors of
   *            the tensor product space
   *          - \f$n_\mu\f$ denotes the dimension of the
   *            \f$\mu\f$-th factor of the image tensor product space
   *          - \f$m_\mu\f$ denotes the dimension of the \f$\mu\f$-th
   *            factor of the preimage tensor product space
   *
   * \tparam T data type of its elements
   */
  template<typename T>
  class TensorBasis
  {
  public:
    TensorBasis(int d, const std::vector<int>& n, const std::vector<int>& m)
      : m_d(d), m_n(n), m_m(m)
    {
      resize_data();
    }

    TensorBasis(const TensorBasis& b)
      : m_d(b.m_d), m_n(b.m_n), m_m(b.m_m), m_data(b.m_data)
    {
      update_weak_vectors();
    }

    void swap(TensorBasis& b)
    {
      using std::swap;
      swap(m_d, b.m_d);
      swap(m_n, b.m_n);
      swap(m_m, b.m_m);
      swap(m_data, b.m_data);
      swap(m_basis_vectors, b.m_basis_vectors); // correctness due to standard!
    }

    TensorBasis& operator= (TensorBasis b)
    {
      swap(b); // doing the swap trick!
      return *this;
    }

    /// Returns the vector \f$(b_{1\mu},\ldots,b_{m_\mu\mu})\f$
    std::vector<T>& getBasis(int mu) {
#ifdef CHECK_RANGES_ON
      return m_data.at(mu);
#else
      return m_data[mu];
#endif
    }

    /// \copybrief getBasis
    const std::vector<T>& getBasis(int mu) const
    {
#ifdef CHECK_RANGES_ON
      return m_data.at(mu);
#else
      return m_data[mu];
#endif
    }

    /// Returns the vector \f$b_{k\mu}\f$
    WeakVector<const T> getBasisVector(int mu, int k) const
    {
#ifdef CHECK_RANGES_ON
      return WeakVector<const T>(m_basis_vectors.at(mu).at(k), m_n[mu]);
#else
      return WeakVector<const T>(m_basis_vectors[mu][k], m_n[mu]);
#endif
    }

    /// \copybrief getBasisVector
    WeakVector<T> getBasisVector(int mu, int k)
    {
#ifdef CHECK_RANGES_ON
      return WeakVector<T>(m_basis_vectors.at(mu).at(k), m_n[mu]);
#else
      return WeakVector<T>(m_basis_vectors[mu][k], m_n[mu]);
#endif
    }

    /// Returns \f$(b_{k\mu})_i\f$
    const T& getBasisVectorEntry(int mu, int k, int i) const
    {
#ifdef CHECK_RANGES_ON
      return m_data.at(mu).at(k*m_n.at(mu)+i);
#else
      return m_data[mu][k*m_n[mu]+i];
#endif
    }

    /// \copybrief getBasisVectorEntry
    T& getBasisVectorEntry(int mu, int k, int i)
    {
#ifdef CHECK_RANGES_ON
      return m_data.at(mu).at(k*m_n.at(mu)+i);
#else
      return m_data[mu][k*m_n[mu]+i];
#endif
    }

    std::vector<T>& operator() (int mu)
    {
      return getBasis(mu);
    }

    const std::vector<T>& operator() (int mu) const
    {
      return getBasis(mu);
    }

    WeakVector<T> operator() (int mu, int k)
    {
      return getBasisVector(mu, k);
    }

    WeakVector<const T> operator() (int mu, int k) const
    {
      return getBasisVector(mu, k);
    }

    T& operator() (int mu, int k, int i)
    {
      return getBasisVectorEntry(mu, k, i);
    }

    const T& operator() (int mu, int k, int i) const
    {
      return getBasisVectorEntry(mu, k, i);
    }

    void resize(int d, const std::vector<int>& n, const std::vector<int>& m)
    {
#ifdef CHECK_RANGES_ON
      if ((n.size() != d) || (m.size() != d))
      {
        throw std::out_of_range("TensorBasis::resize(int,"
                                " const std::vector<int>&,"
                                " const std::vector<int>&)"
                                " : Dimensions mismatch");
      }
#endif
      m_d = d;
      m_n = n;
      m_m = m;
      resize_data();
    }

    /// Returns \f$d\f$
    int getSpatDirs() const { return m_d; }

    /// \copybrief getSpatDirs
    int get_d() const { return m_d; }

    /// Returns the vector \f$n=(n_1,\ldots,n_d)\f$
    const std::vector<int>& getDimsOfSubspaces() const { return m_n; }

    /// \copybrief getDimsOfSubspaces
    const std::vector<int>& get_n() const { return m_n; }

    /// \brief Returns \f$n_\mu\f$
    const int get_n(int mu) const { return m_n.at(mu); }

    /// Returns the vector \f$m=(m_1,\ldots,m_d)\f$
    const std::vector<int>& getDimsOfBases() const { return m_m; }

    /// \copybrief getDimsOfBases
    const std::vector<int>& get_m() const { return m_m; }

    /// \brief Return \f$m_\mu\f$
    const int get_m(int mu) const { return m_m.at(mu); }

  private:

    void resize_data()
    {
      m_data.clear();
      m_data.resize(m_d);
      for (int mu = 0; mu < m_d; mu++) {
        m_data[mu].assign(m_n[mu]*m_m[mu], 0);
      }
      update_weak_vectors();
    }

    void update_weak_vectors()
    {
      m_basis_vectors.clear();
      m_basis_vectors.resize(m_d);

      for (int mu = 0; mu < m_d; ++mu)
      {
        std::vector<T*>& temp = m_basis_vectors[mu];
        temp.reserve(m_m[mu]);

        int offset = 0;
        for (int k = 0; k < m_m[mu]; k++) {
          temp.push_back(&m_data[mu][offset]);
          offset += m_n[mu];
        }
      }
    }

    int m_d;
    std::vector<int> m_n;
    std::vector<int> m_m;

    std::vector<std::vector<T> > m_data;
    std::vector<std::vector<T*> > m_basis_vectors;

  };

} // namespace TensorCalculus

#endif // __TENSOR_BASIS_HPP
