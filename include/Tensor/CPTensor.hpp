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
 * \file CPTensor.hpp
 * \brief Declares the CPTensor class
 * \author Philipp Wähnert
 */

#ifndef __CP_TENSOR_HPP
#define __CP_TENSOR_HPP

#include <vector>          // std::vector<T>
#include <algorithm>       // std::swap
#include <stdexcept>       // std::out_of_range
// #include "Utilities.hpp"
#include <cmath>           // std::abs
#include "Constants.hpp"
#include "Vector/WeakVector.hpp"  // WeakVector<T>

// #define RANGE_CHECKS_ON

namespace TensorCalculus {

  /*!
   * \class CPTensor
   * \brief Class template to represent a tensor by direct Kronecker sums
   *
   * According to the theory a tensor
   * \f$x\in\bigotimes_{\mu=1}^d\mathbf R^{n_\mu}\f$
   * can always represented in the following form
   * \f[
   *   x = \sum_{j=1}^r\bigotimes_{\mu=1}^dx_{j\mu}
   * \f]
   * with \f$x_{j\mu}\in\mathbf R^{n_\mu}\f$. The entries of the tensor are then
   * \f[
   *   x_{\underline i} = \sum_{j=1}^r\prod_{\mu=1}^d(x_{j\mu})_{i_\mu}
   * \f]
   * with \f$\underline i\in\times_{\mu=1}^d\{1,\ldots,n_\mu\}\f$.
   *
   * \remark The names in this class are normally choosen according to the
   *         above notation.
   *          - \f$r\f$ denotes the rang of separation
   *          - \f$d\f$ denotes the number of factors of the tensor product
   *            space
   *          - \f$n_\mu\f$ denotes the dimension of \f$\mu\f$-th factor
   *            of the tensor product space
   *
   * \tparam T data type of its elements, i.e. the \f$x_{j\mu}\f$
   *           are elements of \c T<br>
   *           Must be standard constructable, assignable, if \c x and
   *           \c y are of type \c T then
   *           <tt>x+y</tt> and <tt>x*y</tt> must be convertible to \c T
   *
   * \todo Replace in <tt>std::vector< WeakVector<T> ></tt> the \c vector<T>
   *       with an appropriate container because \c WeakVector<T> isn't standard
   *       constructible nor assignable<br>
   *       <b>Update:</b> Error seems to be fixed. Now m_xjmu is only a vector
   *       of vectors of pointers and the result type of
   *       \c getVectorOfRepresentants is now a copy of a WeakVector filled
   *       right there with the right pointers and lengths.
   */
  template <typename T>
  class CPTensor
  {
  public:

    /// Standard constructor initializes the tensor with \f$r=0\f$ and \f$d=0\f$
    CPTensor()
      : m_r(0), m_d(0), m_n()
    {
      resize_data();
      refresh_weak_vectors();
    }

    /// \brief Constructor initializing the tensor with the
    ///        given \f$r\f$, \f$d\f$,
    ///        constant \f$n_\mu=n\f$ and zeros in the entries
    CPTensor(int r, int d, int n)
      : m_r(r), m_d(d), m_n(d, n)
    {
      resize_data();
      refresh_weak_vectors();
    }

    /// \brief Constructor initializing the tensor with the
    ///        given \f$r\f$, \f$d\f$,
    ///        \f$n=(n_1,\ldots,n_d)\f$ and zeros in the entries
    CPTensor(int r, int d, const std::vector<int>& n)
      : m_r(r), m_d(d), m_n(n)
    {
      resize_data();
      refresh_weak_vectors();
    }

    CPTensor(std::istream& stream) {
      stream.ignore(13, '='); // "Tensor in d ="
      stream >> m_d;
      m_n.resize(m_d);
      for (int mu = 0; mu < m_d; ++mu) {
        stream.ignore(30, '='); // "n[*] ="
        stream >> m_n[mu];
      }
      stream.ignore(8, '='); // "Rank k ="
      stream >> m_r;
      m_data.reserve(m_d);
      for (int mu = 0; mu < m_d; ++mu) {
        m_data.push_back(std::vector<T>(m_n[mu] * m_r));
      }
      for (int j = 0; j < m_r; ++j) {
        for (int mu = 0; mu < m_d; ++mu) {
          const int n_mu = m_n[mu];
          for (int i_mu = 0; i_mu < n_mu; ++i_mu) {
            stream >> m_data[mu][n_mu * j + i_mu];
          }
        }
      }
    }

    /// Copy constructor
    CPTensor(const CPTensor& x)
      : m_r(x.m_r), m_d(x.m_d), m_n(x.m_n), m_data(x.m_data)
    {
      refresh_weak_vectors();
    }

    /// Swapping its contents with these in \c x
    void swap(CPTensor& x)
    {
      using std::swap;
      swap(m_r, x.m_r);
      swap(m_d, x.m_d);
      swap(m_n, x.m_n);
      swap(m_data, x.m_data);
      // refresh_weak_vectors();
      // A refresh isn't necessary because after a swap of
      // a vector all pointers remains valid.
      // Thus we only have to swap the WeakVectors
      swap(m_xjmu, x.m_xjmu);
    }

    /// Assignment operator implemented by the copy and swap
    /// technique to get self assignment safety
    CPTensor& operator= (CPTensor rhs) // no reference!
    {
      swap(rhs); // exception and self assignment safe swap.
      return *this;
    }

    /// Returns the vector \f$(x_{1\mu},\ldots,x_{r\mu})\f$
    std::vector<T>& getVectorOfDimension(int mu)
    {
#ifdef RANGE_CHECKS_ON
      range_check(mu);
#endif
      return m_data[mu];
    }

    /// \copybrief getVectorOfDimension(int)
    const std::vector<T>& getVectorOfDimension(int mu) const
    {
#ifdef RANGE_CHECKS_ON
      range_check(mu);
#endif
      return m_data[mu];
    }

    /// Returns a weak vector to \f$x_{j\mu}\f$
    WeakVector<T> getVectorOfRepresentants(int j, int mu)
    {
#ifdef RANGE_CHECKS_ON
      range_check(j, mu);
#endif
      return WeakVector<T>(m_xjmu[mu][j], m_n[mu]); // [mu*m_r + j];
    }

    /// \copybrief getVectorOfRepresentants(int, int)
    WeakVector<const T> getVectorOfRepresentants(int j, int mu) const
    {
#ifdef RANGE_CHECKS_ON
      range_check(j, mu);
#endif
      return WeakVector<const T>(m_xjmu[mu][j], m_n[mu]); // [mu*m_r + j];
    }

    /// Returns the element \f$(x_{j\mu})_i\f$
    T& getVector(int j, int mu, int i)
    {
#ifdef RANGE_CHECKS_ON
      range_check(j, mu, i);
#endif
      return m_data[mu][j*m_n[mu]+i];
    }

    /// \copybrief getVector(int, int, int)
    const T& getVector(int j, int mu, int i) const
    {
#ifdef RANGE_CHECKS_ON
      range_check(j, mu, i);
#endif
      return  m_data[mu][j*m_n[mu]+i];
    }

    /// Returns the number
    /// \f$\sum_{j=\mathtt{min\_j}+1}^{\mathtt{max\_j-1}}\prod_{\mu=1}^d (x_{j\mu})_{i_\mu}\f$
    const T getTensorEntry(const std::vector<int>& index,
                           int min_j, int max_j) const
    {
#ifdef RANGE_CHECKS_ON
      if ((min_j < 0) || (max_j > m_r) || (max_j < m_r))
      {
        throw std::out_of_range("CPTensor::range_check"
                                " - min_j or max_j out of range");
      }
#endif

      T result = 0;

      for (int j = min_j; j < max_j; j++)
      {
        T product = 1;
        for (int mu = 0; mu < m_d; mu++)
        {
#ifdef RANGE_CHECKS_ON
          product *= getVector(j, mu, index.at(mu));
#else
          product *= getVector(j, mu, index[mu]);
#endif
        }
        result += product;
      }

      return result;
    }

    /// Returns the number \f$\sum_{j=1}^r\prod_{\mu=1}^d (x_{j\mu})_{i_\mu}\f$
    const T getTensorEntry(const std::vector<int>& index) const
    {
      return getTensorEntry(index, 0, m_r);
    }

    /// \copybrief getVectorOfDimension(int)
    std::vector<T>& operator() (int mu)
    {
      return getVectorOfDimension(mu);
    }

    /// \copybrief getVectorOfDimension(int)
    const std::vector<T>& operator() (int mu) const
    {
      return getVectorOfDimension(mu);
    }

    /// \copybrief getVectorOfRepresentants(int, int)
    WeakVector<T> operator() (int j, int mu)
    {
      return getVectorOfRepresentants(j, mu);
    }

    /// \copybrief getVectorOfRepresentants(int, int)
    WeakVector<const T> operator() (int j, int mu) const
    {
      return getVectorOfRepresentants(j, mu);
    }

    /// \copybrief getVector(int, int, int)
    T& operator() (int j, int mu, int i)
    {
      return getVector(j, mu, i);
    }

    /// \copybrief getVector(int, int, int)
    const T& operator() (int j, int mu, int i) const
    {
      return getVector(j, mu, i);
    }

    /// \copybrief getTensorEntry
    const T operator() (const std::vector<int>& index) const
    {
      return getTensorEntry(index, 0, m_r);
    }

    /// \copybrief getTensorEntry
    const T operator() (const std::vector<int>& index,
                        int min_j, int max_j) const
    {
      return getTensorEntry(index, min_j, max_j);
    }

    void range_check(int mu) const {
      if ((mu < 0) || (mu >= m_d))
      {
        throw std::out_of_range("CPTensor::range_check"
                                " - mu out of range");
      }
    }

    void range_check(int j, int mu) const
    {
      range_check(mu);
      if ((j < 0) || (j >= m_r))
      {
        throw std::out_of_range("CPTensor::range_check"
                                " - j out of range");
      }
    }

    void range_check(int j, int mu, int i) const
    {
      range_check(j, mu);
      if ((i < 0) || (i >= m_n[mu]))
      {
        throw std::out_of_range("CPTensor::range_check"
                                "- i out of range");
      }
    }

    /// Resizes the tensor to
    /// \f$r\f$, \f$d\f$ and \f$n=(n_1,\ldots,n_d)\f$ and set
    /// the entries to zero
    void resize(int r, int d, const std::vector<int>& n)
    {
      m_r = r;
      m_d = d;
      m_n = n;

      resize_data();
      refresh_weak_vectors();
    }

    /// Assigns every entry of the tensor to \c value
    void assign(T value)
    {
      for (typename std::vector<std::vector<T> >::iterator i = m_data.begin();
           i != m_data.end(); ++i)
      {
        std::fill((*i).begin(), (*i).end(), value);
      }
    }

    /// Sets the tensor to zero
    void setNull()
    {
      assign(T());
    }

    /// \brief Assigns the tensor with values from the generator \c gen
    /// \remark The order of assignment is
    ///         \f$x_{11}, x_{21}, \ldots, x_{r1}, x_{12},\ldots\f$
    ///
    /// The generator must be something permitting expressions
    /// like <tt>T x = gen()</tt>
    template<typename Generator>
    void assign(Generator gen)
    {
      for (typename std::vector<std::vector<T> >::iterator i = m_data.begin();
           i != m_data.end(); ++i)
      {
        std::generate((*i).begin(), (*i).end(), gen);
      }
    }

    /// Returns \f$r\f$
    const int getSeparationRang() const { return m_r; }

    /// \copybrief getSeparationRang
    const int get_r() const { return m_r; }

    /// Returns \f$d\f$
    const int getSpatialDirections() const { return m_d; }

    /// \copybrief getSpatialDirections
    const int get_d() const { return m_d; }

    /// Returns \f$n_\mu\f$
    const int getDimensionOfSpatDir(int mu) const { return m_n[mu]; }

    /// \copybrief getDimensionOfSpatDir
    const int get_n(int mu) const { return m_n[mu]; }

    //const std::vector<int> getDimensions() { return dims_of_spat_dirs; }
    /// Returns the vector \f$n=(n_1,\ldots,n_d)\f$
    const std::vector<int>& getDimensions() const { return m_n; }

    /// \copybrief getDimensions
    const std::vector<int>& get_n() const { return m_n; }

    void write(std::ostream& stream) const {  
      stream << "Tensor in d = " << m_d << '\n';
      for(int mu = 0; mu < m_d; mu++) {
        stream << "n[" << mu << "] = " << m_n[mu] << '\n';
      }
      stream << "Rank k = " << m_r << '\n';
  
      stream.setf(std::ios::scientific, std::ios::floatfield);
      stream.precision(20);
      
      for(int j = 0; j < m_r; j++) {
        for(int mu = 0; mu < m_d; mu++) {
          const int n_mu = m_n[mu];
          for(int i = 0; i < n_mu; i++) {
            stream << m_data[mu][j * n_mu + i] << '\n';
          }
        }
      }
    }

  private:

    /// Resizes the internal data array to the given lengths
    void resize_data()
    {
      m_data.resize(m_d);
      for (int mu = 0; mu < m_d; ++mu)
      {
        m_data[mu].assign(m_n[mu]*m_r, 0);
      }
    }

    /// Refreshs the weak vectors pointing to the diverse parts
    /// of the internal data
    void refresh_weak_vectors()
    {
      m_xjmu.clear();
      m_xjmu.resize(m_d);

      for (int mu = 0; mu < m_d; ++mu)
      {
        std::vector<T*>& temp = m_xjmu[mu];
        temp.reserve(m_r);

        for (int j = 0; j < m_r; ++j)
        {
          temp.push_back(&m_data[mu][j*m_n[mu]]);
        }
      }

    }

    /// The rang of separation \f$r\f$
    int m_r;

    /// The number of factors of the tensor product space \f$d\f$
    int m_d;

    /// The element \c m_n[mu-1] represents \f$n_\mu\f$ the dimension
    /// of the \f$\mu\f$-th factor of the tensor product
    std::vector<int> m_n;

    /// \brief The element \c m_data[mu-1][j*n_nmu[mu-1]+i-1] represents
    ///        \f$(x_{j\mu})_i\f$ or respectively \c m_data[mu] represents
    ///        the vector \f$(x_{1\mu},\ldots,x_{r\mu})\f$
    std::vector< std::vector<T> > m_data;

    /// The element \c m_xjmu[mu][j] points to the first entry of \f$x_{j\mu}\f$
    std::vector< std::vector<T*> > m_xjmu;

  };

} // namespace TensorCalculus

#endif // __CP_TENSOR_HPP
