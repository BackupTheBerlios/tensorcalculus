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
 * \file CPTensorOperators.hpp
 * \brief Declares the usual operators for the CPTensor class
 * \author Philipp Wähnert
 */

#ifndef __CP_TENSOR_OPERATORS_HPP
#define __CP_TENSOR_OPERATORS_HPP

#include "Vector/WeakVectorOperators.hpp"
#include "Vector/VectorOperators.hpp"
#include "Tensor/CPTensor.hpp"
#include "BlasInterface.hpp"
#include "Utilities/SkippedProducts.hpp"
#include <cmath>              // std::pow
#include <ostream>
//#include <istream>          // not yet needed! operator>> not implemented
#include <iomanip>
#include <stdexcept>          // std::out_of_range
#include <vector>

namespace TensorCalculus {

  /// Returns the <tt>n</tt>-th unit vector of length \c size
  std::vector<double> nth_weight(int size, int n)
  {
    std::vector<double> result(size);
    result[n] = 1;
    return result;
  }

  /// \brief Scales the tensor \c t by the factor \c alpha
  ///        and weights \c weights
  /// \details Assigns \f$\alpha t\f$ to \c t with
  ///   \f$\displaystyle\alpha t = \sum_{j=1}^r\bigotimes_{\mu=1}^d\alpha^{w_\mu}x_{j\mu}\f$
  ///          where \f$w_\mu\f$ are the weights
  ///          satisfying \f$w_1+\ldots+w_d=1\f$.
  template<typename T>
  void scale(CPTensor<T>& t, T alpha,
             const std::vector<double>& weights)
  {
#ifdef CHECK_RANGES_ON
    if (weights.size() != t.getSpatialDirections())
    {
      throw std::out_of_range("scale(CPTensor<T>&,"
                              " T, const std::vector<double>&)"
                              " : Weight's size mismatch");
    }
#endif

    bool negative = alpha < 0;
    if (negative) alpha = -alpha;
    for (int mu = 0; mu < t.getSpatialDirections(); mu++)
    {
#ifdef CHECK_RANGES_ON
      T alpha_weighted = std::pow(alpha, weights.at(mu));
#else
      T alpha_weighted = std::pow(alpha, weights[mu]);
#endif
      using namespace VectorOperators;
      // scale(t(mu), alpha_weighted);
      if (negative && mu == 0) { alpha_weighted = -alpha_weighted; }
      t(mu) *= alpha_weighted;
    }
  }

  /// Scales the tensor \c t by the factor \c alpha and
  /// constant wights \f$\frac1d\f$
  template<typename T>
  void scale(CPTensor<T>& t, T alpha)
  {
    std::vector<double> weights(t.get_d(),
                                1/static_cast<T>(t.get_d()));
    scale(t, alpha, weights); 
  }

  /// Performs \f$r:=x\odot y\f$
  template<typename T>
  void hadamardProduct(const CPTensor<T>& x, const CPTensor<T>& y,
                       CPTensor<T>& r)
  {
    const int d = std::min(x.get_d(), y.get_d());
    std::vector<int> dims(d);

    const std::vector<int> x_dims = x.get_n();
    const std::vector<int> y_dims = y.get_n();

    for (int mu = 0; mu < d; mu++)
    {
      dims[mu] = std::min(x_dims[mu], y_dims[mu]);
    }

    const int rang = x.get_r() * y.get_r();
    int k = 0;
    r.resize(rang, d, dims);

    for (int jx = 0; jx < x.getSeparationRang(); jx++)
    {
      for (int jy = 0; jy < y.getSeparationRang(); jy++)
      {
        for (int mu = 0; mu < d; mu++)
        {
          pointwiseProduct(x(jx, mu), y(jy, mu), r(k, mu));
        }
        k++;
      }
    }

  }

  /// Returns \f$\prod_{\mu=1}^d\|x_{j\mu}\|_2\f$
  template<typename T>
  const T frobeniusNormOfSummand(const CPTensor<T>& x, int j)
  {
    T result = 1;
    for (int mu = 0; mu < x.getSpatialDirections(); mu++)
    {
      result *= l2_norm(x(j, mu));
    }

    return result;
  }

  /// \brief Returns \f$\|x\|_2\f$
  /// \details \f[\|x\|_2^2=\sum_{j=1}^r\sum_{k=1}^r\prod_{\mu=1}^d x_{j\mu}\cdot x_{k\mu}\f]
  template<typename T>
  const T frobeniusNorm(const CPTensor<T>& x)
  {
    T result = 0;

    for (int j = 0; j < x.getSeparationRang(); j++)
    {
      T prod = 1;
      for (int mu = 0; mu < x.getSpatialDirections(); mu++)
      {
        prod *= innerProduct(x(j, mu), x(j, mu));
      }

      result += prod;

      for (int k = j+1; k < x.getSeparationRang(); k++)
      {
        prod = 1;
        for (int mu = 0; mu < x.getSpatialDirections(); mu++)
        {
          prod *= innerProduct(x(j, mu), x(k, mu));
        }

        result += (2*prod);
      }
    }

    return std::sqrt(result);
  }

  /// \brief Returns the inner product \f$x\cdot y\f$
  /// \details \f[\langle x,y\rangle=\sum_{j=1}^r\sum_{k=1}^r\prod_{\mu=1}^dx_{j\mu}\cdot y_{k\mu}\f]
  template<typename T>
  const T innerProduct(const CPTensor<T>& x, const CPTensor<T>& y)
  {
#ifdef RANGE_CHECKS_ON
    if ((x.getSpatialDirections() != y.getSpatialDirections()) ||
        (x.getDimensions() != y.getDimensions()))
    {
      throw std::out_of_range("inner_product(const CPTensor<T>& x,"
                              " const CPTensor<T>& y)"
                              " : x doesn't match y");
    }
#endif

    T result = 0;

    for (int j = 0; j < x.getSeparationRang(); j++)
    {
      for (int k = 0; k < y.getSeparationRang(); k++)
      {
        T prod = 1;

        for (int mu = 0; mu < x.getSpatialDirections(); mu++)
        {
          prod *= innerProduct(x(j, mu), y.getVectorOfRepresentants(k, mu));
        }
        result += prod;
      }
    }

    return result;
  }

  /// \brief Balances the tensor
  /// \details Let \f$x=\sum_j\bigotimes_\mu x_{j\mu}\f$
  ///   then this function alters it to
  ///   \f$x=\sum_j\bigotimes_\mu\left(\prod_{\nu\neq\mu}1/\|x_{j\nu}\|_2\right)x_{j\mu}\f$
  template<typename T>
  void balance(CPTensor<T>& x)
  {
    const int d = x.get_d();

    for (int j = 0; j < x.getSeparationRang(); j++)
    {
      T factor = std::pow(frobeniusNormOfSummand(x, j),
                          1/static_cast<T>(d));

      for (int mu = 0; mu < d; mu++)
      {
        WeakVector<T> a = x.getVectorOfRepresentants(j, mu);
        a *= (factor / l2_norm(a));
      }
    }

  }

  /// Returns \f$\|x-y\|_2\f$
  template<typename T>
  T distance(const CPTensor<T>& x, const CPTensor<T>& y)
  {
    using std::sqrt;
    using std::abs;

    const T a = innerProduct(x, x);
    const T b = innerProduct(x, y);
    const T c = innerProduct(y, y);
    
    const T value = a - 2 * b + c;
    // const T value = sqrt(a) * sqrt(c) * (sqrt(a) / sqrt(c) - 2 * b / (sqrt(a) * sqrt(c)) + sqrt(c) / sqrt(a));

    return sqrt(abs(value));
  }

  /// \brief Performs \f$y:=\alpha x+y\f$
  template<typename T>
  void add(T alpha, const CPTensor<T>& x, CPTensor<T>& y)
  {
    const int rx = x.get_r();
    const int ry = y.get_r();
    const int dx = x.get_d();
    const int dy = y.get_d();

    const std::vector<int> dimsx = x.get_n();

#ifdef RANGE_CHECKS_ON
    if ((dx != dy) || (dimsx != y.get_n()))
    {
      throw std::out_of_range("add(T, const CPTensor<T>&,"
                              "CPTensor<T>&) : Dimensions mismatch");
    }
#endif

    // it might be that x == y! resizing y would resize x too ...
    const CPTensor<T> old_y = y;

    y.resize(rx + ry, dx, dimsx);

    for (int j = 0; j < ry; j++)
    {
      for (int mu = 0; mu < dy; mu++)
      {
        y(j, mu) = old_y(j, mu);
      }
    }

    for (int j = 0; j < rx; j++)
    {
      for (int mu = 0; mu < dx; mu++)
      {
        WeakVector<T> w = y(j + ry, mu);
        w = x(j, mu);
        w *= alpha;
      }
    }

  }

  /// \brief Checks whether \c x and \c y are equal
  /// \details <tt>x==y</tt> iff \f$x_{(i_1,ldots,i_d)}=y_{(i_1,\ldots,i_d)}\f$
  ///          for all \f$i_\mu\in\{1,\ldots,n_\mu\}\f$<br>
  ///          This function tests the equality by testing
  ///          \f$\|x-y\|<\varepsilon\f$
  template<typename T>
  bool equals(const CPTensor<T>& x, const CPTensor<T>& y,
              T epsilon = Constants<T>::epsilon)
  {
    return distance(x, y) < epsilon;
  }

  /// \copybrief CPTensorCalculus::equals(const CPTensor<T>&,const CPTensor<T>&,T)
  /// \details <tt>x==y</tt> iff \f$x_{(i_1,ldots,i_d)}=y_{(i_1,\ldots,i_d)}\f$
  ///             for all \f$i_\mu\in\{1,\ldots,n_\mu\}\f$<br>
  ///          This function tests the equality by testing
  ///          \f$\|x-y\|<\varepsilon\f$
  template<typename T>
  bool operator== (const CPTensor<T>& x, const CPTensor<T>& y)
  {
    return equals(x, y);
  }


  /// \brief Checks whether \c x and \c y aren't equal
  /// \copydetails equals
  template<typename T>
  bool operator!= (const CPTensor<T>& x, const CPTensor<T>& y)
  {
    return !equals(x, y);
  }

  /// \copybrief add
  template<typename T>
  CPTensor<T> operator+(const CPTensor<T>& x, const CPTensor<T>& y)
  {
    CPTensor<T> result(y);
    add(1.0, x, result);
    return result;
  }

  /// \brief Pushes the contents of the CPTensor \c x to the stream \c s
  template<typename T>
  std::ostream& operator << (std::ostream& s, const CPTensor<T>& x)
  {
    const int d = x.get_d();
    const int r = x.get_r();
    s << "(" << d << ", " << r << ", [";

    const std::vector<int> dims = x.get_n();
    int max_dim = 0;

    for (int mu = 0; mu < d; mu++)
    {
      const int dim = dims[mu];
      if (dim > max_dim) max_dim = dim;
      s << dim;
      if (mu < d-1) s << ", ";
    }
    s << "])" << std::endl;

    s << std::setprecision(3)
      << std::scientific
      << std::setiosflags(std::ios::left)
      << std::setw(8);

    for(int i = 0; i < max_dim; i++)
    {
      for(int mu = 0; mu < d; mu++)
      {
        for(int j = 0; j < r; j++)
        {
          if (i < dims[mu])
          {
            s << x(j, mu, i) << " ";
          } else {
            s << "          ";
          }
        }
        s << "|";
      }
      s << std::endl;
    }

    return s;
  }

  /// \brief Squeezes out all summands with at least one vector
  ///        equals to zero and thus reducing the rank of \c x
  template<typename T>
  void compress(CPTensor<T>& x)
  {
    // TODO: Implement it!
  }

  /// \brief Fills in the vector \c x the entries of the CPTensor \c t
  template<typename T>
  void getTensor(const CPTensor<T>& t, std::vector<T>& x)
  {
    const int d = t.get_d();
    const std::vector<int>& n = t.get_n();

    x.clear();

    std::vector<int> index(d, 0);
    while (index[d-1] < n[d-1])
    {
      x.push_back( t.getTensorEntry(index) );
      index[0]++;
      for (int mu = 0; mu < d-1; ++mu)
      {
        if (index[mu] == n[mu])
        {
          index[mu+1]++;
          index[mu] = 0;
        } // else break;
      }
    }

  }
  
  template<typename T>
  void addSquareNormDerivative(const CPTensor<T>& x, const SkippedProducts< std::vector<T> >& skipped_scalar_products, CPTensor<T>& result)
  {
    const int d = x.getSpatialDirections();
    const int r = x.getSeparationRang();
    
    // CPTensor<T> result(r, d, x.getDimensions());
    for (int mu1 = 0; mu1 < d; ++mu1) {
      const int n_mu1 = x.getDimensionOfSpatDir(mu1);
      const std::vector<T>& x_mu1 = x.getVectorOfDimension(mu1);
      const std::vector<T>& scalar_product_skipped = skipped_scalar_products[mu1];
      
      std::vector<T>& dx_mu1 = result.getVectorOfDimension(mu1);
      Blas<T>::gemm('N', 'N', n_mu1, r, r, 1,
                    &x_mu1[0], n_mu1, &scalar_product_skipped[0], r,
                    1, &dx_mu1[0], n_mu1);
    }
  }

  template<typename T>
  std::vector<T> skippedScalarProduct(const CPTensor<T>& x, const int mu1) {
    const int d = x.getSpatialDirections();
    const int r = x.getSeparationRang();
    
    std::vector<T> scalar_product(r * r, 1);
    std::vector<T> scalar_product_update(r * r);
    
    StandardMonoidTraits< std::vector<T> > traits;

    for (int mu = 0; mu < mu1; ++mu) {
      const int n_mu = x.getDimensionOfSpatDir(mu);
      Blas<T>::gemm('T', 'N', r, r, n_mu, 1,
                    &x.getVectorOfDimension(mu)[0], n_mu,
                    &x.getVectorOfDimension(mu)[0], n_mu,
                    0, &scalar_product_update[0], r);
      traits.update(scalar_product_update, scalar_product);
    }
    
    for (int mu = mu1 + 1; mu < d; ++mu) {
      const int n_mu = x.getDimensionOfSpatDir(mu);
      Blas<T>::gemm('T', 'N', r, r, n_mu, 1,
                    &x.getVectorOfDimension(mu)[0], n_mu,
                    &x.getVectorOfDimension(mu)[0], n_mu,
                    0, &scalar_product_update[0], r);
      traits.update(scalar_product_update, scalar_product);
    }

    return scalar_product;
  }

  template<typename T>
  SkippedProducts< std::vector<T> > skippedScalarProducts(const CPTensor<T>& x)
  {
    const int d = x.getSpatialDirections();
    const int r = x.getSeparationRang();
    
    std::vector< std::vector<T> > scalar_products;
    scalar_products.reserve(x.getSpatialDirections());
    for (int mu = 0; mu < d; ++mu) {
      const int n_mu = x.getDimensionOfSpatDir(mu);
      
      std::vector<T> scalar_product(r * r);
      Blas<T>::gemm('T', 'N', r, r, n_mu, 1,
                    &x.getVectorOfDimension(mu)[0], n_mu,
                    &x.getVectorOfDimension(mu)[0], n_mu,
                    0, &scalar_product[0], r);
      scalar_products.push_back(scalar_product);
    }
    return SkippedProducts< std::vector<T> >(scalar_products);
  }

  template<typename T>
  void addSquareNormDerivative(const CPTensor<T>& x, CPTensor<T>& result) {
    addSquareNormDerivative(x, skippedScalarProducts(x), result);
    // return squareNormDerivative(x, skippedScalarProducts(x));
  }
  
  template<typename T>
  void addSecondPenaltyTermDerivatice(const CPTensor<T>& x, const SkippedProducts< std::vector<T> >& skipped_products, const T lambda2, CPTensor<T>& result) {
    const int d = x.get_d();
    const int r = x.get_r();
    
    for (int mu1 = 0; mu1 < d; ++mu1) {
      std::vector<T>& dx_mu1 = result.getVectorOfDimension(mu1);
      const std::vector<T>& x_mu1 = x.getVectorOfDimension(mu1);
      const int n_mu1 = x.getDimensionOfSpatDir(mu1);
      
      for (int j1 = 0; j1 < r; ++j1) {
        const T alpha_j1 = skipped_products[mu1][j1];
        Blas<T>::axpy(n_mu1, lambda2 * alpha_j1, &x_mu1[0], 1, &dx_mu1[0], 1);
      }
    }
  }

  template<typename T>
  void addSecondPenaltyTermDerivatice(const CPTensor<T>& x, const T lambda2, CPTensor<T>& result) {
    addSecondPenaltyTermDerivatice(x, skippedScalarProducts(x), lambda2, result);
  }

  template<typename T>
  struct ParameterVectorSpaceTraits;
  
  template<typename T>
  struct ParameterVectorSpaceTraits< CPTensor<T> > {
    typedef CPTensor<T> Vectors;
    typedef T Scalars;
    
    void update(const Scalars alpha, const Vectors& x, Vectors& y) const
    {
      const int d = x.getSpatialDirections();
      if (d != y.getSpatialDirections() ||
          x.getSeparationRang() != y.getSeparationRang() ||
          x.getDimensions() != y.getDimensions()) {
        y.resize(x.getSeparationRang(), d, x.getDimensions());
      }
      StandardVectorSpaceTraits< std::vector<T> > traits;
      for (int mu = 0; mu < d; ++mu) {
        traits.update(alpha, x.getVectorOfDimension(mu), y.getVectorOfDimension(mu));
      }
    }
    
    void scale(const Scalars alpha, Vectors& x) const
    {
      const int d = x.getSpatialDirections();
      StandardVectorSpaceTraits< std::vector<T> > traits;
      for (int mu = 0; mu < d; ++mu) {
        traits.scale(alpha, x.getVectorOfDimension(mu));
      }
    }
    
    Vectors zero() const {
      return Vectors();
    }
  };

  template<typename T>
  struct ParameterInnerProductSpaceTraits;
  
  template<typename T>
  struct ParameterInnerProductSpaceTraits< CPTensor<T> > : public ParameterVectorSpaceTraits< CPTensor<T> > {
    using ParameterVectorSpaceTraits< CPTensor<T> >::update;
    using ParameterVectorSpaceTraits< CPTensor<T> >::scale;
    using ParameterVectorSpaceTraits< CPTensor<T> >::zero;
    
    typedef CPTensor<T> Vectors;
    typedef T Scalars;
    
    Scalars innerProduct(const Vectors& x, const Vectors& y) const
    {
      const int d = x.getSpatialDirections();
#if defined(_DEBUG) || !defined(NDEBUG) || defined(RANGE_CHECKS_ON)     
      if (d != y.getSpatialDirections() ||
          x.getSeparationRang() != y.getSeparationRang() ||
          x.getDimensions() != y.getDimensions()) {
        throw std::invalid_argument("CP-Tensors don't match");
      }
#endif
      StandardInnerProductSpaceTraits< std::vector<T> > traits;
      Scalars result = 0;
      for (int mu = 0; mu < d; ++mu) {
        result += traits.innerProduct(x.getVectorOfDimension(mu), y.getVectorOfDimension(mu));
      }
      return result;
    }
  };
  
  template<typename T>
  void addCrossTensor(const CPTensor<T>& x, int mu1, const std::vector<int>& i_hat, const T alpha, std::vector<T>& result) {
    const int d = x.get_d();
    const int r = x.get_r();
    
    std::vector<T> prod(r, 1);
    for (int mu = 0; mu < d; ++mu) {
      if (mu != mu1) {
        const std::vector<T>& x_mu = x.getVectorOfDimension(mu);
        const int n_mu = x.getDimensionOfSpatDir(mu);
        
        for (int l = 0; l < r; ++l) {
          prod[l] *= x_mu[i_hat[mu] + n_mu * l];
        }
      }
    }
    
    const int n_mu1 = x.getDimensionOfSpatDir(mu1);
    Blas<T>::gemv('N', n_mu1, r, alpha, &x.getVectorOfDimension(mu1)[0], n_mu1,
                  &prod[0], 1, 1, &result[0], 1);
  }
    
  // f(y) = <x,y>
  // df/d(y_j1mu1)_imu11 = 
  template<typename T>
  void addTargetFunction(const CPTensor<T>& x, const CPTensor<T>& y, const T alpha, CPTensor<T>& result) {
    const int d = x.get_d();
    const int ry = y.get_r();
    const int rx = x.get_r();
    
    SkippedProducts< std::vector<T> > skipped_mixed;
    std::vector<T> sp_mixed(rx * ry); 
    
    // SkippedProducts< std::vector<T> > skipped_nonmixed;
    // std::vector<T> sp_nonmixed(rx * rx);
    
    for (int mu = 0; mu < d; ++mu) {
      const int n_mu = y.get_n(mu);
      const std::vector<T>& x_mu = x.getVectorOfDimension(mu);
      const std::vector<T>& y_mu = y.getVectorOfDimension(mu);
      
      Blas<T>::gemm('T', 'N', rx, ry, n_mu, 1,
                    &x_mu[0], n_mu,
                    &y_mu[0], n_mu,
                    0, &sp_mixed[0], rx);
      skipped_mixed.pushFactor(sp_mixed);
      
      // Blas<T>::gemm('T', 'N', rx, rx, n_mu, 1,
      //               &x_mu[0], n_mu,
      //               &x_mu[0], n_mu,
      //               1, &sp_nonmixed[0], rx);
      // skipped_nonmixed.pushFactor(sp_nonmixed);
    }
    
    for (int mu1 = 0; mu1 < d; ++mu1) {
      const int n_mu1 = y.get_n(mu1);
      const std::vector<T>& x_mu1 = x.getVectorOfDimension(mu1);
      
      Blas<T>::gemm('N', 'N', n_mu1, ry, rx, alpha,
                    &x_mu1[0], n_mu1,
                    &skipped_mixed[mu1][0], rx, 1,
                    &result.getVectorOfDimension(mu1)[0], n_mu1);
      
      //~ Blas<T>::gemm('T', 'N', n_mu1, rx, ry,
                    //~ alpha, &skipped_mixed[mu1][0], ry,
                            //~ &x_mu1[0], n_mu1,
                    //~ 1, &result.getVectorOfDimension(mu1)[0], n_mu1);
      // Blas<T>::gemm('T', 'N', n_mu1, rx, rx,
      //               alpha, &skipped_nonmixed[mu1][0], rx,
      //                      &x_mu1[0], n_mu1,
      //               1, &result.getVectorOfDimension(mu1)[0], n_mu1);
    }
  }

} // namespace TensorCalculus

#endif // __CP_TENSOR_OPERATORS_HPP

