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

#ifndef __SKIPPEDPRODUCTS_HPP
#define __SKIPPEDPRODUCTS_HPP

#include <vector>
#include "StandardTraits.hpp"

namespace TensorCalculus {

  template<typename Monoid,
           template<typename> class MonoidTraits = StandardMonoidTraits>
  class SkippedProducts {
  public:
    SkippedProducts(const std::vector<Monoid>& factors, MonoidTraits<Monoid> traits = MonoidTraits<Monoid>())
      : traits(traits),
        skipped_products(calc_skipped_products(factors)),
        full_product(calc_full_product(factors))
    { }

    explicit SkippedProducts(MonoidTraits<Monoid> traits = MonoidTraits<Monoid>())
      : traits(traits),
        skipped_products(),
        full_product(traits.identityElement())
    { }

    int countSkippedProducts() const { return skipped_products.size(); }
    const Monoid& getSkippedProduct(int i) const { return skipped_products[i]; }
    
    const Monoid& operator [] (int i) const { return skipped_products[i]; }

    const Monoid& getFullProduct() const { return full_product; }

    void setFactors(const std::vector<Monoid>& factors);

    void pushFactor(const Monoid& factor);
    void clear();

  private:
    MonoidTraits<Monoid> traits;

    std::vector<Monoid> skipped_products;
    Monoid full_product;

    std::vector<Monoid> calc_skipped_products(const std::vector<Monoid>& factors) const;

    Monoid calc_full_product(const std::vector<Monoid>& factors) const;
  };

  template<typename Monoid,
           template<typename> class MonoidTraits>
  std::vector<Monoid> SkippedProducts<Monoid, MonoidTraits>::calc_skipped_products(const std::vector<Monoid>& factors) const // , MonoidTraits<Monoid> traits) const
  {
    const int n = factors.size();
    std::vector<Monoid> result;
    if (n == 0) {
      return result;
    }
    result.reserve(n);
    
    Monoid partial_product = traits.identityElement();
    result.push_back(partial_product);
    
    for (int i = 1; i < n; ++i) {
      traits.update(factors[i - 1], partial_product);
      result.push_back(partial_product);
    }
    
    partial_product = traits.identityElement();
    for (int i = n - 2; i >= 0; --i) {
      traits.update(factors[i + 1], partial_product);
      traits.update(partial_product, result[i]);
    }
    
    return result;
  }

  template<typename Monoid,
           template<typename> class MonoidTraits>
  Monoid SkippedProducts<Monoid, MonoidTraits>::calc_full_product(const std::vector<Monoid>& factors) const // , MonoidTraits<Monoid> traits) const
  {
    if (factors.size() == 0) {
      return traits.identityElement();
    } else {
      Monoid result = skipped_products[0];
      traits.update(factors[0], result);
      return result;
    }
  }

  template<typename Monoid,
           template<typename> class MonoidTraits>
  void SkippedProducts<Monoid, MonoidTraits>::setFactors(const std::vector<Monoid>& factors) {
    skipped_products = calc_skipped_products(factors);
    full_product = calc_full_product(factors);
  }
  
  template<typename Monoid,
           template<typename> class MonoidTraits>
  void SkippedProducts<Monoid, MonoidTraits>::pushFactor(const Monoid& factor) {
    if (skipped_products.size() == 0) {
      skipped_products.push_back(traits.identityElement());
      full_product = factor;
    } else {
      for (unsigned int i = 0; i < skipped_products.size(); ++i) {
        traits.update(factor, skipped_products[i]);
      }
      skipped_products.push_back(full_product);
      traits.update(factor, full_product);
    }
  }

  template<typename Monoid,
           template<typename> class MonoidTraits>
  void SkippedProducts<Monoid, MonoidTraits>::clear() {
    skipped_products.clear();
    full_product = traits.identityElement();
  }

} // namespace TensorCalculus

#endif // __SKIPPEDPRODUCTS_HPP
