#ifndef __RANDOM_HPP
#define __RANDOM_HPP

#include <complex>

#include <cstdlib>

namespace TensorCalculus {
  
  template<typename T>
  class Random {
  public:
    Random();
    Random(T min, T max) : min(min), max(max) { }
    
    T operator () () const;
  private:
    T min;
    T max;
  };
  
  template<>
  Random<double>::Random() : min(0.0), max(1.0) { }
  
  template<>
  double Random<double>::operator () () const {
    return (max - min) * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX) + min;
  }

  template<>
  Random<int>::Random() : min(0), max(RAND_MAX) { }
  
  template<>
  int Random<int>::operator () () const {
    return std::rand() % (max - min + 1) + min;
  }

  // Offser no min and max due to the lack of
  // an total ordering on the complex numbers
  template<>
  class Random< std::complex<double> > {
  public:
    std::complex<double> operator () () const {
      return std::complex<double>(static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX),
                                  static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX));
    }
  };

} // namespace TensorCalculus

#endif // __RANDOM_HPP
