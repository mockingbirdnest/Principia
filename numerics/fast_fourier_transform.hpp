
#pragma once

#include <array>
#include <complex>
#include <vector>

#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_fast_fourier_transform {

using quantities::Square;

template<typename Container, typename Scalar, int size_>
class FastFourierTransform {
 public:
  FastFourierTransform(typename Container::const_iterator begin,
                       typename Container::const_iterator end);

  std::array<Square<Scalar>, size_> PowerSpectrum() const;

 private:
  std::array<std::complex<double>> transform_;
};

}  // namespace internal_fast_fourier_transform
}  // namespace numerics
}  // namespace principia

#include "numerics/fast_fourier_transform_body.hpp"