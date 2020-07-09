#pragma once

#include "numerics/fast_fourier_transform.hpp"

#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace internal_fast_fourier_transform {

using quantities::si::Unit;

template<typename Container, typename Scalar, int size_>
FastFourierTransform<Container, Scalar, size_>::FastFourierTransform(
    typename Container::const_iterator begin,
    typename Container::const_iterator end) {
  int i = 0;

  // Type decay and promotion to complex.
  for (auto it = begin; it != end; ++it, ++i) {
    transform_[i] = *it / si::Unit<Scalar>;
  }
  CHECK_EQ(size_, i);


}

template<typename Container, typename Scalar, int size_>
std::array<Square<Scalar>, size_>
FastFourierTransform<Container, Scalar, size_>::PowerSpectrum() const {
}

}  // namespace internal_fast_fourier_transform
}  // namespace numerics
}  // namespace principia
