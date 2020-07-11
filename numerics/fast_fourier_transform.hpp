
#pragma once

#include <array>
#include <complex>
#include <vector>

#include "base/bits.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_fast_fourier_transform {

using base::FloorLog2;
using quantities::Square;

template<typename Container, int size_>
class FastFourierTransform {
 public:
  // The size must be a power of 2.
  static constexpr int size = size_;
  static constexpr int log2_size = FloorLog2(size);
  static_assert(size == 1 << log2_size);

  using Scalar = typename Container::value_type;

  explicit FastFourierTransform(Container const& container);

  FastFourierTransform(typename Container::const_iterator begin,
                       typename Container::const_iterator end);

  std::array<Square<Scalar>, size_> PowerSpectrum() const;

 //private:
  std::array<std::complex<double>, size> transform_;
};

}  // namespace internal_fast_fourier_transform

using internal_fast_fourier_transform::FastFourierTransform;

}  // namespace numerics
}  // namespace principia

#include "numerics/fast_fourier_transform_body.hpp"