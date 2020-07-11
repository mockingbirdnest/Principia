#pragma once

#include "numerics/fast_fourier_transform.hpp"

#include "base/bits.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace internal_fast_fourier_transform {

using base::BitReversedIncrement;
using quantities::si::Unit;

template<typename Container, typename Scalar, int size_>
FastFourierTransform<Container, Scalar, size_>::FastFourierTransform(
    typename Container::const_iterator begin,
    typename Container::const_iterator end) {
  int bit_reversed_index = 0;

  // Type decay, reindexing, and promotion to complex.
  for (auto it = begin;
       it != end;
       ++it,
       bit_reversed_index = BitReversedIncrement(bit_reversed_index,
                                                 log2_size)) {
    transform_[bit_reversed_index] = *it / si::Unit<Scalar>;
  }
  CHECK_EQ(size, bit_reversed_index);


}

template<typename Container, typename Scalar, int size_>
std::array<Square<Scalar>, size_>
FastFourierTransform<Container, Scalar, size_>::PowerSpectrum() const {
}

}  // namespace internal_fast_fourier_transform
}  // namespace numerics
}  // namespace principia
