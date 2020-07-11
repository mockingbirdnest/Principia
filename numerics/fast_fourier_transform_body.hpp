#pragma once

#include "numerics/fast_fourier_transform.hpp"

#include "base/bits.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace internal_fast_fourier_transform {

using base::BitReversedIncrement;
using quantities::Angle;
using quantities::Sin;
using quantities::si::Radian;
using quantities::si::Unit;

template<int size_>
void DanielsonLánczos(
    std::array<std::complex<double>, size_>::iterator const begin) {
  constexpr int size = size_;

  DanielsonLánczos<size / 2>(begin);
  DanielsonLánczos<size / 2>(begin + size / 2);

  Angle const θ = π * Radian / size;
  double const sin_θ = Sin(θ);
  double const cos_2θ_minus_1 = -2 * sin_θ * sin_θ;
  double const sin_2θ = Sin(2 * θ);
  std::complex<double> const e⁻²ⁱᶿ_minus_1(cos_2θ_minus_1, -sin_2θ);
  std::complex<double> const e⁻²ⁱᵏᶿ = 1;
  auto it = begin;
  for (int i = 0; i < size / 2; ++it, ++i) {
    auto t = *it * e⁻²ⁱᵏᶿ;
    *(it + size / 2) = *it - t;
    *it += t;

    e⁻²ⁱᵏᶿ += e⁻²ⁱᵏᶿ * (e⁻²ⁱᶿ_minus_1);
  }
}

template<>
void DanielsonLánczos<1>(
    std::array<std::complex<double>, 1>::iterator const begin) {}

template<typename Container, typename Scalar, int size_>
FastFourierTransform<Container, Scalar, size_>::FastFourierTransform(
    typename Container::const_iterator begin,
    typename Container::const_iterator end) {
  DCHECK_EQ(size, std::distance(begin, end));

  // Type decay, reindexing, and promotion to complex.
  int bit_reversed_index = 0;
  for (auto it = begin;
       it != end;
       ++it,
       bit_reversed_index = BitReversedIncrement(bit_reversed_index,
                                                 log2_size)) {
    transform_[bit_reversed_index] = *it / si::Unit<Scalar>;
  }

  DanielsonLánczos<size>(transform_.begin());
}

template<typename Container, typename Scalar, int size_>
std::array<Square<Scalar>, size_>
FastFourierTransform<Container, Scalar, size_>::PowerSpectrum() const {
}

}  // namespace internal_fast_fourier_transform
}  // namespace numerics
}  // namespace principia
