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
namespace si = quantities::si;

template<int array_size_, int chunk_size_ = array_size_>
class DanielsonLánczos {
 public:
  static void Transform(typename std::array<std::complex<double>,
                                            array_size_>::iterator const begin);
};

template<int array_size_>
class DanielsonLánczos<array_size_, 1> {
 public:
  static void Transform(typename std::array<std::complex<double>,
                                            array_size_>::iterator const begin);
};

template<int array_size_, int chunk_size_>
void DanielsonLánczos<array_size_, chunk_size_>::Transform(
    typename std::array<std::complex<double>, array_size_>::iterator const
        begin) {
  constexpr int N = chunk_size_;

  DanielsonLánczos<array_size_, N / 2>::Transform(begin);
  DanielsonLánczos<array_size_, N / 2>::Transform(begin + N / 2);

  Angle const θ = π * Radian / N;
  double const sin_θ = Sin(θ);
  double const cos_2θ_minus_1 = -2 * sin_θ * sin_θ;
  double const sin_2θ = Sin(2 * θ);
  std::complex<double> const e⁻²ⁱᶿ_minus_1(cos_2θ_minus_1, -sin_2θ);
  std::complex<double> e⁻²ⁱᵏᶿ = 1;
  auto it = begin;
  for (int k = 0; k < N / 2; ++it, ++k) {
    auto t = *it * e⁻²ⁱᵏᶿ;
    *(it + N / 2) = *it - t;
    *it += t;

    e⁻²ⁱᵏᶿ += e⁻²ⁱᵏᶿ * (e⁻²ⁱᶿ_minus_1);
  }
}

template<int array_size_>
void DanielsonLánczos<array_size_, 1>::Transform(
    typename std::array<std::complex<double>, array_size_>::iterator const
        begin) {}

template<typename Container, int size_>
FastFourierTransform<Container, size_>::FastFourierTransform(
    Container const& container)
    : FastFourierTransform(container.cbegin(), container.cend()) {}

template<typename Container, int size_>
FastFourierTransform<Container, size_>::FastFourierTransform(
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

  DanielsonLánczos<size>::Transform(transform_.begin());
}

template<typename Container, int size_>
std::array<Square<typename FastFourierTransform<Container, size_>::Scalar>,
           size_>
FastFourierTransform<Container, size_>::PowerSpectrum() const {
}

}  // namespace internal_fast_fourier_transform
}  // namespace numerics
}  // namespace principia
