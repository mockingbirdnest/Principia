#pragma once

#include "numerics/fast_fourier_transform.hpp"

#include <map>
#include <optional>

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

// Implementation of the Danielson-Lánczos algorithm using templates for
// recursion and template specializations for short FFTs [DL42, Myr07].
template<typename Complex, int array_size_, int chunk_size_ = array_size_>
class DanielsonLánczos {
 public:
  static void Transform(typename std::array<Complex,
                                            array_size_>::iterator const begin);
};

template<typename Complex, int array_size_>
class DanielsonLánczos<Complex, array_size_, 1> {
 public:
  static void Transform(typename std::array<Complex,
                                            array_size_>::iterator const begin);
};

template<typename Complex, int array_size_>
class DanielsonLánczos<Complex, array_size_, 2> {
 public:
  static void Transform(typename std::array<Complex,
                                            array_size_>::iterator const begin);
};

template<typename Complex, int array_size_>
class DanielsonLánczos<Complex, array_size_, 4> {
 public:
  static void Transform(typename std::array<Complex,
                                            array_size_>::iterator const begin);
};

template<typename Complex, int array_size_, int chunk_size_>
void DanielsonLánczos<Complex, array_size_, chunk_size_>::Transform(
    typename std::array<Complex, array_size_>::iterator const begin) {
  constexpr int N = chunk_size_;

  DanielsonLánczos<Complex, array_size_, N / 2>::Transform(begin);
  DanielsonLánczos<Complex, array_size_, N / 2>::Transform(begin + N / 2);

  Angle const θ = π * Radian / N;
  double const sin_θ = Sin(θ);
  double const cos_2θ_minus_1 = -2 * sin_θ * sin_θ;
  double const sin_2θ = Sin(2 * θ);
  // Computing e⁻²ⁱ⁽ᵏ⁺¹⁾ᶿ as e⁻²ⁱᵏᶿ + e⁻²ⁱᵏᶿ (e⁻²ⁱᶿ - 1) rather than e⁻²ⁱᵏᶿe⁻²ⁱᶿ
  // improves accuracy [Myr07].
  Complexification<double> const e⁻²ⁱᶿ_minus_1(cos_2θ_minus_1, -sin_2θ);
  Complexification<double> e⁻²ⁱᵏᶿ = 1;
  auto it = begin;
  for (int k = 0; k < N / 2; ++it, ++k, e⁻²ⁱᵏᶿ += e⁻²ⁱᵏᶿ * (e⁻²ⁱᶿ_minus_1)) {
    auto const t = *(it + N / 2) * e⁻²ⁱᵏᶿ;
    *(it + N / 2) = *it - t;
    *it += t;
  }
}

template<typename Complex, int array_size_>
void DanielsonLánczos<Complex, array_size_, 1>::Transform(
    typename std::array<Complex, array_size_>::iterator const begin) {}

template<typename Complex, int array_size_>
void DanielsonLánczos<Complex, array_size_, 2>::Transform(
    typename std::array<Complex, array_size_>::iterator const begin) {
  auto const t = *(begin + 1);
  *(begin + 1) = *begin - t;
  *begin += t;
}

template<typename Complex, int array_size_>
void DanielsonLánczos<Complex, array_size_, 4>::Transform(
    typename std::array<Complex, array_size_>::iterator const begin) {
  {
    auto const t = *(begin + 1);
    *(begin + 1) = *begin - t;
    *begin += t;
  }
  {
    auto const t = *(begin + 3);
    *(begin + 3) = {(begin + 2)->imaginary_part() - t.imaginary_part(),
                    t.real_part() - (begin + 2)->real_part()};
    *(begin + 2) += t;
  }
  {
    auto const t = *(begin + 2);
    *(begin + 2) = *begin - t;
    *begin += t;
  }
  {
    auto const t = *(begin + 3);
    *(begin + 3) = *(begin + 1) - t;
    *(begin + 1) += t;
  }
}

template<typename Value, typename Argument, std::size_t size_>
template<typename Container, typename>
FastFourierTransform<Value, Argument, size_>::FastFourierTransform(
    Container const& container,
    Difference<Argument> const& Δt)
    : FastFourierTransform(container.cbegin(), container.cend(), Δt) {}

template<typename Value, typename Argument, std::size_t size_>
template<typename Iterator, typename>
FastFourierTransform<Value, Argument, size_>::FastFourierTransform(
    Iterator const begin,
    Iterator const end,
    Difference<Argument> const& Δt)
    : Δt_(Δt),
      Δω_(2 * π * Radian / (size * Δt_)) {
  DCHECK_EQ(size, std::distance(begin, end));

  // Reindexing and promotion to complex.
  int bit_reversed_index = 0;
  for (auto it = begin;
       it != end;
       ++it,
       bit_reversed_index = BitReversedIncrement(bit_reversed_index,
                                                 log2_size)) {
    transform_[bit_reversed_index] = *it;
  }

  DanielsonLánczos<Complexification<Value>, size>::Transform(
      transform_.begin());
}

template<typename Value, typename Argument, std::size_t size_>
FastFourierTransform<Value, Argument, size_>::FastFourierTransform(
    std::array<Value, size> const& container,
    Difference<Argument> const& Δt)
    : FastFourierTransform(container.cbegin(), container.cend(), Δt) {}

template<typename Value, typename Argument, std::size_t size_>
auto FastFourierTransform<Value, Argument, size_>::PowerSpectrum() const
    -> std::map<AngularFrequency, typename Hilbert<Value>::InnerProductType> {
  std::map<AngularFrequency, typename Hilbert<Value>::InnerProductType>
      spectrum;
  int k = 0;
  for (auto const& coefficient : transform_) {
    spectrum.emplace_hint(spectrum.end(), k * Δω_, coefficient.Norm²());
    ++k;
  }
  return spectrum;
}

template<typename Value, typename Argument, std::size_t size_>
auto FastFourierTransform<Value, Argument, size_>::Mode() const
    -> Interval<AngularFrequency> {
  auto const spectrum = PowerSpectrum();
  typename std::map<AngularFrequency,
                    typename Hilbert<Value>::InnerProductType>::const_iterator
      max = spectrum.end();

  // Only look at the first size / 2 + 1 elements because the spectrum is
  // symmetrical.
  auto it = spectrum.begin();
  for (int i = 0; i < size / 2 + 1; ++i, ++it) {
    if (max == spectrum.end() || it->second > max->second) {
      max = it;
    }
  }
  Interval<AngularFrequency> result;
  if (max == spectrum.begin()) {
    result.Include(max->first);
  } else {
    result.Include(std::prev(max)->first);
  }
  result.Include(std::next(max)->first);
  return result;
}

template<typename Value, typename Argument, std::size_t size_>
Complexification<Value> const&
    FastFourierTransform<Value, Argument, size_>::operator[](
        int const s) const {
  return transform_[s];
}

template<typename Value, typename Argument, std::size_t size_>
typename FastFourierTransform<Value, Argument, size_>::AngularFrequency
FastFourierTransform<Value, Argument, size_>::frequency(int const s) const {
  DCHECK_GE(s, 0);
  DCHECK_LT(s, size);
  return s * Δω_;
}

}  // namespace internal_fast_fourier_transform
}  // namespace numerics
}  // namespace principia
