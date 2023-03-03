#pragma once

#include <array>
#include <complex>
#include <map>
#include <type_traits>
#include <vector>

#include "base/bits.hpp"
#include "geometry/complexification.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/interval.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_fast_fourier_transform {

using namespace principia::base::_bits;
using namespace principia::geometry::_complexification;
using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_interval;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// Given (u₀, ..., uₙ₋₁), this class computes the discrete Fourier transform
//   Uₛ = ∑ᵣ uᵣ exp(-2πirs/n),
// corresponding to
//   Fourier[{...}, FourierParameters -> {1, -1}]
// in Mathematica notation (the "signal processing" Fourier transform).
template<typename Value, typename Argument, std::size_t size_>
class FastFourierTransform {
 public:
  // This is only an actual angular frequency if |Argument| is time-like.
  // If |Argument| is an angular frequency, this is a time.
  using AngularFrequency = Derivative<Angle, Argument>;

  // The size must be a power of 2.
  static constexpr int size = size_;
  static constexpr int log2_size = FloorLog2(size);
  static_assert(size == 1 << log2_size);

  // In the constructors, the container must have |size| elements.  For the
  // purpose of expressing the frequencies, the values are assumed to be
  // sampled at intervals of Δt.

  template<typename Container,
           typename = std::enable_if_t<
               std::is_convertible_v<typename Container::value_type, Value>>>
  FastFourierTransform(Container const& container,
                       Difference<Argument> const& Δt);

  template<typename Iterator,
           typename = std::enable_if_t<std::is_convertible_v<
               typename std::iterator_traits<Iterator>::value_type,
               Value>>>
  FastFourierTransform(Iterator begin, Iterator end,
                       Difference<Argument> const& Δt);

  FastFourierTransform(std::array<Value, size> const& container,
                       Difference<Argument> const& Δt);

  std::map<AngularFrequency, typename Hilbert<Value>::Norm²Type>
  PowerSpectrum() const;

  // Returns the interval that contains the largest peak of power in the
  // specifed range.
  Interval<AngularFrequency> Mode(AngularFrequency const& min_ω,
                                  AngularFrequency const& max_ω) const;

  // Given s ∈ [0, size - 1] ∩ ℕ, returns the coefficient Uₛ.
  Complexification<Value> const& operator[](int s) const;

  // Given s ∈ [0, size - 1] ∩ ℕ, returns the frequency corresponding to Uₛ.
  AngularFrequency frequency(int s) const;

 private:
  Difference<Argument> const Δt_;
  AngularFrequency const Δω_;

  // The elements of transform_ are spaced in frequency by ω_.
  std::array<Complexification<Value>, size> transform_;

  friend class FastFourierTransformTest;
};

}  // namespace internal_fast_fourier_transform

using internal_fast_fourier_transform::FastFourierTransform;

}  // namespace numerics
}  // namespace principia

#include "numerics/fast_fourier_transform_body.hpp"
