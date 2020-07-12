
#pragma once

#include <array>
#include <complex>
#include <vector>

#include "base/bits.hpp"
#include "geometry/interval.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_fast_fourier_transform {

using base::FloorLog2;
using geometry::Interval;
using quantities::AngularFrequency;
using quantities::Square;
using quantities::Time;

// This class computes Fourier[{...}, FourierParameters -> {1, -1}] in
// Mathematica notation.  (The "signal processing" Fourier transform.)
template<typename Container, int size_>
class FastFourierTransform {
 public:
  // The size must be a power of 2.
  static constexpr int size = size_;
  static constexpr int log2_size = FloorLog2(size);
  static_assert(size == 1 << log2_size);

  using Scalar = typename Container::value_type;

  // In these constructors the samples are assumed to be separated by Δt.
  FastFourierTransform(Container const& container,
                       Time const& Δt);
  FastFourierTransform(typename Container::const_iterator begin,
                       typename Container::const_iterator end,
                       Time const& Δt);

  std::map<AngularFrequency, Square<Scalar>> PowerSpectrum() const;

  // Return the interval that contains the largest peak of power.
  Interval<AngularFrequency> Mode() const;

 private:
  Time const Δt_;
  AngularFrequency const ω_;

  // The elements of transform_ are in SI units of Scalar.  They are spaced in
  // frequency by ω_.
  std::array<std::complex<double>, size> transform_;

  friend class FastFourierTransformTest;
};

}  // namespace internal_fast_fourier_transform

using internal_fast_fourier_transform::FastFourierTransform;

}  // namespace numerics
}  // namespace principia

#include "numerics/fast_fourier_transform_body.hpp"
