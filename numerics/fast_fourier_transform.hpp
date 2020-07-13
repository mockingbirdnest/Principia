
#pragma once

#include <array>
#include <complex>
#include <type_traits>
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
// Mathematica notation (the "signal processing" Fourier transform).
template<typename Scalar, std::size_t size_>
class FastFourierTransform {
 public:
  // The size must be a power of 2.
  static constexpr int size = size_;
  static constexpr int log2_size = FloorLog2(size);
  static_assert(size == 1 << log2_size);

  // In the constructors, the container must have |size| elements.  The samples
  // are assumed to be separated by Δt.

  template<typename Container,
           typename = std::enable_if_t<
               std::is_convertible_v<typename Container::value_type, Scalar>>>
  FastFourierTransform(Container const& container,
                       Time const& Δt);

  template<typename Iterator,
           typename = std::enable_if_t<std::is_convertible_v<
               typename std::iterator_traits<Iterator>::value_type,
               Scalar>>>
  FastFourierTransform(Iterator begin, Iterator end,
                       Time const& Δt);

  FastFourierTransform(std::array<Scalar, size> const& container,
                       Time const& Δt);

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
