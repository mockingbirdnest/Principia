#pragma once

#include "numerics/fast_fourier_transform.hpp"

namespace principia {
namespace numerics {
namespace internal_fast_fourier_transform {

template<typename Container, typename Scalar, int size_>
FastFourierTransform<Container, Scalar, size_>::FastFourierTransform(
    typename Container::const_iterator begin,
    typename Container::const_iterator end) {}

template<typename Container, typename Scalar, int size_>
std::array<Square<Scalar>, size_>
FastFourierTransform<Container, Scalar, size_>::PowerSpectrum() const {
}

}  // namespace internal_fast_fourier_transform
}  // namespace numerics
}  // namespace principia
