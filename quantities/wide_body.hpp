
#pragma once

#include "quantities/wide.hpp"

namespace principia {
namespace quantities {
namespace internal_wide {

template<typename T>
Wide<T>::Wide(T const x) : wide_(_mm_set1_pd(static_cast<double>(x))) {}

template<typename T>
__m128d Wide<T>::m128d() const {
  return wide_;
}

template<typename D>
Wide<Quantity<D>>::Wide(Quantity<D> const& x)
    : wide_(_mm_set1_pd(x.magnitude_)) {}

template<typename D>
__m128d Wide<Quantity<D>>::m128d() const {
  return wide_;
}

template<typename T>
Wide<Wide<T>>::Wide(Wide<T> const& x) : wide_(x.wide_) {}

template<typename T>
__m128d Wide<Wide<T>>::m128d() const {
  return wide_;
}

}  // namespace internal_wide
}  // namespace quantities
}  // namespace principia
