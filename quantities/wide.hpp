
#pragma once

#include <nmmintrin.h>

#include <type_traits>

#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
namespace internal_wide {

template<typename T>
class Wide final {
  static_assert(std::is_arithmetic<T>::value, "Nonarithmetic type");
 public:
  explicit Wide(T x);

  __m128d m128d() const;

 private:
  __m128d wide_;
};

template<typename D>
class Wide<Quantity<D>> final {
 public:
  explicit Wide(Quantity<D> const& x);

  __m128d m128d() const;

 private:
  __m128d wide_;
};

template<typename T>
class Wide<Wide<T>> final {
 public:
  explicit Wide(Wide<T> const& x);

  __m128d m128d() const;

 private:
  __m128d wide_;
};

}  // namespace internal_wide

using internal_wide::Wide;

}  // namespace quantities
}  // namespace principia

#include "quantities/wide_body.hpp"
