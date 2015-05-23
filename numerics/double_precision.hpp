#pragma once

#include "quantities/quantities.hpp"

namespace principia {

using quantities::Difference;

namespace numerics {

// A simple container for accumulating a value using compensated summation.  The
// type of the value must be an affine space.  The value constructor is not
// explicit to make it easy to construct an object with no error.
template<typename T>
struct DoublePrecision {
  DoublePrecision() = default;
  DoublePrecision(T const& value);  // NOLINT(runtime/explicit)

  void Increment(Difference<T> const& increment);

  T value;
  Difference<T> error;
};

}  // namespace numerics
}  // namespace principia

#include "numerics/double_precision_body.hpp"
