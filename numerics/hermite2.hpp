
#pragma once

#include <functional>
#include <utility>
#include <vector>

#include "base/array.hpp"
#include "geometry/hilbert.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_hermite2 {

using quantities::Derivative;
using quantities::Difference;
using geometry::Hilbert;

// A 2nd degree Hermite polynomial defined by its values at the bounds of some
// interval and one of its derivatives.
template<typename Value, typename Argument>
class Hermite2 final {
 public:
  using Derivative1 = Derivative<Value, Argument>;

  // TODO(phl): Add support for providing the derivative at the upper bound.
  Hermite2(std::pair<Argument, Argument> arguments,
           std::pair<Value, Value> const& values,
           Derivative1 const& derivative_lo);

  Value Evaluate(Argument const& argument) const;
  Derivative1 EvaluateDerivative(Argument const& argument) const;

  Argument FindExtremum() const;

 private:
  using Derivative2 = Derivative<Derivative1, Argument>;

  std::pair<Argument, Argument> const arguments_;
  Value a0_;
  Derivative1 a1_;
  Derivative2 a2_;
};

}  // namespace internal_hermite2

using internal_hermite2::Hermite2;

}  // namespace numerics
}  // namespace principia

#include "numerics/hermite2_body.hpp"
