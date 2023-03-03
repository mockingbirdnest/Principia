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

using namespace principia::geometry::_hilbert;
using namespace principia::quantities::_named_quantities;

// A 2nd degree Hermite polynomial defined by its values at the bounds of some
// interval and one of its derivatives.
template<typename Value, typename Argument>
class Hermite2 final {
 public:
  using Derivative1 = Derivative<Value, Argument>;

  // The derivative corresponds to the first argument.  Note that we don't
  // require any ordering of the arguments, so this constructor is usable if the
  // derivative is known at the upper bound of an interval.
  Hermite2(std::pair<Argument, Argument> arguments,
           std::pair<Value, Value> const& values,
           Derivative1 const& derivative_first);

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
