
#pragma once

#include <utility>
#include <vector>

#include "base/array.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_hermite3 {

using base::BoundedArray;
using quantities::Derivative;

// A 3rd degree Hermite polynomial defined by its values and derivatives at the
// bounds of some interval.
template<typename Argument, typename Value>
class Hermite3 final {
 public:
  using Derivative1 = Derivative<Value, Argument>;

  Hermite3(std::pair<Argument, Argument> const& arguments,
           std::pair<Value, Value> const& values,
           std::pair<Derivative1, Derivative1> const& derivatives);

  // NOTE(egg): this does not appear to use Casteljau's algorithm; perhaps it
  // should?
  Value Evaluate(Argument const& argument) const;
  Derivative1 EvaluateDerivative(Argument const& argument) const;

  // The result is sorted.
  BoundedArray<Argument, 2> FindExtrema() const;

 private:
  using Derivative2 = Derivative<Derivative1, Argument>;
  using Derivative3 = Derivative<Derivative2, Argument>;

  std::pair<Argument, Argument> const arguments_;
  Value a0_;
  Derivative1 a1_;
  Derivative2 a2_;
  Derivative3 a3_;
};

}  // namespace internal_hermite3

using internal_hermite3::Hermite3;

}  // namespace numerics
}  // namespace principia

#include "numerics/hermite3_body.hpp"
