#pragma once

#include <set>
#include <utility>

#include "quantities/quantities.hpp"

namespace principia {

using quantities::Difference;
using quantities::Exponentiation;
using quantities::Quotient;

namespace numerics {

// A 3rd-degree Hermite polynomial defined by its values and derivatives at the
// bounds of some interval.
template<typename Argument, typename Value>
class Hermite3 {
 public:
  using Derivative = Quotient<Value, Difference<Argument>>;

  Hermite3(std::pair<Argument, Argument> const& arguments,
           std::pair<Value, Value> const& values,
           std::pair<Derivative, Derivative> const& derivatives);

  Value Evaluate(Argument const& argument) const;
  Derivative EvaluateDerivative(Argument const& argument) const;

  std::set<Argument> FindExtrema() const;

 private:
  using Derivative2 = Quotient<Derivative, Difference<Argument>>;
  using Derivative3 = Quotient<Derivative2, Difference<Argument>>;

  std::pair<Argument, Argument> const arguments_;
  Value a0_;
  Derivative a1_;
  Derivative2 a2_;
  Derivative3 a3_;
};

}  // namespace numerics
}  // namespace principia

#include "numerics/hermite3_body.hpp"
