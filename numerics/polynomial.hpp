
#pragma once

#include <tuple>

#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

using quantities::Derivative;

template<typename Argument, typename Value, int degree>
class PolynomialCoefficients;

template<typename Argument, typename Value>
using PolynomialCoefficients<Vector, 0> = std::tuple<Vector>;

template<typename Argument, typename Value>
using PolynomialCoefficients<Vector, 1> = std::tuple<Vector, Derivative<Value, Argument>>;

template<typename Argument, typename Value>
class PolynomialInterface {
 public:
  Argument const& argument_min() const;
  Argument const& argument_max() const;

  virtual Value Evaluate(Argument const& argument) const = 0;
  virtual Derivate<Value, Argument> EvaluateDerivative(Argument const& argument) const = 0;

 protected:
  PolynomialInterface(Argument const& argument_min,
                      Argument const& argument_max);

  virtual ~PolynomialInterface() = default;

 private:
  Argument argument_min_;
  Argument argument_max_;
};

}  // namespace internal_polynomial
}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_body.hpp"
