
#pragma once

#include <tuple>
#include <utility>

#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

using quantities::Derivative;
using quantities::Difference;

//TODO(phl): move elsewhere?
template<typename Value, typename Argument, int order>
struct NthDerivativeGenerator {
  using Type = Derivative<
      typename NthDerivativeGenerator<Value, Argument, order - 1>::Type,
      Argument>;
};
template<typename Value, typename Argument>
struct NthDerivativeGenerator<Value, Argument, 0> {
  using Type = Value;
};

template<typename Value, typename Argument, int order>
using NthDerivative =
    typename NthDerivativeGenerator<Value, Argument, order>::Type;

template<typename Value, typename Argument, typename>
struct NthDerivativesGenerator;
template<typename Value, typename Argument, int... orders>
struct NthDerivativesGenerator<Value,
                               Argument,
                               std::integer_sequence<int, orders...>> {
  using Type = std::tuple<typename NthDerivative<Value, Argument, orders>...>;
};

template<typename Value, typename Argument, typename Sequence>
using NthDerivatives =
    typename NthDerivativesGenerator<Value, Argument, Sequence>::Type;

template<typename Value, typename Argument>
class Polynomial {
 public:
  Argument const& argument_min() const;
  Argument const& argument_max() const;

  virtual Value Evaluate(Argument const& argument) const = 0;
  virtual Derivative<Value, Argument> EvaluateDerivative(
      Argument const& argument) const = 0;

 protected:
  Polynomial(Argument const& argument_min,
             Argument const& argument_max);

  virtual ~Polynomial() = default;

 private:
  Argument argument_min_;
  Argument argument_max_;
};

template<typename Value, typename Argument, int degree>
class PolynomialInMonomialBasis : public Polynomial<Value, Argument> {
 public:
  // Equivalent to:
  //   std::tuple<Value,
  //              Derivative<Value, Argument>,
  //              Derivative<Derivative<Value, Argument>>...>
  using Coefficients =
      typename NthDerivatives<Value,
                              Argument,
                              std::make_integer_sequence<int, degree + 1>>;

  PolynomialInMonomialBasis(Coefficients const& coefficients,
                            Argument const& argument_min,
                            Argument const& argument_max);

  Value Evaluate(Argument const& argument) const override;
  Derivative<Value, Argument> EvaluateDerivative(
      Argument const& argument) const override;

 private:
  Coefficients coefficients_;
};

}  // namespace internal_polynomial

using internal_polynomial::Polynomial;
using internal_polynomial::PolynomialInMonomialBasis;

}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_body.hpp"
