
#pragma once

#include <tuple>
#include <utility>

#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

using quantities::Derivative;

//TODO(phl): move elsewhere?
template<typename Value, typename Argument, int order>
struct NthDerivative {
  using type =
      Derivative<typename NthDerivative<Value, Argument, order - 1>::type,
                 Argument>;
};
template<typename Value, typename Argument>
struct NthDerivative<Value, Argument, 0> {
  using type = Value;
};

template<typename Value, typename Argument, typename>
struct NthDerivatives;
template<typename Value, typename Argument, int... orders>
struct NthDerivatives<Value, Argument, std::integer_sequence<int, orders...>> {
  using type =
      std::tuple<typename NthDerivative<Value, Argument, orders>::type...>;
};

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
      typename NthDerivatives<
          Value,
          Argument,
          std::make_integer_sequence<int, degree + 1>>::type;

  PolynomialInMonomialBasis(Coefficients const& coeffients,
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
