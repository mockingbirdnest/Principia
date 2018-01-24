
#pragma once

#include <tuple>
#include <utility>

#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

using quantities::Derivative;

// TODO(phl): We would like to define NthDerivative in named_quantities.hpp
// thus:
//
//   template<typename Value, typename Argument, int order>
//   using NthDerivative = typename std::conditional_t<
//       order == 0,
//       Value,
//       Quotient<Difference<Value>,
//                Exponentiation<Difference<Argument>, order>>>;
//
//   template<typename Value, typename Argument>
//   using Derivative = NthDerivative<Value, Argument, 1>;
//
// Unfortunately VS2015 is buggy and this interacts poorly with the
// std::integer_sequence below (we get the wrong types).  Revisit once MSFT has
// fixed their bugs.

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
  using Type = std::tuple<NthDerivative<Value, Argument, orders>...>;
};

template<typename Value, typename Argument, typename Sequence>
using NthDerivatives =
    typename NthDerivativesGenerator<Value, Argument, Sequence>::Type;

// |Value| must belong to an affine space.  |Argument| must belong to a vector
// space.
template<typename Value, typename Argument>
class Polynomial {
 public:
  virtual Value Evaluate(Argument const& argument) const = 0;
  virtual Derivative<Value, Argument> EvaluateDerivative(
      Argument const& argument) const = 0;

  // Only useful for benchmarking or analyzing performance.  Do not use in real
  // code.
  virtual int degree() const = 0;

 protected:
  virtual ~Polynomial() = default;
};

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
class PolynomialInMonomialBasis : public Polynomial<Value, Argument> {
 public:
  // Equivalent to:
  //   std::tuple<Value,
  //              Derivative<Value, Argument>,
  //              Derivative<Derivative<Value, Argument>>...>
  using Coefficients =
      NthDerivatives<Value,
                     Argument,
                     std::make_integer_sequence<int, degree_ + 1>>;

  explicit PolynomialInMonomialBasis(Coefficients const& coefficients);

  FORCE_INLINE(inline) Value Evaluate(Argument const& argument) const override;
  FORCE_INLINE(inline) Derivative<Value, Argument> EvaluateDerivative(
      Argument const& argument) const override;

  constexpr int degree() const override;

 private:
  Coefficients coefficients_;
};

}  // namespace internal_polynomial

using internal_polynomial::NthDerivative;
using internal_polynomial::Polynomial;
using internal_polynomial::PolynomialInMonomialBasis;

}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_body.hpp"
