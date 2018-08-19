#pragma once

#include "numerics/polynomial.hpp"

#include <tuple>
#include <utility>

#include "base/not_constructible.hpp"
#include "geometry/serialization.hpp"
#include "quantities/tuples.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

using base::make_not_null_unique;
using base::not_constructible;
using geometry::DoubleOrQuantityOrMultivectorSerializer;
using quantities::Apply;

template<typename Scalar,
         typename Tuple,
         typename = std::make_integer_sequence<int, std::tuple_size_v<Tuple>>>
struct TupleArithmetic;

template<typename Scalar, typename Tuple, int... indices>
struct TupleArithmetic<Scalar, Tuple, std::integer_sequence<int, indices...>>
    : not_constructible {
  template<typename T>
  using ScalarLeftMultiplier = Product<Scalar, T>;
  template<typename T>
  using ScalarRightMultiplier = Product<T, Scalar>;
  template<typename T>
  using ScalarRightDivider = Quotient<T, Scalar>;

  static constexpr Tuple Add(Tuple const& left, Tuple const& right);
  static constexpr Tuple Subtract(Tuple const& left, Tuple const& right);

  static constexpr typename Apply<ScalarLeftMultiplier, Tuple>
  Multiply(Scalar const& left, Tuple const& right);
  static constexpr typename Apply<ScalarRightMultiplier, Tuple>
  Multiply(Tuple const& left, Scalar const& right);

  static constexpr typename Apply<ScalarRightDivider, Tuple>
  Divide(Tuple const& left, Scalar const& right);
};

template<typename Scalar, typename Tuple, int... indices>
constexpr Tuple
TupleArithmetic<Scalar, Tuple, std::integer_sequence<int, indices...>>::Add(
    Tuple const& left, Tuple const& right) {
  return {std::get<indices>(left) + std::get<indices>(right)...};
}

template<typename Scalar, typename Tuple, int... indices>
constexpr Tuple
TupleArithmetic<Scalar, Tuple, std::integer_sequence<int, indices...>>::
    Subtract(Tuple const& left, Tuple const& right) {
  return {std::get<indices>(left) - std::get<indices>(right)...};
}

template<typename Scalar, typename Tuple, int... indices>
constexpr auto
TupleArithmetic<Scalar, Tuple, std::integer_sequence<int, indices...>>::
    Multiply(Scalar const& left, Tuple const& right) ->
    typename Apply<ScalarLeftMultiplier, Tuple> {
  return {left * std::get<indices>(right)...};
}

template<typename Scalar, typename Tuple, int... indices>
constexpr auto
TupleArithmetic<Scalar, Tuple, std::integer_sequence<int, indices...>>::
    Multiply(Tuple const& left, Scalar const& right) ->
    typename Apply<ScalarRightMultiplier, Tuple> {
  return {std::get<indices>(left) * right...};
}

template<typename Scalar, typename Tuple, int... indices>
constexpr auto
TupleArithmetic<Scalar, Tuple, std::integer_sequence<int, indices...>>::
    Divide(Tuple const& left, Scalar const& right) ->
    typename Apply<ScalarRightDivider, Tuple> {
  return {std::get<indices>(left) / right...};
}

template<typename Tuple, int k, int size = std::tuple_size_v<Tuple>>
struct TupleSerializer : not_constructible {
  static void WriteToMessage(
      Tuple const& tuple,
      not_null<serialization::PolynomialInMonomialBasis*> message);
  static void FillFromMessage(
      serialization::PolynomialInMonomialBasis const& message,
      Tuple& tuple);
};

template<typename Tuple, int size>
struct TupleSerializer<Tuple, size, size> : not_constructible {
  static void WriteToMessage(
      Tuple const& tuple,
      not_null<serialization::PolynomialInMonomialBasis*> message);
  static void FillFromMessage(
      serialization::PolynomialInMonomialBasis const& message,
      Tuple& tuple);
};

template<typename Tuple, int k, int size>
void TupleSerializer<Tuple, k, size>::WriteToMessage(
    Tuple const& tuple,
    not_null<serialization::PolynomialInMonomialBasis*> message) {
  DoubleOrQuantityOrMultivectorSerializer<
      std::tuple_element_t<k, Tuple>,
      serialization::PolynomialInMonomialBasis::Coefficient>::
      WriteToMessage(std::get<k>(tuple), message->add_coefficient());
  TupleSerializer<Tuple, k + 1, size>::WriteToMessage(tuple, message);
}

template<typename Tuple, int k, int size>
void TupleSerializer<Tuple, k, size>::FillFromMessage(
    serialization::PolynomialInMonomialBasis const& message,
    Tuple& tuple) {
  std::get<k>(tuple) =
      DoubleOrQuantityOrMultivectorSerializer<
          std::tuple_element_t<k, Tuple>,
          serialization::PolynomialInMonomialBasis::Coefficient>::
          ReadFromMessage(message.coefficient(k));
  TupleSerializer<Tuple, k + 1, size>::FillFromMessage(message, tuple);
}

template<typename Tuple, int size>
void TupleSerializer<Tuple, size, size>::WriteToMessage(
    Tuple const& tuple,
    not_null<serialization::PolynomialInMonomialBasis*> message) {}

template<typename Tuple, int size>
void TupleSerializer<Tuple, size, size>::FillFromMessage(
    serialization::PolynomialInMonomialBasis const& message,
    Tuple& tuple) {}


#define PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(value)                  \
  case value:                                                          \
    return make_not_null_unique<                                       \
        PolynomialInMonomialBasis<Value, Argument, value, Evaluator>>( \
        PolynomialInMonomialBasis<Value, Argument, value, Evaluator>:: \
            ReadFromMessage(message))

template<typename Value, typename Argument>
template<template<typename, typename, int> class Evaluator>
not_null<std::unique_ptr<Polynomial<Value, Argument>>>
Polynomial<Value, Argument>::ReadFromMessage(
    serialization::Polynomial const& message) {
  // 24 is the largest exponent that we can serialize for Quantity.
  switch (message.degree()) {
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(1);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(2);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(3);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(4);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(5);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(6);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(7);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(8);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(9);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(10);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(11);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(12);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(13);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(14);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(15);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(16);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(17);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(18);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(19);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(20);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(21);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(22);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(23);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(24);
    default:
      LOG(FATAL) << "Unexpected degree " << message.degree();
      break;
  }
}

#undef PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
PolynomialInMonomialBasis(Coefficients const& coefficients)
    : coefficients_(coefficients) {}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
Value PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
Evaluate(Argument const& argument) const {
  return Evaluator<Value, Argument, degree_>::Evaluate(coefficients_, argument);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
Derivative<Value, Argument>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
EvaluateDerivative(Argument const& argument) const {
  return Evaluator<Value, Argument, degree_>::EvaluateDerivative(
      coefficients_, argument);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
constexpr int
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::degree() const {
  return degree_;
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
void PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
    WriteToMessage(not_null<serialization::Polynomial*> message) const {
  message->set_degree(degree_);
  auto* const extension =
      message->MutableExtension(
          serialization::PolynomialInMonomialBasis::extension);
  TupleSerializer<Coefficients, 0>::WriteToMessage(coefficients_, extension);
  // No |origin|.
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::ReadFromMessage(
    serialization::Polynomial const& message) {
  CHECK_EQ(degree_, message.degree()) << message.DebugString();
  CHECK(message.HasExtension(
           serialization::PolynomialInMonomialBasis::extension))
      << message.DebugString();
  auto const& extension =
      message.GetExtension(
          serialization::PolynomialInMonomialBasis::extension);
  Coefficients coefficients;
  TupleSerializer<Coefficients, 0>::FillFromMessage(extension, coefficients);
  CHECK(!extension.has_origin()) << message.DebugString();
  return PolynomialInMonomialBasis(coefficients);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
PolynomialInMonomialBasis(Coefficients const& coefficients,
                          Point<Argument> const& origin)
    : coefficients_(coefficients),
      origin_(origin) {}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
Value PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
Evaluate(Point<Argument> const& argument) const {
  return Evaluator<Value, Argument, degree_>::Evaluate(
      coefficients_, argument - origin_);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
Derivative<Value, Argument>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
EvaluateDerivative(Point<Argument> const& argument) const {
  return Evaluator<Value, Argument, degree_>::EvaluateDerivative(
      coefficients_, argument - origin_);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
constexpr int
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
degree() const {
  return degree_;
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
void PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
    WriteToMessage(not_null<serialization::Polynomial*> message) const {
  message->set_degree(degree_);
  auto* const extension =
      message->MutableExtension(
          serialization::PolynomialInMonomialBasis::extension);
  TupleSerializer<Coefficients, 0>::WriteToMessage(coefficients_, extension);
  origin_.WriteToMessage(extension->mutable_origin());
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
ReadFromMessage(serialization::Polynomial const& message) {
  CHECK_EQ(degree_, message.degree()) << message.DebugString();
  CHECK(message.HasExtension(
           serialization::PolynomialInMonomialBasis::extension))
      << message.DebugString();
  auto const& extension =
      message.GetExtension(
          serialization::PolynomialInMonomialBasis::extension);
  Coefficients coefficients;
  TupleSerializer<Coefficients, 0>::FillFromMessage(extension, coefficients);
  auto const origin = Point<Argument>::ReadFromMessage(extension.origin());
  return PolynomialInMonomialBasis(coefficients, origin);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> operator+(
    PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const& left,
    PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
        right) {
  return PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>(
      TupleArithmetic<double,
                      typename PolynomialInMonomialBasis<
                                   Value, Argument, degree_, Evaluator>::
                          Coefficients>::Add(left.coefficients_,
                                             right.coefficients_));
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> operator-(
    PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const& left,
    PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
        right) {
  return PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>(
      TupleArithmetic<double,
                      typename PolynomialInMonomialBasis<
                                   Value, Argument, degree_, Evaluator>::
                          Coefficients>::Subtract(left.coefficients_,
                                                  right.coefficients_));
}

template<typename Scalar,
         typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Product<Scalar, Value>, Argument, degree_, Evaluator>
operator*(Scalar const& left,
          PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
              right) {
  return PolynomialInMonomialBasis<
             Product<Scalar, Value>, Argument, degree_, Evaluator>(
      TupleArithmetic<Scalar,
                      typename PolynomialInMonomialBasis<
                                   Value, Argument, degree_, Evaluator>::
                          Coefficients>::Multiply(left, right.coefficients_));
}

template<typename Scalar,
         typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Product<Value, Scalar>, Argument, degree_, Evaluator>
operator*(PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
              left,
          Scalar const& right) {
  return PolynomialInMonomialBasis<
             Product<Value, Scalar>, Argument, degree_, Evaluator>(
      TupleArithmetic<Scalar,
                      typename PolynomialInMonomialBasis<
                                   Value, Argument, degree_, Evaluator>::
                          Coefficients>::Multiply(left.coefficients_, right));
}

template<typename Scalar,
         typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Quotient<Scalar, Value>, Argument, degree_, Evaluator>
operator/(PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
              left,
          Scalar const& right) {
  return PolynomialInMonomialBasis<
             Quotient<Scalar, Value>, Argument, degree_, Evaluator>(
      TupleArithmetic<Scalar,
                      typename PolynomialInMonomialBasis<
                                   Value, Argument, degree_, Evaluator>::
                          Coefficients>::Divide(left.coefficients_, right));
}

}  // namespace internal_polynomial
}  // namespace numerics
}  // namespace principia
