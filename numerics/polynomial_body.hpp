#pragma once

#include "numerics/polynomial.hpp"

#include <algorithm>
#include <tuple>
#include <utility>

#include "base/not_constructible.hpp"
#include "geometry/cartesian_product.hpp"
#include "geometry/serialization.hpp"
#include "numerics/combinatorics.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

using base::check_not_null;
using base::make_not_null_unique;
using base::not_constructible;
using geometry::DoubleOrQuantityOrMultivectorSerializer;
using geometry::cartesian_product::operator+;
using geometry::cartesian_product::operator-;
using geometry::cartesian_product::operator*;
using geometry::cartesian_product::operator/;
using geometry::polynomial_ring::operator*;
using quantities::Apply;

template<typename Tuple, int order,
         typename = std::make_index_sequence<std::tuple_size_v<Tuple> - order>>
struct TupleDerivation;

template<typename Tuple, int order, std::size_t... indices>
struct TupleDerivation<Tuple, order, std::index_sequence<indices...>> {
  static constexpr auto Derive(Tuple const& tuple);
};

template<typename Tuple, int order, std::size_t... indices>
constexpr auto
TupleDerivation<Tuple, order, std::index_sequence<indices...>>::Derive(
    Tuple const& tuple) {
  return std::make_tuple(FallingFactorial(order + indices, order) *
                         std::get<order + indices>(tuple)...);
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

#define PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(value)                      \
  case value:                                                              \
    return check_not_null(                                                 \
        google::protobuf::Arena::Create<                                   \
            PolynomialInMonomialBasis<Value, Argument, value, Evaluator>>( \
            &arena,                                                        \
            PolynomialInMonomialBasis<Value, Argument, value, Evaluator>:: \
                ReadFromMessage(message)))

template<typename Value, typename Argument>
template<template<typename, typename, int> class Evaluator>
not_null<Polynomial<Value, Argument>*>
Polynomial<Value, Argument>::ReadFromMessage(
    serialization::Polynomial const& message,
    google::protobuf::Arena& arena) {
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
constexpr PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
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
template<int order>
PolynomialInMonomialBasis<
    Derivative<Value, Argument, order>, Argument, degree_ - order, Evaluator>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
Derivative() const {
  return PolynomialInMonomialBasis<
             quantities::Derivative<Value, Argument, order>, Argument,
             degree_ - order, Evaluator>(
             TupleDerivation<Coefficients, order>::Derive(coefficients_));
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
constexpr
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

template<typename Value, typename Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Value, Argument,
                          std::max(ldegree_, rdegree_), Evaluator>
operator+(
    PolynomialInMonomialBasis<Value, Argument, ldegree_, Evaluator> const& left,
    PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator> const&
        right) {
  return PolynomialInMonomialBasis<Value, Argument,
                                   std::max(ldegree_, rdegree_), Evaluator>(
      left.coefficients_ + right.coefficients_);
}

template<typename Value, typename Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Value, Argument,
                          std::max(ldegree_, rdegree_), Evaluator>
operator-(
    PolynomialInMonomialBasis<Value, Argument, ldegree_, Evaluator> const& left,
    PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator> const&
        right) {
  return PolynomialInMonomialBasis<Value, Argument,
                                   std::max(ldegree_, rdegree_), Evaluator>(
      left.coefficients_ - right.coefficients_);
}

template<typename Scalar,
         typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Product<Scalar, Value>, Argument,
                          degree_, Evaluator>
operator*(Scalar const& left,
          PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
              right) {
  return PolynomialInMonomialBasis<Product<Scalar, Value>, Argument, degree_,
                                   Evaluator>(left * right.coefficients_);
}

template<typename Scalar,
         typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Product<Value, Scalar>, Argument,
                          degree_, Evaluator>
operator*(PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
              left,
          Scalar const& right) {
  return PolynomialInMonomialBasis<Product<Value, Scalar>, Argument, degree_,
                                   Evaluator>(left.coefficients_ * right);
}

template<typename Scalar,
         typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Quotient<Value, Scalar>, Argument,
                          degree_, Evaluator>
operator/(PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
              left,
          Scalar const& right) {
  return PolynomialInMonomialBasis<Quotient<Value, Scalar>, Argument, degree_,
                                   Evaluator>(left.coefficients_ / right);
}

template<typename LValue, typename RValue,
         typename Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Product<LValue, RValue>, Argument,
                          ldegree_ + rdegree_, Evaluator>
operator*(
    PolynomialInMonomialBasis<LValue, Argument, ldegree_, Evaluator> const&
        left,
    PolynomialInMonomialBasis<RValue, Argument, rdegree_, Evaluator> const&
        right) {
  return PolynomialInMonomialBasis<Product<LValue, RValue>, Argument,
                                   ldegree_ + rdegree_, Evaluator>(
             left.coefficients_ * right.coefficients_);
}

}  // namespace internal_polynomial
}  // namespace numerics
}  // namespace principia
