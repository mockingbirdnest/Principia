#pragma once

#include "numerics/polynomial.hpp"

#include <algorithm>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "absl/strings/str_join.h"
#include "absl/strings/str_cat.h"
#include "base/not_constructible.hpp"
#include "base/traits.hpp"
#include "geometry/cartesian_product.hpp"
#include "geometry/serialization.hpp"
#include "numerics/combinatorics.hpp"
#include "polynomial.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

using base::is_instance_of_v;
using base::make_not_null_unique;
using base::not_constructible;
using geometry::DoubleOrQuantityOrMultivectorSerializer;
using geometry::cartesian_product::operator+;
using geometry::cartesian_product::operator-;
using geometry::cartesian_product::operator*;
using geometry::cartesian_product::operator/;
using geometry::polynomial_ring::operator*;
using quantities::Apply;
using quantities::Exponentiation;
using quantities::Pow;
using quantities::Time;

// Index-by-index assignment of RTuple to LTuple, which must have at least as
// many elements (and the types must match).
template<typename LTuple, typename RTuple,
         typename = std::make_index_sequence<std::tuple_size_v<RTuple>>>
struct TupleAssigner;

template<typename LTuple, typename RTuple, std::size_t... indices>
struct TupleAssigner<LTuple, RTuple, std::index_sequence<indices...>> {
  static void Assign(LTuple& left_tuple, RTuple const& right_tuple);
};

template<typename LTuple, typename RTuple, std::size_t... indices>
void TupleAssigner<LTuple, RTuple, std::index_sequence<indices...>>::Assign(
    LTuple& left_tuple,
    RTuple const& right_tuple) {
  // This fold expression effectively implements repeated assignments.
  ((std::get<indices>(left_tuple) = std::get<indices>(right_tuple)), ...);
}

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

template<typename Argument, typename Tuple,
         typename = std::make_index_sequence<std::tuple_size_v<Tuple>>>
struct TupleIntegration;

template<typename Argument, typename Tuple, std::size_t... indices>
struct TupleIntegration<Argument, Tuple, std::index_sequence<indices...>> {
  static constexpr auto Integrate(Tuple const& tuple);
};

template<typename Argument, typename Tuple, std::size_t... indices>
constexpr auto
TupleIntegration<Argument, Tuple, std::index_sequence<indices...>>::Integrate(
    Tuple const& tuple) {
  constexpr auto zero = std::tuple_element_t<0, Tuple>{} * Argument{};
  return std::make_tuple(
      zero, std::get<indices>(tuple) / static_cast<double>(indices + 1)...);
}

template<typename Tuple, int k, int size = std::tuple_size_v<Tuple>>
struct TupleSerializer : not_constructible {
  static void WriteToMessage(
      Tuple const& tuple,
      not_null<serialization::PolynomialInMonomialBasis*> message);
  static void FillFromMessage(
      serialization::PolynomialInMonomialBasis const& message,
      Tuple& tuple);
  static std::vector<std::string> DebugString(Tuple const& tuple,
                                              std::string const& argument);
};

template<typename Tuple, int size>
struct TupleSerializer<Tuple, size, size> : not_constructible {
  static void WriteToMessage(
      Tuple const& tuple,
      not_null<serialization::PolynomialInMonomialBasis*> message);
  static void FillFromMessage(
      serialization::PolynomialInMonomialBasis const& message,
      Tuple& tuple);
  static std::vector<std::string> DebugString(Tuple const& tuple,
                                              std::string const& argument);
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

template<typename Tuple, int k, int size>
std::vector<std::string> TupleSerializer<Tuple, k, size>::DebugString(
    Tuple const& tuple,
    std::string const& argument) {
  auto tail =
      TupleSerializer<Tuple, k + 1, size>::DebugString(tuple, argument);
  auto const coefficient = std::get<k>(tuple);
  if (coefficient == std::tuple_element_t<k, Tuple>{}) {
    return tail;
  }
  std::string head;
  switch (k) {
    case 0:
      head = quantities::DebugString(coefficient);
      break;
    case 1:
      head = absl::StrCat(
        quantities::DebugString(coefficient), " * ", argument);
      break;
    default:
      head = absl::StrCat(
          quantities::DebugString(coefficient), " * ", argument, "^", k);
      break;
  }
  tail.insert(tail.begin(), head);
  return tail;
}

template<typename Tuple, int size>
void TupleSerializer<Tuple, size, size>::WriteToMessage(
    Tuple const& tuple,
    not_null<serialization::PolynomialInMonomialBasis*> message) {}

template<typename Tuple, int size>
void TupleSerializer<Tuple, size, size>::FillFromMessage(
    serialization::PolynomialInMonomialBasis const& message,
    Tuple& tuple) {}

template<typename Tuple, int size>
std::vector<std::string> TupleSerializer<Tuple, size, size>::DebugString(
    Tuple const& tuple,
    std::string const& argument) {
  return {};
}


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
constexpr PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
PolynomialInMonomialBasis(Coefficients coefficients)
    : coefficients_(std::move(coefficients)) {}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
template<int higher_degree_,
         template<typename, typename, int> class HigherEvaluator>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
operator PolynomialInMonomialBasis<Value, Argument, higher_degree_,
                                   HigherEvaluator>() const {
  static_assert(degree_ <= higher_degree_);
  using Result = PolynomialInMonomialBasis<Value, Argument, higher_degree_,
                                           HigherEvaluator>;
  typename Result::Coefficients higher_coefficients;
  TupleAssigner<typename Result::Coefficients, Coefficients>::Assign(
      higher_coefficients, coefficients_);
  return Result(higher_coefficients);
}

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
bool PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::is_zero()
    const {
  return coefficients_ == Coefficients{};
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
PolynomialInMonomialBasis<
    Primitive<Value, Argument>, Argument, degree_ + 1, Evaluator>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
Primitive() const {
  return PolynomialInMonomialBasis<
             quantities::Primitive<Value, Argument>, Argument,
             degree_ + 1, Evaluator>(
             TupleIntegration<Argument, Coefficients>::
                Integrate(coefficients_));
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>&
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::operator+=(
    PolynomialInMonomialBasis const& right) {
  *this = *this + right;
  return *this;
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>&
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::operator-=(
    PolynomialInMonomialBasis const& right) {
  *this = *this - right;
  return *this;
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
PolynomialInMonomialBasis(Coefficients coefficients,
                          Point<Argument> const& origin)
    : coefficients_(std::move(coefficients)),
      origin_(origin) {}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
template<int higher_degree_,
         template<typename, typename, int> class HigherEvaluator>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
operator PolynomialInMonomialBasis<Value, Point<Argument>, higher_degree_,
                                   HigherEvaluator>() const {
  static_assert(degree_ <= higher_degree_);
  using Result = PolynomialInMonomialBasis<
                     Value, Point<Argument>, higher_degree_, HigherEvaluator>;
  typename Result::Coefficients higher_coefficients;
  TupleAssigner<typename Result::Coefficients, Coefficients>::Assign(
      higher_coefficients, coefficients_);
  return Result(higher_coefficients, origin_);
}

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
bool PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
is_zero() const {
  return coefficients_ == Coefficients{};
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
Point<Argument> const&
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
origin() const {
  return origin_;
}

template<typename Value, typename Argument, int degree, int n,
         template<typename, typename, int> class Evaluator,
         typename = std::make_index_sequence<degree + 1>>
struct Bino;

template<typename Value, typename Argument, int degree, int n,
         template<typename, typename, int> class Evaluator,
         std::size_t... k>
struct Bino<Value, Argument, degree, n, Evaluator, std::index_sequence<k...>> {
  using Coefficients =
      typename
      PolynomialInMonomialBasis<Value, Point<Argument>, degree, Evaluator>::
          Coefficients;
  static typename Coefficients Make(
      std::tuple_element_t<n, Coefficients> const& coefficient,
      Argument const& shift) {
    return {(k <= n ? coefficient * Binomial(n, k) *
                          Pow<static_cast<int>(n - k)>(shift)
                    : std::tuple_element_t<k, Coefficients>{})...};
  }
};

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator,
         typename = std::make_index_sequence<degree_ + 1>>
struct DuDdu;

template<typename Value, typename Argument, int degree,
         template<typename, typename, int> class Evaluator,
         std::size_t... indices>
struct DuDdu<Value, Argument, degree, Evaluator,
             std::index_sequence<indices...>> {
  using Polynomial =
      PolynomialInMonomialBasis<Value, Point<Argument>, degree, Evaluator>;
  static Polynomial Make(typename Polynomial::Coefficients const& coefficients,
                         Point<Argument> const& from_origin,
                         Point<Argument> const& to_origin) {
    Argument const shift = to_origin - from_origin;
    return Polynomial(
        typename Polynomial::Coefficients{(
            Bino<Value, Argument, degree, indices, Evaluator>::Make(
                std::get<indices>(coefficients), shift) +
            ...)},
        to_origin);
  }
};

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::AtOrigin(
    Point<Argument> const& origin) const {
  // (x - x0)^n = (x - x1 + x1 - x0)^n = Sum (n|k)(x - x1)^k(x1 - x0)^(n - k)
  // Time^n
  return DuDdu<Value, Argument, degree_, Evaluator>::Make(
      coefficients_, origin_, origin);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
template<int order>
PolynomialInMonomialBasis<
    Derivative<Value, Argument, order>, Point<Argument>, degree_ - order,
    Evaluator>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
Derivative() const {
  return PolynomialInMonomialBasis<
             quantities::Derivative<Value, Argument, order>, Point<Argument>,
             degree_ - order, Evaluator>(
             TupleDerivation<Coefficients, order>::Derive(coefficients_),
             origin_);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<
    Primitive<Value, Argument>, Point<Argument>, degree_ + 1, Evaluator>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
Primitive() const {
  return PolynomialInMonomialBasis<
             quantities::Primitive<Value, Argument>, Point<Argument>,
             degree_ + 1, Evaluator>(
             TupleIntegration<Argument, Coefficients>::
                Integrate(coefficients_),
             origin_);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>&
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
operator+=(PolynomialInMonomialBasis const& right) {
  *this = *this + right;
  return *this;
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>&
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
operator-=(PolynomialInMonomialBasis const& right) {
  *this = *this - right;
  return *this;
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

template<typename Value, typename Argument, int rdegree_,
         template<typename, typename, int> class Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator>
operator+(PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator> const&
              right) {
  return right;
}

template<typename Value, typename Argument, int rdegree_,
         template<typename, typename, int> class Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator>
operator-(PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator> const&
              right) {
  if constexpr (is_instance_of_v<Point, Argument>) {
    return PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator>(
        -right.coefficients_,
        right.origin_);
  } else {
    return PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator>(
        -right.coefficients_);
  }
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
  if constexpr (is_instance_of_v<Point, Argument>) {
    CONSTEXPR_CHECK(left.origin_ == right.origin_);
    return PolynomialInMonomialBasis<Value, Argument,
                                     std::max(ldegree_, rdegree_), Evaluator>(
        left.coefficients_ + right.coefficients_,
        left.origin_);
  } else {
    return PolynomialInMonomialBasis<Value, Argument,
                                     std::max(ldegree_, rdegree_), Evaluator>(
        left.coefficients_ + right.coefficients_);
  }
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
  if constexpr (is_instance_of_v<Point, Argument>) {
    CONSTEXPR_CHECK(left.origin_ == right.origin_);
    return PolynomialInMonomialBasis<Value, Argument,
                                     std::max(ldegree_, rdegree_), Evaluator>(
        left.coefficients_ - right.coefficients_,
        left.origin_);
  } else {
    return PolynomialInMonomialBasis<Value, Argument,
                                     std::max(ldegree_, rdegree_), Evaluator>(
        left.coefficients_ - right.coefficients_);
  }
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
  if constexpr (is_instance_of_v<Point, Argument>) {
    return PolynomialInMonomialBasis<Product<Scalar, Value>, Argument, degree_,
                                     Evaluator>(left * right.coefficients_,
                                                right.origin_);
  } else {
    return PolynomialInMonomialBasis<Product<Scalar, Value>, Argument, degree_,
                                     Evaluator>(left * right.coefficients_);
  }
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
  if constexpr (is_instance_of_v<Point, Argument>) {
    return PolynomialInMonomialBasis<Product<Value, Scalar>, Argument, degree_,
                                     Evaluator>(left.coefficients_ * right,
                                                left.origin_);
  } else {
    return PolynomialInMonomialBasis<Product<Value, Scalar>, Argument, degree_,
                                     Evaluator>(left.coefficients_ * right);
  }
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
  if constexpr (is_instance_of_v<Point, Argument>) {
    return PolynomialInMonomialBasis<Quotient<Value, Scalar>, Argument, degree_,
                                     Evaluator>(left.coefficients_ / right,
                                                left.origin_);
  } else {
    return PolynomialInMonomialBasis<Quotient<Value, Scalar>, Argument, degree_,
                                     Evaluator>(left.coefficients_ / right);
  }
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
  if constexpr (is_instance_of_v<Point, Argument>) {
    CONSTEXPR_CHECK(left.origin_ == right.origin_);
    return PolynomialInMonomialBasis<Product<LValue, RValue>, Argument,
                                     ldegree_ + rdegree_, Evaluator>(
               left.coefficients_ * right.coefficients_,
               left.origin_);
  } else {
    return PolynomialInMonomialBasis<Product<LValue, RValue>, Argument,
                                     ldegree_ + rdegree_, Evaluator>(
               left.coefficients_ * right.coefficients_);
  }
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
std::ostream& operator<<(
    std::ostream& out,
    PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
        polynomial) {
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
          Coefficients;
  std::vector<std::string> debug_string;
  if constexpr (is_instance_of_v<Point, Argument>) {
    debug_string = TupleSerializer<Coefficients, 0>::DebugString(
        polynomial.coefficients_,
        "(T - " + DebugString(polynomial.origin_) + ")");
  } else {
    debug_string = TupleSerializer<Coefficients, 0>::DebugString(
        polynomial.coefficients_, "T");
  }
  if (debug_string.empty()) {
    out << quantities::DebugString(Value{});
  } else {
    out << absl::StrJoin(debug_string, " + ");
  }
  return out;
}

}  // namespace internal_polynomial
}  // namespace numerics
}  // namespace principia
