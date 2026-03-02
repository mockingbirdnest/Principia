#pragma once

#include "numerics/polynomial_in_monomial_basis.hpp"

#include <algorithm>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "base/concepts.hpp"
#include "base/for_all_of.hpp"
#include "base/not_constructible.hpp"
#include "base/tags.hpp"
#include "geometry/direct_sum.hpp"
#include "geometry/serialization.hpp"
#include "google/protobuf/util/message_differencer.h"
#include "numerics/combinatorics.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/quadrature.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _polynomial_in_monomial_basis {
namespace internal {

using namespace principia::base::_concepts;
using namespace principia::base::_for_all_of;
using namespace principia::base::_not_constructible;
using namespace principia::base::_tags;
using namespace principia::geometry::_direct_sum;
using namespace principia::geometry::_serialization;
using namespace principia::numerics::_combinatorics;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_quadrature;
using namespace principia::quantities::_quantities;

// A helper for changing the origin of a monomial (x - x₁)ⁿ.  It computes the
// coefficients of the same monomial as a polynomial of (x - x₂), i.e.:
//  cₙ(x - x₁)ⁿ = cₙ((x - x₂) + (x₁ - x₂))ⁿ =
//                Σ cₙ(n k)(x - x₂)ᵏ(x₁ - x₂)ⁿ⁻ᵏ
// where (n k) is the binomial coefficient.  The coefficients are for a
// polynomial of the given degree, with zeros for the unneeded high-degree
// terms.
template<typename Value, typename Argument, int degree, int n,
         typename = std::make_index_sequence<degree + 1>>
struct MonomialAtOrigin;

template<typename Value, typename Argument, int degree, int n, std::size_t... k>
struct MonomialAtOrigin<Value, Argument, degree, n, std::index_sequence<k...>> {
  using Coefficients =
      Derivatives<Value, Argument, degree + 1>;

  // The parameter coefficient is the coefficient of the monomial.  The
  // parameter shift is x₁ - x₂, computed only once by the caller.
  static Coefficients MakeCoefficients(
      std::tuple_element_t<n, Coefficients> const& coefficient,
      Difference<Argument> const& shift);
};

template<typename Value, typename Argument, int degree, int n, std::size_t... k>
auto MonomialAtOrigin<Value, Argument, degree, n,
                      std::index_sequence<k...>>::MakeCoefficients(
    std::tuple_element_t<n, Coefficients> const& coefficient,
    Difference<Argument> const& shift) -> Coefficients {
  return {(k <= n ? coefficient * Binomial(n, k) *
                        Pow<static_cast<int>(n - k)>(shift)
                  : std::tuple_element_t<k, Coefficients>{})...};
}

// A helper for changing the origin of an entire polynomial, by repeatedly
// using `MonomialAtOrigin`.  We need two helpers because changing the origin is
// a quadratic operation in terms of the degree.
template<typename Value, typename Argument, int degree,
         template<typename, typename, int> typename Evaluator_,
         typename = std::make_index_sequence<degree + 1>>
struct PolynomialAtOrigin;

template<typename Value, typename Argument, int degree,
         template<typename, typename, int> typename Evaluator_,
         std::size_t... indices>
struct PolynomialAtOrigin<Value, Argument, degree, Evaluator_,
                          std::index_sequence<indices...>> {
  using Polynomial =
      PolynomialInMonomialBasis<Value, Argument, degree, Evaluator_>;

  static Polynomial MakePolynomial(
      typename Polynomial::Coefficients const& coefficients,
      Argument const& from_origin,
      Argument const& to_origin);
};

template<typename Value, typename Argument, int degree,
         template<typename, typename, int> typename Evaluator_,
         std::size_t ...indices>
auto PolynomialAtOrigin<Value, Argument, degree, Evaluator_,
                        std::index_sequence<indices...>>::
MakePolynomial(typename Polynomial::Coefficients const& coefficients,
               Argument const& from_origin,
               Argument const& to_origin) -> Polynomial {
  Difference<Argument> const shift = to_origin - from_origin;
  std::array<typename Polynomial::Coefficients, degree + 1> const
      all_coefficients{
          MonomialAtOrigin<Value, Argument, degree, indices>::
              MakeCoefficients(get<indices>(coefficients), shift)...};

  // It would be nicer to compute the sum using a fold expression, but Clang
  // refuses to find the operator + in that context.  Fold expressions, the
  // final frontier...
  typename Polynomial::Coefficients sum_coefficients;
  for (auto const& coefficients : all_coefficients) {
    sum_coefficients = sum_coefficients + coefficients;
  }
  return Polynomial(sum_coefficients, to_origin);
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
  return DirectSum{FallingFactorial(order + indices, order) *
                   get<order + indices>(tuple)...};
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
  constexpr auto zero = Primitive<std::tuple_element_t<0, Tuple>, Argument>{};
  return DirectSum{
      zero, get<indices>(tuple) / static_cast<double>(indices + 1)...};
}

// Helper for composition.
template<int exponent, typename P>
constexpr auto Pow(P const& polynomial) {
  static_assert(exponent > 0, "Cannot raise a polynomial to the zero-th power");
  if constexpr (exponent == 1) {
    return polynomial;
  } else if constexpr (exponent % 2 == 0) {
    return Pow<exponent / 2>(polynomial * polynomial);
  } else {
    return Pow<exponent / 2>(polynomial * polynomial) * polynomial;
  }
}

template<typename Tuple, int k, int size = std::tuple_size_v<Tuple>>
struct TupleSerializer : not_constructible {
  static void WriteToMessage(
      Tuple const& tuple,
      not_null<serialization::PolynomialInMonomialBasis*> message);
  static void FillFromMessage(
      serialization::PolynomialInMonomialBasis const& message,
      Tuple& tuple);
  // Cannot be called DebugString because of ADL in its body.
  static std::vector<std::string> TupleDebugString(Tuple const& tuple,
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
  // Cannot be called DebugString because of ADL in its body.
  static std::vector<std::string> TupleDebugString(Tuple const& tuple,
                                                   std::string const& argument);
};

template<typename Tuple, int k, int size>
void TupleSerializer<Tuple, k, size>::WriteToMessage(
    Tuple const& tuple,
    not_null<serialization::PolynomialInMonomialBasis*> message) {
  DoubleOrQuantityOrPointOrMultivectorSerializer<
      std::tuple_element_t<k, Tuple>,
      serialization::PolynomialInMonomialBasis::Coefficient>::
      WriteToMessage(get<k>(tuple), message->add_coefficient());
  TupleSerializer<Tuple, k + 1, size>::WriteToMessage(tuple, message);
}

template<typename Tuple, int k, int size>
void TupleSerializer<Tuple, k, size>::FillFromMessage(
    serialization::PolynomialInMonomialBasis const& message,
    Tuple& tuple) {
  get<k>(tuple) =
      DoubleOrQuantityOrPointOrMultivectorSerializer<
          std::tuple_element_t<k, Tuple>,
          serialization::PolynomialInMonomialBasis::Coefficient>::
          ReadFromMessage(message.coefficient(k));
  TupleSerializer<Tuple, k + 1, size>::FillFromMessage(message, tuple);
}

template<typename Tuple, int k, int size>
std::vector<std::string> TupleSerializer<Tuple, k, size>::TupleDebugString(
    Tuple const& tuple,
    std::string const& argument) {
  auto tail =
      TupleSerializer<Tuple, k + 1, size>::TupleDebugString(tuple, argument);
  auto const coefficient = get<k>(tuple);
  if (coefficient == std::tuple_element_t<k, Tuple>{}) {
    return tail;
  }
  std::string head;
  switch (k) {
    case 0:
      head = DebugString(coefficient);
      break;
    case 1:
      head = absl::StrCat(DebugString(coefficient), " * ", argument);
      break;
    default:
      head = absl::StrCat(DebugString(coefficient), " * ", argument, "^", k);
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
std::vector<std::string> TupleSerializer<Tuple, size, size>::TupleDebugString(
    Tuple const& tuple,
    std::string const& argument) {
  return {};
}


template<typename Value, typename Argument, int degree,
         template<typename, typename, int> typename Evaluator_>
not_null<std::unique_ptr<Polynomial<Value, Argument>>>
Policy::WithEvaluator(
    PolynomialInMonomialBasis<Value, Argument, degree, Evaluator_>&& polynomial)
    const {
  switch (kind_) {
    case serialization::PolynomialInMonomialBasis::Policy::
        ALWAYS_ESTRIN_WITHOUT_FMA: {
      auto result = static_cast<
          PolynomialInMonomialBasis<Value, Argument, degree, EstrinWithoutFMA>>(
          std::move(polynomial));
      return make_not_null_unique<decltype(result)>(std::move(result));
    }
    case serialization::PolynomialInMonomialBasis::Policy::ALWAYS_ESTRIN: {
      auto result = static_cast<
          PolynomialInMonomialBasis<Value, Argument, degree, Estrin>>(
          std::move(polynomial));
      return make_not_null_unique<decltype(result)>(std::move(result));
    }
  }
  LOG(FATAL) << "Unexpected policy " << kind_;
}

inline constexpr Policy Policy::AlwaysEstrin() {
  return Policy(
      serialization::PolynomialInMonomialBasis::Policy::ALWAYS_ESTRIN);
}

inline constexpr Policy Policy::AlwaysEstrinWithoutFMA() {
  return Policy(
      serialization::PolynomialInMonomialBasis::Policy::
                    ALWAYS_ESTRIN_WITHOUT_FMA);
}

inline void Policy::WriteToMessage(
    not_null<serialization::PolynomialInMonomialBasis::Policy*> message) const {
  message->set_kind(kind_);
}

inline Policy Policy::ReadFromMessage(
    serialization::PolynomialInMonomialBasis::Policy const& message) {
  return Policy(message.kind());
}

inline constexpr Policy::Policy(
    serialization::PolynomialInMonomialBasis::Policy::Kind const kind)
    : kind_(kind) {}


template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
PolynomialInMonomialBasis(Coefficients coefficients,
                          Argument const& origin)
    : coefficients_(std::move(coefficients)),
      origin_(origin) {}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
PolynomialInMonomialBasis(Coefficients coefficients)
  requires additive_group<Argument>
    : coefficients_(std::move(coefficients)),
      origin_(Argument{}) {}
template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
PolynomialInMonomialBasis()
  requires additive_group<Value>
    : coefficients_(Coefficients{}),
      origin_(Argument{}) {}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
template<int higher_degree_, typename Self>
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
operator PolynomialInMonomialBasis<Value_, Argument_, higher_degree_,
                                   Evaluator_>(this Self&& self) {
  static_assert(degree_ <= higher_degree_);
  using Result =
      PolynomialInMonomialBasis<Value, Argument, higher_degree_, Evaluator_>;
  typename Result::Coefficients higher_coefficients(uninitialized);
  for_all_of(higher_coefficients).loop_indexed([&]<int i>(auto& coefficient) {
    if constexpr (i <= degree_) {
      coefficient = std::move(get<i>(std::forward<Self>(self).coefficients_));
    } else {
      coefficient = {};
    }
  });
  return Result(higher_coefficients,
                std::move(std::forward<Self>(self).origin_));
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
template<template<typename, typename, int> typename OtherEvaluator,
         typename Self>
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
operator PolynomialInMonomialBasis<Value_, Argument_, degree_,
                                   OtherEvaluator>(this Self&& self) {
  return PolynomialInMonomialBasis<Value, Argument, degree_, OtherEvaluator>(
      std::move(std::forward<Self>(self).coefficients_),
      std::move(std::forward<Self>(self).origin_));
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>&
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
operator+=(PolynomialInMonomialBasis const& right) {
  *this = *this + right;
  return *this;
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>&
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
operator-=(PolynomialInMonomialBasis const& right) {
  *this = *this - right;
  return *this;
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
template<typename S>
constexpr PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>&
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::operator*=(
    S const& right) requires module<Value, S> {
  return *this = *this * right;
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
template<typename S>
constexpr PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>&
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::operator/=(
    S const& right) requires vector_space<Value, S> {
  return *this = *this / right;
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
Value_ PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
operator()(Argument const argument) const {
  return Evaluator<Value_, Difference<Argument_>, degree_>::Evaluate(
      coefficients_, argument - origin_);
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
Derivative<Value_, Argument_>
 PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
EvaluateDerivative(Argument const argument) const {
  return Evaluator<Value_, Difference<Argument_>, degree_>::EvaluateDerivative(
      coefficients_, argument - origin_);
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
void PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
EvaluateWithDerivative(Argument const argument,
                       Value& value,
                       base::_algebra::Derivative<Value, Argument>&
                           derivative) const {
  Evaluator<Value_, Difference<Argument_>, degree_>::EvaluateWithDerivative(
      coefficients_, argument - origin_, value, derivative);
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
constexpr
int PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
degree() const {
  return degree_;
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
bool PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
is_zero() const {
  return coefficients_ == Coefficients{};
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
typename PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
Coefficients const&
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
coefficients() const {
  return coefficients_;
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
Argument_ const&
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
origin() const {
  return origin_;
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
AtOrigin(Argument const& origin) const {
  return PolynomialAtOrigin<Value, Argument, degree_, Evaluator_>::
      MakePolynomial(coefficients_,
                     /*from_origin=*/origin_,
                     /*to_origin=*/origin);
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
template<int order>
PolynomialInMonomialBasis<
    Derivative<Value_, Argument_, order>, Argument_, degree_ - order,
    Evaluator_>
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
Derivative() const {
  return PolynomialInMonomialBasis<
             base::_algebra::Derivative<Value, Argument, order>,
             Argument,
             degree_ - order, Evaluator_>(
             TupleDerivation<Coefficients, order>::Derive(coefficients_),
             origin_);
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
PolynomialInMonomialBasis<Primitive<Value_, Argument_>, Argument_, degree_ + 1,
                          Evaluator_>
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
Primitive() const
  requires additive_group<Value> {
  return PolynomialInMonomialBasis<
             base::_algebra::Primitive<Value, Argument>,
             Argument,
             degree_ + 1, Evaluator_>(
             TupleIntegration<Argument, Coefficients>::Integrate(coefficients_),
             origin_);
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
Primitive<Value_, Argument_>
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
Integrate(Argument const& argument1,
          Argument const& argument2) const
  requires additive_group<Value> {
  // + 2 is to take into account the truncation resulting from integer division.
  return _quadrature::GaussLegendre<(degree_ + 2) / 2>(*this,
                                                       argument1, argument2);
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
void PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
    WriteToMessage(not_null<serialization::Polynomial*> message) const {
  using ArgumentSerializer = DoubleOrQuantityOrPointOrMultivectorSerializer<
      Argument,
      serialization::PolynomialInMonomialBasis>;
  if constexpr (serializable<ArgumentSerializer>) {
    message->set_degree(degree_);
    auto* const extension = message->MutableExtension(
        serialization::PolynomialInMonomialBasis::extension);
    TupleSerializer<Coefficients, 0>::WriteToMessage(coefficients_, extension);
    ArgumentSerializer::WriteToMessage(origin_, extension);
    Evaluator<Value_, Difference<Argument_>, degree_>::WriteToMessage(
        extension->mutable_evaluator());
  } else {
    LOG(FATAL) << "Non serializable PolynomialInMonomialBasis";
  }
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
ReadFromMessage(serialization::Polynomial const& message) {
  CHECK_EQ(degree_, message.degree()) << message.DebugString();
  CHECK(message.HasExtension(
           serialization::PolynomialInMonomialBasis::extension))
      << message.DebugString();
  auto const& extension =
      message.GetExtension(
          serialization::PolynomialInMonomialBasis::extension);

  bool const is_pre_gröbner = extension.origin_case() ==
    serialization::PolynomialInMonomialBasis::ORIGIN_NOT_SET;
  bool const is_pre_καραθεοδωρή = !extension.has_evaluator();
  LOG_IF(WARNING, is_pre_καραθεοδωρή)
      << "Reading pre-"
      << (is_pre_gröbner ? "Gröbner" : "Καραθεοδωρή")
      << " PolynomialInMonomialBasis";

  if (!is_pre_καραθεοδωρή) {
    serialization::PolynomialInMonomialBasis::Evaluator evaluator_message;
    Evaluator<Value, Argument, degree_>::WriteToMessage(&evaluator_message);
    CHECK(::google::protobuf::util::MessageDifferencer::Equals(
        extension.evaluator(), evaluator_message))
        << "Evaluator mismatch for post-Καραθεοδωρή deserialization"
        << extension.DebugString();
  }

  Coefficients coefficients;
  TupleSerializer<Coefficients, 0>::FillFromMessage(extension, coefficients);

  auto const origin = is_pre_gröbner
                          ? Argument{}
                          : DoubleOrQuantityOrPointOrMultivectorSerializer<
                                Argument,
                                serialization::PolynomialInMonomialBasis>::
                                ReadFromMessage(extension);
  return PolynomialInMonomialBasis(coefficients, origin);
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
Value_ PRINCIPIA_VECTORCALL
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
VirtualEvaluate(Argument argument) const {
  return (*this)(argument);
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
Derivative<Value_, Argument_> PRINCIPIA_VECTORCALL
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
VirtualEvaluateDerivative(Argument argument) const {
  return EvaluateDerivative(argument);
}

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_>
void PRINCIPIA_VECTORCALL
PolynomialInMonomialBasis<Value_, Argument_, degree_, Evaluator_>::
VirtualEvaluateWithDerivative(Argument argument,
                              Value& value,
                              D<Value, Argument>& derivative) const {
  return EvaluateWithDerivative(argument, value, derivative);
}

template<additive_group Value, affine Argument, int rdegree,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree, Evaluator>
operator+(
    PolynomialInMonomialBasis<Value, Argument, rdegree, Evaluator> const&
        right) {
  return right;
}

template<additive_group Value, affine Argument, int rdegree,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree, Evaluator>
operator-(PolynomialInMonomialBasis<Value, Argument, rdegree, Evaluator> const&
              right) {
  return PolynomialInMonomialBasis<Value, Argument, rdegree, Evaluator>(
      -right.coefficients_,
      right.origin_);
}

template<additive_group Value, affine Argument, int ldegree, int rdegree,
         template<typename, typename, int> typename Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Value, Argument, std::max(ldegree, rdegree),
                          Evaluator>
operator+(
    PolynomialInMonomialBasis<Value, Argument, ldegree, Evaluator> const& left,
    PolynomialInMonomialBasis<Value, Argument, rdegree, Evaluator> const&
        right) {
  CONSTEXPR_CHECK(left.origin_ == right.origin_);
  typename PolynomialInMonomialBasis<Value, Argument,
                                     std::max(ldegree, rdegree),
                                     Evaluator>::Coefficients
  result_coefficients(uninitialized);
  for_integer_range<0, std::max(ldegree, rdegree) + 1>::loop([&]<int i> {
    if constexpr (i <= std::min(ldegree, rdegree)) {
      get<i>(result_coefficients) =
          get<i>(left.coefficients_) + get<i>(right.coefficients_);
    } else if constexpr (i <= ldegree) {
      get<i>(result_coefficients) = get<i>(left.coefficients_);
    } else {
      get<i>(result_coefficients) = get<i>(right.coefficients_);
    }
  });
  return PolynomialInMonomialBasis<Value, Argument,
                                   std::max(ldegree, rdegree), Evaluator>(
      std::move(result_coefficients),
      left.origin_);
}

template<affine Value, affine Argument, int ldegree, int rdegree,
         template<typename, typename, int> typename Evaluator>
  requires (!additive_group<Value>)  // NOLINT
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Value, Argument, std::max(ldegree, rdegree),
                          Evaluator>
operator+(PolynomialInMonomialBasis<Difference<Value>, Argument, ldegree,
                                    Evaluator> const& left,
          PolynomialInMonomialBasis<Value, Argument, rdegree,
                                    Evaluator> const& right) {
  CONSTEXPR_CHECK(left.origin_ == right.origin_);
  typename PolynomialInMonomialBasis<Value, Argument,
                                     std::max(ldegree, rdegree),
                                     Evaluator>::Coefficients
  result_coefficients(uninitialized);
  for_integer_range<0, std::max(ldegree, rdegree) + 1>::loop([&]<int i> {
    if constexpr (i <= std::min(ldegree, rdegree)) {
      get<i>(result_coefficients) =
          get<i>(left.coefficients_) + get<i>(right.coefficients_);
    } else if constexpr (i <= ldegree) {
      get<i>(result_coefficients) = get<i>(left.coefficients_);
    } else {
      get<i>(result_coefficients) = get<i>(right.coefficients_);
    }
  });
  return PolynomialInMonomialBasis<Value, Argument,
                                   std::max(ldegree, rdegree), Evaluator>(
      std::move(result_coefficients),
      left.origin_);
}

template<affine Value, affine Argument, int ldegree, int rdegree,
         template<typename, typename, int> typename Evaluator>
  requires (!additive_group<Value>)  // NOLINT
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Value, Argument, std::max(ldegree, rdegree),
                          Evaluator>
operator+(PolynomialInMonomialBasis<Value, Argument, ldegree,
                                    Evaluator> const& left,
          PolynomialInMonomialBasis<Difference<Value>, Argument, rdegree,
                                    Evaluator> const& right) {
  CONSTEXPR_CHECK(left.origin_ == right.origin_);
  typename PolynomialInMonomialBasis<Value, Argument,
                                     std::max(ldegree, rdegree),
                                     Evaluator>::Coefficients
  result_coefficients(uninitialized);
  for_integer_range<0, std::max(ldegree, rdegree) + 1>::loop([&]<int i> {
    if constexpr (i <= std::min(ldegree, rdegree)) {
      get<i>(result_coefficients) =
          get<i>(left.coefficients_) + get<i>(right.coefficients_);
    } else if constexpr (i <= ldegree) {
      get<i>(result_coefficients) = get<i>(left.coefficients_);
    } else {
      get<i>(result_coefficients) = get<i>(right.coefficients_);
    }
  });
  return PolynomialInMonomialBasis<Value, Argument,
                                   std::max(ldegree, rdegree), Evaluator>(
      std::move(result_coefficients),
      left.origin_);
}

template<affine LValue, affine RValue, affine Argument,
         int ldegree, int rdegree,
         template<typename, typename, int> typename Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Difference<LValue, RValue>, Argument,
                          std::max(ldegree, rdegree),
                          Evaluator>
operator-(
    PolynomialInMonomialBasis<LValue, Argument, ldegree, Evaluator> const& left,
    PolynomialInMonomialBasis<RValue, Argument, rdegree, Evaluator> const&
        right) {
  CONSTEXPR_CHECK(left.origin_ == right.origin_);
  typename PolynomialInMonomialBasis<Difference<LValue, RValue>, Argument,
                                     std::max(ldegree, rdegree),
                                     Evaluator>::Coefficients
  result_coefficients(uninitialized);
  for_integer_range<0, std::max(ldegree, rdegree) + 1>::loop([&]<int i> {
    if constexpr (i <= std::min(ldegree, rdegree)) {
      get<i>(result_coefficients) =
          get<i>(left.coefficients_) - get<i>(right.coefficients_);
    } else if constexpr (i <= ldegree) {
      get<i>(result_coefficients) = get<i>(left.coefficients_);
    } else {
      get<i>(result_coefficients) = -get<i>(right.coefficients_);
    }
  });
  return PolynomialInMonomialBasis<Difference<LValue, RValue>, Argument,
                                   std::max(ldegree, rdegree), Evaluator>(
      std::move(result_coefficients),
      left.origin_);
}

template<typename Scalar,
         typename Value, affine Argument, int degree,
         template<typename, typename, int> typename Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Product<Scalar, Value>, Argument, degree, Evaluator>
operator*(Scalar const& left,
          PolynomialInMonomialBasis<Value, Argument, degree, Evaluator> const&
              right) {
  return PolynomialInMonomialBasis<Product<Scalar, Value>, Argument, degree,
                                   Evaluator>(left * right.coefficients_,
                                              right.origin_);
}

template<typename Scalar,
         typename Value, affine Argument, int degree,
         template<typename, typename, int> typename Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Product<Value, Scalar>, Argument, degree, Evaluator>
operator*(
    PolynomialInMonomialBasis<Value, Argument, degree, Evaluator> const& left,
    Scalar const& right) {
  return PolynomialInMonomialBasis<Product<Value, Scalar>, Argument, degree,
                                   Evaluator>(left.coefficients_ * right,
                                              left.origin_);
}

template<typename Scalar,
         typename Value, affine Argument, int degree,
         template<typename, typename, int> typename Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Quotient<Value, Scalar>, Argument, degree, Evaluator>
operator/(
    PolynomialInMonomialBasis<Value, Argument, degree, Evaluator> const& left,
    Scalar const& right) {
  return PolynomialInMonomialBasis<Quotient<Value, Scalar>, Argument, degree,
                                   Evaluator>(left.coefficients_ / right,
                                              left.origin_);
}

template<typename LValue, typename RValue,
         affine Argument, int ldegree, int rdegree,
         template<typename, typename, int> typename Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<Product<LValue, RValue>, Argument, ldegree + rdegree,
                          Evaluator>
operator*(
    PolynomialInMonomialBasis<LValue, Argument, ldegree, Evaluator> const& left,
    PolynomialInMonomialBasis<RValue, Argument, rdegree, Evaluator> const&
        right) {
  CONSTEXPR_CHECK(left.origin_ == right.origin_);
  typename PolynomialInMonomialBasis<Product<LValue, RValue>, Argument,
                                     ldegree + rdegree,
                                     Evaluator>::Coefficients
      result_coefficients;
  for_all_of(left.coefficients_)
      .loop_indexed([&]<int i>(auto const& l) {
        for_all_of(right.coefficients_)
            .loop_indexed([&]<int j>(auto const& r) {
              get<i + j>(result_coefficients) += l * r;
            });
      });
  return PolynomialInMonomialBasis<Product<LValue, RValue>, Argument,
                                   ldegree + rdegree,
                                   Evaluator>(std::move(result_coefficients),
                                              left.origin_);
}

#if PRINCIPIA_COMPILER_MSVC_HANDLES_POLYNOMIAL_OPERATORS
template<typename Value, affine Argument, int ldegree,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, ldegree, Evaluator>
operator+(PolynomialInMonomialBasis<Difference<Value>, Argument,
                                    ldegree,
                                    Evaluator> const& left,
          Value const& right) {
#else
template<typename Value,
         std::same_as<Difference<Value>> ValueDifference, affine Argument,
         int ldegree,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, ldegree, Evaluator>
operator+(PolynomialInMonomialBasis<ValueDifference, Argument,
                                    ldegree,
                                    Evaluator> const& left,
          Value const& right) {
#endif
  typename PolynomialInMonomialBasis<Value, Argument, ldegree, Evaluator>::
      Coefficients result_coefficients(uninitialized);
  for_all_of(left.coefficients_).loop_indexed([&]<int i>(auto const& l) {
    if constexpr (i == 0) {
      get<i>(result_coefficients) = l + right;
    } else {
      get<i>(result_coefficients) = l;
    }
  });
  return PolynomialInMonomialBasis<Value, Argument, ldegree, Evaluator>(
      std::move(result_coefficients),
      left.origin_);
}

template<typename Value, affine Argument, int rdegree,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree, Evaluator>
operator+(Value const& left,
          PolynomialInMonomialBasis<Difference<Value>, Argument, rdegree,
                                    Evaluator> const& right) {
  typename PolynomialInMonomialBasis<Value, Argument, rdegree, Evaluator>::
      Coefficients result_coefficients(uninitialized);
  for_all_of(right.coefficients_).loop_indexed([&]<int i>(auto const& r) {
    if constexpr (i == 0) {
      get<i>(result_coefficients) = left + r;
    } else {
      get<i>(result_coefficients) = r;
    }
  });
  return PolynomialInMonomialBasis<Value, Argument, rdegree, Evaluator>(
      std::move(result_coefficients),
      right.origin_);
}

template<typename Value, affine Argument, int ldegree,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Difference<Value>, Argument, ldegree,
                                    Evaluator>
operator-(
    PolynomialInMonomialBasis<Value, Argument, ldegree, Evaluator> const& left,
    Value const& right) {
  typename PolynomialInMonomialBasis<Difference<Value>, Argument,
                                     ldegree,
                                     Evaluator>::Coefficients
      result_coefficients(uninitialized);
  for_all_of(left.coefficients_).loop_indexed([&]<int i>(auto const& l) {
    if constexpr (i == 0) {
      get<i>(result_coefficients) = l - right;
    } else {
      get<i>(result_coefficients) = l;
    }
  });
  return PolynomialInMonomialBasis<Difference<Value>, Argument, ldegree,
                                   Evaluator>(
      std::move(result_coefficients),
      left.origin_);
}

template<typename Value, affine Argument, int rdegree,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Difference<Value>, Argument, rdegree,
                                    Evaluator>
operator-(Value const& left,
          PolynomialInMonomialBasis<Value, Argument, rdegree, Evaluator> const&
              right) {
  typename PolynomialInMonomialBasis<Difference<Value>, Argument,
                                     rdegree,
                                     Evaluator>::Coefficients
      result_coefficients(uninitialized);
  for_all_of(right.coefficients_).loop_indexed([&]<int i>(auto const& r) {
    if constexpr (i == 0) {
      get<i>(result_coefficients) = left - r;
    } else {
      get<i>(result_coefficients) = r;
    }
  });
  return PolynomialInMonomialBasis<Difference<Value>, Argument, rdegree,
                                   Evaluator>(
      std::move(result_coefficients),
      right.origin_);
}

template<typename LValue, typename RValue,
         typename RArgument, int ldegree, int rdegree,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<LValue, RArgument, ldegree * rdegree,
                                    Evaluator>
Compose(
    PolynomialInMonomialBasis<LValue, RValue, ldegree, Evaluator> const& left,
    PolynomialInMonomialBasis<RValue, RArgument, rdegree, Evaluator> const&
        right) {
  typename PolynomialInMonomialBasis<LValue, RArgument,
                                     ldegree * rdegree,
                                     Evaluator>::Coefficients
      result_coefficients;
  for_all_of(left.coefficients_).loop_indexed([&]<int i>(auto const& l) {
    if constexpr (i == 0) {
      get<i>(result_coefficients) = l;
    } else {
      // NOTE(egg):
      // `for_all_of((l * Pow<i>(right - left.origin_)).coefficients_)...`
      // does not compile.  The temporary polynomial would outlive the loop
      // so it should be fine, but presumably we would need to be clever
      // about value categories somewhere in `for_all_of`.
      auto const left_monomial = l * Pow<i>(right - left.origin_);
      for_all_of(left_monomial.coefficients_)
          .loop_indexed([&]<int j>(auto const& c) {
            get<j>(result_coefficients) += c;
          });
    }
  });
  return PolynomialInMonomialBasis<LValue, RArgument, ldegree * rdegree,
                                   Evaluator>(
      std::move(result_coefficients),
      right.origin_);
}

template<typename LValue, typename RValue,
         affine Argument, int ldegree, int rdegree,
         template<typename, typename, int> typename Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<
    InnerProductType<LValue, RValue>, Argument,
    ldegree + rdegree, Evaluator>
PointwiseInnerProduct(
    PolynomialInMonomialBasis<LValue, Argument, ldegree,
                              Evaluator> const& left,
    PolynomialInMonomialBasis<RValue, Argument, rdegree,
                              Evaluator> const& right) {
  typename PolynomialInMonomialBasis<InnerProductType<LValue, RValue>, Argument,
                                     ldegree + rdegree,
                                     Evaluator>::Coefficients
      result_coefficients;
  CONSTEXPR_CHECK(left.origin_ == right.origin_);
  for_all_of(left.coefficients_).loop_indexed([&]<int i>(auto l) {
    for_all_of(right.coefficients_).loop_indexed([&]<int j>(auto r) {
      get<i + j>(result_coefficients) += InnerProduct(l, r);
    });
  });
  return PolynomialInMonomialBasis<InnerProductType<LValue, RValue>, Argument,
                                   ldegree + rdegree,
                                   Evaluator>(std::move(result_coefficients),
                                              left.origin_);
}

template<typename Value, affine Argument, int degree,
         template<typename, typename, int> typename Evaluator>
std::ostream& operator<<(
    std::ostream& out,
    PolynomialInMonomialBasis<Value, Argument, degree,
                              Evaluator> const& polynomial) {
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree, Evaluator>::
          Coefficients;
  std::vector<std::string> debug_string;
  debug_string = TupleSerializer<Coefficients, 0>::TupleDebugString(
      polynomial.coefficients_,
      absl::StrCat("(T - ", DebugString(polynomial.origin_), ")"));
  if (debug_string.empty()) {
    out << DebugString(Value{});
  } else {
    out << absl::StrJoin(debug_string, " + ");
  }
  return out;
}

}  // namespace internal
}  // namespace _polynomial_in_monomial_basis
}  // namespace numerics
}  // namespace principia
