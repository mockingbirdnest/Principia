#ifndef PRINCIPIA_NUMERICS_POLYNOMIAL_HPP_
#define PRINCIPIA_NUMERICS_POLYNOMIAL_HPP_

#include <memory>
#include <string>

#include "base/algebra.hpp"
#include "base/macros.hpp"  // 🧙 For forward declarations.
#include "base/not_null.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
namespace _polynomial {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::base::_not_null;

template<affine Value_, affine Argument_>
  requires homogeneous_affine_module<Value_, Difference<Argument_>>
class Polynomial {
 public:
  using Argument = Argument_;
  using Value = Value_;

  // This virtual destructor makes this class and its subclasses non-literal, so
  // constexpr-ness is a bit of a lie for polynomials.
  // TODO(phl): Consider providing an explicit deleter function that would allow
  // making the destructor protected and nonvirtual.
  virtual ~Polynomial() = default;

  friend constexpr bool operator==(Polynomial const& left,
                                   Polynomial const& right) = default;
  friend constexpr bool operator!=(Polynomial const& left,
                                   Polynomial const& right) = default;

  Value PRINCIPIA_VECTORCALL operator()(Argument argument) const;
  Derivative<Value, Argument> PRINCIPIA_VECTORCALL EvaluateDerivative(
      Argument argument) const;
  void PRINCIPIA_VECTORCALL EvaluateWithDerivative(
      Argument argument,
      Value& value,
      Derivative<Value, Argument>& derivative) const;

  // Only useful for benchmarking, analyzing performance or for downcasting.  Do
  // not use in other circumstances.
  virtual int degree() const = 0;

  // Only useful for logging.  Do not use in real code.
  virtual bool is_zero() const = 0;

  virtual void WriteToMessage(
      not_null<serialization::Polynomial*> message) const = 0;

  static not_null<std::unique_ptr<Polynomial>> ReadFromMessage(
      serialization::Polynomial const& message);
  // Compatibility deserialization, when the evaluator is not present in
  // `message`.
  template<template<typename, typename, int> typename Evaluator>
  static not_null<std::unique_ptr<Polynomial>> ReadFromMessage(
      serialization::Polynomial const& message);

 protected:
  virtual Value PRINCIPIA_VECTORCALL
  VirtualEvaluate(Argument argument) const = 0;
  virtual Derivative<Value, Argument> PRINCIPIA_VECTORCALL
  VirtualEvaluateDerivative(Argument argument) const = 0;
  virtual void PRINCIPIA_VECTORCALL VirtualEvaluateWithDerivative(
      Argument argument,
      Value& value,
      Derivative<Value, Argument>& derivative) const = 0;
};

}  // namespace internal

using internal::Polynomial;

}  // namespace _polynomial
}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_body.hpp"

#endif  // PRINCIPIA_NUMERICS_POLYNOMIAL_HPP_
