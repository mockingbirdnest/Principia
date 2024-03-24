#ifndef PRINCIPIA_NUMERICS_POLYNOMIAL_HPP_
#define PRINCIPIA_NUMERICS_POLYNOMIAL_HPP_

#include <memory>
#include <string>

#include "base/macros.hpp"  // 🧙 For forward declarations.
#include "base/not_null.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
namespace _polynomial {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::quantities::_named_quantities;

// |Value_| must belong to an affine space.  |Argument_| must belong to a ring
// or to Point based on a ring.
template<typename Value_, typename Argument_>
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

  virtual Value operator()(Argument const& argument) const = 0;
  virtual Derivative<Value, Argument> EvaluateDerivative(
      Argument const& argument) const = 0;

  // Only useful for benchmarking, analyzing performance or for downcasting.  Do
  // not use in other circumstances.
  virtual int degree() const = 0;

  // Only useful for logging.  Do not use in real code.
  virtual bool is_zero() const = 0;

  virtual void WriteToMessage(
      not_null<serialization::Polynomial*> message) const = 0;

  static not_null<std::unique_ptr<Polynomial>> ReadFromMessage(
      serialization::Polynomial const& message);
};

}  // namespace internal

using internal::Polynomial;

}  // namespace _polynomial
}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_body.hpp"

#endif  // PRINCIPIA_NUMERICS_POLYNOMIAL_HPP_
