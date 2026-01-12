#pragma once

#include <functional>
#include <utility>
#include <vector>

#include "base/algebra.hpp"
#include "base/array.hpp"
#include "geometry/hilbert.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"

namespace principia {
namespace numerics {
namespace _hermite3 {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::base::_array;
using namespace principia::geometry::_hilbert;
using namespace principia::numerics::_polynomial_in_monomial_basis;

// A 3rd degree Hermite polynomial defined by its values and derivatives at the
// bounds of some interval.
template<affine Value_, affine Argument_>
class Hermite3 final {
  using NormType = typename Hilbert<Difference<Value_>>::NormType;

 public:
  using Argument = Argument_;
  using Value = Value_;
  using Derivative1 = Derivative<Value, Argument>;

  Hermite3(std::pair<Argument, Argument> const& arguments,
           std::pair<Value, Value> const& values,
           std::pair<Derivative1, Derivative1> const& derivatives);

  // NOTE(egg): this does not appear to use Casteljau's algorithm; perhaps it
  // should?
  Value Evaluate(Argument const& argument) const;
  Derivative1 EvaluateDerivative(Argument const& argument) const;
  void EvaluateWithDerivative(Argument const& argument,
                              Value& value,
                              Derivative1& derivative) const;

  // The result is sorted.
  BoundedArray<Argument, 2> FindExtrema() const;

  NormType LInfinityL₁Norm() const;

  // `samples` must be a container; `get_argument` and `get_value` on the
  // its elements must return `Argument` and `Value` respectively.  If `h` is
  // this polynomial, `tᵢ` an argument from `samples` and `xᵢ` the corresponding
  // value, this function returns `maxᵢ(‖h(tᵢ) - xᵢ‖₂)`.  The complexity is
  // linear in the size of `samples`.
  template<typename Samples>
  NormType LInfinityL₂Error(
      Samples const& samples,
      std::function<Argument const&(typename Samples::value_type const&)> const&
          get_argument,
      std::function<Value const&(typename Samples::value_type const&)> const&
          get_value) const;

  // Returns true if the `LInfinityL₂Error` is less than `tolerance`.  More
  // efficient than the above function in the case where it returns false.
  template<typename Samples>
  bool LInfinityL₂ErrorIsWithin(
      Samples const& samples,
      std::function<Argument const&(typename Samples::value_type const&)> const&
          get_argument,
      std::function<Value const&(typename Samples::value_type const&)> const&
          get_value,
      NormType const& tolerance) const;

 private:
  static PolynomialInMonomialBasis<Value, Argument, 3> MakePolynomial(
      std::pair<Argument, Argument> const& arguments,
      std::pair<Value, Value> const& values,
      std::pair<Derivative1, Derivative1> const& derivatives);

  PolynomialInMonomialBasis<Value, Argument, 3> p_;
  PolynomialInMonomialBasis<Derivative1, Argument, 2> pʹ_;
};

}  // namespace internal

using internal::Hermite3;

}  // namespace _hermite3
}  // namespace numerics
}  // namespace principia

#include "numerics/hermite3_body.hpp"
