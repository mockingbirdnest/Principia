#pragma once

#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "numerics/polynomial_in_чебышёв_basis.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _newhall {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::numerics::_polynomial;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::numerics::_polynomial_in_чебышёв_basis;
using namespace principia::quantities::_named_quantities;

// Computes a Newhall approximation of the given `degree` in the Чебышёв basis.
// `q` and `v` are the positions and velocities over a constant division of
// [t_min, t_max].  `error_estimate` gives an estimate of the error between the
// approximation the input data.  The client probably wants to compute some
// norm of that estimate.
template<typename Value, int degree>
PolynomialInЧебышёвBasis<Value, Instant, degree>
NewhallApproximationInЧебышёвBasis(std::vector<Value> const& q,
                                   std::vector<Variation<Value>> const& v,
                                   Instant const& t_min,
                                   Instant const& t_max,
                                   Value& error_estimate);

// Same as above but the `degree` is not a constant expression.
template<typename Value>
not_null<std::unique_ptr<PolynomialInЧебышёвBasis<Value, Instant>>>
NewhallApproximationInЧебышёвBasis(int degree,
                                   std::vector<Value> const& q,
                                   std::vector<Variation<Value>> const& v,
                                   Instant const& t_min,
                                   Instant const& t_max,
                                   Value& error_estimate);

// Computes a Newhall approximation of the given `degree` in the monomial basis.
// The parameters have the same meaning as in the preceding function.  The
// result is a polynomial of `Time` valid around `(t_min + t_max) / 2` with
// an argument in the range [(t_min - t_max) / 2, (t_max - t_min) / 2].  The
// polynomial is created using the given evaluator policy.
template<typename Value, int degree>
PolynomialInMonomialBasis<Value, Instant, degree>
NewhallApproximationInMonomialBasis(std::vector<Value> const& q,
                                    std::vector<Variation<Value>> const& v,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Policy const& policy,
                                    Difference<Value>& error_estimate);

// Same as above but the `degree` is not a constant expression.
template<typename Value>
not_null<std::unique_ptr<Polynomial<Value, Instant>>>
NewhallApproximationInMonomialBasis(int degree,
                                    std::vector<Value> const& q,
                                    std::vector<Variation<Value>> const& v,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Policy const& policy,
                                    Difference<Value>& error_estimate);

}  // namespace internal

using internal::NewhallApproximationInЧебышёвBasis;
using internal::NewhallApproximationInMonomialBasis;

}  // namespace _newhall
}  // namespace numerics
}  // namespace principia

#include "numerics/newhall_body.hpp"
