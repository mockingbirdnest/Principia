#pragma once

#include <memory>
#include <vector>

#include "base/box.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/чебышёв_series.hpp"
#include "numerics/polynomial.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_newhall {

using base::Box;
using geometry::Instant;
using quantities::Variation;

// Computes a Newhall approximation of the given |degree| in the Чебышёв basis.
// |q| and |v| are the positions and velocities over a constant division of
// [t_min, t_max].  |error_estimate| gives an estimate of the error between the
// approximation the input data.  The client probably wants to compute some
// norm of that estimate.
template<typename Vector>
ЧебышёвSeries<Vector>
NewhallApproximationInЧебышёвBasis(int degree,
                                   std::vector<Vector> const& q,
                                   std::vector<Variation<Vector>> const& v,
                                   Instant const& t_min,
                                   Instant const& t_max,
                                   Vector& error_estimate);

// Computes a Newhall approximation of the given |degree| in the monomial basis.
// The parameters have the same meaning as in the preceding function.  The
// result is a polynomial of |Time| valid around |(t_min + t_max) / 2| with
// an argument in the range [(t_min - t_max) / 2, (t_max - t_min) / 2].
template<typename Vector, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Vector, Instant, degree, Evaluator>
NewhallApproximationInMonomialBasis(std::vector<Vector> const& q,
                                    std::vector<Variation<Vector>> const& v,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Vector& error_estimate);

// Same as above but the |degree| is not a constant expression.
template<typename Vector,
         template<typename, typename, int> class Evaluator>
Box<Polynomial<Vector, Instant>>
NewhallApproximationInMonomialBasis(int degree,
                                    std::vector<Vector> const& q,
                                    std::vector<Variation<Vector>> const& v,
                                    Instant const& t_min,
                                    Instant const& t_max,
                                    Vector& error_estimate);

}  // namespace internal_newhall

using internal_newhall::NewhallApproximationInЧебышёвBasis;
using internal_newhall::NewhallApproximationInMonomialBasis;

}  // namespace numerics
}  // namespace principia

#include "numerics/newhall_body.hpp"
