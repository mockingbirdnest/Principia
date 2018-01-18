#pragma once

#include <vector>

#include "geometry/named_quantities.hpp"
#include "numerics/чебышёв_series.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_newhall {

using geometry::Instant;
using quantities::Variation;

// Computes a Newhall approximation of the given |degree| in the Чебышёв basis.
// |q| and |v| are the positions and velocities over a constant division of
// [t_min, t_max].
template<typename Vector>
ЧебышёвSeries<Vector> NewhallApproximationInЧебышёвBasis(
    int degree,
    std::vector<Vector> const& q,
    std::vector<Variation<Vector>> const& v,
    Instant const& t_min,
    Instant const& t_max);

}  // namespace internal_newhall

using internal_newhall::NewhallApproximationInЧебышёвBasis;

}  // namespace numerics
}  // namespace principia

#include "numerics/newhall_body.hpp"
