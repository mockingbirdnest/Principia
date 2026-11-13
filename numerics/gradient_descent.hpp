#pragma once

#include <concepts>
#include <functional>

#include "absl/status/statusor.h"
#include "base/algebra.hpp"
#include "geometry/differential_geometry.hpp"
#include "geometry/hilbert.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _gradient_descent {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::geometry::_differential_geometry;
using namespace principia::geometry::_hilbert;
using namespace principia::quantities::_quantities;

// Statuses returned by the gradient descent.
namespace termination_condition {
constexpr absl::StatusCode Done = absl::StatusCode::kOk;
constexpr absl::StatusCode NoMinimum = absl::StatusCode::kNotFound;
}  // namespace termination_condition

// Stops when the search displacement is smaller than `tolerance`.  Returns
// `NoMinimum` if no minimum is found within distance `radius` of
// `start_argument`.  The first step size defaults to the tolerance.  The
// assumption is that, if the caller provides a reasonable value then (1) we
// won't miss "interesting features" of f; (2) the finite differences won't
// underflow or have other unpleasant properties.
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
  requires(std::totally_ordered<Scalar>)
absl::StatusOr<Argument> BroydenFletcherGoldfarbShanno(
    Argument const& start_argument,
    VectorField<Scalar, Argument> const& f,
    GradientField<Scalar, Argument> const& grad_f,
    NormType<Difference<Argument>> const& tolerance,
    NormType<Difference<Argument>> const& radius =
        Infinity<NormType<Difference<Argument>>>,
    std::optional<NormType<Difference<Argument>>> const&
        first_step = std::nullopt);

// Same as above, but the Gateaux derivative of f is passed in addition to its
// gradient.  Useful when the Gateaux derivative is significantly less expensive
// to compute than the full gradient.
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
  requires(std::totally_ordered<Scalar>)
absl::StatusOr<Argument> BroydenFletcherGoldfarbShanno(
    Argument const& start_argument,
    VectorField<Scalar, Argument> const& f,
    GradientField<Scalar, Argument> const& grad_f,
    GateauxDerivative<Scalar, Argument> const& gateaux_derivative_f,
    NormType<Difference<Argument>> const& tolerance,
    NormType<Difference<Argument>> const& radius =
        Infinity<NormType<Difference<Argument>>>,
    std::optional<NormType<Difference<Argument>>> const&
        first_step = std::nullopt);

}  // namespace internal

using internal::BroydenFletcherGoldfarbShanno;
namespace termination_condition = internal::termination_condition;

}  // namespace _gradient_descent
}  // namespace numerics
}  // namespace principia

#include "numerics/gradient_descent_body.hpp"
