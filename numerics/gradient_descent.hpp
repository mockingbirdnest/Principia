#pragma once

#include <concepts>
#include <functional>

#include "absl/status/statusor.h"
#include "base/algebra.hpp"
#include "geometry/hilbert.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _gradient_descent {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::geometry::_hilbert;
using namespace principia::quantities::_quantities;

// A helper for generating declarations derived from `Scalar` and `Argument`.
// It must declare the type of the gradient and a function returning the inner
// product form thus:
//   using Gradient = ...;
//   static ... InnerProductForm();
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
struct Generator;

template<typename Value, euclidean Argument>
using Field = std::function<Value(Argument const&)>;

template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
using Gradient = typename Generator<Scalar, Argument>::Gradient;

// Must return InnerProduct(grad_f(argument), direction);
template<homogeneous_field Scalar, homogeneous_euclidean_space<Scalar> Argument>
using GateauxDerivative =
    std::function<Scalar(Argument const& argument,
                         Difference<Argument> const& direction)>;

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
    Field<Scalar, Argument> const& f,
    Field<Gradient<Scalar, Argument>, Argument> const& grad_f,
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
    Field<Scalar, Argument> const& f,
    Field<Gradient<Scalar, Argument>, Argument> const& grad_f,
    GateauxDerivative<Scalar, Argument> const& gateaux_derivative_f,
    NormType<Difference<Argument>> const& tolerance,
    NormType<Difference<Argument>> const& radius =
        Infinity<NormType<Difference<Argument>>>,
    std::optional<NormType<Difference<Argument>>> const&
        first_step = std::nullopt);

}  // namespace internal

using internal::BroydenFletcherGoldfarbShanno;
using internal::Gradient;
using internal::Field;
namespace termination_condition = internal::termination_condition;

}  // namespace _gradient_descent
}  // namespace numerics
}  // namespace principia

#include "numerics/gradient_descent_body.hpp"
