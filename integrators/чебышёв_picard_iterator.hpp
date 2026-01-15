// The files containing the tree of of child classes of `Integrator` must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#include "integrators/integrators.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_ЧЕБЫШЁВ_PICARD_ITERATOR_HPP_
#define PRINCIPIA_INTEGRATORS_ЧЕБЫШЁВ_PICARD_ITERATOR_HPP_

#include <memory>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/instant.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/quantities.hpp"
#include "serialization/integrators.pb.h"

namespace principia {
namespace integrators {
namespace _чебышёв_picard_iterator {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::integrators::_integrators;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_quantities;

struct ЧебышёвPicardIterationParams {
  // Controls the number of nodes at which the function will be evaluated.
  //
  // Note that this is the highest node _index_ rather than the number of nodes;
  // the actual number of nodes is M + 1.
  //
  // Must be at least 1.
  int M;

  // The order of the Чебышёв series used to approximate the system state.
  //
  // Must be at least 1 (if you want to approximate your system with a constant,
  // use some other method).
  int N;

  // The maximum allowed number of Picard iterations per step. If iteration has
  // not stopped (according to the stopping criterion) by the final step, the
  // iteration will be considered to have diverged.
  int max_iterations;

  // If the maximum absolute difference between successive state approximations
  // is less than this for two Picard iterations in a row, iteration will be
  // considered to have converged.
  double stopping_criterion;
};

// This class solves ordinary differential equations of the form x′ = f(x, t)
// using Чебышёв-Picard iteration.
//
// Чебышёв-Picard iteration combines Чебышёв interpolation with Picard
// iteration; it was first proposed in [CN63]. It works as follows:
//
// * Start with some initial approximation x⁰ (the constant function
//   x⁰(t) = x(t₀) is typically used).
// * Repeatedly construct xⁱ⁺¹ from xⁱ by
//   * Approximating f(t, xⁱ) using Чебышёв interpolation.
//   * Approximating x₀ + ∫_t₀^t f(t, xⁱ(t)) dt by integrating the above
//     approximation. This will be xⁱ⁺¹.
// * The sequence xⁱ should converge to (a Чебышёв approximation of) x within
//   some interval of t₀.
//
// Due to the properties of Чебышёв polynomials, the iteration may peformed
// simply using linear algeba operations (in addition to the required evaluation
// of f at Чебышёв nodes) [FN83].
//
// This code uses the formulation from [Mac15].

template <typename ODE_>
class ЧебышёвPicardIterator : public FixedStepSizeIntegrator<ODE_> {
 public:
  using ODE = ODE_;
  using AppendState = typename Integrator<ODE>::AppendState;

  class Instance : public FixedStepSizeIntegrator<ODE>::Instance {
   public:
    absl::Status Solve(ODE::IndependentVariable const& t_final) override;

    ЧебышёвPicardIterator const& integrator() const override;

    not_null<std::unique_ptr<typename Integrator<ODE>::Instance> > Clone()
        const override;

   private:
    Instance(InitialValueProblem<ODE> const& problem,
             AppendState const& append_state, Time const& step,
             ЧебышёвPicardIterator const& integrator);

    ЧебышёвPicardIterator const& integrator_;

    // Working variables which are stored here so we don't need to reallocate
    // them on each Solve call.

    // Stores the nodes rescaled to the current step.
    std::vector<typename ODE::IndependentVariable> t_;

    // Controls the boundary condition.
    UnboundedMatrix<double> CₓX₀_;

    // Xⁱ is an (M + 1)×n matrix containing the values of the dependent
    // variables at each node.
    UnboundedMatrix<double> Xⁱ_;
    UnboundedMatrix<double> Xⁱ⁺¹_;

    // The computed derivative (at each node, for the current iteration).
    UnboundedMatrix<double> yʹ_;

    friend class ЧебышёвPicardIterator;
  };

  // Constructs a ЧебышёвPicardIterator with the given parameters.
  ЧебышёвPicardIterator(ЧебышёвPicardIterationParams const& params);

  ЧебышёвPicardIterationParams const& params() const;

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance> > NewInstance(
      InitialValueProblem<ODE> const& problem, AppendState const& append_state,
      Time const& step) const override;

  void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> message) const override;

 private:
  ЧебышёвPicardIterationParams params_;

  // The nodes used for function evaluation.
  //
  // These are Чебышёв nodes of the second kind.
  UnboundedVector<double> nodes_;

  // The product of 1.31a and 1.31b from [Mac15].
  UnboundedMatrix<double> CₓCα_;
};

}  // namespace internal

using internal::ЧебышёвPicardIterationParams;
using internal::ЧебышёвPicardIterator;

}  // namespace _чебышёв_picard_iterator
}  // namespace integrators
}  // namespace principia

#include "integrators/чебышёв_picard_iterator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_ЧЕБЫШЁВ_PICARD_ITERATOR_HPP_
#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
