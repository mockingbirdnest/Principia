// The files containing the tree of of child classes of `Integrator` must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#include "integrators/integrators.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_ЧЕБЫШЁВ_PICARD_INTEGRATOR_HPP_
#define PRINCIPIA_INTEGRATORS_ЧЕБЫШЁВ_PICARD_INTEGRATOR_HPP_

#include <memory>
#include <vector>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "geometry/cartesian_product.hpp"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/fixed_arrays.hpp"
#include "quantities/quantities.hpp"
#include "serialization/integrators.pb.h"

namespace principia {
namespace integrators {
namespace _чебышёв_picard_integrator {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_cartesian_product;
using namespace principia::integrators::_integrators;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::quantities::_quantities;

struct ЧебышёвPicardIterationParams {
  // The maximum allowed number of Picard iterations per step. If iteration has
  // not stopped (according to the stopping criterion) by the final step, the
  // iteration will be considered to have diverged.
  std::int64_t max_iterations;

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
//   * Approximating x₀ + ∫ₜ₀ᵗ f(t, xⁱ(t)) dt by integrating the above
//     approximation. This will be xⁱ⁺¹.
// * The sequence xⁱ should converge to (a Чебышёв approximation of) x within
//   some interval of t₀.
//
// Due to the properties of Чебышёв polynomials, the iteration may be peformed
// simply using linear algebra operations (in addition to the required
// evaluation of f at Чебышёв nodes) [FN83].
//
// This code uses the formulation from [Mac15].

template<ЧебышёвPicardMethod Method, typename ODE_>
class ЧебышёвPicardIntegrator : public FixedStepSizeIntegrator<ODE_> {
 public:
  using ODE = ODE_;
  using AppendState = typename Integrator<ODE>::AppendState;

  class Instance : public FixedStepSizeIntegrator<ODE>::Instance {
   public:
    absl::Status Solve(ODE::IndependentVariable const& t_final) override;

    ЧебышёвPicardIntegrator const& integrator() const override;

    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

   private:
    Instance(InitialValueProblem<ODE> const& problem,
             AppendState const& append_state,
             Time const& step,
             ЧебышёвPicardIntegrator const& integrator,
             ЧебышёвPicardIterationParams const& params);

    static constexpr std::int64_t M = Method::M;

    // The dimension of the ODE.
    static constexpr std::int64_t n =
        std::tuple_size_v<typename ODE::DependentVariables>;

    ЧебышёвPicardIntegrator const& integrator_;
    ЧебышёвPicardIterationParams params_;

    // Working variables which are stored here so we don't need to reallocate
    // them on each Solve call.

    // Stores the nodes rescaled to the current step.
    std::vector<typename ODE::IndependentVariable> t_;

    // Controls the boundary condition.
    FixedVector<direct_sum_t<typename ODE::DependentVariables>,
                M + 1,
                /*use_heap=*/true>
        CₓX₀_;

    // Xⁱ is an (M + 1)×n matrix containing the values of the dependent
    // variables at each node.
    FixedVector<direct_sum_t<typename ODE::DependentVariables>,
                M + 1,
                /*use_heap=*/true>
        Xⁱ_;
    FixedVector<direct_sum_t<typename ODE::DependentVariables>,
                M + 1,
                /*use_heap=*/true>
        Xⁱ⁺¹_;

    // The computed derivative (at each node, for the current iteration).
    FixedVector<direct_sum_t<typename ODE::DependentVariableDerivatives>,
                M + 1,
                /*use_heap=*/true>
        yʹ_;

    friend class ЧебышёвPicardIntegrator;
  };

  // Constructs a ЧебышёвPicardIntegrator.
  ЧебышёвPicardIntegrator();

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      InitialValueProblem<ODE> const& problem,
      AppendState const& append_state,
      Time const& step) const override;

  // Temporary NewInstance overload taking parameters. This will be removed once
  // instance parameters are added to FixedStepSizeIntegrator.
  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      InitialValueProblem<ODE> const& problem,
      AppendState const& append_state,
      Time const& step,
      ЧебышёвPicardIterationParams const& params) const;

  void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> message) const override;

 private:
  static constexpr std::int64_t M = Method::M;
  static constexpr std::int64_t N = Method::N;

  // The nodes used for function evaluation.
  //
  // These are Чебышёв nodes of the second kind.
  FixedVector<double, M + 1> nodes_;

  // The product of 1.31a and 1.31b from [Mac15].
  FixedMatrix<double, M + 1, M + 1, /*use_heap=*/true> CₓCα_;
};

}  // namespace internal

using internal::ЧебышёвPicardIntegrator;
using internal::ЧебышёвPicardIterationParams;

}  // namespace _чебышёв_picard_integrator
}  // namespace integrators
}  // namespace principia

#include "integrators/чебышёв_picard_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_ЧЕБЫШЁВ_PICARD_INTEGRATOR_HPP_
#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
