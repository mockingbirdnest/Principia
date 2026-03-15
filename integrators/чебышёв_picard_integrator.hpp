// The files containing the tree of of child classes of `Integrator` must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#include "integrators/integrators.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_ЧЕБЫШЁВ_PICARD_INTEGRATOR_HPP_
#define PRINCIPIA_INTEGRATORS_ЧЕБЫШЁВ_PICARD_INTEGRATOR_HPP_

#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "base/algebra.hpp"
#include "base/not_null.hpp"
#include "geometry/direct_sum.hpp"
#include "geometry/instant.hpp"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/quantities.hpp"
#include "serialization/integrators.pb.h"

namespace principia {
namespace integrators {
namespace _чебышёв_picard_integrator {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::base::_not_null;
using namespace principia::geometry::_direct_sum;
using namespace principia::geometry::_instant;
using namespace principia::integrators::_integrators;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_quantities;

template<typename ODE>
struct ЧебышёвPicardIterationParams {
  // The maximum allowed number of Picard iterations per step. If iteration has
  // not stopped (according to the stopping criterion) by the final step, the
  // iteration will be considered to have diverged.
  std::int64_t max_iterations;

  // This function will be called on the iteration deltas for each node. If it
  // returns true for all nodes for two iterations in a row, iteration will be
  // considered to have converged.
  std::function<bool(typename ODE::State::Error const&)>
      stopping_criterion;
};

// Matrices used in Чебышёв-Picard iteration (specialized by ODE order).
template<ЧебышёвPicardMethod Method, std::int64_t order>
struct ЧебышёвPicardMatrices {};

// Matrices for first-order ODEs.
template<ЧебышёвPicardMethod Method>
struct ЧебышёвPicardMatrices<Method, 1> {
  explicit ЧебышёвPicardMatrices(
      FixedVector<double, Method::M + 1> const& nodes);

  static constexpr std::int64_t M = Method::M;
  static constexpr std::int64_t N = Method::N;

  // The product of 1.31a and 1.31b from [Mac15].
  FixedMatrix<double, M + 1, M + 1, /*use_heap=*/true> CₓCα;
};

// Matrices for second-order ODEs.
template<ЧебышёвPicardMethod Method>
struct ЧебышёвPicardMatrices<Method, 2> {
  explicit ЧебышёвPicardMatrices(
      FixedVector<double, Method::M + 1> const& nodes);

  static constexpr std::int64_t M = Method::M;
  static constexpr std::int64_t N = Method::N;

  // The product of 1.53b (v-type) and 1.53a (β-type) from [Mac15].
  FixedMatrix<double, M + 1, M + 1, /*use_heap=*/true> vCₓᵝCα;
  // The product of 1.53b (x-type), 1.53c (α-type), and 1.53a (β-type) from
  // [Mac15].
  FixedMatrix<double, M + 1, M + 1, /*use_heap=*/true> xCₓᵅCᵧᵝCα;
};

// Helper struct to select iteration state types.
template<std::int64_t M, typename ODE>
struct ЧебышёвPicardIterationState {};

template<std::int64_t M,
         typename IndependentVariable,
         typename... DependentVariable>
struct ЧебышёвPicardIterationState<
    M,
    ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable,
                                                   DependentVariable...>> {
  using ODE =
      ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable,
                                                     DependentVariable...>;
  using DependentVariableMatrix =
      FixedVector<DirectSum<DependentVariable...>, M + 1, true>;
  using RightHandSideMatrix = FixedVector<
      DirectSum<Derivative<DependentVariable, IndependentVariable>...>,
      M + 1,
      true>;

  static DependentVariableMatrix UninitializedDependentVariableMatrix(
      InitialValueProblem<ODE> const& problem);
  static RightHandSideMatrix UninitializedRightHandSideMatrix(
      InitialValueProblem<ODE> const& problem);

  static IndependentVariable GetIndependentVariable(
      typename ODE::State const& state);
};

template<std::int64_t M, typename DependentVariable>
struct ЧебышёвPicardIterationState<
    M,
    ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>> {
  using ODE =
      ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>;
  // TODO(rnlahaye): these should use some sort of "partially-fixed matrix".
  using DependentVariableMatrix =
      std::pair<UnboundedMatrix<DependentVariable>,
                UnboundedMatrix<Derivative<DependentVariable, Instant>>>;
  using RightHandSideMatrix =
      UnboundedMatrix<Derivative<DependentVariable, Instant, 2>>;

  static DependentVariableMatrix UninitializedDependentVariableMatrix(
      InitialValueProblem<ODE> const& problem);
  static RightHandSideMatrix UninitializedRightHandSideMatrix(
      InitialValueProblem<ODE> const& problem);

  static Instant GetIndependentVariable(ODE::State const& state);
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
             ЧебышёвPicardIterationParams<ODE> const& params);

    static constexpr std::int64_t M = Method::M;

    ЧебышёвPicardIntegrator const& integrator_;
    ЧебышёвPicardIterationParams<ODE> params_;

    // Working variables which are stored here so we don't need to reallocate
    // them on each Solve call.

    // Stores the nodes rescaled to the current step.
    std::vector<typename ODE::IndependentVariable> t_;

    // Controls the boundary condition.
    ЧебышёвPicardIterationState<M, ODE>::DependentVariableMatrix boundary_;

    // Xⁱ is an (M + 1)×n matrix containing the values of the dependent
    // variables at each node.
    ЧебышёвPicardIterationState<M, ODE>::DependentVariableMatrix Xⁱ_;
    ЧебышёвPicardIterationState<M, ODE>::DependentVariableMatrix Xⁱ⁺¹_;

    // The computed right-hand-side of the ODE (at each node, for the current
    // iteration).
    //
    // For a first-order ODE, this contains first derivatives. For a
    // second-order ODE, it contains second derivatives.
    ЧебышёвPicardIterationState<M, ODE>::RightHandSideMatrix f_;

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
      ЧебышёвPicardIterationParams<ODE> const& params) const;

  void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> message) const override;

 private:
  static constexpr std::int64_t M = Method::M;
  static constexpr std::int64_t N = Method::N;

  // The nodes used for function evaluation.
  //
  // These are Чебышёв nodes of the second kind.
  FixedVector<double, M + 1> nodes_;

  // The matrices used for iteration (which depend on the ODE order).
  ЧебышёвPicardMatrices<Method, ODE::order> matrices_;
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
