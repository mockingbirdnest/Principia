// The files containing the tree of of child classes of `Integrator` must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#include "integrators/integrators.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_CHEBYSHEV_PICARD_ITERATOR_HPP_
#define PRINCIPIA_INTEGRATORS_CHEBYSHEV_PICARD_ITERATOR_HPP_

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
namespace _chebyshev_picard_iterator {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::integrators::_integrators;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_quantities;

struct ChebyshevPicardIterationParams {
  // Controls the number of nodes at which the function will be evaluated.
  //
  // Note that this is the highest node _index_ rather than the number of nodes;
  // the actual number of nodes is M + 1.
  //
  // Must be at least 1.
  int M;

  // The order of the Chebyshev sequence used to approximate the system state.
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

template <typename ODE_>
class ChebyshevPicardIterator : public FixedStepSizeIntegrator<ODE_> {
 public:
  using ODE = ODE_;
  using AppendState = typename Integrator<ODE>::AppendState;

  class Instance : public FixedStepSizeIntegrator<ODE>::Instance {
   public:
    absl::Status Solve(ODE::IndependentVariable const& t_final) override;

    ChebyshevPicardIterator const& integrator() const override;

    not_null<std::unique_ptr<typename Integrator<ODE>::Instance> > Clone()
        const override;

   private:
    Instance(InitialValueProblem<ODE> const& problem,
             AppendState const& append_state, Time const& step,
             ChebyshevPicardIterator const& integrator);

    ChebyshevPicardIterator const& integrator_;

    friend class ChebyshevPicardIterator;
  };

  // Constructs a ChebyshevPicardIterator with the given parameters.
  ChebyshevPicardIterator(ChebyshevPicardIterationParams const& params);

  ChebyshevPicardIterationParams const& params() const;

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance> > NewInstance(
      InitialValueProblem<ODE> const& problem, AppendState const& append_state,
      Time const& step) const override;

  void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> message) const override;

 private:
  ChebyshevPicardIterationParams params_;

  // The nodes used for function evaluation.
  //
  // These are Chebyshev nodes of the second kind.
  UnboundedVector<double> nodes_;

  // 1.31b from Macomber's thesis.
  UnboundedMatrix<double> cx_;
  // The product of 1.31a and 1.31b from Macomber's thesis.
  UnboundedMatrix<double> cx_cα_;
};

}  // namespace internal

using internal::ChebyshevPicardIterationParams;
using internal::ChebyshevPicardIterator;

}  // namespace _chebyshev_picard_iterator
}  // namespace integrators
}  // namespace principia

#include "integrators/chebyshev_picard_iterator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_CHEBYSHEV_PICARD_ITERATOR_HPP_
#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
