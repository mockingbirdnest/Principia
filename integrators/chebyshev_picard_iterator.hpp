#pragma once

#include <memory>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/instant.hpp"
#include "integrators/integrators.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace _chebyshev_picard_iterator {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::integrators::_integrators;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_quantities;

template <typename ODE_>
class ChebyshevPicardIterator : public FixedStepSizeIntegrator<ODE_> {
 public:
  using ODE = ODE_;
  using AppendState = typename Integrator<ODE>::AppendState;

  class Instance : public FixedStepSizeIntegrator<ODE>::Instance {
   public:
    absl::Status Solve(ODE::IndependentVariable const& t_final) override;

    ChebyshevPicardIterator const& integrator() const override;

    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

   private:
    Instance(InitialValueProblem<ODE> const& problem,
             AppendState const& append_state, Time const& step,
             ChebyshevPicardIterator const& integrator);

    ChebyshevPicardIterator const& integrator_;

    friend class ChebyshevPicardIterator;
  };

  // Constructs a ChebyshevPicardIterator the given Chebyshev order and number
  // of sample points.
  //
  // The value of sample_points must be at least order.
  ChebyshevPicardIterator(int order, int sample_points, int max_iterations,
                          double stopping_criterion);

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      InitialValueProblem<ODE> const& problem, AppendState const& append_state,
      Time const& step) const override;

  void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> message) const override;

 private:
  // The maximum number of iterations to do per step. If we exceed this an error
  // is returned.
  int max_iterations_;

  // If the max relative difference between successive iterations of the
  // dependent variables are less than this value twice in a row, the sequence
  // will be considered converged.
  double stopping_criterion_;

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

using internal::ChebyshevPicardIterator;

}  // namespace _chebyshev_picard_iterator
}  // namespace integrators
}  // namespace principia

#include "integrators/chebyshev_picard_iterator_body.hpp"
