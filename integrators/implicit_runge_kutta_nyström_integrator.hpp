#pragma once

#include <functional>
#include <vector>

#include "base/not_null.hpp"
#include "integrators/motion_integrator.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using base::not_null;
using quantities::Variation;

namespace integrators {

class IRKNIntegrator : public MotionIntegrator {
 public:
  IRKNIntegrator(std::vector<std::vector<double>> const& a,
                 std::vector<double> const& b,
                 std::vector<double> const& c);

  virtual ~IRKNIntegrator() = default;

  IRKNIntegrator() = delete;
  IRKNIntegrator(IRKNIntegrator const&) = delete;
  IRKNIntegrator(IRKNIntegrator&&) = delete;
  IRKNIntegrator& operator=(IRKNIntegrator const&) = delete;
  IRKNIntegrator& operator=(IRKNIntegrator&&) = delete;

  template<typename Position>
  using IRKNRightHandSideComputation =
      std::function<
          void(Time const& t,
               std::vector<Position> const&,
               not_null<std::vector<Variation<Variation<Position>>>*> const)>;

  template<typename Position>
  void SolveTrivialKineticEnergyIncrement(
      IRKNRightHandSideComputation<Position> compute_acceleration,
      Parameters<Position, Variation<Position>> const& parameters,
      not_null<Solution<Position, Variation<Position>>*> const solution) const;

 protected:
  int stages_;
  // The Runge-Kutta matrix.
  std::vector<std::vector<double>> a_;
  // The weights.
  std::vector<double> b_;
  // The nodes.
  std::vector<double> c_;
  // dᵀ = bᵀA⁻¹.
  std::vector<double> d_;

}

}  // namespace integrators
}  // namespace principia

#include "integrators/implicit_runge_kutta_nyström_integrator_body.hpp"
