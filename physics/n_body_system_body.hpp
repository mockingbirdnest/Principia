#include "physics/n_body_system.hpp"

#include <cmath>
#include <functional>
#include <vector>

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "physics/frame.hpp"
#include "quantities/quantities.hpp"

using principia::integrators::SPRKIntegrator;
using principia::integrators::SymplecticIntegrator;
using principia::quantities::Length;
using principia::quantities::Momentum;
using principia::si::Metre;
using principia::si::Second;

namespace principia {
namespace physics {

NBodySystem::NBodySystem(std::vector<Body<InertialFrame>*> const* bodies)
    : bodies_(bodies) {}

NBodySystem::~NBodySystem() {
  for (Body<InertialFrame>* body : *bodies_) {
    delete body;
  }
}

void NBodySystem::Integrate(SymplecticIntegrator const& integrator,
                            Time const& tmax,
                            Time const& Δt,
                            int const sampling_period) {
  SymplecticIntegrator::Parameters parameters;
  SymplecticIntegrator::Solution solution;

  Vector<Length, InertialFrame> position;
  Vector<Momentum, InertialFrame> momentum;
  Time time;
  for (const Body<InertialFrame>* body : *bodies_) {
    body->GetLast(&position, &momentum, &time);
    parameters.q0.push_back(position);
    parameters.p0.push_back(momentum);
    parameters.t0.push_back(time);
  }

  parameters.tmax = (tmax / (1 * Time::SIUnit())).value();
  parameters.Δt = (Δt / (1 * Time::SIUnit())).value();
  parameters.sampling_period = sampling_period;
  dynamic_cast<const SPRKIntegrator*>(&integrator)->Solve(
      std::bind(&NBodySystem::ComputeGravitationalForces, this,
                std::placeholders::_1,
                std::placeholders::_2,
                std::placeholders::_3),
      &ComputeGravitationalVelocities,
      parameters, &solution);
}

void NBodySystem::ComputeGravitationalForces(
    double const t,
    std::vector<double> const& q,
    std::vector<double>* result) const {
  static auto dimension_factor = Metre.Pow<-3>() * Second.Pow<2>();
  // TODO(phl): Used to deal with proper accelerations here.
  for (size_t b1 = 0; b1 < bodies_->size(); ++b1) {
    for (size_t b2 = b1 + 1; b2 < bodies_->size(); ++b2) {
      if (!(*bodies_)[b1]->is_massless() || !(*bodies_)[b2]->is_massless()) {
        double const Δq0 = q[3 * b1] - q[3 * b2];
        double const Δq1 = q[3 * b1 + 1] - q[3 * b2 + 1];
        double const Δq2 = q[3 * b1 + 2] - q[3 * b2 + 2];

        double const squared_distance = Δq0 * Δq0 + Δq1 * Δq1 + Δq2 * Δq2;
        double const denominator = squared_distance * sqrt(squared_distance);

        if (!(*bodies_)[b2]->is_massless()) {
          double const μ2OverRSquared =
              (((*bodies_)[b2]->gravitational_parameter() / denominator) *
               dimension_factor).value();
          (*result)[3 * b1] -= Δq0 * μ2OverRSquared;
          (*result)[3 * b1 + 1] -= Δq1 * μ2OverRSquared;
          (*result)[3 * b1 + 2] -= Δq2 * μ2OverRSquared;
        }
        // Lex. III. Actioni contrariam semper & æqualem esse reactionem:
        // sive corporum duorum actiones in se mutuo semper esse æquales &
        // in partes contrarias dirigi.
        if (!(*bodies_)[b1]->is_massless()) {
          double const μ1OverRSquared =
              (((*bodies_)[b1]->gravitational_parameter() / denominator) *
               dimension_factor).value();
          (*result)[3 * b2] += Δq0 * μ1OverRSquared;
          (*result)[3 * b2 + 1] += Δq1 * μ1OverRSquared;
          (*result)[3 * b2 + 2] += Δq2 * μ1OverRSquared;
        }
      }
    }
  }
}

void NBodySystem::ComputeGravitationalVelocities(std::vector<double> const& p,
                                                 std::vector<double>* result) {
  *result = p;
}

}  // namespace physics
}  // namespace principia
