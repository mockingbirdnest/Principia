#include "physics/n_body_system.hpp"

#include <cmath>
#include <functional>
#include <vector>

#include "geometry/r3_element.hpp"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "quantities/quantities.hpp"

// TODO(phl): This is a header file, you're polluting the root namespace!
// Put that in the function bodies.
using principia::geometry::R3Element;
using principia::integrators::SPRKIntegrator;
using principia::integrators::SymplecticIntegrator;
using principia::quantities::Length;
using principia::quantities::Speed;
using principia::quantities::Quantity;

namespace principia {
namespace physics {

namespace {

template<typename Scalar>
Scalar FromDouble(double const quantity) {
  return quantity * quantities::SIUnit<Scalar>();
}

template<typename Scalar, typename Frame>
Vector<Scalar, Frame> FromDouble(double const x,
                                 double const y,
                                 double const z) {
  R3Element<Scalar> coordinates;
  coordinates.x = FromDouble<Scalar>(x);
  coordinates.y = FromDouble<Scalar>(y);
  coordinates.z = FromDouble<Scalar>(z);
  return Vector<Scalar, Frame>(coordinates);
}

template<typename Scalar>
double ToDouble(Scalar const& quantity) {
  return quantity / quantities::SIUnit<Scalar>();
}

template<typename Scalar, typename Frame>
std::vector<double> ToDouble(Vector<Scalar, Frame> const& vector) {
  R3Element<Scalar> const& coordinates = vector.coordinates();
  return {ToDouble(coordinates.x),
          ToDouble(coordinates.y),
          ToDouble(coordinates.z)};
}

}  // namespace

template<typename InertialFrame>
NBodySystem<InertialFrame>::NBodySystem(
    std::vector<Body<InertialFrame>*> const* bodies) : bodies_(bodies) {}

template<typename InertialFrame>
NBodySystem<InertialFrame>::~NBodySystem() {
  for (Body<InertialFrame>* body : *bodies_) {
    delete body;
  }
}

template<typename InertialFrame>
std::vector<Body<InertialFrame>*> const& NBodySystem<InertialFrame>::bodies() {
  return *bodies_;
}


template<typename InertialFrame>
void NBodySystem<InertialFrame>::Integrate(
    SymplecticIntegrator const& integrator,
    Time const& tmax,
    Time const& Δt,
    int const sampling_period) {
  using quantities::SIUnit;
  SymplecticIntegrator::Parameters parameters;
  SymplecticIntegrator::Solution solution;

  // Prepare the input data.
  std::unique_ptr<Time> reference_time;
  for (Body<InertialFrame> const* body : *bodies_) {
    Vector<Length, InertialFrame> const& position = body->positions().back();
    Vector<Speed, InertialFrame> const& velocity = body->velocities().back();
    Time const& time = body->times().back();
    for (double const q : ToDouble(position)) {
      parameters.q0.push_back(q);
    }
    for (double const p : ToDouble(velocity)) {
      parameters.p0.push_back(p);
    }
    parameters.t0 = ToDouble(time);
    // All the positions/velocities must be for the same time.
    if (reference_time == nullptr) {
      reference_time.reset(new Time(time));
    } else {
      CHECK_EQ(*reference_time, time);
    }
  }

  parameters.tmax = tmax / SIUnit<Time>();
  parameters.Δt = Δt / SIUnit<Time>();
  parameters.sampling_period = sampling_period;
  dynamic_cast<const SPRKIntegrator*>(&integrator)->Solve(
      std::bind(&NBodySystem::ComputeGravitationalAccelerations, this,
                std::placeholders::_1,
                std::placeholders::_2,
                std::placeholders::_3),
      &ComputeGravitationalVelocities,
      parameters, &solution);

  // TODO(phl): Ignoring errors for now.
  CHECK_EQ(solution.position.size(), solution.momentum.size());
  std::vector<double> const& t = solution.time.quantities;
  // Loop over the bodies.
  // TODO(phl): It looks like we are transposing in the integrator and then
  // transposing here again.
  for (size_t i = 0; i < solution.position.size(); i += 3) {
    Body<InertialFrame>* body = (*bodies_)[i / 3];
    std::vector<double> const& q0 = solution.position[i + 0].quantities;
    std::vector<double> const& q1 = solution.position[i + 1].quantities;
    std::vector<double> const& q2 = solution.position[i + 2].quantities;
    std::vector<double> const& p0 = solution.momentum[i + 0].quantities;
    std::vector<double> const& p1 = solution.momentum[i + 1].quantities;
    std::vector<double> const& p2 = solution.momentum[i + 2].quantities;
    CHECK_EQ(t.size(), q0.size());
    CHECK_EQ(t.size(), q1.size());
    CHECK_EQ(t.size(), q2.size());
    CHECK_EQ(t.size(), p0.size());
    CHECK_EQ(t.size(), p1.size());
    CHECK_EQ(t.size(), p2.size());
    for (size_t j = 0; j < t.size(); ++j) {
      Vector<Length, InertialFrame> const position =
          FromDouble<Length, InertialFrame>(q0[j], q1[j], q2[j]);
      Vector<Speed, InertialFrame> const velocity =
          FromDouble<Speed, InertialFrame>(p0[j], p1[j], p2[j]);
      Time const time = FromDouble<Time>(t[j]);
      body->AppendToTrajectory(position, velocity, time);
    }
  }
}

template<typename InertialFrame>
void NBodySystem<InertialFrame>::ComputeGravitationalAccelerations(
    double const t,
    std::vector<double> const& q,
    std::vector<double>* result) const {
  using quantities::SIUnit;
  static auto const dimension_factor = 1 / SIUnit<GravitationalParameter>();
  result->assign(result->size(), 0);

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
              ((*bodies_)[b2]->gravitational_parameter() / denominator) *
                  dimension_factor;
          (*result)[3 * b1] -= Δq0 * μ2OverRSquared;
          (*result)[3 * b1 + 1] -= Δq1 * μ2OverRSquared;
          (*result)[3 * b1 + 2] -= Δq2 * μ2OverRSquared;
        }
        // Lex. III. Actioni contrariam semper & æqualem esse reactionem:
        // sive corporum duorum actiones in se mutuo semper esse æquales &
        // in partes contrarias dirigi.
        if (!(*bodies_)[b1]->is_massless()) {
          double const μ1OverRSquared =
              ((*bodies_)[b1]->gravitational_parameter() / denominator) *
              dimension_factor;
          (*result)[3 * b2] += Δq0 * μ1OverRSquared;
          (*result)[3 * b2 + 1] += Δq1 * μ1OverRSquared;
          (*result)[3 * b2 + 2] += Δq2 * μ1OverRSquared;
        }
      }
    }
  }
}

template<typename InertialFrame>
void NBodySystem<InertialFrame>::ComputeGravitationalVelocities(
    std::vector<double> const& p,
    std::vector<double>* result) {
  *result = p;
}

}  // namespace physics
}  // namespace principia
