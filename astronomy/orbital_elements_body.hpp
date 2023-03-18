#pragma once

#include "astronomy/orbital_elements.hpp"

#include <algorithm>
#include <tuple>
#include <vector>

#include "absl/strings/str_cat.h"
#include "base/jthread.hpp"
#include "base/status_utilities.hpp"
#include "numerics/quadrature.hpp"
#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"
#include "integrators/methods.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/kepler_orbit.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace astronomy {
namespace _orbital_elements {
namespace internal {

using namespace principia::base::_jthread;
using namespace principia::geometry::_named_quantities;
using namespace principia::integrators::_embedded_explicit_runge_kutta_integrator;  // NOLINT
using namespace principia::integrators::_integrators;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::numerics::_quadrature;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_kepler_orbit;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

constexpr int osculating_equinoctial_elements_per_sidereal_period = 256;
constexpr double max_clenshaw_curtis_relative_error = 1.0e-6;
constexpr int max_clenshaw_curtis_points = 2000;
// Carefully tuned based on MercuryOrbiter test.
constexpr double max_clenshaw_curtis_relative_error_for_initial_integration =
    1.0e-8;
constexpr Length eerk_a_tolerance = 10 * Milli(Metre);

template<typename Inertial, typename PrimaryCentred>
absl::StatusOr<OrbitalElements> OrbitalElements::ForTrajectory(
    Trajectory<Inertial> const& secondary_trajectory,
    DynamicFrame<Inertial, PrimaryCentred> const& primary_centred,
    MassiveBody const& primary,
    Body const& secondary,
    bool fill_osculating_equinoctial_elements) {
  return ForRelativeDegreesOfFreedom<PrimaryCentred>(
      [&primary_centred, &secondary_trajectory](Instant const& t) {
        return primary_centred.ToThisFrameAtTime(t)(
                   secondary_trajectory.EvaluateDegreesOfFreedom(t)) -
               DegreesOfFreedom<PrimaryCentred>{PrimaryCentred::origin,
                                                PrimaryCentred::unmoving};
      },
      std::max(primary_centred.t_min(), secondary_trajectory.t_min()),
      std::min(primary_centred.t_max(), secondary_trajectory.t_max()),
      primary,
      secondary,
      fill_osculating_equinoctial_elements);
}

template<typename PrimaryCentred>
absl::StatusOr<OrbitalElements> OrbitalElements::ForTrajectory(
    Trajectory<PrimaryCentred> const& trajectory,
    MassiveBody const& primary,
    Body const& secondary,
    bool const fill_osculating_equinoctial_elements) {
  return ForRelativeDegreesOfFreedom<PrimaryCentred>(
      [&trajectory](Instant const& t) {
        return trajectory.EvaluateDegreesOfFreedom(t) -
               DegreesOfFreedom<PrimaryCentred>{PrimaryCentred::origin,
                                                PrimaryCentred::unmoving};
      },
      trajectory.t_min(),
      trajectory.t_max(),
      primary,
      secondary,
      fill_osculating_equinoctial_elements);
}

inline std::vector<OrbitalElements::ClassicalElements> const&
OrbitalElements::mean_elements() const {
  return mean_classical_elements_;
}

inline Time OrbitalElements::sidereal_period() const {
  return sidereal_period_;
}

inline Time OrbitalElements::nodal_period() const {
  return nodal_period_;
}

inline Time OrbitalElements::anomalistic_period() const {
  return anomalistic_period_;
}

inline AngularFrequency OrbitalElements::nodal_precession() const {
  return nodal_precession_;
}

inline Interval<Length> OrbitalElements::mean_semimajor_axis_interval() const {
  return mean_semimajor_axis_interval_;
}

inline Interval<double> OrbitalElements::mean_eccentricity_interval() const {
  return mean_eccentricity_interval_;
}

inline Interval<Angle> OrbitalElements::mean_inclination_interval() const {
  return mean_inclination_interval_;
}

inline Interval<Angle>
OrbitalElements::mean_longitude_of_ascending_node_interval() const {
  return mean_longitude_of_ascending_node_interval_;
}

inline Interval<Angle> OrbitalElements::mean_argument_of_periapsis_interval()
    const {
  return mean_argument_of_periapsis_interval_;
}

inline Interval<Length> OrbitalElements::mean_periapsis_distance_interval()
    const {
  return mean_periapsis_distance_interval_;
}

inline Interval<Length> OrbitalElements::mean_apoapsis_distance_interval()
    const {
  return mean_apoapsis_distance_interval_;
}

inline std::vector<OrbitalElements::EquinoctialElements> const&
OrbitalElements::osculating_equinoctial_elements() const {
  return osculating_equinoctial_elements_;
}

inline std::vector<OrbitalElements::EquinoctialElements> const&
OrbitalElements::mean_equinoctial_elements() const {
  return mean_equinoctial_elements_;
}

template<typename Frame, typename RelativeDegreesOfFreedomComputation>
absl::StatusOr<OrbitalElements> OrbitalElements::ForRelativeDegreesOfFreedom(
    RelativeDegreesOfFreedomComputation const&
        relative_degrees_of_freedom_at_time,
    Instant const& t_min,
    Instant const& t_max,
    MassiveBody const& primary,
    Body const& secondary,
    bool fill_osculating_equinoctial_elements) {
  OrbitalElements orbital_elements;
  if (t_min >= t_max) {
    return absl::InvalidArgumentError(
        absl::StrCat("trajectory has min time ",
                     DebugString(t_min),
                     " and max time ",
                     DebugString(t_max)));
  }

  auto const osculating_elements =
      [&primary, &secondary, &relative_degrees_of_freedom_at_time](
          Instant const& time) -> KeplerianElements<Frame> {
    return KeplerOrbit<Frame>(primary,
                              secondary,
                              relative_degrees_of_freedom_at_time(time),
                              time)
        .elements_at_epoch();
  };

  auto const wound_osculating_λ =
      [&osculating_elements](Instant const& time) -> Angle {
    auto const elements = osculating_elements(time);
    Angle const& ϖ = *elements.longitude_of_periapsis;
    Angle const& M = *elements.mean_anomaly;
    return ϖ + M;
  };

  KeplerianElements<Frame> const initial_osculating_elements =
      osculating_elements(t_min);
  Time const estimated_period = *initial_osculating_elements.period;
  if (!IsFinite(estimated_period) || estimated_period <= Time{}) {
    return absl::OutOfRangeError("estimated period is " +
                                 DebugString(estimated_period));
  }

  std::vector<Angle> unwound_λs;
  // 3 is greater than 2 to make sure that we start in the right direction.
  Time const third_of_estimated_period = estimated_period / 3;
  unwound_λs.reserve(std::floor((t_max - t_min) / third_of_estimated_period) +
                     1);
  for (Instant t = t_min; t <= t_max; t += third_of_estimated_period) {
    Angle const λ = wound_osculating_λ(t);
    unwound_λs.push_back(unwound_λs.empty() ? λ
                                            : UnwindFrom(unwound_λs.back(), λ));
  }

  auto const osculating_equinoctial_elements =
      [&osculating_elements, t_min, third_of_estimated_period, &unwound_λs](
          Instant const& time) -> EquinoctialElements {
    auto const elements = osculating_elements(time);
    double const& e = *elements.eccentricity;
    Angle const& ϖ = *elements.longitude_of_periapsis;
    Angle const& Ω = elements.longitude_of_ascending_node;
    Angle const& M = *elements.mean_anomaly;
    Angle const& i = elements.inclination;
    double const tg_½i = Tan(i / 2);
    double const cotg_½i = 1 / tg_½i;
    double const sin_Ω = Sin(Ω);
    double const cos_Ω = Cos(Ω);
    return {.t = time,
            .a = *elements.semimajor_axis,
            .h = e * Sin(ϖ),
            .k = e * Cos(ϖ),
            .λ = UnwindFrom(
                unwound_λs[(time - t_min) / third_of_estimated_period], ϖ + M),
            .p = tg_½i * sin_Ω,
            .q = tg_½i * cos_Ω,
            .pʹ = cotg_½i * sin_Ω,
            .qʹ = cotg_½i * cos_Ω};
  };

  auto const sidereal_period =
      SiderealPeriod(osculating_equinoctial_elements, t_min, t_max);
  RETURN_IF_ERROR(sidereal_period);
  orbital_elements.sidereal_period_ = sidereal_period.value();
  if (!IsFinite(orbital_elements.sidereal_period_) ||
      orbital_elements.sidereal_period_ <= Time{}) {
    // Guard against NaN sidereal periods (from hyperbolic orbits) or negative
    // sidereal periods (from aberrant trajectories, see #2811).
    return absl::OutOfRangeError(
        "sidereal period is " + DebugString(orbital_elements.sidereal_period_));
  }

  auto mean_equinoctial_elements =
      MeanEquinoctialElements(osculating_equinoctial_elements,
                              t_min, t_max,
                              orbital_elements.sidereal_period_);
  RETURN_IF_ERROR(mean_equinoctial_elements);
  orbital_elements.mean_equinoctial_elements_ =
      std::move(mean_equinoctial_elements).value();

  if (fill_osculating_equinoctial_elements) {
    for (Instant t = t_min;
         t <= t_max;
         t += orbital_elements.sidereal_period_ /
              osculating_equinoctial_elements_per_sidereal_period) {
      orbital_elements.osculating_equinoctial_elements_.push_back(
          osculating_equinoctial_elements(t));
    }
  }

  if (orbital_elements.mean_equinoctial_elements_.size() < 2) {
    return absl::OutOfRangeError(
        "trajectory does not span one sidereal period: sidereal period is " +
        DebugString(orbital_elements.sidereal_period_) + ", trajectory spans " +
        DebugString(t_max - t_min));
  }

  auto mean_classical_elements =
      ToClassicalElements(orbital_elements.mean_equinoctial_elements_);
  RETURN_IF_ERROR(mean_classical_elements);
  orbital_elements.mean_classical_elements_ =
      std::move(mean_classical_elements).value();
  RETURN_IF_ERROR(orbital_elements.ComputePeriodsAndPrecession());
  RETURN_IF_ERROR(orbital_elements.ComputeIntervals());
  return orbital_elements;
}

template<typename EquinoctialElementsComputation>
absl::StatusOr<Time> OrbitalElements::SiderealPeriod(
    EquinoctialElementsComputation const& equinoctial_elements,
    Instant const& t_min,
    Instant const& t_max) {
  Time const Δt = t_max - t_min;
  Instant const t0 = t_min + Δt / 2;
  Product<Angle, Square<Time>> const ʃ_λt_dt = AutomaticClenshawCurtis(
      [&equinoctial_elements, &t0](
          Instant const& t) -> Product<Angle, Time> {
        // TODO(egg): Consider computing only λ.
        return equinoctial_elements(t).λ * (t - t0);
      },
      t_min,
      t_max,
      max_clenshaw_curtis_relative_error,
      max_clenshaw_curtis_points);
  return 2 * π * Radian * Pow<3>(Δt) / (12 * ʃ_λt_dt);
}

template<typename EquinoctialElementsComputation>
absl::StatusOr<std::vector<OrbitalElements::EquinoctialElements>>
OrbitalElements::MeanEquinoctialElements(
    EquinoctialElementsComputation const& equinoctial_elements,
    Instant const& t_min,
    Instant const& t_max,
    Time const& period) {
  // This function averages the elements in |osculating| over |period|.  For
  // each |mean_elements| in the result, for all э in the set of equinoctial
  // elements {a, h, k, λ, p, q, pʹ, qʹ}, |mean_elements.э| is the integral of
  // the osculating э from |mean_elements.t - period / 2| to
  // |mean_elements.t + period / 2|, divided by |period|.

  // We integrate the function (э(t + period / 2) - э(t - period / 2)) / period
  // using as the initial value an integral obtained by Clenshaw-Curtis.

  using ODE =
      ExplicitFirstOrderOrdinaryDifferentialEquation<Instant,
                                                     Length,
                                                     double,
                                                     double,
                                                     Angle,
                                                     double,
                                                     double,
                                                     double,
                                                     double>;
  ODE const equation = {
      .compute_derivative = [&equinoctial_elements, period](
                                Instant const& t,
                                ODE::DependentVariables const& y,
                                ODE::DependentVariableDerivatives& yʹ) {
        auto const [_1, a₋, h₋, k₋, λ₋, p₋, q₋, pʹ₋, qʹ₋] =
            equinoctial_elements(t - period / 2);
        auto const [_2, a₊, h₊, k₊, λ₊, p₊, q₊, pʹ₊, qʹ₊] =
            equinoctial_elements(t + period / 2);
        yʹ = {(a₊ - a₋) / period,
              (h₊ - h₋) / period,
              (k₊ - k₋) / period,
              (λ₊ - λ₋) / period,
              (p₊ - p₋) / period,
              (q₊ - q₋) / period,
              (pʹ₊ - pʹ₋) / period,
              (qʹ₊ - qʹ₋) / period};
        return absl::OkStatus();
      }};

  std::vector<EquinoctialElements> mean_elements;
  auto const append_state = [&mean_elements](ODE::State const& state) {
    Instant const& t = state.s.value;
    auto const& [a, h, k, λ, p, q, pʹ, qʹ] = state.y;
    mean_elements.push_back(EquinoctialElements{.t = t,
                                                .a = a.value,
                                                .h = h.value,
                                                .k = k.value,
                                                .λ = λ.value,
                                                .p = p.value,
                                                .q = q.value,
                                                .pʹ = pʹ.value,
                                                .qʹ = qʹ.value});
  };

  auto const tolerance_to_error_ratio =
      [](Time const& step,
         ODE::State const& state,
         ODE::State::Error const& error) -> double {
    auto const& [Δa, Δh, Δk, Δλ, Δp, Δq, Δpʹ, Δqʹ] = error;
    return eerk_a_tolerance / Abs(Δa);
  };

  auto const initial_integration =
      [&equinoctial_elements, period, t_min, t_max](auto const element) {
        return AutomaticClenshawCurtis(
                   [element, &equinoctial_elements](Instant const& t) {
                     return equinoctial_elements(t).*element;
                   },
                   t_min,
                   t_min + period,
                   max_clenshaw_curtis_relative_error_for_initial_integration,
                   /*max_points=*/max_clenshaw_curtis_points) /
               period;
      };

  // Ensure that Clenshaw-Curtis will not go out of the bounds of the
  // trajectory.
  if (t_max < t_min + period) {
    return mean_elements;
  }

  ODE::DependentVariables const initial_mean_elements{
      initial_integration(&EquinoctialElements::a),
      initial_integration(&EquinoctialElements::h),
      initial_integration(&EquinoctialElements::k),
      initial_integration(&EquinoctialElements::λ),
      initial_integration(&EquinoctialElements::p),
      initial_integration(&EquinoctialElements::q),
      initial_integration(&EquinoctialElements::pʹ),
      initial_integration(&EquinoctialElements::qʹ)};

  // Compute bounds that make sure that the ODE integrator never evaluate the
  // trajectory outside of its bounds.
  Instant t₁ = t_min + period / 2;
  if (t₁ - period / 2 < t_min) {
    t₁ = NextUp(t₁);
  }
  Instant t₂ = t_max - period / 2;
  if (t₂ + period / 2 > t_max) {
    t₂ = NextDown(t₂);
  }

  InitialValueProblem<ODE> const problem = {
      .equation = equation,
      .initial_state = ODE::State(t₁, initial_mean_elements)};
  append_state(problem.initial_state);

  Time first_step = t₂ - t₁;
  if (t₁ + first_step > t₂) {
    first_step = NextDown(first_step);
  }
  if (first_step <= Time{}) {
    return mean_elements;
  }

  auto const instance =
      EmbeddedExplicitRungeKuttaIntegrator<
          integrators::methods::DormandPrince1986RK547FC,
          ODE>()
          .NewInstance(problem,
                       append_state,
                       tolerance_to_error_ratio,
                       AdaptiveStepSizeIntegrator<ODE>::Parameters(
                           /*first_step=*/first_step,
                           /*safety_factor=*/0.9));
  RETURN_IF_ERROR(instance->Solve(t₂));

  return mean_elements;
}

inline absl::StatusOr<std::vector<OrbitalElements::ClassicalElements>>
OrbitalElements::ToClassicalElements(
    std::vector<EquinoctialElements> const& equinoctial_elements) {
  std::vector<ClassicalElements> classical_elements;
  classical_elements.reserve(equinoctial_elements.size());
  for (auto const& equinoctial : equinoctial_elements) {
    RETURN_IF_STOPPED;
    double const tg_½i = Sqrt(Pow<2>(equinoctial.p) + Pow<2>(equinoctial.q));
    double const cotg_½i =
        Sqrt(Pow<2>(equinoctial.pʹ) + Pow<2>(equinoctial.qʹ));
    Angle const i =
        cotg_½i > tg_½i ? 2 * ArcTan(tg_½i) : 2 * ArcTan(1 / cotg_½i);
    Angle const Ω = cotg_½i > tg_½i ? ArcTan(equinoctial.p, equinoctial.q)
                                    : ArcTan(equinoctial.pʹ, equinoctial.qʹ);
    double const e = Sqrt(Pow<2>(equinoctial.h) + Pow<2>(equinoctial.k));
    Angle const ϖ = ArcTan(equinoctial.h, equinoctial.k);
    Angle const ω = ϖ - Ω;
    Angle const M = equinoctial.λ - ϖ;
    classical_elements.push_back(
        {.time = equinoctial.t,
         .semimajor_axis = equinoctial.a,
         .eccentricity = e,
         .inclination = i,
         .longitude_of_ascending_node = classical_elements.empty()
             ? Mod(Ω, 2 * π * Radian)
             : UnwindFrom(classical_elements.back().longitude_of_ascending_node,
                          Ω),
         .argument_of_periapsis = classical_elements.empty()
             ? Mod(ω, 2 * π * Radian)
             : UnwindFrom(classical_elements.back().argument_of_periapsis, ω),
         .mean_anomaly = classical_elements.empty()
             ? Mod(M, 2 * π * Radian)
             : UnwindFrom(classical_elements.back().mean_anomaly, M),
         .periapsis_distance = (1 - e) * equinoctial.a,
         .apoapsis_distance = (1 + e) * equinoctial.a});
  }
  return classical_elements;
}

inline absl::Status OrbitalElements::ComputePeriodsAndPrecession() {
  Time const Δt = mean_classical_elements_.back().time -
                  mean_classical_elements_.front().time;
  auto const Δt³ = Pow<3>(Δt);
  // We compute the mean rate (slope) of the mean anomaly M(t), the mean
  // argument of latitude u(t), and the longitude of the ascending node Ω(t).
  // On an interval [t̄ - Δt/2, t̄ + Δt/2], the slope of э is computed as
  //   ∫ (э(t) - э̄) (t - t̄) dt / ∫ (t - t̄)² dt;
  // this is the continuous analogue of a simple linear regression.
  // With ∫ (t - t̄)² dt = Δt³ / 12 and ∫ э̄ (t - t̄) = э̄ ∫ (t - t̄) dt = 0,
  // this simplifies to
  //   12 ∫ э(t) (t - t̄) dt / Δt³.
  // We first compute ∫ э(t) (t - t̄) dt for the three elements of interest.

  auto const interpolate_function_of_mean_classical_element =
      [this](auto const f, Instant const& t) {
        CHECK_LE(t, mean_classical_elements_.back().time);
        auto const it =
            std::partition_point(mean_classical_elements_.begin(),
                                 mean_classical_elements_.end(),
                                 [&t](ClassicalElements const& elements) {
                                   return elements.time < t;
                                 });
        ClassicalElements const& high = *it;
        if (it == mean_classical_elements_.begin()) {
          return f(high);
        } else {
          ClassicalElements const& low = *std::prev(it);
          double const α = (t - low.time) / (high.time - low.time);
          return f(low) + α * (f(high) - f(low));
        }
      };

  Instant const t̄ = mean_classical_elements_.front().time + Δt / 2;

  Product<Angle, Square<Time>> const ʃ_Mt_dt = AutomaticClenshawCurtis(
      [&interpolate_function_of_mean_classical_element, &t̄](Instant const& t) {
        return interpolate_function_of_mean_classical_element(
            [&t, &t̄](ClassicalElements const& elements) {
              return elements.mean_anomaly * (t - t̄);
            },
            t);
      },
      mean_classical_elements_.front().time,
      mean_classical_elements_.back().time,
      max_clenshaw_curtis_relative_error,
      /*max_points=*/mean_classical_elements_.size());
  Product<Angle, Square<Time>> const ʃ_ut_dt = AutomaticClenshawCurtis(
      [&interpolate_function_of_mean_classical_element, &t̄](Instant const& t) {
        return interpolate_function_of_mean_classical_element(
            [&t, &t̄](ClassicalElements const& elements) {
              return (elements.argument_of_periapsis + elements.mean_anomaly) *
                     (t - t̄);
            },
            t);
      },
      mean_classical_elements_.front().time,
      mean_classical_elements_.back().time,
      max_clenshaw_curtis_relative_error,
      /*max_points=*/mean_classical_elements_.size());
  Product<Angle, Square<Time>> const ʃ_Ωt_dt = AutomaticClenshawCurtis(
      [&interpolate_function_of_mean_classical_element, &t̄](Instant const& t) {
        return interpolate_function_of_mean_classical_element(
            [&t, &t̄](ClassicalElements const& elements) {
              return elements.longitude_of_ascending_node * (t - t̄);
            },
            t);
      },
      mean_classical_elements_.front().time,
      mean_classical_elements_.back().time,
      max_clenshaw_curtis_relative_error,
      /*max_points=*/mean_classical_elements_.size());

  // The periods are 2π over the mean rate of the relevant element; the nodal
  // precession is the mean rate of Ω.

  anomalistic_period_ = 2 * π * Radian * Δt³ / (12 * ʃ_Mt_dt);
  nodal_period_ = 2 * π * Radian * Δt³ / (12 * ʃ_ut_dt);
  nodal_precession_ = 12 * ʃ_Ωt_dt / Δt³;
  return absl::OkStatus();
}

inline absl::Status OrbitalElements::ComputeIntervals() {
  for (auto const& elements : mean_classical_elements_) {
    RETURN_IF_STOPPED;
    mean_semimajor_axis_interval_.Include(elements.semimajor_axis);
    mean_eccentricity_interval_.Include(elements.eccentricity);
    mean_inclination_interval_.Include(elements.inclination);
    mean_longitude_of_ascending_node_interval_.Include(
        elements.longitude_of_ascending_node);
    mean_argument_of_periapsis_interval_.Include(
        elements.argument_of_periapsis);
    mean_periapsis_distance_interval_.Include(elements.periapsis_distance);
    mean_apoapsis_distance_interval_.Include(elements.apoapsis_distance);
  }
  return absl::OkStatus();
}


}  // namespace internal
}  // namespace _orbital_elements
}  // namespace astronomy
}  // namespace principia
