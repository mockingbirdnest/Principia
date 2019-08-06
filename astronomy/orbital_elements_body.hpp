#pragma once

#include "astronomy/orbital_elements.hpp"

#include <algorithm>
#include <vector>

#include "quantities/elementary_functions.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbital_elements {

using base::Error;
using base::Status;
using quantities::ArcTan;
using quantities::Cos;
using quantities::NaN;
using quantities::Pow;
using quantities::Product;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Square;
using quantities::Tan;
using quantities::si::Radian;

// Returns the element of {α + 2nπ| n ∈ ℤ} which is closest to |previous_angle|.
inline Angle UnwindFrom(Angle const& previous_angle, Angle const& α) {
  return α + std::nearbyint((previous_angle - α) / (2 * π * Radian)) * 2 * π *
                 Radian;
}

template<typename T>
Difference<T> OrbitalElements::Interval<T>::measure() const {
  return max >= min ? max - min : Difference<T>{};
}

template<typename T>
T OrbitalElements::Interval<T>::midpoint() const {
  return max >= min ? min + measure() / 2 : NaN<T>();
}

template<typename T>
void OrbitalElements::Interval<T>::Include(T const& x) {
  min = std::min(min, x);
  max = std::max(max, x);
}

template<typename PrimaryCentred>
StatusOr<OrbitalElements> OrbitalElements::ForTrajectory(
    DiscreteTrajectory<PrimaryCentred> const& trajectory,
    MassiveBody const& primary,
    Body const& secondary) {
  OrbitalElements result;
  if (trajectory.Size() < 2) {
    return Status(Error::INVALID_ARGUMENT,
                  trajectory.Empty() ? "empty trajectory"
                                     : "insufficient trajectory size");
  }
  result.osculating_equinoctial_elements_ =
      OsculatingEquinoctialElements(trajectory, primary, secondary);
  result.sidereal_period_ =
      SiderealPeriod(result.osculating_equinoctial_elements_);
  result.mean_equinoctial_elements_ = MeanEquinoctialElements(
      result.osculating_equinoctial_elements_, result.sidereal_period_);
  if (result.mean_equinoctial_elements_.size() < 2) {
    return Status(
        Error::OUT_OF_RANGE,
        "trajectory does not span one sidereal period: sidereal period is " +
            DebugString(result.sidereal_period_) + ", trajectory spans " +
            DebugString(trajectory.last().time() - trajectory.Begin().time()));
  }
  result.mean_classical_elements_ =
      ToClassicalElements(result.mean_equinoctial_elements_);
  result.ComputePeriodsAndPrecession();
  result.ComputeMeanElementRanges();
  return result;
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

inline OrbitalElements::Interval<Length>
OrbitalElements::mean_semimajor_axis_range() const {
  return mean_semimajor_axis_range_;
}

inline OrbitalElements::Interval<double>
OrbitalElements::mean_eccentricity_range() const {
  return mean_eccentricity_range_;
}

inline OrbitalElements::Interval<Angle>
OrbitalElements::mean_inclination_range() const {
  return mean_inclination_range_;
}

inline OrbitalElements::Interval<Angle>
OrbitalElements::mean_longitude_of_ascending_node_range() const {
  return mean_longitude_of_ascending_node_range_;
}

inline OrbitalElements::Interval<Angle>
OrbitalElements::mean_argument_of_periapsis_range() const {
  return mean_argument_of_periapsis_range_;
}

template<typename PrimaryCentred>
std::vector<OrbitalElements::EquinoctialElements>
OrbitalElements::OsculatingEquinoctialElements(
    DiscreteTrajectory<PrimaryCentred> const& trajectory,
    MassiveBody const& primary,
    Body const& secondary) {
  DegreesOfFreedom<PrimaryCentred> const primary_dof{
      PrimaryCentred::origin, Velocity<PrimaryCentred>{}};
  std::vector<EquinoctialElements> result;
  for (auto it = trajectory.Begin(); it != trajectory.End(); ++it) {
    auto const osculating_elements =
        KeplerOrbit<PrimaryCentred>(primary,
                                    secondary,
                                    it.degrees_of_freedom() - primary_dof,
                                    it.time())
            .elements_at_epoch();
    double const& e = *osculating_elements.eccentricity;
    Angle const& ϖ = *osculating_elements.longitude_of_periapsis;
    Angle const& Ω = osculating_elements.longitude_of_ascending_node;
    Angle const& M = *osculating_elements.mean_anomaly;
    Angle const& i = osculating_elements.inclination;
    double const tg_½i = Tan(i / 2);
    double const cotg_½i = 1 / tg_½i;
    result.push_back(
        {/*.t = */ it.time(),
         /*.a = */ *osculating_elements.semimajor_axis,
         /*.h = */ e * Sin(ϖ),
         /*.k = */ e * Cos(ϖ),
         /*.λ = */ result.empty() ? ϖ + M : UnwindFrom(result.back().λ, ϖ + M),
         /*.p = */ tg_½i * Sin(Ω),
         /*.q = */ tg_½i * Cos(Ω),
         /*.pʹ = */ cotg_½i * Sin(Ω),
         /*.qʹ = */ cotg_½i * Cos(Ω)});
  }
  return result;
}

inline std::vector<OrbitalElements::EquinoctialElements> const&
OrbitalElements::osculating_equinoctial_elements() const {
  return osculating_equinoctial_elements_;
}

inline std::vector<OrbitalElements::EquinoctialElements> const&
OrbitalElements::mean_equinoctial_elements() const {
  return mean_equinoctial_elements_;
}

inline Time OrbitalElements::SiderealPeriod(
    std::vector<EquinoctialElements> const& equinoctial_elements) {
  Time const Δt =
      equinoctial_elements.back().t - equinoctial_elements.front().t;
  Instant const t0 = equinoctial_elements.front().t + Δt / 2;
  Product<Angle, Square<Time>> ſ_λt_dt;

  auto const first = equinoctial_elements.begin();
  Time const first_t = first->t - t0;
  Product<Angle, Time> previous_λt = first->λ * first_t;
  for (auto previous = first, it = first + 1; it != equinoctial_elements.end();
       previous = it, ++it) {
    Time const t = it->t - t0;
    auto const λt = it->λ * t;
    Time const dt = it->t - previous->t;
    ſ_λt_dt += (λt + previous_λt) / 2 * dt;
    previous_λt = λt;
  }
  return 2 * π * Radian * Pow<3>(Δt) / (12 * ſ_λt_dt);
}

inline std::vector<OrbitalElements::EquinoctialElements>
OrbitalElements::MeanEquinoctialElements(
    std::vector<EquinoctialElements> const& osculating,
    Time const& period) {
  Instant const& t_min = osculating.front().t;
  struct IntegratedEquinoctialElements {
    Instant t_max;
    // The integrals are from t_min to t_max.
    Product<Length, Time> ſ_a_dt;
    Time ſ_h_dt;
    Time ſ_k_dt;
    Product<Angle, Time> ſ_λ_dt;
    Time ſ_p_dt;
    Time ſ_q_dt;
    Time ſ_pʹ_dt;
    Time ſ_qʹ_dt;
  };
  std::vector<IntegratedEquinoctialElements> integrals;
  integrals.push_back({t_min});
  for (auto previous = osculating.begin(), it = osculating.begin() + 1;
       it != osculating.end();
       previous = it, ++it) {
    integrals.push_back(integrals.back());
    integrals.back().t_max = it->t;
    Time const dt = it->t - previous->t;
    integrals.back().ſ_a_dt += (it->a + previous->a) / 2 * dt;
    integrals.back().ſ_h_dt += (it->h + previous->h) / 2 * dt;
    integrals.back().ſ_k_dt += (it->k + previous->k) / 2 * dt;
    integrals.back().ſ_λ_dt += (it->λ + previous->λ) / 2 * dt;
    integrals.back().ſ_p_dt += (it->p + previous->p) / 2 * dt;
    integrals.back().ſ_q_dt += (it->q + previous->q) / 2 * dt;
    integrals.back().ſ_pʹ_dt += (it->pʹ + previous->pʹ) / 2 * dt;
    integrals.back().ſ_qʹ_dt += (it->qʹ + previous->qʹ) / 2 * dt;
  }
  std::vector<EquinoctialElements> result;
  int i = 0;
  for (auto first = integrals.begin(); first != integrals.end(); ++first) {
    Instant const t_first = first->t_max;
    while (integrals[i].t_max - t_first < period) {
      if (++i == integrals.size()) {
        return result;
      }
    }
    auto const& last = integrals[i - 1];
    // |element| should be a pointer to a member of |EquinoctialElements|;
    // Integrates from |last.t_max| to |t_first + period|.
    auto ſ = [i, &period, &t_first, &osculating](auto element) {
      auto const& next_osculating = osculating[i];
      auto const& last_osculating = osculating[i - 1];
      Time const Δt = next_osculating.t - last_osculating.t;
      Time const dt = t_first + period - last_osculating.t;
      auto const element_at_end =
          last_osculating.*element +
          (next_osculating.*element - last_osculating.*element) * (dt / Δt);
      return (last_osculating.*element + element_at_end) / 2 * dt;
    };
    result.emplace_back();
    result.back().t = t_first + period / 2;
    result.back().a =
        (last.ſ_a_dt - first->ſ_a_dt + ſ(&EquinoctialElements::a)) / period;
    result.back().h =
        (last.ſ_h_dt - first->ſ_h_dt + ſ(&EquinoctialElements::h)) / period;
    result.back().k =
        (last.ſ_k_dt - first->ſ_k_dt + ſ(&EquinoctialElements::k)) / period;
    result.back().λ =
        (last.ſ_λ_dt - first->ſ_λ_dt + ſ(&EquinoctialElements::λ)) / period;
    result.back().p =
        (last.ſ_p_dt - first->ſ_p_dt + ſ(&EquinoctialElements::p)) / period;
    result.back().q =
        (last.ſ_q_dt - first->ſ_q_dt + ſ(&EquinoctialElements::q)) / period;
    result.back().pʹ =
        (last.ſ_pʹ_dt - first->ſ_pʹ_dt + ſ(&EquinoctialElements::pʹ)) / period;
    result.back().qʹ =
        (last.ſ_qʹ_dt - first->ſ_qʹ_dt + ſ(&EquinoctialElements::qʹ)) / period;
  }
  return result;
}

inline std::vector<OrbitalElements::ClassicalElements>
OrbitalElements::ToClassicalElements(
    std::vector<EquinoctialElements> const& equinoctial_elements) {
  std::vector<ClassicalElements> result;
  for (auto const& equinoctial : equinoctial_elements) {
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
    result.push_back(
        {equinoctial.t,
         equinoctial.a,
         e,
         i,
         result.empty()
             ? Ω
             : UnwindFrom(result.back().longitude_of_ascending_node, Ω),
         result.empty() ? ω
                        : UnwindFrom(result.back().argument_of_periapsis, ω),
         result.empty() ? M : UnwindFrom(result.back().mean_anomaly, M)});
  }
  return result;
}

inline void OrbitalElements::ComputePeriodsAndPrecession() {
  Time const Δt = mean_classical_elements_.back().time -
                  mean_classical_elements_.front().time;
  Instant const t0 = mean_classical_elements_.front().time + Δt / 2;
  Product<Angle, Square<Time>> ſ_Mt_dt;
  Product<Angle, Square<Time>> ſ_ut_dt;
  Product<Angle, Square<Time>> ſ_Ωt_dt;

  auto const first = mean_classical_elements_.begin();
  Time first_t = first->time - t0;
  Product<Angle, Time> previous_Mt = first->mean_anomaly * first_t;
  Product<Angle, Time> previous_ut =
      (first->argument_of_periapsis + first->mean_anomaly) * first_t;
  Product<Angle, Time> previous_Ωt =
      first->longitude_of_ascending_node * first_t;
  for (auto previous = first, it = first + 1;
       it != mean_classical_elements_.end();
       previous = it, ++it) {
    Angle const u = it->argument_of_periapsis + it->mean_anomaly;
    Angle const& M = it->mean_anomaly;
    Angle const& Ω = it->longitude_of_ascending_node;
    Time const t = it->time - t0;
    auto const Mt = M * t;
    auto const ut = u * t;
    auto const Ωt = Ω * t;
    Time const dt = it->time - previous->time;
    ſ_Mt_dt += (Mt + previous_Mt) / 2 * dt;
    ſ_ut_dt += (ut + previous_ut) / 2 * dt;
    ſ_Ωt_dt += (Ωt + previous_Ωt) / 2 * dt;
    previous_Mt = Mt;
    previous_ut = ut;
    previous_Ωt = Ωt;
  }
  anomalistic_period_ = 2 * π * Radian * Pow<3>(Δt) / (12 * ſ_Mt_dt);
  nodal_period_ = 2 * π * Radian * Pow<3>(Δt) / (12 * ſ_ut_dt);
  nodal_precession_ = 12 * ſ_Ωt_dt / Pow<3>(Δt);
}

inline void OrbitalElements::ComputeMeanElementRanges() {
  for (auto const& elements : mean_classical_elements_) {
    mean_semimajor_axis_range_.Include(elements.semimajor_axis);
    mean_eccentricity_range_.Include(elements.eccentricity);
    mean_inclination_range_.Include(elements.inclination);
    mean_longitude_of_ascending_node_range_.Include(
        elements.longitude_of_ascending_node);
    mean_argument_of_periapsis_range_.Include(elements.argument_of_periapsis);
  }
}


}  // namespace internal_orbital_elements
}  // namespace astronomy
}  // namespace principia
