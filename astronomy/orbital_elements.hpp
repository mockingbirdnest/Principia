#pragma once

#include <vector>

#include "base/status_or.hpp"
#include "geometry/interval.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/body.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbital_elements {

using base::StatusOr;
using geometry::Instant;
using geometry::Interval;
using physics::Body;
using physics::DiscreteTrajectory;
using physics::MassiveBody;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Difference;
using quantities::Infinity;
using quantities::Length;
using quantities::Time;

class OrbitalElements {
 public:
  template<typename PrimaryCentred>
  static StatusOr<OrbitalElements> ForTrajectory(
      DiscreteTrajectory<PrimaryCentred> const& trajectory,
      MassiveBody const& primary,
      Body const& secondary);

  // The classical Keplerian elements (a, e, i, Ω, ω, M),
  // together with an epoch.
  struct ClassicalElements {
    Instant time;
    Length semimajor_axis;
    double eccentricity;
    Angle inclination;
    Angle longitude_of_ascending_node;
    Angle argument_of_periapsis;
    Angle mean_anomaly;
  };

  // Mean element time series.  These elements are free of short-period
  // variations, i.e., variations whose period is the orbital period.
  std::vector<ClassicalElements> const& mean_elements() const;

  // The period of the (osculating) mean longitude λ = Ω + ω + M.
  // Note that since our mean elements are filtered by integration over this
  // period, it does not make much sense to recompute it based on our mean
  // elements.
  Time sidereal_period() const;
  // The period of the (mean) mean argument of latitude u = ω + M.
  Time nodal_period() const;
  // The period of the (mean) mean anomaly M.
  Time anomalistic_period() const;

  // The rate of precession of Ω.
  AngularFrequency nodal_precession() const;

  // NOTE(egg): The argument of periapsis ω typically precesses as well.
  // However, long-period variations tend to be comparatively large, so that a
  // precession rate computed over a few orbits would be highly inaccurate.
  // More importantly, whereas the actual value of Ω′ is relevant to, e.g.,
  // orbit recurrence computation or sun-synchronicity, one typically cares
  // about ω′ only when requiring that ω′ be 0 (in a frozen orbit), in which
  // case the more relevant requirement is that ω stays close to some reference
  // value.

  // Of the mean classical elements (a, e, i, Ω, ω, M), under the influence of
  // gravitational forces,
  // — M always exhibits a fast secular variation (anomalistic mean motion);
  // — Ω often exhibits a secular variation (nodal precession); there are
  //   however rare cases where it is kept constant (so-called inertial orbits
  //   that achieve Ω′ = 0 by being polar, e.g., CoRoT or Gravity Probe B); in
  //   that case, the frozen value may occasionally be relevant: for CoRoT, it
  //   determines the region of the sky that may be observed.
  // — ω exhibits a secular variation, except for frozen orbits or orbits at the
  //   critical inclination; For frozen orbits (type II frozen orbits in the
  //   terminology of Ulrich Walter (2018), Astronautics: The Physics of Space
  //   Flight), its constant value must be either 90° or 270°; for orbits at the
  //   critical inclination (type I frozen orbits), ω is arbitrary; in highly
  //   eccentric cases, it is often chosen to be 270° so that the apogee is at
  //   high latitudes (Молния, みちびき, etc.).
  // — a, e, i exhibit no secular variation.
  // However, the elements that exhibit no secular variation still have
  // long-period variations; instead of trying to characterize these complex
  // effects, we provide the interval of values taken by these elements over the
  // trajectory being analysed.

  Interval<Length> mean_semimajor_axis_interval() const;
  Interval<double> mean_eccentricity_interval() const;
  Interval<Angle> mean_inclination_interval() const;
  Interval<Angle> mean_longitude_of_ascending_node_interval() const;
  Interval<Angle> mean_argument_of_periapsis_interval() const;

  // The equinoctial elements, and in particular the osculating equinoctial
  // elements, are not directly interesting; anything that could be derived from
  // them should be directly computed by this class instead.  They are however
  // useful for experimentation in Mathematica, to see whether the
  // transformation from osculating to mean elements is well-behaved, whether
  // the mean elements are stable, and what useful quantities can be derived
  // from the mean elements.

  // The equinoctial elements, together with an epoch.
  // See Broucke and Cefola (1972), On the equinoctial orbit elements.
  struct EquinoctialElements {
    Instant t;  // The epoch of the elements.
    Length a;   // The semimajor axis.
    double h;   // e sin ϖ = e sin (Ω + ω).
    double k;   // e cos ϖ = e cos (Ω + ω).
    Angle λ;    // The mean longitude ϖ + M = Ω + ω + M.
    double p;   // tg i/2 sin Ω.
    double q;   // tg i/2 cos Ω.
    // pʹ and qʹ use the cotangent of the half-inclination instead of its
    // tangent; they are better suited to retrograde orbits.
    double pʹ;  // cotg i/2 sin Ω.
    double qʹ;  // cotg i/2 cos Ω.
  };

  std::vector<EquinoctialElements> const& osculating_equinoctial_elements()
      const;
  std::vector<EquinoctialElements> const& mean_equinoctial_elements() const;

  class TentativeElements {
   public:
    template<typename PrimaryCentred>
    static StatusOr<TentativeElements> ForTrajectory(
        DiscreteTrajectory<PrimaryCentred> const& trajectory,
        MassiveBody const& primary,
        Body const& secondary);

    Time sidereal_period() const;

   private:
    std::vector<EquinoctialElements> osculating_equinoctial_elements_;
    Time sidereal_period_;
    friend class OrbitalElements;
  };

 private:
  OrbitalElements() = default;

  // Returns the osculating equinoctial elements for each point in |trajectory|.
  // If at any point the osculating orbit is hyperbolic (resulting in NaN
  // elements), that NaN is appended to the elements, and further elements are
  // omitted.
  // Thus, only the last element of the result may be NaN; if it is not NaN, the
  // size of the result equals |trajectory.size()|.
  template<typename PrimaryCentred>
  static std::vector<EquinoctialElements> OsculatingEquinoctialElements(
      DiscreteTrajectory<PrimaryCentred> const& trajectory,
      MassiveBody const& primary,
      Body const& secondary);

  // |equinoctial_elements| must contain at least 2 elements.
  static Time SiderealPeriod(
      std::vector<EquinoctialElements> const& equinoctial_elements);

  // |osculating| must contain at least 2 elements.
  // The resulting elements are averaged over one period, centred on
  // their |EquinoctialElements::t|.
  static std::vector<EquinoctialElements> MeanEquinoctialElements(
      std::vector<EquinoctialElements> const& osculating,
      Time const& period);

  static std::vector<ClassicalElements> ToClassicalElements(
      std::vector<EquinoctialElements> const& equinoctial_elements);

  // |mean_classical_elements_| must have been computed; sets
  // |anomalistic_period_|, |nodal_period_|, and |nodal_precession_|
  // accordingly. Note that this does not compute |sidereal_period_| (our mean
  // element computation is based on it, so it gets computed earlier).
  void ComputePeriodsAndPrecession();

  // |mean_classical_elements_| must have been computed; sets
  // |mean_*_interval_| accordingly.
  void ComputeMeanElementIntervals();

  std::vector<EquinoctialElements> osculating_equinoctial_elements_;
  Time sidereal_period_;
  std::vector<EquinoctialElements> mean_equinoctial_elements_;
  std::vector<ClassicalElements> mean_classical_elements_;
  Time anomalistic_period_;
  Time nodal_period_;
  AngularFrequency nodal_precession_;

  Interval<Length> mean_semimajor_axis_interval_;
  Interval<double> mean_eccentricity_interval_;
  Interval<Angle> mean_inclination_interval_;
  Interval<Angle> mean_longitude_of_ascending_node_interval_;
  Interval<Angle> mean_argument_of_periapsis_interval_;
};

}  // namespace internal_orbital_elements

using internal_orbital_elements::OrbitalElements;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/orbital_elements_body.hpp"
