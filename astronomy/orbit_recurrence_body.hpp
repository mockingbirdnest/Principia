#pragma once

#include "astronomy/orbit_recurrence.hpp"

#include <limits>
#include <numeric>

#include "base/mod.hpp"
#include "geometry/sign.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbit_recurrence {

using base::mod;
using geometry::Sign;
using quantities::si::Radian;
using quantities::Abs;

inline OrbitRecurrence::OrbitRecurrence(int const νₒ,
                                        int const Dᴛₒ,
                                        int const Cᴛₒ)
    : νₒ_(νₒ), Dᴛₒ_(Dᴛₒ), Cᴛₒ_(Cᴛₒ) {
  CHECK_NE(Cᴛₒ, 0) << *this;
  if (νₒ != 0) {
    CHECK_EQ(Sign(νₒ), Sign(Cᴛₒ)) << *this;
  }
  CHECK_LE(Abs(2 * Dᴛₒ), Abs(Cᴛₒ)) << *this;
  CHECK_EQ(std::gcd(Dᴛₒ, Cᴛₒ), 1) << *this;

  int const sign_Cᴛₒ = Sign(Cᴛₒ) * 1;
  int& Eᴛₒ = subcycle_;
  if (Dᴛₒ == 0) {
    Eᴛₒ = 0;
    return;
  }
  Eᴛₒ = sign_Cᴛₒ;
  // See 11.5.3; the termination condition is (11.25).
  // By trying the values in ascending order, we get the smallest solution Eᴛₒ*
  // for Eᴛₒ.
  while (
      !(mod(Eᴛₒ * Dᴛₒ, Cᴛₒ) == sign_Cᴛₒ || mod(Eᴛₒ * Dᴛₒ, -Cᴛₒ) == -sign_Cᴛₒ)) {
    Eᴛₒ += sign_Cᴛₒ;
  }
}

template<typename Frame>
OrbitRecurrence OrbitRecurrence::ClosestRecurrence(
    Time const& nodal_period,
    AngularFrequency const& nodal_precession,
    RotatingBody<Frame> const& primary,
    int max_abs_Cᴛₒ) {
  AngularFrequency const& Ωʹ = nodal_precession;
  AngularFrequency const& Ωʹᴛ = primary.angular_frequency();

  // Nodal mean motion.
  AngularFrequency const& nd = 2 * π * Radian / nodal_period;
  // Daily recurrence frequency, see (7.41).
  double const κ = nd / (Ωʹᴛ - Ωʹ);
  // Look for the closest rational approximation Nᴛₒ / Cᴛₒ to κ whose
  // denominator is at most max_Cᴛₒ.
  // The notation follows section 11.7.2.
  int Cᴛₒ;
  double min_frac_abs_κ_J = std::numeric_limits<double>::infinity();
  for (int J = 1; J <= max_abs_Cᴛₒ; ++J) {
    double const abs_κ_J = Abs(κ * J);
    double const frac_abs_κ_J = Abs(abs_κ_J - std::nearbyint(abs_κ_J));
    if (frac_abs_κ_J < min_frac_abs_κ_J) {
      min_frac_abs_κ_J = frac_abs_κ_J;
      Cᴛₒ = Sign(κ) * J;
    }
  }

  int const νₒ = std::nearbyint(κ);
  int const Dᴛₒ = std::nearbyint((κ - νₒ) * Cᴛₒ);
  return OrbitRecurrence(νₒ, Dᴛₒ, Cᴛₒ);
}

inline int OrbitRecurrence::νₒ() const {
  return νₒ_;
}

inline int OrbitRecurrence::Dᴛₒ() const {
  return Dᴛₒ_;
}

inline int OrbitRecurrence::Cᴛₒ() const {
  return Cᴛₒ_;
}

inline int OrbitRecurrence::number_of_revolutions() const {
  // See (11.13).
  return νₒ_ * Cᴛₒ_ + Dᴛₒ_;
}

inline Angle OrbitRecurrence::equatorial_shift() const {
  double const Nᴛₒ = number_of_revolutions();
  double const ⅟κ = Cᴛₒ_ / Nᴛₒ;
  // See (8.24).
  return -2 * π * Radian * ⅟κ;
}

inline Angle OrbitRecurrence::base_interval() const {
  return Abs(equatorial_shift());
}

inline Angle OrbitRecurrence::grid_interval() const {
  int const Nᴛₒ = number_of_revolutions();
  // See (11.20).
  return 2 * π * Radian / Nᴛₒ;
}

inline int OrbitRecurrence::subcycle() const {
  return subcycle_;
}

inline bool operator==(OrbitRecurrence const& left,
                       OrbitRecurrence const& right) {
  return left.νₒ() == right.νₒ() &&
         left.Dᴛₒ() == right.Dᴛₒ() &&
         left.Cᴛₒ() == right.Cᴛₒ();
}

inline bool operator!=(OrbitRecurrence const& left,
                       OrbitRecurrence const& right) {
  return left.νₒ() != right.νₒ() ||
         left.Dᴛₒ() != right.Dᴛₒ() ||
         left.Cᴛₒ() != right.Cᴛₒ();
}

inline std::ostream& operator<<(std::ostream& out,
                                OrbitRecurrence const& recurrence) {
  return out << "[" << recurrence.νₒ() << "; " << recurrence.Dᴛₒ() << "; "
             << recurrence.Cᴛₒ() << "]";
}

}  // namespace internal_orbit_recurrence
}  // namespace astronomy
}  // namespace principia
