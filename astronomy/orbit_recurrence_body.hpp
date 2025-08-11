#pragma once

#include "astronomy/orbit_recurrence.hpp"

#include <cmath>
#include <limits>
#include <numeric>

#include "base/mod.hpp"
#include "base/status_utilities.hpp"  // 🧙 For RETURN_IF_ERROR.
#include "geometry/sign.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace astronomy {
namespace _orbit_recurrence {
namespace internal {

using namespace principia::base::_mod;
using namespace principia::geometry::_sign;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_si;

inline absl::StatusOr<int> SafeNearbyInt(double const x) {
  if (std::isfinite(x)) {
    if (x < std::numeric_limits<int>::lowest() ||
        x > std::numeric_limits<int>::max()) {
      return absl::InvalidArgumentError(
          "Floating-point value too large or too small");
    } else {
      return std::nearbyint(x);
    }
  } else {
    return absl::InvalidArgumentError("Nonfinite floating-point value");
  }
}

inline OrbitRecurrence::OrbitRecurrence(int const νₒ,
                                        int const Dᴛₒ,
                                        int const Cᴛₒ)
    : νₒ_(νₒ), Dᴛₒ_(Dᴛₒ), Cᴛₒ_(Cᴛₒ) {
  CHECK_NE(Cᴛₒ, 0) << *this;
  Sign const sign_Cᴛₒ = Sign::OfNonZero(Cᴛₒ);
  if (νₒ != 0) {
    CHECK_EQ(Sign::OfNonZero(νₒ), sign_Cᴛₒ) << *this;
  }
  CHECK_LE(Abs(2 * Dᴛₒ), Abs(Cᴛₒ)) << *this;
  CHECK_EQ(std::gcd(Dᴛₒ, Cᴛₒ), 1) << *this;

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
absl::StatusOr<OrbitRecurrence> OrbitRecurrence::ClosestRecurrence(
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
    auto const status_or_int_abs_κ_J = SafeNearbyInt(abs_κ_J);
    RETURN_IF_ERROR(status_or_int_abs_κ_J);
    double const frac_abs_κ_J = Abs(abs_κ_J - status_or_int_abs_κ_J.value());
    if (frac_abs_κ_J < min_frac_abs_κ_J) {
      min_frac_abs_κ_J = frac_abs_κ_J;
      Cᴛₒ = Sign(κ) * J;
    }
  }

  auto const status_or_νₒ = SafeNearbyInt(κ);
  RETURN_IF_ERROR(status_or_νₒ);
  int const νₒ = status_or_νₒ.value();
  auto const status_or_Dᴛₒ = SafeNearbyInt((κ - νₒ) * Cᴛₒ);
  RETURN_IF_ERROR(status_or_Dᴛₒ);
  int const Dᴛₒ = status_or_Dᴛₒ.value();
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
  double const κ⁻¹ = Cᴛₒ_ / Nᴛₒ;
  // See (8.24).
  return -2 * π * Radian * κ⁻¹;
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

}  // namespace internal
}  // namespace _orbit_recurrence
}  // namespace astronomy
}  // namespace principia
