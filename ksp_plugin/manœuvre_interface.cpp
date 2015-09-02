#include "ksp_plugin/manœuvre_interface.hpp"

#include "quantities/constants.hpp"
#include "quantities/si.hpp"

namespace principia {

using constants::StandardGravity;
using si::Degree;
using si::Kilogram;
using si::Metre;
using si::Newton;
using si::Second;

namespace ksp_plugin {

namespace {

// TODO(egg): bad code replication! bad!

XYZ ToXYZ(R3Element<double> const& r3_element) {
  return {r3_element.x, r3_element.y, r3_element.z};
}

// Returns a unit vector pointing in the direction defined by |right_ascension|
// and |declination|, assuming the reference system is as follows:
// xy-plane: plane of the Earth's mean equator at the reference epoch
// x-axis  : out along ascending node of instantaneous plane of the Earth's
//           orbit and the Earth's mean equator at the reference epoch
// z-axis  : along the Earth mean north pole at the reference epoch
Vector<double, Barycentric> Direction(Angle const& right_ascension,
                                      Angle const& declination) {
  // Positive angles map {1, 0, 0} to the positive z hemisphere, which is north.
  // An angle of 0 keeps {1, 0, 0} on the equator.
  auto const decline = Rotation<Barycentric, Barycentric>(
                           declination,
                           Bivector<double, Barycentric>({0, -1, 0}));
  // Rotate counterclockwise around {0, 0, 1} (north), i.e., eastward.
  auto const ascend = Rotation<Barycentric, Barycentric>(
                          right_ascension,
                          Bivector<double, Barycentric>({0, 0, 1}));
  return ascend(decline(Vector<double, Barycentric>({1, 0, 0})));
}

}  // namespace

Manœuvre<Barycentric>* principia__NewManœuvreIspByWeight(
    double thrust,
    double initial_mass,
    double specific_impulse_by_weight,
    double right_ascension,
    double declination) {
  return std::make_unique<Manœuvre<Barycentric>>(
      thrust * Newton,
      initial_mass * Kilogram,
      specific_impulse_by_weight * Second * StandardGravity,
      Direction(right_ascension * Degree, declination * Degree)).release();
}

double principia__thrust(Manœuvre<Barycentric> const* manœuvre) {
  return CHECK_NOTNULL(manœuvre)->thrust() / Newton;
}
double principia__initial_mass(Manœuvre<Barycentric> const* manœuvre) {
  return CHECK_NOTNULL(manœuvre)->initial_mass() / Kilogram;
}
double principia__specific_impulse_by_weight(
    Manœuvre<Barycentric> const* manœuvre) {
  return (CHECK_NOTNULL(manœuvre)->effective_exhaust_velocity() /
          StandardGravity) /
         Second;
}

XYZ principia__Δv(Manœuvre<Barycentric> const* manœuvre) {
  CHECK_NOTNULL(manœuvre);
  return ToXYZ((manœuvre->direction() * manœuvre->Δv() /
                (Metre / Second)).coordinates());
}

void principia__set_duration(Manœuvre<Barycentric>* manœuvre,
                             double duration) {
  CHECK_NOTNULL(manœuvre)->set_duration(duration * Second);
}

void principia__set_Δv(Manœuvre<Barycentric>* manœuvre, double Δv) {
  CHECK_NOTNULL(manœuvre)->set_Δv(Δv * (Metre / Second));
}

double principia__initial_time(Manœuvre<Barycentric> const* manœuvre) {
  return (CHECK_NOTNULL(manœuvre)->initial_time() - Instant()) / Second;
}

void principia__set_initial_time(Manœuvre<Barycentric>* manœuvre,
                                 double initial_time) {
  return CHECK_NOTNULL(manœuvre)->
      set_initial_time(Instant(initial_time * Second));
}

double principia__time_of_half_Δv(Manœuvre<Barycentric> const* manœuvre) {
  return (CHECK_NOTNULL(manœuvre)->time_of_half_Δv() - Instant()) / Second;
}

void principia__set_time_of_half_Δv(Manœuvre<Barycentric>* manœuvre,
                                    double time_of_half_Δv) {
  return CHECK_NOTNULL(manœuvre)->
      set_time_of_half_Δv(Instant(time_of_half_Δv * Second));
}

}  // namespace ksp_plugin
}  // namespace principia
