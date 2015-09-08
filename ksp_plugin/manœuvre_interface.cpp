#include "ksp_plugin/manœuvre_interface.hpp"

#include "quantities/constants.hpp"
#include "quantities/si.hpp"

namespace principia {

using constants::StandardGravity;
using geometry::RadiusLatitudeLongitude;
using si::Degree;
using si::Kilogram;
using si::Metre;
using si::Newton;
using si::Second;

namespace ksp_plugin {

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
      Vector<double, Barycentric>(
          RadiusLatitudeLongitude(
              1.0,
              declination * Degree,
              right_ascension * Degree).ToCartesian())).release();
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

double principia__right_ascension(Manœuvre<Barycentric> const* manœuvre) {
  return CHECK_NOTNULL(manœuvre)->
             direction().coordinates().ToSpherical().longitude / Degree;
}

double principia__declination(Manœuvre<Barycentric> const* manœuvre) {
  return CHECK_NOTNULL(manœuvre)->
             direction().coordinates().ToSpherical().latitude / Degree;
}

void principia__set_duration(Manœuvre<Barycentric>* manœuvre,
                             double duration) {
  CHECK_NOTNULL(manœuvre)->set_duration(duration * Second);
}

double principia__duration(Manœuvre<Barycentric> const* manœuvre) {
  return CHECK_NOTNULL(manœuvre)->duration() / Second;
}

void principia__set_Δv(Manœuvre<Barycentric>* manœuvre, double Δv) {
  CHECK_NOTNULL(manœuvre)->set_Δv(Δv * (Metre / Second));
}

double principia__Δv(Manœuvre<Barycentric> const* manœuvre) {
  return CHECK_NOTNULL(manœuvre)->Δv() / (Metre / Second);
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

double principia__final_time(Manœuvre<Barycentric> const* manœuvre) {
  return (CHECK_NOTNULL(manœuvre)->final_time() - Instant()) / Second;
}

double principia__final_mass(Manœuvre<Barycentric> const* manœuvre) {
  return CHECK_NOTNULL(manœuvre)->final_mass() / Kilogram;
}

}  // namespace ksp_plugin
}  // namespace principia
