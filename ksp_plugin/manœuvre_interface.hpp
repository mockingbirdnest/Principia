#pragma once

#include "ksp_plugin/interface.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "ksp_plugin/plugin.hpp"

namespace principia {

using quantities::SpecificImpulse;

namespace ksp_plugin {

struct FlightPlanUI {};

// The underlying persisted state of the UI for editing a single manoeuvre.
struct ManœuvreEditor {
  std::string name;

  bool show_engine_config;
  Force thrust;
  std::string thrust_text_;
  SpecificImpulse specific_impulse;
  std::string specific_impulse_text_;

  bool show_engine_mass;
  bool use_final_mass_of_previous_burn;
  Mass vehicle_mass;
  std::string vehicle_mass_text_;

  Instant initial_time;
  // These manœuvres get moved when this one moves, keeping a constant time
  // offset.
  std::list<ManœuvreUI> dependent_manœuvres;

  // Here we should have a transform together with a choice of frame field.

  Velocity<Barycentric> Δv;
  bool spherical;
  std::string Δv_x_text;
  std::string vehicle_mass_text_;
  std::string vehicle_mass_text_;


};

// TODO(egg): that constructor is going to get annoying, set everything with
// single-use mutators.
extern "C" DLLEXPORT
Manœuvre<Barycentric>* principia__NewManœuvreIspByWeight(
    double thrust,
    double initial_mass,
    double specific_impulse_by_weight,
    double right_ascension,
    double declination);

extern "C" DLLEXPORT
double principia__thrust(Manœuvre<Barycentric> const* manœuvre);
extern "C" DLLEXPORT
double principia__initial_mass(Manœuvre<Barycentric> const* manœuvre);
extern "C" DLLEXPORT
double principia__specific_impulse_by_weight(
    Manœuvre<Barycentric> const* manœuvre);

extern "C" DLLEXPORT
double principia__right_ascension(Manœuvre<Barycentric> const* manœuvre);
extern "C" DLLEXPORT
double principia__declination(Manœuvre<Barycentric> const* manœuvre);

extern "C" DLLEXPORT
void principia__set_duration(Manœuvre<Barycentric>* manœuvre,
                             double duration);
extern "C" DLLEXPORT
double principia__duration(Manœuvre<Barycentric> const* manœuvre);
extern "C" DLLEXPORT
void principia__set_Δv(Manœuvre<Barycentric>* manœuvre, double Δv);
extern "C" DLLEXPORT
double principia__Δv(Manœuvre<Barycentric> const* manœuvre);

extern "C" DLLEXPORT
double principia__initial_time(Manœuvre<Barycentric> const* manœuvre);
extern "C" DLLEXPORT
void principia__set_initial_time(Manœuvre<Barycentric>* manœuvre,
                                 double initial_time);
extern "C" DLLEXPORT
double principia__time_of_half_Δv(Manœuvre<Barycentric> const* manœuvre);
extern "C" DLLEXPORT
void principia__set_time_of_half_Δv(Manœuvre<Barycentric>* manœuvre,
                                    double time_of_half_Δv);

extern "C" DLLEXPORT
double principia__final_time(Manœuvre<Barycentric> const* manœuvre);
extern "C" DLLEXPORT
double principia__final_mass(Manœuvre<Barycentric> const* manœuvre);

}  // namespace ksp_plugin
}  // namespace principia
