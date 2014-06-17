#include "ksp/principia.hpp"

#include <map>
#include <string>
#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/permutation.hpp"
#include "geometry/rotation.hpp"
#include "ksp/utilities.hpp"
#include "physics/body.hpp"
#include "physics/n_body_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

using principia::geometry::DebugString;
using principia::geometry::Permutation;
using principia::geometry::Quaternion;
using principia::geometry::Rotation;
using principia::geometry::Vector;
using principia::physics::Body;
using principia::physics::NBodySystem;
using principia::quantities::GravitationalParameter;
using principia::quantities::Length;
using principia::quantities::SIUnit;
using principia::quantities::DebugString;
using principia::si::Hour;
using principia::si::Second;

#define LOG_UNITY(message) UnityEngine::Debug::Log(                            \
    gcnew System::String(                                                      \
        (std::string(__FILE__) + ":" + std::to_string(__LINE__) + ", in " +    \
         __FUNCSIG__ + "]" + (message)).c_str()));

namespace principia {
namespace ksp {

Principia::~Principia() {
  FreeAllOwnedPointers();
}

Principia::!Principia() {
  FreeAllOwnedPointers();
}

void Principia::Start() {
  gui_style = gcnew UnityEngine::GUIStyle(UnityEngine::GUI::skin->button);

  gui_style->normal->textColor    = UnityEngine::Color::white;
  gui_style->focused->textColor   = UnityEngine::Color::white;
  gui_style->hover->textColor     = UnityEngine::Color::yellow;
  gui_style->active->textColor    = UnityEngine::Color::yellow;
  gui_style->onNormal->textColor  = UnityEngine::Color::green;
  gui_style->onFocused->textColor = UnityEngine::Color::green;
  gui_style->onHover->textColor   = UnityEngine::Color::green;
  gui_style->onActive->textColor  = UnityEngine::Color::green;
  gui_style->padding              = gcnew UnityEngine::RectOffset(8, 8, 8, 8);

  RenderingManager::AddToPostDrawQueue(
      3,                                           // queueSpot
      gcnew Callback(this, &Principia::DrawGUI));  // drawFunction
  main_window_position_ = UnityEngine::Rect(
      UnityEngine::Screen::width / 2.0f,   // left
      UnityEngine::Screen::height / 2.0f,  // top
      10.0f,                               // width
      10.0f);                              // height
  reference_frame_window_position_ = UnityEngine::Rect(
      UnityEngine::Screen::width / 3.0f,   // left
      UnityEngine::Screen::height / 3.0f,  // top
      10.0f,                               // width
      10.0f);                              // height
  integrator_->Initialize(integrator_->Order5Optimal());
}

void Principia::FixedUpdate() {
  if (simulating_) {
    if (system_->bodies().front()->times().back() - UniversalTime() >
            *Δt_ * sampling_period_) {
      system_->Integrate(*integrator_, UniversalTime(), *Δt_, sampling_period_);
    }
  }
}

void Principia::Update() {
  // TODO(egg): draw trajectories.
}

void Principia::OnDestroy() {
  FreeAllOwnedPointers();
}

void Principia::DrawGUI() {
  UnityEngine::GUI::skin = HighLogic::Skin;
  main_window_position_ = UnityEngine::GUILayout::Window(
      1,  // id
      main_window_position_,
      gcnew UnityEngine::GUI::WindowFunction(this, &Principia::DrawMainWindow),
      "Traces of Various Descriptions",
      UnityEngine::GUILayout::MinWidth(500));
  if (simulating_) {
    reference_frame_window_position_ = UnityEngine::GUILayout::Window(
        2,  // id
        reference_frame_window_position_,
        gcnew UnityEngine::GUI::WindowFunction(
            this,
            &Principia::DrawReferenceFrameWindow),
        "Reference Frame",
        UnityEngine::GUILayout::MinWidth(500));
  }
}

void Principia::DrawMainWindow(int window_id) {
  UnityEngine::GUILayout::BeginVertical();

  if (UnityEngine::GUILayout::Button(
          simulating_ ? "Switch to Keplerian" : "Switch to Newtonian"),
          gui_style,
          UnityEngine::GUILayout::ExpandWidth(true)) {
    simulating_ = !simulating_;
    if (simulating_) {
      SetUpSystem();
    }
  }

  UnityEngine::GUILayout::BeginHorizontal();
  UnityEngine::GUILayout::Label("prediction_length_ =");
  System::String^ prediction_length_string = UnityEngine::GUILayout::TextArea(
      (*prediction_length_ / Hour).ToString(),  // text
      UnityEngine::GUILayout::Width(250));      // options
  double prediction_length_hours;
  if (System::Double::TryParse(prediction_length_string,
                               prediction_length_hours)) {
    *prediction_length_ = prediction_length_hours * Hour;
  }
  UnityEngine::GUILayout::Label("h");
  UnityEngine::GUILayout::EndHorizontal();

  UnityEngine::GUILayout::BeginHorizontal();
  UnityEngine::GUILayout::Label(L"Δt_ =");
  System::String^ timestep_string = UnityEngine::GUILayout::TextArea(
      (*Δt_ / Second).ToString(),           // text
      UnityEngine::GUILayout::Width(250));  // options
  double timestep_seconds;
  if (System::Double::TryParse(timestep_string, timestep_seconds)) {
    *Δt_ = timestep_seconds * Second;
  }
  UnityEngine::GUILayout::Label("s");
  UnityEngine::GUILayout::EndHorizontal();

  UnityEngine::GUILayout::BeginHorizontal();
  UnityEngine::GUILayout::Label("sampling_period_ =");
  System::String^ sampling_period_string = UnityEngine::GUILayout::TextArea(
      sampling_period_.ToString(),          // text
      UnityEngine::GUILayout::Width(250));  // options
  int new_sampling_period;
  if (System::Int32::TryParse(sampling_period_string, new_sampling_period)) {
    sampling_period_ = new_sampling_period;
  }
  UnityEngine::GUILayout::EndHorizontal();

  UnityEngine::GUILayout::EndVertical();
  UnityEngine::GUI::DragWindow(UnityEngine::Rect(0.0f, 0.0f, 10000.0f, 20.0f));
}

// Layout:
// o Sun        o Surface
// o Planet 1   o Centric
// o Planet 2   o Barycentric rotating
// o Planet 3
// ...
void Principia::DrawReferenceFrameWindow(int window_id) {
  UnityEngine::GUILayout::BeginHorizontal();
  UnityEngine::GUILayout::BeginVertical(UnityEngine::GUILayout::Width(250));
  for each (CelestialBody^ body in FlightGlobals::Bodies) {
    if (UnityEngine::GUILayout::Toggle(
            rendering_reference_body_ == body,  // value
            body->name) &&                     // text
        rendering_reference_body_ != body) {
      rendering_reference_body_ = body;
    }
  }
  UnityEngine::GUILayout::EndVertical();
  UnityEngine::GUILayout::BeginVertical(
      UnityEngine::GUILayout::ExpandWidth(true));
  for (auto const& type : kRenderingFrameTypes) {
    // Do not show the barycentric option for the sun, since the sun has no
    // primary.
    if ((rendering_reference_body_ != sun_ || type != kBarycentric) &&
        UnityEngine::GUILayout::Toggle(
            rendering_frame_type_ == type,        // value
            Manage(kFrameTypeNames.at(type))) &&  // text
        rendering_frame_type_ != type) {
      rendering_frame_type_ = type;
    }
  }
  UnityEngine::GUILayout::EndVertical();
  UnityEngine::GUILayout::EndHorizontal();
  UnityEngine::GUI::DragWindow(UnityEngine::Rect(0.0f, 0.0f, 10000.0f, 20.0f));
}

Displacement<IntegrationFrame> IntegrationPosition(CelestialBody^ const body);
VelocityChange<IntegrationFrame> IntegrationVelocity(CelestialBody^ const body);

Displacement<IntegrationFrame> IntegrationPosition(Vessel^ const body);
VelocityChange<IntegrationFrame> IntegrationVelocity(Vessel^ const body);

void Principia::SetUpSystem() {
  *celestials_ = std::map<std::string, Body<IntegrationFrame>*>();
  *vessels_    = std::map<std::string, Body<IntegrationFrame>*>();
  Time const universal_time = UniversalTime();
  std::vector<Body<IntegrationFrame>*>* const bodies =
      new std::vector<Body<IntegrationFrame>*>;
  sun_ = FlightGlobals::Bodies[0];
  rendering_reference_body_ = sun_;
  for each (CelestialBody^ body in FlightGlobals::Bodies) {
    GravitationalParameter const gravitational_parameter =
        body->gravParameter * SIUnit<GravitationalParameter>();
    bodies->push_back(new Body<IntegrationFrame>(gravitational_parameter));
    Displacement<IntegrationFrame> const   position =
      IntegrationPosition(body) - IntegrationPosition(sun_);
    VelocityChange<IntegrationFrame> const velocity =
      IntegrationVelocity(body) - IntegrationVelocity(sun_);
    bodies->back()->AppendToTrajectory({position}, {velocity}, universal_time);
    celestials_->insert({Unmanage(body->name),
                         bodies->back()});
    LOG_UNITY(std::string("\nAdded CelestialBody ") + Unmanage(body->name) +
              "\nGM = " + DebugString(gravitational_parameter) +
              "\nq  = " + DebugString(position) +
              "\nv  = " + DebugString(velocity));
  }
  for each (Vessel^ body in FlightGlobals::Vessels) {
    bodies->push_back(new Body<IntegrationFrame>(GravitationalParameter()));
    Displacement<IntegrationFrame> const   position =
        IntegrationPosition(body) - IntegrationPosition(sun_);
    VelocityChange<IntegrationFrame> const velocity =
        IntegrationVelocity(body) - IntegrationVelocity(sun_);
    bodies->back()->AppendToTrajectory({position}, {velocity}, universal_time);
    LOG_UNITY(std::string("\nAdded Vessel ") + Unmanage(body->name) +
              "\nq  = " + DebugString(position) +
              "\nv  = " + DebugString(velocity));
  }
  Reset(system_, new NBodySystem<IntegrationFrame>(bodies));
}

void Principia::FreeAllOwnedPointers() {
  Reset(prediction_length_, nullptr);
  Reset(Δt_, nullptr);
  Reset(system_, nullptr);
  Reset(celestials_, nullptr);
  Reset(vessels_, nullptr);
}

struct World;
struct AliceWorld;

Permutation<AliceWorld, World> const LookingGlass =
    Permutation<AliceWorld, World>(Permutation<AliceWorld, World>::XZY);

Rotation<IntegrationFrame, World> PlanetariumRotation() {
  UnityEngine::QuaternionD planetarium_rotation = Planetarium::Rotation;
  return Rotation<IntegrationFrame, World>(
      Quaternion(planetarium_rotation.w,
                 {planetarium_rotation.x,
                  planetarium_rotation.y,
                  planetarium_rotation.z}));
}

Displacement<World> WorldPosition(CelestialBody^ const body) {
  Vector3d const position = body->position;
  return Displacement<World>({
      position.x * SIUnit<Length>(),
      position.y * SIUnit<Length>(),
      position.z * SIUnit<Length>()});
}

Displacement<IntegrationFrame> IntegrationPosition(CelestialBody^ const body) {
  return PlanetariumRotation().Inverse()(WorldPosition(body));
}

VelocityChange<World> WorldVelocity(CelestialBody^ const body) {
  Vector3d const frame_velocity = body->GetFrameVel();
  return VelocityChange<World>({
      frame_velocity.x * SIUnit<Speed>(),
      frame_velocity.y * SIUnit<Speed>(),
      frame_velocity.z * SIUnit<Speed>()});
}

VelocityChange<IntegrationFrame> IntegrationVelocity(
    CelestialBody^ const body) {
  return PlanetariumRotation().Inverse()(WorldVelocity(body));
}

Displacement<AliceWorld> OffsetFromReferenceBody(Vessel^ const body) {
  Vector3d const result = body->orbitDriver->pos;
    return Displacement<AliceWorld>({
      result.x * SIUnit<Length>(),
      result.y * SIUnit<Length>(),
      result.z * SIUnit<Length>()});
}

Displacement<World> WorldPosition(Vessel^ const body) {
  return WorldPosition(body->orbitDriver->orbit->referenceBody) +
      LookingGlass(OffsetFromReferenceBody(body));
}

Displacement<IntegrationFrame> IntegrationPosition(Vessel^ const body) {
  return PlanetariumRotation().Inverse()(WorldPosition(body));
}

VelocityChange<World> WorldVelocity(Vessel^ const body) {
  Vector3d const frame_velocity = body->orbitDriver->orbit->GetFrameVel();
  return VelocityChange<World>({
      frame_velocity.x * SIUnit<Speed>(),
      frame_velocity.y * SIUnit<Speed>(),
      frame_velocity.z * SIUnit<Speed>()});
}

VelocityChange<IntegrationFrame> IntegrationVelocity(Vessel^ const body) {
  return PlanetariumRotation().Inverse()(WorldVelocity(body));
}

}  // namespace ksp
}  // namespace principia
