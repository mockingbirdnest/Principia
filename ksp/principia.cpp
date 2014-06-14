#include "ksp/principia.hpp"

#include <string>
#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/rotation.hpp"
#include <msclr/marshal.h>
#include "physics/body.hpp"
#include "physics/n_body_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

using principia::geometry::DebugString;
using principia::geometry::Quaternion;
using principia::geometry::Rotation;
using principia::geometry::Vector;
using principia::physics::Body;
using principia::physics::NBodySystem;
using principia::quantities::GravitationalParameter;
using principia::quantities::Length;
using principia::quantities::SIUnit;

#define LOG_UNITY(message) UnityEngine::Debug::Log(                            \
    gcnew System::String(                                                      \
        (std::string(__FILE__) + ":" + std::to_string(__LINE__) + ", in " +    \
         __FUNCSIG__ + "]" + (message)).c_str()));

namespace principia {
namespace ksp {

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
}

void Principia::FixedUpdate() {

}

void Principia::Update() {
  // TODO(egg): draw trajectories.
}

void Principia::OnPreCull() {

}

void Principia::OnGUI() {

}

void Principia::OnDestroy() {
  delete system_;
}

void Foo(int, ... array<int^>^ a) {}
void Bar(System::String const^ bar);
/*
void Principia::DrawGUI() {
  Foo(1, "bar");
  Bar("bar");
}
*/

void Principia::DrawGUI() {
  UnityEngine::GUI::skin = HighLogic::Skin;
  main_window_position_ = UnityEngine::GUILayout::Window(
      1,  // id
      main_window_position_,
      gcnew UnityEngine::GUI::WindowFunction(this, &Principia::DrawMainWindow),
      "Traces of Various Descriptions",
      UnityEngine::GUILayout::MinWidth(500));
  reference_frame_window_position_ = UnityEngine::GUILayout::Window(
      2,  // id
      reference_frame_window_position_,
      gcnew UnityEngine::GUI::WindowFunction(
          this,
          &Principia::DrawReferenceFrameWindow),
      "Reference Frame",
      UnityEngine::GUILayout::MinWidth(500));
}

void Principia::DrawMainWindow(int window_id) {
  UnityEngine::GUILayout::BeginVertical();
  UnityEngine::GUILayout::TextArea("Planetarium rotation");
  if (UnityEngine::GUILayout::Button(
          simulating_ ? "Switch to Keplerian" : "Switch to Newtonian")) {
    simulating_ = !simulating_;
    if (simulating_) {
      SetUpSystem();
    }
  }
}

void Principia::DrawReferenceFrameWindow(int window_id) {

}

struct World;

template <typename Frame>
using Displacement = Vector<Length, Frame>;
template <typename Frame>
using VelocityChange = Vector<Speed, Frame>;

Displacement<IntegrationFrame> IntegrationPosition(CelestialBody^ const body);
VelocityChange<IntegrationFrame> IntegrationVelocity(CelestialBody^ const body);

void Principia::SetUpSystem() {
  Time const universal_time = Planetarium::GetUniversalTime() * SIUnit<Time>();
  std::vector<Body<IntegrationFrame>*>* const bodies =
      new std::vector<Body<IntegrationFrame>*>;
  sun_ = FlightGlobals::Bodies[0];
  for each (CelestialBody^ body in FlightGlobals::Bodies) {
    GravitationalParameter const gravitational_parameter =
        body->gravParameter * SIUnit<GravitationalParameter>();
    bodies->push_back(new Body<IntegrationFrame>(gravitational_parameter));
    Displacement<IntegrationFrame>   const position = IntegrationPosition(body);
    VelocityChange<IntegrationFrame> const velocity = IntegrationVelocity(body);
    bodies->back()->AppendToTrajectory({position}, {velocity}, universal_time);
    msclr::interop::marshal_context^ context =
      gcnew msclr::interop::marshal_context();
    LOG_UNITY(std::string("\nAdded CelestialBody ") +
              context->marshal_as<const char*>(body->name) +
              "\nGM = " + DebugString(gravitational_parameter) +
              "\nq  = " + DebugString(position) +
              "\nv  = " + DebugString(velocity));
    delete context;
  }
  system_ = new NBodySystem<IntegrationFrame>(bodies);
}

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

}  // namespace ksp
}  // namespace principia
