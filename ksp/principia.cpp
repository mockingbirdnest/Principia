#include "ksp/principia.hpp"

#include <vector>

#include "geometry/grassmann.hpp"
#include "physics/body.hpp"
#include "physics/n_body_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"

using principia::geometry::Vector;
using principia::physics::Body;
using principia::physics::NBodySystem;
using principia::quantities::GravitationalParameter;
using principia::quantities::Length;
using principia::quantities::SIUnit;

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

void Principia::SetUpSystem() {
  double universal_time = Planetarium::GetUniversalTime();
  std::vector<Body<IntegrationFrame>*>* const bodies =
      new std::vector<Body<IntegrationFrame>*>;
  for each (CelestialBody^ body in FlightGlobals::Bodies) {
    bodies->push_back(new Body<IntegrationFrame>(
      body->gravParameter * SIUnit<GravitationalParameter>()));
    bodies->back()->AppendToTrajectory(
        {Vector<Length, IntegrationFrame>()},
        {Vector<Speed, IntegrationFrame>()},
        universal_time * SIUnit<Time>());
  }
  system_ = new NBodySystem<IntegrationFrame>(bodies);
}

struct World;

Vector<Length, IntegrationFrame> IntegrationPosition(CelestialBody^ body);
Vector<Length, World> WorldPosition(CelestialBody^ body) {
  return Vector<Length, World>({
      body->position.x * SIUnit<Length>(),
      body->position.y * SIUnit<Length>(),
      body->position.z * SIUnit<Length>()});
}

}  // namespace ksp
}  // namespace principia
