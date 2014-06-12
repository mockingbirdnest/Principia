#include "ksp/principia.hpp"

#include "physics/body.hpp"
#include "physics/n_body_system.hpp"

using principia::physics::NBodySystem;

namespace principia {
namespace ksp {

void Principia::Start() {
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

}

void Principia::DrawReferenceFrameWindow(int window_id) {

}

}  // namespace ksp
}  // namespace principia
