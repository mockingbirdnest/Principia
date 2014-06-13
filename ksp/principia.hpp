#pragma once

#include "physics/body.hpp"
#include "physics/n_body_system.hpp"

namespace principia {
namespace ksp {

// The reference frame in which the integration is performed.
struct IntegrationFrame;

[KSPAddon(KSPAddon::Startup::Flight, false)]
public ref class Principia : public UnityEngine::MonoBehaviour {
 private:

  // Unity event functions, accessed through reflection.
  // See http://docs.unity3d.com/Manual/ExecutionOrder.html.

  // Called before the first frame update.
  void Start();
  // Unity physics step.
  void FixedUpdate();
  // Called once per frame.
  void Update();
  // Rendering (called before camera culling).
  void OnPreCull();
  // GUI rendering.
  void OnGUI();
  // Called after all frame updates for the last frame of existence.
  void OnDestroy();

  void DrawGUI();
  void DrawMainWindow(int window_id);
  void DrawReferenceFrameWindow(int window_id);

  void SetUpSystem();

  UnityEngine::Rect main_window_position_;
  UnityEngine::Rect reference_frame_window_position_;

  UnityEngine::GUIStyle^ gui_style;

  bool simulating_;

  // We own this pointer, but we cannot use |std::unique_ptr| because that's
  // an unmanaged type.
  physics::NBodySystem<IntegrationFrame>* system_ = nullptr;
};

}  // namespace ksp
}  // namespace principia
