#pragma once

#include <map>

#include "physics/body.hpp"
#include "physics/n_body_system.hpp"

namespace principia {
namespace ksp {

// The reference frame in which the integration is performed.
struct IntegrationFrame;

// The type of frame in which the orbits are rendered.
// Examples given for Earth as the rendering reference body.
enum class RenderingFrameType {
  kSurface,     // e.g., Earth surface.
  kCentric,     // e.g., Geocentric.
  kBarycentric  // e.g., Sun-Earth barycentre.
};

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

  CelestialBody^ sun_;
  CelestialBody^ renderer_reference_body_;

  // We own these pointers, but we cannot use |std::unique_ptr| because that's
  // an unmanaged type.
  physics::NBodySystem<IntegrationFrame>* system_ = nullptr;
  std::map<std::string, physics::Body<IntegrationFrame>*>* celestials_;
  std::map<std::string, physics::Body<IntegrationFrame>*>* vessels_;
  std::map<RenderingFrameType, std::string> const* const frame_type_names_ =
      new std::map<RenderingFrameType, std::string>{
        {RenderingFrameType::kSurface, "Surface"},
        {RenderingFrameType::kCentric, "Centric"},
        {RenderingFrameType::kBarycentric, "Barycentre"},
      };
};

}  // namespace ksp
}  // namespace principia
