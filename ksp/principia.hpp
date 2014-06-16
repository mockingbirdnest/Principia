#pragma once

#include <map>
#include <vector>

#include "physics/body.hpp"
#include "physics/n_body_system.hpp"

namespace principia {
namespace ksp {

// The reference frame in which the integration is performed.
struct IntegrationFrame;

// The reference frame in which the orbits are rendered.
struct RenderingFrame;

// The type of frame in which the orbits are rendered.
// Examples given for Earth as the rendering reference body.
enum RenderingFrameType {
  kSurface,           // e.g., Earth surface.
  kCentric,           // e.g., Geocentric.
  kBarycentric,       // e.g., Sun-Earth barycentre.
};

std::vector<RenderingFrameType> const kRenderingFrameTypes = {kSurface,
                                                              kCentric,
                                                              kBarycentric};

// User-friendly names for the reference frame types.
std::map<RenderingFrameType, std::string> const kFrameTypeNames =
    std::map<RenderingFrameType, std::string>{
        {kSurface, "Surface"},
        {kCentric, "Centric"},
        {kBarycentric, "Barycentric rotating"}};

template <typename Frame>
using Displacement = Vector<Length, Frame>;
template <typename Frame>
using VelocityChange = Vector<Speed, Frame>;

[KSPAddon(KSPAddon::Startup::Flight, false)]
public ref class Principia : public UnityEngine::MonoBehaviour {
 private:
   ~Principia();
   !Principia();

  // Unity event functions, accessed through reflection.
  // See http://docs.unity3d.com/Manual/ExecutionOrder.html.

  // Called before the first frame update.
  void Start();
  // Unity physics step.
  void FixedUpdate();
  // Called once per frame.
  void Update();
  // GUI rendering.
  void OnGUI();
  // Called after all frame updates for the last frame of existence.
  void OnDestroy();

  void DrawGUI();
  void DrawMainWindow(int window_id);
  void DrawReferenceFrameWindow(int window_id);

  void SetUpSystem();

  void FreeAllOwnedPointers();

  UnityEngine::Rect main_window_position_;
  UnityEngine::Rect reference_frame_window_position_;

  UnityEngine::GUIStyle^ gui_style;

  bool simulating_;

  // We own these pointers, but we cannot use |std::unique_ptr| because that's
  // an unmanaged type.
  quantities::Time* Δt_                = new quantities::Time(60 * si::Second);
  quantities::Time* prediction_length_ = new quantities::Time(12 * si::Hour);
  int sampling_period_                 = 10;

  CelestialBody^ sun_;
  CelestialBody^ rendering_reference_body_;
  RenderingFrameType rendering_frame_type_;

  // We own these pointers, but we cannot use |std::unique_ptr| because that's
  // an unmanaged type.
  physics::NBodySystem<IntegrationFrame>* system_ = nullptr;
  std::map<std::string, physics::Body<IntegrationFrame>*>* celestials_ =
      new std::map<std::string, physics::Body<IntegrationFrame>*>();
  std::map<std::string, physics::Body<IntegrationFrame>*>* vessels_ =
      new std::map<std::string, physics::Body<IntegrationFrame>*>();
  integrators::SPRKIntegrator<Length, Speed>* integrator_ =
      new integrators::SPRKIntegrator<Length, Speed>();
};

}  // namespace ksp
}  // namespace principia
