#pragma once

namespace principia {
namespace ksp {

// A singleton that draws on the map view.
[KSPAddon(KSPAddon::Startup::EveryScene, false)]
public ref class MapRenderer : public UnityEngine::MonoBehaviour {
 public:
 private:
  void Start();
  void Update();
  static MapRenderer^ instance_ = nullptr;
  UnityEngine::Mesh^  mesh_     = nullptr;
};

}  // namespace ksp
}  // namespace principia