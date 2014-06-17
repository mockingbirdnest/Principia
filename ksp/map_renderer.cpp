
#include "ksp/map_renderer.hpp"

namespace principia {
namespace ksp {

void MapRenderer::Start() {
  if (instance_ != nullptr) {
    Destroy(instance_->gameObject);
  }
  instance_ = this;

  if (HighLogic::LoadedScene == GameScenes::FLIGHT ||
      HighLogic::LoadedScene == GameScenes::TRACKSTATION) {
    gameObject->layer = 10;
    mesh_ = gameObject->AddComponent<UnityEngine::MeshFilter^>()->mesh;
    UnityEngine::MeshRenderer^ =
        gameObject->AddComponent(UnityEngine::MeshRenderer^);
  }
}

void MapRenderer::Update() {
}

}  // namespace ksp
}  // namespace principia
