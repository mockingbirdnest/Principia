
#include "ksp/map_renderer.hpp"

namespace principia {
namespace ksp {

void MapRenderer::OnPreCull() {
  draw_();
}

MapRenderer^ CreateMapRenderer(void (*draw)()) {
  MapRenderer^ renderer;// =
      //MapView::MapCamera->gameObject->GetComponent<MapRenderer^>();
  if (renderer != nullptr) {
    UnityEngine::Object::Destroy(renderer);
  }
  PlanetariumCamera^ x = MapView::MapCamera;
  UnityEngine::MonoBehaviour^ y = static_cast<UnityEngine::MonoBehaviour^>(x);
  UnityEngine::Component^ x;// = MapView::MapCamera;
  x->gameObject;
  //renderer = MapView::MapCamera->gameObject->AddComponent<MapRenderer^>()
  renderer->draw_ = draw;
  return renderer;
}

}  // namespace ksp
}  // namespace principia
