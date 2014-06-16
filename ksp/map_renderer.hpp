#pragma once

namespace principia {
namespace ksp {

public ref class MapRenderer : public UnityEngine::MonoBehaviour {
 public:
  void (*draw_)();
 private:
  // Unity event functions, accessed through reflection.
  // See http://docs.unity3d.com/Manual/ExecutionOrder.html.
  // Rendering (called before camera culling).
  void OnPreCull();
};

// Returns a MapRenderer which will call |*draw| before the map camera culls.
MapRenderer^ CreateMapRenderer(void (*draw)());

}  // namespace ksp
}  // namespace principia