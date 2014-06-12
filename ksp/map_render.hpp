#pragma once

namespace principia {
namespace ksp {

public ref class MapRenderer : UnityEngine::MonoBehaviour {
 public:
  static MapRenderer CreateAndAttach(void (*draw)());
 private:
  void (*draw)();
};

}  // namespace ksp
}  // namespace principia