using System.Collections.Generic;

namespace principia {
namespace ksp_plugin_adapter {

internal class CollisionReporter : UnityEngine.MonoBehaviour {
  private void FixedUpdate() {
    collisions.Clear();
  }

  private void OnCollisionEnter(UnityEngine.Collision collision) {
    collisions.Add(collision);
  }

  private void OnCollisionStay(UnityEngine.Collision collision) {
    collisions.Add(collision);
  }

  public List<UnityEngine.Collision> collisions =
      new List<UnityEngine.Collision>();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
