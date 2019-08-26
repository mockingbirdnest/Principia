using System;

namespace principia {
namespace ksp_plugin_adapter {

class RenderingActions : UnityEngine.MonoBehaviour {
  internal Action? post_render {
    private get;
    set;
  }

  private void OnPostRender() {
    post_render?.Invoke();
  }

  internal Action? pre_cull {
    private get;
    set;
  }

  private void OnPreCull() {
    pre_cull?.Invoke();
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
