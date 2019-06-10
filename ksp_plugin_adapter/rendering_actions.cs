using System;

namespace principia {
namespace ksp_plugin_adapter {

class RenderingActions : UnityEngine.MonoBehaviour {
  internal Action post_render {
    private get;
    set;
  }

  private void OnPostRender() {
    if (post_render != null) {
      post_render();
    }
  }

  internal Action pre_cull {
    private get;
    set;
  }

  private void OnPreCull() {
    if (pre_cull != null) {
      pre_cull();
    }
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
