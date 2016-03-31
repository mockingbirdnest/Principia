using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

class MapRenderer : UnityEngine.MonoBehaviour {
  internal Action post_render {
    private get;
    set;
  }

  private void OnPostRender() {
    if (post_render != null) {
      post_render();
    }
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
