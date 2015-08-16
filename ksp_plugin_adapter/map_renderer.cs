using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

class MapRenderer : UnityEngine.MonoBehaviour {
  internal delegate void Task();

  internal Task render_on_pre_cull {
    private get;
    set;
  }

  private void OnPreCull() {
    if (render_on_pre_cull != null) {
      render_on_pre_cull();
    }
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
