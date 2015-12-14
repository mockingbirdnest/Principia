using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal abstract class WindowRenderer {
  public interface ManagerInterface {
    event Action render_windows;
  }

  public WindowRenderer(ManagerInterface manager) {
    manager_ = manager;
    manager_.render_windows += RenderWindow;
  }

  public void Deregister() {
    manager_.render_windows -= RenderWindow;
  }

  abstract protected void RenderWindow();

  private ManagerInterface manager_;
}

internal struct Controlled<T> where T : WindowRenderer {
  public T renderer {
    get { return renderer_; }
    set {
      if (renderer_ != null) {
        renderer_.Deregister();
      }
      renderer_ = value;
    }
  }

  private T renderer_;
}

}  // namespace principia
}  // namespace ksp_plugin_adapter
