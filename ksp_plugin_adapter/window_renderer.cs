using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal abstract class WindowRenderer : IDisposable {
  public interface ManagerInterface {
    event Action render_windows;
  }

  public WindowRenderer(ManagerInterface manager) {
    manager_ = manager;
    manager_.render_windows += RenderWindow;
  }

  ~WindowRenderer() {
    manager_.render_windows -= RenderWindow;
  }

  public void Dispose() {
    manager_.render_windows -= RenderWindow;
    GC.SuppressFinalize(this);
  }

  abstract protected void RenderWindow();

  private ManagerInterface manager_;
}

internal struct Controlled<T> where T : IDisposable {
  public T all {
    get {
      return all_;
    }
    set {
      if (all_ != null) {
        all_.Dispose();
      }
      all_ = value;
    }
  }

  private T all_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
