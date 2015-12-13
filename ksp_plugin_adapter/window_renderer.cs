using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

// A class for delegating calls for window rendering without leakage.
// The references are as follows, where => stands for strong references and ->
// for weak references.
//   ManagerInterface <=> Handle -> WindowRenderer
//                               <=
// This does not include any strong refererences to the |WindowRenderer|, so
// that when the last strong reference to the |WindowRenderer| disappears, it
// may be collected.  When that happens, its finalizer breaks the strong
// reference |ManagerInterface => Handle|, leaving the references
//   UIManagerInterface <= Handle -> null
// so that the |Handle| may in turn be collected, leaving no dangling weak
// reference.
// Note that |Handle| is private, preventing any other references from foiling
// the last part of the scheme.
// In a language with unique ownership we would just do
//   manager => function pointers <-> renderer
//           =======================>
// and our lives would be simpler.
internal abstract class WindowRenderer {
  public interface ManagerInterface {
    event Action render_windows;
  }

  public WindowRenderer(ManagerInterface manager) {
    handle_ = new Handle(manager, this);
  }

  ~WindowRenderer() {
    handle_.Deregister();
  }

  protected abstract void RenderWindow();

  private sealed class Handle {
    public Handle(ManagerInterface manager,
                  WindowRenderer window_renderer) {
      manager_ = manager;
      manager_.render_windows += RenderWindow;
      window_renderer_ = new WeakReference(window_renderer);
    }

    private void RenderWindow() {
      WindowRenderer renderer = (WindowRenderer)window_renderer_.Target;
      if (renderer != null) {
        renderer.RenderWindow();
      } else {
        Deregister();
      }
    }

    public void Deregister() {
      manager_.render_windows -= RenderWindow;
    }

    private ManagerInterface manager_;
    private WeakReference window_renderer_;
  }

  private Handle handle_;
}

}  // namespace principia
}  // namespace ksp_plugin_adapter
