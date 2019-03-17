using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

// TODO(egg): eventually |WindowRenderer| should own a rectangle and this should
// not be static.
internal static class WindowUtilities {
  public static void EnsureOnScreen(ref UnityEngine.Rect window_rectangle) {
    const float min_width_on_screen = 50;
    const float min_height_on_screen = 50;
    window_rectangle.x =
        UnityEngine.Mathf.Clamp(
            window_rectangle.x,
            -window_rectangle.width + min_width_on_screen,
            UnityEngine.Screen.width - min_width_on_screen);
    window_rectangle.y =
        UnityEngine.Mathf.Clamp(
            window_rectangle.y,
            -window_rectangle.height + min_height_on_screen,
            UnityEngine.Screen.height - min_height_on_screen);
  }

  public static void InputLock(this UnityEngine.Rect window_rectangle,
                               object window_owner) {
    string name = window_owner.GetType().ToString() + ":lock:" +
                  window_owner.GetHashCode();
    if (window_rectangle.ContainsMouse()) {
      InputLockManager.SetControlLock(PrincipiaLock, name);
    } else {
      InputLockManager.RemoveControlLock(name);
    }
  }

  public static void ClearLock(object window_owner) {
    string name = window_owner.GetType().ToString() + ":lock:" +
                  window_owner.GetHashCode();
    InputLockManager.RemoveControlLock(name);
  }

  private static bool ContainsMouse(this UnityEngine.Rect window_rectangle) {
    UnityEngine.Vector3 mouse = UnityEngine.Input.mousePosition;
    mouse.y = UnityEngine.Screen.height - mouse.y;
    return window_rectangle.Contains(mouse);
  }

  private static readonly ControlTypes PrincipiaLock =
      ControlTypes.ALLBUTCAMERAS &
      ~ControlTypes.ALL_SHIP_CONTROLS;
};

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
    WindowUtilities.ClearLock(this);
  }

  public void Dispose() {
    manager_.render_windows -= RenderWindow;
    WindowUtilities.ClearLock(this);
    GC.SuppressFinalize(this);
  }

  abstract protected void RenderWindow();

  private ManagerInterface manager_;
}

internal struct Controlled<T> where T : class, IDisposable {
  public T get() {
    return all_;
  }

  public void reset(T value = null) {
    if (all_ != null) {
      all_.Dispose();
    }
    all_ = value;
  }

  private T all_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
