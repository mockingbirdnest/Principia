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

internal abstract class BaseWindowRenderer : IConfigNode {
  protected BaseWindowRenderer() {
    lock_name_ = GetType().ToString() + ":lock:" + GetHashCode();
  }

  // Locking.

  protected void ClearLock() {
    InputLockManager.RemoveControlLock(lock_name_);
  }

  private void InputLock() {
    UnityEngine.Vector3 mouse = UnityEngine.Input.mousePosition;
    mouse.y = UnityEngine.Screen.height - mouse.y;
    if (rectangle_.Contains(mouse)) {
      InputLockManager.SetControlLock(PrincipiaLock, lock_name_);
    } else {
      InputLockManager.RemoveControlLock(lock_name_);
    }
  }

  // Rendering.

  protected void Window(UnityEngine.GUI.WindowFunction func,
                        String text) {
    rectangle_ = UnityEngine.GUILayout.Window(
                     id         : this.GetHashCode(),
                     screenRect : rectangle_,
                     func       : func,
                     text       : text,
                     options    : UnityEngine.GUILayout.MinWidth(min_width_));
    EnsureOnScreen();
    InputLock();
  }

  private void EnsureOnScreen() {
    rectangle_.x = UnityEngine.Mathf.Clamp(
                       rectangle_.x,
                       -rectangle_.width + min_width_on_screen_,
                       UnityEngine.Screen.width - min_width_on_screen_);
    rectangle_.y = UnityEngine.Mathf.Clamp(
                       rectangle_.y,
                       -rectangle_.height + min_height_on_screen_,
                       UnityEngine.Screen.height - min_height_on_screen_);
  }

  protected void Shrink() {
    rectangle_.height = 0.0f;
    rectangle_.width = 0.0f;
  }

  // Persistence.

  public void Load(ConfigNode node) {
    String x_value = node.GetAtMostOneValue("x");
    if (x_value != null) {
      rectangle_.x = System.Convert.ToSingle(x_value);
    }
    String y_value = node.GetAtMostOneValue("y");
    if (y_value != null) {
      rectangle_.y = System.Convert.ToSingle(y_value);
    }
  }

  public void Save(ConfigNode node) {
    node.SetValue("x", rectangle_.x, createIfNotFound : true);
    node.SetValue("y", rectangle_.y, createIfNotFound : true);
  }

  abstract protected void RenderWindow();

  private static readonly ControlTypes PrincipiaLock =
      ControlTypes.ALLBUTCAMERAS &
      ~ControlTypes.ALL_SHIP_CONTROLS;

  private static readonly float min_height_on_screen_ = 50;
  private static readonly float min_width_on_screen_ = 50;

  private static readonly float min_width_ = 500;

  private String lock_name_;
  private UnityEngine.Rect rectangle_ =
      new UnityEngine.Rect(x      : (UnityEngine.Screen.width - min_width_) / 2,
                           y      : UnityEngine.Screen.height / 3,
                           width  : min_width_,
                           height : 0);
}

internal abstract class SupervisedWindowRenderer : BaseWindowRenderer {
  public interface ISupervisor {
    event Action dispose_windows;
    event Action render_windows;
  }

  public SupervisedWindowRenderer(ISupervisor supervisor) : base() {
    supervisor_ = supervisor;
    supervisor_.dispose_windows += DisposeWindow;
    supervisor_.render_windows += RenderWindow;
  }

  public void DisposeWindow() {
    supervisor_.dispose_windows -= DisposeWindow;
    supervisor_.render_windows -= RenderWindow;
    ClearLock();
  }

  private ISupervisor supervisor_;
}

internal abstract class UnsupervisedWindowRenderer : BaseWindowRenderer {}

}  // namespace ksp_plugin_adapter
}  // namespace principia
