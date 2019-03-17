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
    lock_name_ = GetType().ToString() + ":lock:" + GetHashCode();
  }

  ~WindowRenderer() {
    manager_.render_windows -= RenderWindow;
    ClearLock();
  }

  public void Dispose() {
    manager_.render_windows -= RenderWindow;
    ClearLock();
    GC.SuppressFinalize(this);
  }

  protected void EnsureOnScreen() {
    rectangle_.x = UnityEngine.Mathf.Clamp(
                       rectangle_.x,
                       -rectangle_.width + min_width_on_screen_,
                       UnityEngine.Screen.width - min_width_on_screen_);
    rectangle_.y = UnityEngine.Mathf.Clamp(
                       rectangle_.y,
                       -rectangle_.height + min_height_on_screen_,
                       UnityEngine.Screen.height - min_height_on_screen_);
  }

  protected void ClearLock() {
    InputLockManager.RemoveControlLock(lock_name_);
  }

  protected void InputLock() {
    UnityEngine.Vector3 mouse = UnityEngine.Input.mousePosition;
    mouse.y = UnityEngine.Screen.height - mouse.y;
    if (rectangle_.Contains(mouse)) {
      InputLockManager.SetControlLock(PrincipiaLock, lock_name_);
    } else {
      InputLockManager.RemoveControlLock(lock_name_);
    }
  }

  protected void Shrink() {
    rectangle_.height = 0.0f;
    rectangle_.width = 0.0f;
  }

  protected void Window(UnityEngine.GUI.WindowFunction func,
                        String text) {
    rectangle_ = UnityEngine.GUILayout.Window(
                     id         : this.GetHashCode(),
                     screenRect : rectangle_,
                     func       : func,
                     text       : text,
                     options    : UnityEngine.GUILayout.MinWidth(min_width_));
  }

  abstract public void RenderWindow();

  private static readonly ControlTypes PrincipiaLock =
      ControlTypes.ALLBUTCAMERAS &
      ~ControlTypes.ALL_SHIP_CONTROLS;

  private static readonly float min_height_on_screen_ = 50;
  private static readonly float min_width_on_screen_ = 50;

  private static readonly float min_width_ = 500;

  private ManagerInterface manager_;
  private String lock_name_;
  private UnityEngine.Rect rectangle_ =
      new UnityEngine.Rect(x      : (UnityEngine.Screen.width - min_width_) / 2,
                           y      : UnityEngine.Screen.height / 3,
                           width  : min_width_,
                           height : 0);
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
