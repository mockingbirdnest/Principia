using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal abstract class BaseWindowRenderer : IConfigNode {
  protected BaseWindowRenderer(UnityEngine.GUILayoutOption[] options) {
    options_ = options.Length == 0 ? default_options_ : options;
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

  public void RenderWindow() {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    if (show_) {
      rectangle_ = UnityEngine.GUILayout.Window(
                       id         : this.GetHashCode(),
                       screenRect : rectangle_,
                       func       : RenderWindow,
                       text       : Title,
                       options    : options_);

      // The first time a window is shown, we have a moral duty to place it at
      // the centre of the screen.  This is tricky because we don't know its
      // width.  So we go through this method until the width has been
      // determined, and we update the window based on the width.  Note that
      // |must_centre_| has to be persisted because there can be multiple
      // scene changes because a window is shown for the first time.
      if (must_centre_ && rectangle_.width > 0) {
        rectangle_ =
            UnityEngine.GUILayout.Window(
                id         : this.GetHashCode(),
                screenRect : new UnityEngine.Rect(
                    x      : (UnityEngine.Screen.width - rectangle_.width) / 2,
                    y      : UnityEngine.Screen.height / 3,
                    width  : 0,
                    height : 0),
                func       : RenderWindow,
                text       : Title,
                options    : options_);
        must_centre_ = false;
      }
      EnsureOnScreen();
      InputLock();
    } else {
      ClearLock();
    }
    UnityEngine.GUI.skin = old_skin;
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

  public void Shrink() {
    rectangle_.height = 0.0f;
    rectangle_.width = 0.0f;
  }

  // Visibility.

  public void Hide() {
    show_ = false;
  }

  public void Show() {
    show_ = true;
  }

  public bool Shown() {
    return show_;
  }

  public void Toggle() {
    show_ = !show_;
  }

  // Persistence.

  public void Load(ConfigNode node) {
    String must_centre_value = node.GetAtMostOneValue("must_centre");
    if (must_centre_value != null) {
      must_centre_ = Convert.ToBoolean(must_centre_value);
    }
    String show_value = node.GetAtMostOneValue("show");
    if (show_value != null) {
      show_ = Convert.ToBoolean(show_value);
    }
    String x_value = node.GetAtMostOneValue("x");
    if (x_value != null) {
      rectangle_.x = Convert.ToSingle(x_value);
    }
    String y_value = node.GetAtMostOneValue("y");
    if (y_value != null) {
      rectangle_.y = Convert.ToSingle(y_value);
    }
  }

  public void Save(ConfigNode node) {
    node.SetValue("must_centre", must_centre_, createIfNotFound : true);
    node.SetValue("show", show_, createIfNotFound : true);
    node.SetValue("x", rectangle_.x, createIfNotFound : true);
    node.SetValue("y", rectangle_.y, createIfNotFound : true);
  }

  abstract protected String Title { get; }
  abstract protected void RenderWindow(int window_id);

  private static readonly ControlTypes PrincipiaLock =
      ControlTypes.ALLBUTCAMERAS &
      ~ControlTypes.ALL_SHIP_CONTROLS;

  private const float min_height_on_screen_ = 50;
  private const float min_width_on_screen_ = 50;
  private static readonly UnityEngine.GUILayoutOption[] default_options_ =
        {UnityEngine.GUILayout.MinWidth(500)};

  private readonly UnityEngine.GUILayoutOption[] options_;
  private readonly String lock_name_;
  private bool must_centre_ = true;
  private bool show_ = false;
  private UnityEngine.Rect rectangle_ = UnityEngine.Rect.zero;
}

internal abstract class SupervisedWindowRenderer : BaseWindowRenderer {
  public interface ISupervisor {
    event Action clear_locks;
    event Action dispose_windows;
    event Action render_windows;
  }

  public SupervisedWindowRenderer(ISupervisor supervisor,
                                  params UnityEngine.GUILayoutOption[] options)
        : base(options) {
    supervisor_ = supervisor;
    supervisor_.clear_locks += ClearLock;
    supervisor_.dispose_windows += DisposeWindow;
    supervisor_.render_windows += RenderWindow;
  }

  public void DisposeWindow() {
    supervisor_.clear_locks -= ClearLock;
    supervisor_.dispose_windows -= DisposeWindow;
    supervisor_.render_windows -= RenderWindow;
  }

  private ISupervisor supervisor_;
}

internal abstract class UnsupervisedWindowRenderer : BaseWindowRenderer {
  public UnsupervisedWindowRenderer(
      params UnityEngine.GUILayoutOption[] options) : base(options) {}
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
