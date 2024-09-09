using System;

namespace principia {
namespace ksp_plugin_adapter {

// A helper class for scaling the UI.  Unlike the WindowRenderer below, it can
// be used for elements that are not windows.
internal class ScalingRenderer {
  protected ScalingRenderer() {
    // All dimensions are expressed in "units".  By default a unit is 25 pixels
    // but it can scale up or down based on the KSP UI scale.
    scale_ = GameSettings.UI_SCALE * GameSettings.UI_SCALE_APPS;
    unit_ = 25 * scale_;
  }

  protected UnityEngine.GUILayoutOption GUILayoutHeight(int units) {
    return UnityEngine.GUILayout.Height(unit_ * units);
  }

  protected UnityEngine.GUILayoutOption GUILayoutMinWidth(int units) {
    return UnityEngine.GUILayout.MinWidth(Width(units));
  }

  protected UnityEngine.GUILayoutOption GUILayoutWidth(int units) {
    return UnityEngine.GUILayout.Width(Width(units));
  }

  protected UnityEngine.GUILayoutOption GUILayoutWidth(float units) {
    return UnityEngine.GUILayout.Width(Width(units));
  }

  protected float Height(float units) {
    return unit_ * units;
  }

  protected float Width(float units) {
    return unit_ * units;
  }

  protected UnityEngine.GUISkin MakeSkin(UnityEngine.GUISkin template) {
    UnityEngine.GUISkin skin = UnityEngine.Object.Instantiate(template);

    // Creating a dynamic font as is done below results in Unity producing
    // incorrect character bounds and everything looks ugly.  They even
    // "document" it in their source code, see
    // https://github.com/Unity-Technologies/UnityCsReference/blob/57f723ec72ca50427e5d17cad0ec123be2372f67/Modules/GraphViewEditor/Views/GraphView.cs#L262.
    // So here I am, sizing a pangram to get an idea of the shape of things and
    // nudging pixels by hand.  It's the 90's, go for it!
    var pangram = new UnityEngine.GUIContent(
        "Portez ce vieux whisky au juge blond qui fume.");
    float button_height = skin.button.CalcHeight(pangram, width : 1000);
    float label_height = skin.label.CalcHeight(pangram, width : 1000);
    float text_area_height = skin.textArea.CalcHeight(pangram, width : 1000);
    float text_field_height = skin.textField.CalcHeight(pangram, width : 1000);
    float toggle_height = skin.toggle.CalcHeight(pangram, width : 1000);

    skin.font = UnityEngine.Font.CreateDynamicFontFromOSFont(
        skin.font.fontNames,
        (int)(skin.font.fontSize * scale_));

    skin.button.alignment = UnityEngine.TextAnchor.MiddleCenter;
    skin.button.contentOffset =
        new UnityEngine.Vector2(0, -button_height * scale_ / 10);
    skin.button.fixedHeight = button_height * scale_;
    skin.horizontalSlider.fixedHeight = 21 * scale_;
    skin.horizontalSliderThumb.fixedHeight = 21 * scale_;
    skin.horizontalSliderThumb.fixedWidth = 12 * scale_;
    skin.label.alignment = UnityEngine.TextAnchor.MiddleLeft;
    skin.label.contentOffset =
        new UnityEngine.Vector2(0, -label_height * scale_ / 20);
    skin.label.fixedHeight = label_height * scale_;
    skin.textArea.alignment = UnityEngine.TextAnchor.MiddleLeft;
    skin.textArea.contentOffset =
        new UnityEngine.Vector2(0, -text_area_height * scale_ / 20);
    skin.textArea.fixedHeight = text_area_height * scale_;
    skin.textField.alignment = UnityEngine.TextAnchor.MiddleLeft;
    skin.textField.contentOffset =
        new UnityEngine.Vector2(0, -text_area_height * scale_ / 20);
    skin.textField.fixedHeight = text_field_height * scale_;
    skin.toggle.fixedHeight = toggle_height * scale_;
    skin.toggle.contentOffset =
        new UnityEngine.Vector2(0, -toggle_height * (scale_ - 1) / 3);
    skin.toggle.alignment = UnityEngine.TextAnchor.UpperLeft;
    skin.toggle.margin = new UnityEngine.RectOffset(
        (int)(skin.toggle.margin.left * scale_),
        skin.toggle.margin.right,
        (int)(skin.toggle.margin.top * 1.7 * scale_),
        skin.toggle.margin.bottom);
    return skin;
  }

  private readonly float scale_;
  private readonly float unit_;
}

// A class that gather all the mechanisms for rendering entire windows.  It
// deals with skins, input locking, sizing, placement, and hiding.  It also
// persists the position and size of the window.
internal abstract class BaseWindowRenderer : ScalingRenderer, IConfigNode {
  protected abstract string Title { get; }
  protected abstract void RenderWindow(int window_id);

  protected BaseWindowRenderer(UnityEngine.GUILayoutOption[] options) {
    UnityEngine.GUILayoutOption[] default_options = {GUILayoutMinWidth(20)};
    options_ = options.Length == 0 ? default_options : options;
    lock_name_ = $"{GetType()}:lock:{GetHashCode()}";
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
    if (skin_ == null) {
      UnityEngine.GUI.skin = null;
      skin_ = MakeSkin(UnityEngine.GUI.skin);
    }
    UnityEngine.GUI.skin = skin_;
    if (show_) {
      // NOTE: Calling `Shrink` here (in a Layout, before drawing the window)
      // satisfies the conditions noted in its doc comment and is safe.
      Shrink();
      rectangle_ = UnityEngine.GUILayout.Window(
          id         : this.GetHashCode(),
          screenRect : rectangle_,
          func       : RenderWindowAndRecordTooltip,
          text       : Title,
          options    : options_);

      // The first time a window is shown, we have a moral duty to place it at
      // the centre of the screen.  This is tricky because we don't know its
      // width.  So we go through this method until the width has been
      // determined, and we update the window based on the width.  Note that
      // `must_centre_` has to be persisted because there can be multiple
      // scene changes before a window is shown for the first time.  Also note
      // the conversion to float to avoid losing decimals and to avoid double
      // rounding.
      if (must_centre_ && rectangle_.width > 0) {
        rectangle_ = UnityEngine.GUILayout.Window(
            id         : GetHashCode(),
            screenRect : new UnityEngine.Rect(
                x      : (UnityEngine.Screen.width - rectangle_.width) / 2,
                y      : (float)UnityEngine.Screen.height / 3f,
                width  : 0,
                height : 0),
            func       : RenderWindowAndRecordTooltip,
            text       : Title,
            options    : options_);
        must_centre_ = false;
      }
      ShowTooltip();
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
        -rectangle_.width + min_width_on_screen,
        UnityEngine.Screen.width - min_width_on_screen);
    rectangle_.y = UnityEngine.Mathf.Clamp(
        rectangle_.y,
        -rectangle_.height + min_height_on_screen,
        UnityEngine.Screen.height - min_height_on_screen);
  }

  private void RenderWindowAndRecordTooltip(int window_id) {
    RenderWindow(window_id);
    if (UnityEngine.Event.current.type == UnityEngine.EventType.Repaint &&
        tooltip_ != UnityEngine.GUI.tooltip) {
      if (tooltip_ == "") {
        tooltip_begin_ = DateTime.UtcNow;
      }
      tooltip_ = UnityEngine.GUI.tooltip;
      var height = Style.Multiline(UnityEngine.GUI.skin.textArea).CalcHeight(
          new UnityEngine.GUIContent(tooltip_), Width(8));
      tooltip_rectangle_ = new UnityEngine.Rect(
          UnityEngine.Input.mousePosition.x + Width(1) / 2,
          UnityEngine.Screen.height -
          (UnityEngine.Input.mousePosition.y - Width(1) / 2),
          Width(8), height);
    }
  }

  private void ShowTooltip() {
    if (tooltip_ != "" &&
        (DateTime.UtcNow - tooltip_begin_).TotalMilliseconds > 500) {
      var tooltip_style = Style.Multiline(UnityEngine.GUI.skin.textArea);
      tooltip_style.font = UnityEngine.GUI.skin.font;
      UnityEngine.GUI.Window(
          tooltip_.GetHashCode(),
          tooltip_rectangle_,
          (int window_id) => {},
          tooltip_,
          tooltip_style);
      UnityEngine.GUI.BringWindowToFront(tooltip_.GetHashCode());
    }
  }

  public void ScheduleShrink() {
    shrink_scheduled_ = true;
  }

  // NOTE(al2me6): Empirical observation led to the following conclusions:
  // In each frame, `OnGUI` (from whence this method is ultimately called)
  // is called at least twice, with the following order of event types:
  // 1. Layout
  // 2. For each interaction during that frame (e.g., mouse event):
  //   a. Interaction event
  //   b. Layout
  // 3. Repaint
  // cf. <https://docs.unity3d.com/ScriptReference/Event.html>.
  // Furthermore, this function is not safe to call between the call to
  // `GUILayout.Window` in the last Layout of the frame and the subsequent
  // Repaint; doing so may cause the window to become blank for one frame.
  private void Shrink() {
    if (shrink_scheduled_
        && UnityEngine.Event.current.type == UnityEngine.EventType.Layout) {
      rectangle_.height = 0.0f;
      rectangle_.width = 0.0f;
      shrink_scheduled_ = false;
    }
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

  public virtual void Load(ConfigNode node) {
    string must_centre_value = node.GetAtMostOneValue("must_centre");
    if (must_centre_value != null) {
      must_centre_ = Convert.ToBoolean(must_centre_value);
    }
    string show_value = node.GetAtMostOneValue("show");
    if (show_value != null) {
      show_ = Convert.ToBoolean(show_value);
    }
    string x_value = node.GetAtMostOneValue("x");
    if (x_value != null) {
      rectangle_.x = Convert.ToSingle(x_value);
    }
    string y_value = node.GetAtMostOneValue("y");
    if (y_value != null) {
      rectangle_.y = Convert.ToSingle(y_value);
    }
  }

  public virtual void Save(ConfigNode node) {
    node.SetValue("must_centre", must_centre_, createIfNotFound : true);
    node.SetValue("show", show_, createIfNotFound : true);
    node.SetValue("x", rectangle_.x, createIfNotFound : true);
    node.SetValue("y", rectangle_.y, createIfNotFound : true);
  }

  // Lock everything but camera mode switching and ship controls.
  private static readonly ControlTypes PrincipiaLock =
      (ControlTypes.ALLBUTCAMERAS | ControlTypes.CAMERACONTROLS) &
          ~ControlTypes.ALL_SHIP_CONTROLS;

  private const float min_height_on_screen = 50;
  private const float min_width_on_screen = 50;

  private readonly UnityEngine.GUILayoutOption[] options_;
  private readonly string lock_name_;
  private UnityEngine.GUISkin skin_;
  private bool must_centre_ = true;
  private bool shrink_scheduled_ = false;
  private bool show_ = false;
  protected UnityEngine.Rect rectangle_;
  private DateTime tooltip_begin_;
  private string tooltip_ = "";
  private UnityEngine.Rect tooltip_rectangle_;
}

// The supervisor of a window decides when to clear input locks, when to render
// the window and when to delete it.  A supervised window has a “Hide” button,
// in addition to the external toggle handled by the supervisor.
internal abstract class SupervisedWindowRenderer : BaseWindowRenderer {
  public interface ISupervisor {
    event Action LockClearing;
    event Action WindowsDisposal;
    event Action WindowsRendering;
  }

  protected SupervisedWindowRenderer(ISupervisor supervisor,
                                     params UnityEngine.GUILayoutOption[]
                                         options) : base(options) {
    supervisor_ = supervisor;
    supervisor_.LockClearing += ClearLock;
    supervisor_.WindowsDisposal += DisposeWindow;
    supervisor_.WindowsRendering += RenderWindow;
  }

  protected abstract void RenderWindowContents(int window_id);

  protected sealed override void RenderWindow(int window_id) {
    if (UnityEngine.GUI.Button(new UnityEngine.Rect(
            x: rectangle_.width - Width(1),
            y: 0,
            width: Width(1),
            height: Width(1)),
            "×")) {
      Hide();
    }
    RenderWindowContents(window_id);
  }

  public void DisposeWindow() {
    supervisor_.LockClearing -= ClearLock;
    supervisor_.WindowsDisposal -= DisposeWindow;
    supervisor_.WindowsRendering -= RenderWindow;
  }

  private readonly ISupervisor supervisor_;
}

// A supervised window that displays properties of a vessel, given by a
// delegate.  The delegate must return a nonnull value if and only if there is
// a currently selected/active vessel *and* it is known to the plugin.
// Subclasses must access the vessel through the predicted_vessel property and
// should not cache it.
internal abstract class
    VesselSupervisedWindowRenderer : SupervisedWindowRenderer {
  public delegate Vessel PredictedVessel();

  protected VesselSupervisedWindowRenderer(ISupervisor supervisor,
                                           PredictedVessel predicted_vessel,
                                           params UnityEngine.GUILayoutOption[]
                                               options) : base(
      supervisor,
      options) {
    predicted_vessel_ = predicted_vessel;
  }

  // A helper for implementing the RenderButton() method of the subclasses.
  protected void RenderButton(string text,
                              params UnityEngine.GUILayoutOption[] options) {
    if (UnityEngine.GUILayout.Button(text, options)) {
      Toggle();
    }
    // Override the state of the toggle if there is no predicted vessel.
    if (predicted_vessel == null) {
      Hide();
    }
  }

  protected Vessel predicted_vessel => predicted_vessel_();

  private readonly PredictedVessel predicted_vessel_;
}

// A window without a supervisor is effectively modal.
internal abstract class UnsupervisedWindowRenderer : BaseWindowRenderer {
  protected UnsupervisedWindowRenderer(
      params UnityEngine.GUILayoutOption[] options) : base(options) {}
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
