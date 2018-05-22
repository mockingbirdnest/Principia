#undef HAS_SURFACE
#undef HAS_BODY_CENTRED_ALIGNED_WITH_PARENT

using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

static class CelestialExtensions {
  public static bool is_leaf(this CelestialBody celestial) {
    return celestial.orbitingBodies.Count == 0;
  }

  public static bool is_root(this CelestialBody celestial) {
    return celestial.orbit == null;
  }
}

class ReferenceFrameSelector : WindowRenderer {
  public enum FrameType {
    BARYCENTRIC_ROTATING = 6001,
    BODY_CENTRED_NON_ROTATING = 6000,
    BODY_CENTRED_PARENT_DIRECTION = 6002,
    BODY_SURFACE = 6003,
  }

  public delegate void Callback(NavigationFrameParameters frame_parameters);

  public ReferenceFrameSelector(
      ManagerInterface manager,
      IntPtr plugin,
      Callback on_change,
      string name,
      bool is_plotting_frame_selector = true) : base(manager) {
    plugin_ = plugin;
    on_change_ = on_change;
    name_ = name;
    is_plotting_frame_selector_ = is_plotting_frame_selector;
    frame_type = FrameType.BODY_CENTRED_NON_ROTATING;
    expanded_ = new Dictionary<CelestialBody, bool>();
    foreach (CelestialBody celestial in FlightGlobals.Bodies) {
      if (!celestial.is_leaf() && !celestial.is_root()) {
        expanded_.Add(celestial, false);
      }
    }
    selected_celestial =
        FlightGlobals.currentMainBody ?? FlightGlobals.GetHomeBody();
    for (CelestialBody celestial = selected_celestial;
         celestial.orbit != null;
         celestial = celestial.orbit.referenceBody) {
      if (!celestial.is_leaf()) {
        expanded_[celestial] = true;
      }
    }
    on_change_(FrameParameters());
    window_rectangle_.x = UnityEngine.Screen.width / 2;
    window_rectangle_.y = UnityEngine.Screen.height / 3;

    InitializeToolbarToggle();
  }

  public FrameType frame_type { get; private set; }
  public CelestialBody selected_celestial { get; private set; }
  public Vessel target_override { get; set; }

  // Sets the |frame_type| to |type| unless this would be invalid for the
  // |selected_celestial|, in which case |frame_type| is set to
  // |BODY_CENTRED_NON_ROTATING|.
  public void SetFrameType(FrameType type) {
    if (selected_celestial.is_root() &&
        (type == FrameType.BARYCENTRIC_ROTATING ||
         type == FrameType.BODY_CENTRED_PARENT_DIRECTION)) {
      frame_type = FrameType.BODY_CENTRED_NON_ROTATING;
    } else {
      frame_type = type;
    }
    on_change_(FrameParameters());
  }

  // Sets the |frame_type| to |BODY_SURFACE| and sets |selected_celestial| to
  // the given |celestial|.
  public void SetToSurfaceFrameOf(CelestialBody celestial) {
    frame_type = FrameType.BODY_SURFACE;
    selected_celestial = celestial;
    on_change_(FrameParameters());
  }

  public static String Name(FrameType type,
                            CelestialBody selected,
                            Vessel target_override) {
   if (target_override) {
     return "Target Local Vert./Horiz. at " + selected.NameWithArticle();
   }
   switch (type) {
     case FrameType.BODY_CENTRED_NON_ROTATING:
       return selected.name + "-Centred Inertial";
     case FrameType.BARYCENTRIC_ROTATING:
        if (selected.is_root()) {
          throw Log.Fatal("Naming barycentric rotating frame of root body");
        } else {
          return selected.referenceBody.name + "-" + selected.name +
                 " Barycentric";
        }
     case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        if (selected.is_root()) {
          throw Log.Fatal(
              "Naming parent-direction rotating frame of root body");
        } else {
          return selected.name + "-Centred " + selected.referenceBody.name +
                 "-Aligned";
        }
     case FrameType.BODY_SURFACE:
       return selected.name + "-Centred " + selected.name + "-Fixed";
     default:
       throw Log.Fatal("Unexpected type " + type.ToString());
   }
  }

  public static String ShortName(FrameType type,
                                 CelestialBody selected,
                                 Vessel target_override) {
    if (target_override) {
      return "Tgt LVLH@" + selected.name[0];
    }
    switch (type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
        return selected.name[0] + "CI";
      case FrameType.BARYCENTRIC_ROTATING:
        if (selected.is_root()) {
          throw Log.Fatal("Naming barycentric rotating frame of root body");
        } else {
          return selected.referenceBody.name[0] + (selected.name[0] + "B");
        }
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        if (selected.is_root()) {
          throw Log.Fatal(
              "Naming parent-direction rotating frame of root body");
        } else {
          return selected.name[0] + "C" + selected.referenceBody.name[0] +
                 "A";
        }
      case FrameType.BODY_SURFACE:
        return selected.name[0] + "C" + selected.name[0] + "F";
      default:
        throw Log.Fatal("Unexpected type " + type.ToString());
    }
  }

  public static String Description(FrameType type,
                                   CelestialBody selected,
                                   Vessel target_override) {
    if (target_override) {
      return "Reference frame fixing the target vessel (" +
             target_override.vesselName + "), the plane of its orbit around " +
             selected.NameWithArticle() + ", and the line between them";
    }
    switch (type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
        return "Non-rotating reference frame fixing the centre of " +
               selected.NameWithArticle();
      case FrameType.BARYCENTRIC_ROTATING:
        if (selected.is_root()) {
          throw Log.Fatal("Describing barycentric rotating frame of root body");
        } else {
          return "Reference frame fixing the barycentre of " +
                 selected.NameWithArticle() + " and " +
                 selected.referenceBody.NameWithArticle() +
                 ", the plane in which they move about the barycentre, and" +
                 " the line between them";
        }
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        if (selected.is_root()) {
          throw Log.Fatal(
              "Describing parent-direction rotating frame of root body");
        } else {
          return "Reference frame fixing the centre of " +
                 selected.NameWithArticle() +
                 ", the plane of its orbit around " +
                 selected.referenceBody.NameWithArticle() +
                 ", and the line between them";
        }
      case FrameType.BODY_SURFACE:
        return "Reference frame fixing the surface of " +
               selected.NameWithArticle();
      default:
        throw Log.Fatal("Unexpected type " + type.ToString());
    }
  }

  public String Name() {
    return Name(frame_type, selected_celestial, target_override);
  }

  public String ShortName() {
    return ShortName(frame_type, selected_celestial, target_override);
  }

  public String Description() {
    return Description(frame_type, selected_celestial, target_override);
  }

  public CelestialBody[] FixedBodies() {
    if (target_override) {
      return new CelestialBody[]{};
    }
    switch (frame_type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
      case FrameType.BODY_SURFACE:
        return new CelestialBody[]{selected_celestial};
      case FrameType.BARYCENTRIC_ROTATING:
        return new CelestialBody[]{};
      default:
        throw Log.Fatal("Unexpected frame_type " + frame_type.ToString());
    }
  }

  public NavigationFrameParameters FrameParameters() {
    switch (frame_type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
      case FrameType.BODY_SURFACE:
        return new NavigationFrameParameters{
            extension = (int)frame_type,
            centre_index = selected_celestial.flightGlobalsIndex};
      case FrameType.BARYCENTRIC_ROTATING:
        return new NavigationFrameParameters{
            extension = (int)frame_type,
            primary_index =
                selected_celestial.referenceBody.flightGlobalsIndex,
            secondary_index = selected_celestial.flightGlobalsIndex};
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        // We put the primary body as secondary, because the one we want fixed
        // is the secondary body (which means it has to be the primary in the
        // terminology of |BodyCentredBodyDirection|).
        return new NavigationFrameParameters{
            extension = (int)frame_type,
            primary_index = selected_celestial.flightGlobalsIndex,
            secondary_index =
                selected_celestial.referenceBody.flightGlobalsIndex};
      default:
        throw Log.Fatal("Unexpected frame_type " + frame_type.ToString());
    }
  }

  public void Reset(NavigationFrameParameters parameters) {
    frame_type = (FrameType)parameters.extension;
    switch (frame_type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
      case FrameType.BODY_SURFACE:
        selected_celestial = FlightGlobals.Bodies[parameters.centre_index];
        break;
      case FrameType.BARYCENTRIC_ROTATING:
        selected_celestial = FlightGlobals.Bodies[parameters.secondary_index];
        break;
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        selected_celestial = FlightGlobals.Bodies[parameters.primary_index];
        break;
    }
  }

  public void Hide() {
    ShowSelector = false;
  }

  // Toolbar toggle
  private static ToolbarToggle show_selector_toggle_;

  bool ShowSelector
  {
    get => is_plotting_frame_selector_ ? show_selector_toggle_.Value : show_selector_;
    set {
      if (is_plotting_frame_selector_) {
        show_selector_toggle_.Value = value;
      } else {
        show_selector_ = value;
      }
    }
  }

  private void InitializeToolbarToggle() {
    if (show_selector_toggle_ == null) {
      show_selector_toggle_ = new ToolbarToggle(
          id: "PlottingFrame",
          tooltip: "Plotting Frame Selection",
          visibility: new[] { GameScenes.FLIGHT, GameScenes.TRACKSTATION });
    }
  }

  public void RenderButton() {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    if (UnityEngine.GUILayout.Button(name_ + " selection (" + Name() +
                                     ")...")) {
      ShowSelector = !ShowSelector;
    }
    UnityEngine.GUI.skin = old_skin;
  }

  protected override void RenderWindow() {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    if (ShowSelector) {
      window_rectangle_ = UnityEngine.GUILayout.Window(
                              id         : this.GetHashCode(),
                              screenRect : window_rectangle_,
                              func       : RenderSelector,
                              text       : name_ + " selection (" + Name() +
                                           ")");
      WindowUtilities.EnsureOnScreen(ref window_rectangle_);
      window_rectangle_.InputLock(this);
    } else {
      WindowUtilities.ClearLock(this);
    }
    UnityEngine.GUI.skin = old_skin;
  }

  protected override void SetToolbarButtonState(bool is_enabled) {
    if (show_selector_toggle_ != null) {
      show_selector_toggle_.IsEnabled = is_enabled;
    }
  }

  private void RenderSelector(int window_id) {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;

    using (new HorizontalLayout()) {
      // Left-hand side: tree view for celestial selection.
      using (new VerticalLayout(UnityEngine.GUILayout.Width(200))) {
        RenderSubtree(celestial : Planetarium.fetch.Sun, depth : 0);
      }

      // Right-hand side: toggles for reference frame type selection.
      using (new VerticalLayout()) {
        if (target_override) {
          UnityEngine.GUILayout.Label(
              "Using target-centred frame selected on navball speed display",
              UnityEngine.GUILayout.Width(150));
          UnityEngine.GUILayout.Label(
              Description(frame_type, selected_celestial, target_override),
              UnityEngine.GUILayout.Width(150));
        } else {
          TypeSelector(FrameType.BODY_SURFACE);
          TypeSelector(FrameType.BODY_CENTRED_NON_ROTATING);
          if (!selected_celestial.is_root()) {
            CelestialBody parent = selected_celestial.orbit.referenceBody;
            TypeSelector(FrameType.BARYCENTRIC_ROTATING);
            TypeSelector(FrameType.BODY_CENTRED_PARENT_DIRECTION);
          }
        }
      }
    }

    UnityEngine.GUI.DragWindow(
        position : new UnityEngine.Rect(x      : 0f,
                                        y      : 0f,
                                        width  : 10000f,
                                        height : 10000f));

    UnityEngine.GUI.skin = old_skin;
  }

  private void RenderSubtree(CelestialBody celestial, int depth) {
    // Horizontal offset between a node and its children.
    const int offset = 20;
    using (new HorizontalLayout()) {
      if (!celestial.is_root()) {
        UnityEngine.GUILayout.Label(
            "",
            UnityEngine.GUILayout.Width(offset * (depth - 1)));
        string button_text;
        if (celestial.is_leaf()) {
          button_text = "";
        } else if (expanded_[celestial]) {
          button_text = "-";
        } else {
          button_text = "+";
        }
        if (UnityEngine.GUILayout.Button(button_text,
                                         UnityEngine.GUILayout.Width(offset))) {
          Shrink();
          if (!celestial.is_leaf()) {
            expanded_[celestial] = !expanded_[celestial];
          }
        }
      }
      if (UnityEngine.GUILayout.Toggle(selected_celestial == celestial,
                                       celestial.name)) {
        if (selected_celestial != celestial) {
          selected_celestial = celestial;
          if (celestial.is_root() && frame_type != FrameType.BODY_SURFACE) {
            frame_type = FrameType.BODY_CENTRED_NON_ROTATING;
          }
          on_change_(FrameParameters());
        }
      }
    }
    if (celestial.is_root() || (!celestial.is_leaf() && expanded_[celestial])) {
      foreach (CelestialBody child in celestial.orbitingBodies) {
        RenderSubtree(child, depth + 1);
      }
    }
  }

  private void TypeSelector(FrameType value) {
   bool old_wrap = UnityEngine.GUI.skin.toggle.wordWrap;
   UnityEngine.GUI.skin.toggle.wordWrap = true;
   if (UnityEngine.GUILayout.Toggle(
           frame_type == value,
           Description(value, selected_celestial, target_override),
           UnityEngine.GUILayout.Width(150),
           UnityEngine.GUILayout.Height(120))) {
     if (frame_type != value) {
       frame_type = value;
       on_change_(FrameParameters());
     }
    }
    UnityEngine.GUI.skin.toggle.wordWrap = old_wrap;
  }

  private void Shrink() {
    window_rectangle_.height = 0.0f;
    window_rectangle_.width = 0.0f;
  }

  private Callback on_change_;
  // Not owned.
  private IntPtr plugin_;
  private bool show_selector_;
  private UnityEngine.Rect window_rectangle_;
  private Dictionary<CelestialBody, bool> expanded_;
  private readonly string name_;
  private readonly bool is_plotting_frame_selector_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
