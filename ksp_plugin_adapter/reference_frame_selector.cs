﻿#undef HAS_SURFACE
#undef HAS_BODY_CENTRED_ALIGNED_WITH_PARENT

using System;
using System.Collections.Generic;

namespace principia {
namespace ksp_plugin_adapter {

internal static class CelestialExtensions {
  public static bool is_leaf(this CelestialBody celestial) {
    return celestial.orbitingBodies.Count == 0;
  }

  public static bool is_root(this CelestialBody celestial) {
    return celestial.orbit == null;
  }
}

internal class ReferenceFrameSelector : SupervisedWindowRenderer {
  public enum FrameType {
    BARYCENTRIC_ROTATING = 6001,
    BODY_CENTRED_NON_ROTATING = 6000,
    BODY_CENTRED_PARENT_DIRECTION = 6002,
    BODY_SURFACE = 6003,
  }

  public delegate void Callback(NavigationFrameParameters frame_parameters);

  public ReferenceFrameSelector(ISupervisor supervisor,
                                Callback on_change,
                                string name)
      : base(supervisor, UnityEngine.GUILayout.MinWidth(0)) {
    on_change_ = on_change;
    name_ = name;

    // TODO(phl): Bogus initialization.  Find a way to get these data from the
    // C++ side (we do for flight planning).
    frame_type = FrameType.BODY_CENTRED_NON_ROTATING;
    selected_celestial = FlightGlobals.GetHomeBody();
    is_freshly_constructed_ = true;

    expanded_ = new Dictionary<CelestialBody, bool>();
    foreach (CelestialBody celestial in FlightGlobals.Bodies) {
      if (!celestial.is_leaf() && !celestial.is_root()) {
        expanded_.Add(celestial, false);
      }
    }
  }

  public void UpdateMainBody() {
    EffectChange(() => {
      frame_type = FrameType.BODY_CENTRED_NON_ROTATING;
      selected_celestial =
          FlightGlobals.currentMainBody ?? FlightGlobals.GetHomeBody();
      for (CelestialBody celestial = selected_celestial;
           celestial.orbit != null;
           celestial = celestial.orbit.referenceBody) {
        if (!celestial.is_leaf()) {
          expanded_[celestial] = true;
        }
      }
    });
  }

  public void SetFrameParameters(NavigationFrameParameters parameters) {
    EffectChange(() => {
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
    });
  }

  // Sets the |frame_type| to |type| unless this would be invalid for the
  // |selected_celestial|, in which case |frame_type| is set to
  // |BODY_CENTRED_NON_ROTATING|.
  public void SetFrameType(FrameType type) {
    EffectChange(() => {
      if (selected_celestial.is_root() &&
          (type == FrameType.BARYCENTRIC_ROTATING ||
           type == FrameType.BODY_CENTRED_PARENT_DIRECTION)) {
        frame_type = FrameType.BODY_CENTRED_NON_ROTATING;
      } else {
        frame_type = type;
      }
    });
  }

  // Sets the |frame_type| to |BODY_SURFACE| and sets |selected_celestial| to
  // the given |celestial|.
  public void SetToSurfaceFrameOf(CelestialBody celestial) {
    EffectChange(() => {
      frame_type = FrameType.BODY_SURFACE;
      selected_celestial = celestial;
    });
  }

  public static string Name(FrameType type,
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

  public static string ShortName(FrameType type,
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

  public static string Description(FrameType type,
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

  public string Name() {
    return Name(frame_type, selected_celestial, target_override);
  }

  public string ShortName() {
    return ShortName(frame_type, selected_celestial, target_override);
  }

  public string ReferencePlaneDescription() {
    if (!target_override &&
        (frame_type == FrameType.BODY_CENTRED_NON_ROTATING ||
         frame_type == FrameType.BODY_SURFACE)) {
      return $"equator of {selected_celestial.NameWithArticle()}";
    }
    string secondary = target_override != null
        ? "the target"
        : selected_celestial.NameWithArticle();
    string primary = target_override != null
        ? selected_celestial.NameWithArticle()
        : selected_celestial.referenceBody.NameWithArticle();
    return $"orbit of {secondary} around {primary}";
  }

  public CelestialBody[] FixedBodies() {
    if (target_override) {
      return new CelestialBody[]{};
    }
    switch (frame_type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
      case FrameType.BODY_SURFACE:
        return new []{selected_celestial};
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

  public void RenderButton() {
    if (UnityEngine.GUILayout.Button(name_ + " selection (" + Name() +
                                     ")...")) {
      Toggle();
    }
  }

  public FrameType frame_type { get; private set; }
  public CelestialBody selected_celestial { get; private set; }
  public Vessel? target_override { get; set; }

  protected override string Title => name_ + " selection (" + Name() + ")";

  protected override void RenderWindow(int window_id) {
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      // Left-hand side: tree view for celestial selection.
      using (new UnityEngine.GUILayout.VerticalScope(GUILayoutWidth(8))) {
        RenderSubtree(celestial : Planetarium.fetch.Sun, depth : 0);
      }

      // Right-hand side: toggles for reference frame type selection.
      using (new UnityEngine.GUILayout.VerticalScope()) {
        if (target_override) {
          UnityEngine.GUILayout.Label(
              "Using target-centred frame selected on navball speed " +
              "display\n\n" +
              Description(frame_type, selected_celestial, target_override),
              Style.Multiline(UnityEngine.GUI.skin.label),
              GUILayoutWidth(6));
        } else {
          TypeSelector(FrameType.BODY_SURFACE);
          TypeSelector(FrameType.BODY_CENTRED_NON_ROTATING);
          if (!selected_celestial.is_root()) {
            TypeSelector(FrameType.BARYCENTRIC_ROTATING);
            TypeSelector(FrameType.BODY_CENTRED_PARENT_DIRECTION);
          }
        }
      }
    }
    UnityEngine.GUI.DragWindow();
  }

  private void RenderSubtree(CelestialBody celestial, int depth) {
    // Horizontal offset between a node and its children.
    const int offset = 1;
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      if (!celestial.is_root()) {
        UnityEngine.GUILayout.Label("",
                                    GUILayoutWidth(offset * (depth - 1)));
        string button_text;
        if (celestial.is_leaf()) {
          button_text = "";
        } else if (expanded_[celestial]) {
          button_text = "-";
        } else {
          button_text = "+";
        }
        if (UnityEngine.GUILayout.Button(button_text,
                                         GUILayoutWidth(offset))) {
          Shrink();
          if (!celestial.is_leaf()) {
            expanded_[celestial] = !expanded_[celestial];
          }
        }
      }
      if (UnityEngine.GUILayout.Toggle(selected_celestial == celestial,
                                       celestial.name)) {
        EffectChange(() => {
          selected_celestial = celestial;
          if (celestial.is_root() && frame_type != FrameType.BODY_SURFACE) {
            frame_type = FrameType.BODY_CENTRED_NON_ROTATING;
          }
        });
      }
    }
    if (celestial.is_root() || (!celestial.is_leaf() && expanded_[celestial])) {
      foreach (CelestialBody child in celestial.orbitingBodies) {
        RenderSubtree(child, depth + 1);
      }
    }
  }

  private void TypeSelector(FrameType value) {
    var style = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.toggle);
    style.fixedHeight = 0;
    style.wordWrap = true;
    if (UnityEngine.GUILayout.Toggle(
            frame_type == value,
            Description(value, selected_celestial, target_override),
            style,
            GUILayoutWidth(6),
            GUILayoutHeight(5))) {
      EffectChange(() => {
        frame_type = value;
      });
    }
  }

  // Runs an action that may change the frame and run on_change_ if there
  // actually was a change.
  private void EffectChange(Action action) {
    var old_frame_type = frame_type;
    var old_selected_celestial = selected_celestial;
    action();
    if (is_freshly_constructed_ ||
        frame_type != old_frame_type ||
        selected_celestial != old_selected_celestial) {
      on_change_(FrameParameters());
      is_freshly_constructed_ = false;
    }
  }

  private readonly Callback on_change_;
  private readonly string name_;
  private readonly Dictionary<CelestialBody, bool> expanded_;
  private bool is_freshly_constructed_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
