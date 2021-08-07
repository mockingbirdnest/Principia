#undef HAS_SURFACE
#undef HAS_BODY_CENTRED_ALIGNED_WITH_PARENT

using System;
using System.Collections.Generic;
using KSP.Localization;

namespace principia {
namespace ksp_plugin_adapter {

internal static class CelestialExtensions {
  public static bool is_leaf(this CelestialBody celestial, Vessel target) {
    return celestial.orbitingBodies.Count == 0 &&
           target?.orbit.referenceBody != celestial;
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

  public delegate void Callback(NavigationFrameParameters? frame_parameters,
                                Vessel target_vessel);

  public ReferenceFrameSelector(ISupervisor supervisor,
                                Callback on_change,
                                string name) : base(
      supervisor,
      UnityEngine.GUILayout.MinWidth(0)) {
    on_change_ = on_change;
    name_ = name;

    // TODO(phl): Bogus initialization.  Find a way to get these data from the
    // C++ side (we do for flight planning).
    frame_type = FrameType.BODY_CENTRED_NON_ROTATING;
    selected_celestial = FlightGlobals.GetHomeBody();
    is_freshly_constructed_ = true;

    expanded_ = new Dictionary<CelestialBody, bool>();
    foreach (CelestialBody celestial in FlightGlobals.Bodies) {
      if (!celestial.is_root()) {
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
        if (!celestial.is_leaf(target)) {
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

  // Sets the frame type to the surface frame of some body; this is used when the
  // speed display is switched to surface mode.
  public void SetToSurfaceFrame() {
    EffectChange(() => {
      target_frame_selected = false;
      frame_type = FrameType.BODY_SURFACE;
    });
  }

  public bool IsSurfaceFrame() {
    return frame_type == FrameType.BODY_SURFACE;
  }

  public void SetToOrbitalFrame() {
    EffectChange(() => {
      target_frame_selected = false;
      if (selected_celestial.is_root() &&
          (last_orbital_type_ == FrameType.BARYCENTRIC_ROTATING ||
            last_orbital_type_ == FrameType.BODY_CENTRED_PARENT_DIRECTION)) {
        frame_type = FrameType.BODY_CENTRED_NON_ROTATING;
      } else {
        frame_type = last_orbital_type_;
      }
    });
  }

  public void SetTargetFrame() {
    EffectChange(() => {
      target_frame_selected = true;
    });
  }

  public void UnsetTargetFrame() {
    EffectChange(() => {
      target_frame_selected = false;
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

  private static string TargetFrameName(Vessel target) {
    return Localizer.Format("#Principia_ReferenceFrameSelector_Name_Target",
                            target.orbit.referenceBody.name);
  }

  private static string Name(FrameType type,
                             CelestialBody selected) {
    switch (type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
        return Localizer.Format(
            "#Principia_ReferenceFrameSelector_Name_BodyCentredNonRotating",
            selected.name);
      case FrameType.BARYCENTRIC_ROTATING:
        if (selected.is_root()) {
          throw Log.Fatal("Naming barycentric rotating frame of root body");
        } else {
          return Localizer.Format(
              "#Principia_ReferenceFrameSelector_Name_BarycentricRotating",
              selected.referenceBody.name,
              selected.name);
        }
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        if (selected.is_root()) {
          throw Log.Fatal(
              "Naming parent-direction rotating frame of root body");
        } else {
          return Localizer.Format(
              "#Principia_ReferenceFrameSelector_Name_BodyCentredParentDirection",
              selected.name,
              selected.referenceBody.name);
        }
      case FrameType.BODY_SURFACE:
        return Localizer.Format(
            "#Principia_ReferenceFrameSelector_Name_BodySurface",
            selected.name);
      default:
        throw Log.Fatal("Unexpected type " + type.ToString());
    }
  }

  private static string TargetFrameShortName(Vessel target) {
    return Localizer.Format(
        "#Principia_ReferenceFrameSelector_ShortName_Target",
        target.orbit.referenceBody.name[0]);
  }

  private static string ShortName(FrameType type,
                                  CelestialBody selected) {
    switch (type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
        return Localizer.Format(
            "#Principia_ReferenceFrameSelector_ShortName_BodyCentredNonRotating",
            selected.name[0]);
      case FrameType.BARYCENTRIC_ROTATING:
        if (selected.is_root()) {
          throw Log.Fatal("Naming barycentric rotating frame of root body");
        } else {
          return Localizer.Format(
              "#Principia_ReferenceFrameSelector_ShortName_BarycentricRotating",
              selected.referenceBody.name[0],
              selected.name[0]);
        }
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        if (selected.is_root()) {
          throw Log.Fatal(
              "Naming parent-direction rotating frame of root body");
        } else {
          return Localizer.Format(
              "#Principia_ReferenceFrameSelector_ShortName_BodyCentredParentDirection",
              selected.name[0],
              selected.referenceBody.name[0]);
        }
      case FrameType.BODY_SURFACE:
        return Localizer.Format(
            "#Principia_ReferenceFrameSelector_ShortName_BodySurface",
            selected.name[0]);
      default:
        throw Log.Fatal("Unexpected type " + type.ToString());
    }
  }

  private static string TargetFrameDescription(Vessel target) {
    return Localizer.Format(
        "#Principia_ReferenceFrameSelector_Description_Target",
        target.vesselName,
        target.orbit.referenceBody.NameWithArticle());
  }

  private static string Description(FrameType type,
                                    CelestialBody selected) {
    switch (type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
        return Localizer.Format(
            "#Principia_ReferenceFrameSelector_Description_BodyCentredNonRotating",
            selected.NameWithArticle());
      case FrameType.BARYCENTRIC_ROTATING:
        if (selected.is_root()) {
          throw Log.Fatal("Describing barycentric rotating frame of root body");
        } else {
          return Localizer.Format(
              "#Principia_ReferenceFrameSelector_Description_BarycentricRotating",
              selected.NameWithArticle(),
              selected.referenceBody.NameWithArticle());
        }
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        if (selected.is_root()) {
          throw Log.Fatal(
              "Describing parent-direction rotating frame of root body");
        } else {
          return Localizer.Format(
              "#Principia_ReferenceFrameSelector_Description_BodyCentredParentDirection",
              selected.NameWithArticle(),
              selected.referenceBody.NameWithArticle());
        }
      case FrameType.BODY_SURFACE:
        return Localizer.Format(
            "#Principia_ReferenceFrameSelector_Description_BodySurface",
            selected.NameWithArticle());
      default:
        throw Log.Fatal("Unexpected type " + type.ToString());
    }
  }

  public string Name() {
    return target_frame_selected ? TargetFrameName(target)
                                 : Name(frame_type, selected_celestial);
  }

  public string ShortName() {
    return target_frame_selected ? TargetFrameShortName(target)
                                 : ShortName(frame_type, selected_celestial);
  }

  public string ReferencePlaneDescription() {
    if (!target_frame_selected &&
        (frame_type == FrameType.BODY_CENTRED_NON_ROTATING ||
         frame_type == FrameType.BODY_SURFACE)) {
      return Localizer.Format(
          "#Principia_ReferenceFrameSelector_ReferencePlane_Centred",
          selected_celestial.NameWithArticle());
    }
    string secondary =
        target_frame_selected
            ? Localizer.Format(
                "#Principia_ReferenceFrameSelector_ReferencePlane_Secondary_Target")
            : selected_celestial.NameWithArticle();
    string primary = target_frame_selected
                         ? selected_celestial.NameWithArticle()
                         : selected_celestial.referenceBody.NameWithArticle();
    return Localizer.Format("#Principia_ReferenceFrameSelector_ReferencePlane",
                            secondary,
                            primary);
  }

  // If the reference frames is defined by two bodies, |OrientingBody()| is the
  // one that is not fixed, but instead defines the orientation.  If the
  // reference frame is defined from a single body, |OrientingBody()| is null.
  public CelestialBody OrientingBody() {
    if (target_frame_selected) {
      return target.orbit.referenceBody;
    }
    switch (frame_type) {
      case FrameType.BARYCENTRIC_ROTATING:
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        return selected_celestial.referenceBody;
      default:
        return null;
    }
  }

  public bool FixesBody(CelestialBody celestial) {
    // TODO(egg): When we have the rotating-pulsating frame, this should return
    // true for both bodies.
    return celestial == Centre();
  }

  public CelestialBody Centre() {
    if (target_frame_selected) {
      return null;
    }
    switch (frame_type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
      case FrameType.BODY_SURFACE:
        return selected_celestial;
      case FrameType.BARYCENTRIC_ROTATING:
        return null;
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
            centre_index = selected_celestial.flightGlobalsIndex
        };
      case FrameType.BARYCENTRIC_ROTATING:
        return new NavigationFrameParameters{
            extension = (int)frame_type,
            primary_index = selected_celestial.referenceBody.flightGlobalsIndex,
            secondary_index = selected_celestial.flightGlobalsIndex
        };
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        // We put the primary body as secondary, because the one we want fixed
        // is the secondary body (which means it has to be the primary in the
        // terminology of |BodyCentredBodyDirection|).
        return new NavigationFrameParameters{
            extension = (int)frame_type,
            primary_index = selected_celestial.flightGlobalsIndex,
            secondary_index =
                selected_celestial.referenceBody.flightGlobalsIndex
        };
      default:
        throw Log.Fatal("Unexpected frame_type " + frame_type.ToString());
    }
  }

  public void RenderButton() {
    if (UnityEngine.GUILayout.Button(
        Localizer.Format("#Principia_ReferenceFrameSelector_ToggleButton",
                         name_,
                         Name()))) {
      Toggle();
    }
  }

  public FrameType frame_type { get; private set; }
  private CelestialBody selected_celestial { get; set; }
  public Vessel target { get; set; }
  public bool target_frame_selected { get; private set; }

  protected override string Title =>
      Localizer.Format("#Principia_ReferenceFrameSelector_Title",
                       name_,
                       Name());

  protected override void RenderWindow(int window_id) {
    using (new UnityEngine.GUILayout.VerticalScope()) {
      UnityEngine.GUILayout.Label(
          Localizer.Format(
              "#Principia_ReferenceFrameSelector_Description_Heading",
              Name(),
              ShortName()));
      UnityEngine.GUILayout.Label(
          target_frame_selected ? TargetFrameDescription(target)
                                : Description(frame_type, selected_celestial),
          Style.Multiline(UnityEngine.GUI.skin.label),
          GUILayoutHeight(4));
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        // Left-hand side: tree view of celestials.
        using (new UnityEngine.GUILayout.VerticalScope(GUILayoutWidth(8))) {
          RenderSubtree(celestial : Planetarium.fetch.Sun, depth : 0);
        }

        // Right-hand side: toggles for reference frame type selection.
        using (new UnityEngine.GUILayout.VerticalScope()) {
          RenderSubtreeToggleGrid(Planetarium.fetch.Sun);
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
        UnityEngine.GUILayout.Space(Width(offset * (depth - 1)));
        if (celestial.is_leaf(target)) {
          UnityEngine.GUILayout.Button(
              "", UnityEngine.GUI.skin.label, GUILayoutWidth(offset));
        } else {
          string button_text = expanded_[celestial] ? "−" : "+";
          if (UnityEngine.GUILayout.Button(
                  button_text, GUILayoutWidth(offset))) {
            Shrink();
            expanded_[celestial] = !expanded_[celestial];
          }
        }
      }
      UnityEngine.GUILayout.Label(celestial.name);
    }
    if (celestial.is_root() ||
        (!celestial.is_leaf(target) && expanded_[celestial])) {
      if (target?.orbit.referenceBody == celestial) {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Space(Width(offset * depth));
          UnityEngine.GUILayout.Button("", GUILayoutWidth(offset));
          UnityEngine.GUILayout.Label(
              Localizer.Format("#Principia_ReferenceFrameSelector_Target",
                               target.vesselName));
        }
      }
      foreach (CelestialBody child in celestial.orbitingBodies) {
        RenderSubtree(child, depth + 1);
      }
    }
  }

  private void RenderSubtreeToggleGrid(CelestialBody celestial) {
    int column_width = 2;
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      if (ToggleButton(
              SelectedFrameIs(celestial, FrameType.BODY_CENTRED_NON_ROTATING),
              Localizer.Format("#Principia_ReferenceFrameSelector_Centre"),
              GUILayoutWidth(column_width))) {
        EffectChange(() => {
          target_frame_selected = false;
          selected_celestial = celestial;
          frame_type = FrameType.BODY_CENTRED_NON_ROTATING;
        });
      }
      if (ToggleButton(
              SelectedFrameIs(celestial, FrameType.BODY_SURFACE),
              Localizer.Format("#Principia_ReferenceFrameSelector_Surface"),
              GUILayoutWidth(column_width))) {
        EffectChange(() => {
          target_frame_selected = false;
          selected_celestial = celestial;
          frame_type = FrameType.BODY_SURFACE;
        });
      }
      if (!celestial.is_root() &&
          ToggleButton(
              SelectedFrameIs(celestial,
                              FrameType.BODY_CENTRED_PARENT_DIRECTION),
              Localizer.Format("#Principia_ReferenceFrameSelector_Orbit"),
              GUILayoutWidth(column_width))) {
        EffectChange(() => {
          target_frame_selected = false;
          selected_celestial = celestial;
          frame_type = FrameType.BODY_CENTRED_PARENT_DIRECTION;
        });
      }
    }
    if (celestial.is_root() ||
        (!celestial.is_leaf(target) && expanded_[celestial])) {
      if (target?.orbit.referenceBody == celestial) {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Button("", UnityEngine.GUI.skin.label, GUILayoutWidth(column_width));
          UnityEngine.GUILayout.Button("", UnityEngine.GUI.skin.label, GUILayoutWidth(column_width));
          if (ToggleButton(
                  target_frame_selected,
                  Localizer.Format("#Principia_ReferenceFrameSelector_Orbit"),
                  GUILayoutWidth(column_width))) {
            EffectChange(() => {
              target_frame_selected = true;
            });
          }
        }
      }
      foreach (CelestialBody child in celestial.orbitingBodies) {
        RenderSubtreeToggleGrid(child);
      }
    }
  }

  private bool SelectedFrameIs(CelestialBody celestial, FrameType type) {
    return !target_frame_selected &&
        selected_celestial == celestial && frame_type == type;
  }

  // Runs an action that may change the frame and run on_change_ if there
  // actually was a change.
  private void EffectChange(Action action) {
    var old_frame_type = frame_type;
    var old_selected_celestial = selected_celestial;
    var target_frame_was_selected = target_frame_selected;
    action();
    if (is_freshly_constructed_ ||
        frame_type != old_frame_type ||
        selected_celestial != old_selected_celestial ||
        target_frame_selected != target_frame_was_selected) {
      on_change_(
          target_frame_selected ? null
                                : (NavigationFrameParameters?)FrameParameters(),
          target_frame_selected ? target : null);
      is_freshly_constructed_ = false;
    }
    if (frame_type != FrameType.BODY_SURFACE) {
      last_orbital_type_ = frame_type;
    }
  }

  // Functionally like UnityEngine.GuiLayout.Toggle, but as a button rather than
  // a checkbox.
  private bool ToggleButton(bool value, string text,
                            params UnityEngine.GUILayoutOption[] options) {
    return UnityEngine.GUILayout.Toolbar(selected: value ? 0 : -1,
                                         new string[1] {text},
                                         Style.LitToggleButton(),
                                         options) == 0;
  }

  private readonly Callback on_change_;
  private readonly string name_;
  private readonly Dictionary<CelestialBody, bool> expanded_;
  private bool is_freshly_constructed_;
  private ReferenceFrameSelector.FrameType last_orbital_type_ =
      ReferenceFrameSelector.FrameType.BODY_CENTRED_NON_ROTATING;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
