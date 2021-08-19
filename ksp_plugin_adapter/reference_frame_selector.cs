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
    pinned_ = new Dictionary<CelestialBody, bool>();
    foreach (CelestialBody celestial in FlightGlobals.Bodies) {
      expanded_.Add(celestial, false);
      pinned_.Add(celestial, false);
    }
  }

  public void UpdateMainBody() {
    EffectChange(() => {
      frame_type = FrameType.BODY_CENTRED_NON_ROTATING;
      selected_celestial =
          FlightGlobals.currentMainBody ?? FlightGlobals.GetHomeBody();
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
    return !target_frame_selected && frame_type == FrameType.BODY_SURFACE;
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
    return L10N.CelestialString(
        "#Principia_ReferenceFrameSelector_Name_Target",
        new []{target.orbit.referenceBody});
  }

  private static string Name(FrameType type,
                             CelestialBody selected) {
    switch (type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
        return L10N.CelestialString(
            "#Principia_ReferenceFrameSelector_Name_BodyCentredNonRotating",
            new[]{selected});
      case FrameType.BARYCENTRIC_ROTATING:
        if (selected.is_root()) {
          throw Log.Fatal("Naming barycentric rotating frame of root body");
        } else {
          return "DEPRECATED";
        }
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        if (selected.is_root()) {
          throw Log.Fatal(
              "Naming parent-direction rotating frame of root body");
        } else {
          return L10N.ZWSPToHyphenBetweenNonCJK(L10N.CelestialString(
              "#Principia_ReferenceFrameSelector_Name_BodyCentredParentDirection",
              new[]{selected, selected.referenceBody}));
        }
      case FrameType.BODY_SURFACE:
        return L10N.CelestialString(
            "#Principia_ReferenceFrameSelector_Name_BodySurface",
            new[]{selected});
      default:
        throw Log.Fatal("Unexpected type " + type.ToString());
    }
  }

  private static string TargetFrameAbbreviation(Vessel target) {
    return L10N.CelestialStringOrNull(
        "#Principia_ReferenceFrameSelector_Abbreviation_Target",
        new[]{target.orbit.referenceBody});
  }

  private static string TargetFrameNavballName(Vessel target) {
    return TargetFrameAbbreviation(target) ?? L10N.CelestialStringOrNull(
        "#Principia_ReferenceFrameSelector_NavballName_Target",
        new[]{target.orbit.referenceBody}) ?? TargetFrameName(target);
  }

  private static string TargetFrameSelectorText(Vessel target) {
    return TargetFrameAbbreviation(target) ?? Localizer.Format(
        "#Principia_ReferenceFrameSelector_SelectorText_Target");
  }

  private static string TargetFrameSelectorTooltip(Vessel target) {
    string name = TargetFrameName(target);
    return Localizer.Format(
        "#Principia_ReferenceFrameSelector_Tooltip_Target",
        name);
  }

  private static string Abbreviation(FrameType type, CelestialBody selected) {
    switch (type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
        return L10N.CelestialStringOrNull(
            "#Principia_ReferenceFrameSelector_Abbreviation_BodyCentredNonRotating",
            new[]{selected});
      case FrameType.BARYCENTRIC_ROTATING:
        return "DEPRECATED";
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        if (selected.is_root()) {
          throw Log.Fatal(
              "Naming parent-direction rotating frame of root body");
        } else {
          return L10N.CelestialStringOrNull(
              "#Principia_ReferenceFrameSelector_Abbreviation_BodyCentredParentDirection",
              new[]{selected, selected.referenceBody});
        }
      case FrameType.BODY_SURFACE:
        return L10N.CelestialStringOrNull(
            "#Principia_ReferenceFrameSelector_Abbreviation_BodySurface",
            new[]{selected});
      default:
        throw Log.Fatal("Unexpected type " + type.ToString());
    }
  }

  private static string NavballName(FrameType type,
                                    CelestialBody selected) {
    string result = Abbreviation(type, selected);
    if (result != null) {
      return result;
    }
    switch (type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
        result = L10N.CelestialStringOrNull(
            "#Principia_ReferenceFrameSelector_NavballName_BodyCentredNonRotating",
            new[]{selected});
        break;
      case FrameType.BARYCENTRIC_ROTATING:
        result = "DEPRECATED";
        break;
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        if (selected.is_root()) {
          throw Log.Fatal(
              "Naming parent-direction rotating frame of root body");
        } else {
          result = L10N.CelestialStringOrNull(
              "#Principia_ReferenceFrameSelector_NavballName_BodyCentredParentDirection",
              new[]{selected, selected.referenceBody});
        }
        break;
      case FrameType.BODY_SURFACE:
        result = L10N.CelestialStringOrNull(
            "#Principia_ReferenceFrameSelector_NavballName_BodySurface",
            new[]{selected});
        break;
      default:
        throw Log.Fatal("Unexpected type " + type.ToString());
    }
    if (result != null) {
      return result;
    }
    return Name(type, selected);
  }

  private static string SelectorText(FrameType type,
                                     CelestialBody selected) {
    string abbreviation = Abbreviation(type, selected);
    if (abbreviation != null) {
      return abbreviation;
    }
    switch (type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
        return Localizer.Format(
            "#Principia_ReferenceFrameSelector_SelectorText_BodyCentredNonRotating");
      case FrameType.BARYCENTRIC_ROTATING:
        return "DEPRECATED";
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        return Localizer.Format(
            "#Principia_ReferenceFrameSelector_SelectorText_BodyCentredParentDirection");
      case FrameType.BODY_SURFACE:
        return Localizer.Format(
            "#Principia_ReferenceFrameSelector_SelectorText_BodySurface");
      default:
        throw Log.Fatal("Unexpected type " + type.ToString());
    }
  }

  private static string SelectorTooltip(FrameType type,
                                        CelestialBody selected) {
    string name = Name(type, selected);
    switch (type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
        return Localizer.Format(
            "#Principia_ReferenceFrameSelector_Tooltip_BodyCentredNonRotating",
            name,
            selected.Name());
      case FrameType.BARYCENTRIC_ROTATING:
        return "DEPRECATED";
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        return Localizer.Format(
            "#Principia_ReferenceFrameSelector_Tooltip_BodyCentredParentDirection",
            name,
            selected.Name());
      case FrameType.BODY_SURFACE:
        return Localizer.Format(
            "#Principia_ReferenceFrameSelector_Tooltip_BodySurface",
            name,
            selected.Name());
      default:
        throw Log.Fatal("Unexpected type " + type.ToString());
    }
    
  }

  private static string TargetFrameDescription(Vessel target) {
    return Localizer.Format(
        "#Principia_ReferenceFrameSelector_Description_Target",
        target.vesselName,
        target.orbit.referenceBody.Name());
  }

  private static string Description(FrameType type,
                                    CelestialBody selected) {
    switch (type) {
      case FrameType.BODY_CENTRED_NON_ROTATING:
        return L10N.CelestialString(
            "#Principia_ReferenceFrameSelector_Description_BodyCentredNonRotating",
            new[]{selected});
      case FrameType.BARYCENTRIC_ROTATING:
        return "DEPRECATED";
      case FrameType.BODY_CENTRED_PARENT_DIRECTION:
        if (selected.is_root()) {
          throw Log.Fatal(
              "Describing parent-direction rotating frame of root body");
        } else {
          return L10N.CelestialString(
              "#Principia_ReferenceFrameSelector_Description_BodyCentredParentDirection",
              new[]{selected, selected.referenceBody});
        }
      case FrameType.BODY_SURFACE:
        return L10N.CelestialString(
            "#Principia_ReferenceFrameSelector_Description_BodySurface",
            new[]{selected});
      default:
        throw Log.Fatal("Unexpected type " + type.ToString());
    }
  }

  private string Abbreviation() {
    return target_frame_selected ? TargetFrameAbbreviation(target)
                                 : Abbreviation(frame_type, selected_celestial);
  }

  public string Name() {
    return target_frame_selected ? TargetFrameName(target)
                                 : Name(frame_type, selected_celestial);
  }

  public string NavballName() {
    return target_frame_selected ? TargetFrameNavballName(target)
                                 : NavballName(frame_type, selected_celestial);
  }

  public string ReferencePlaneDescription() {
    if (!target_frame_selected &&
        (frame_type == FrameType.BODY_CENTRED_NON_ROTATING ||
         frame_type == FrameType.BODY_SURFACE)) {
      return Localizer.Format(
          "#Principia_ReferenceFrameSelector_ReferencePlane_Centred",
          selected_celestial.Name());
    }
    string secondary =
        target_frame_selected
            ? Localizer.Format(
                "#Principia_ReferenceFrameSelector_ReferencePlane_Secondary_Target")
            : selected_celestial.Name();
    string primary = target_frame_selected
                         ? selected_celestial.Name()
                         : selected_celestial.referenceBody.Name();
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
      case FrameType.BODY_CENTRED_NON_ROTATING:
      case FrameType.BODY_SURFACE:
        return null;
      default:
        throw Log.Fatal("Unexpected frame_type " + frame_type.ToString());
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
              Abbreviation()));
      UnityEngine.GUILayout.Label(
          target_frame_selected ? TargetFrameDescription(target)
                                : Description(frame_type, selected_celestial),
          Style.Multiline(UnityEngine.GUI.skin.label),
          GUILayoutHeight(3));
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

  private bool AnyDescendantPinned(CelestialBody celestial) {
    if (pinned_[celestial]) {
      return true;
    }
    if (target_pinned_ && target?.orbit.referenceBody == celestial) {
      return true;
    }
    foreach (CelestialBody body in celestial.orbitingBodies) {
      if (AnyDescendantPinned(body)) {
        return true;
      }
    }
    return false;
  }

  private void RenderSubtree(CelestialBody celestial, int depth) {
    // Horizontal offset between a node and its children.
    const int offset = 1;
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Space(Width(offset * depth));
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
      UnityEngine.GUILayout.Label(celestial.StandaloneName());
      UnityEngine.GUILayout.FlexibleSpace();
      if (celestial.is_root()) {
        UnityEngine.GUILayout.Label(
            Localizer.Format("#Principia_ReferenceFrameSelector_Pin"));
      } else if (UnityEngine.GUILayout.Toggle(pinned_[celestial], "") !=
                 pinned_[celestial]) {
        pinned_[celestial] = !pinned_[celestial];
        Shrink();
      }
    }
    if (!celestial.is_leaf(target)) {
      if ((expanded_[celestial] || target_pinned_) &&
          target?.orbit.referenceBody == celestial) {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Space(Width(offset * (depth + 1)));
          UnityEngine.GUILayout.Button(
              "", UnityEngine.GUI.skin.label, GUILayoutWidth(offset));
          UnityEngine.GUILayout.Label(
              Localizer.Format("#Principia_ReferenceFrameSelector_Target",
                               target.vesselName));
          UnityEngine.GUILayout.FlexibleSpace();
          if (UnityEngine.GUILayout.Toggle(target_pinned_, "") !=
              target_pinned_) {
            target_pinned_ = !target_pinned_;
            Shrink();
          }
        }
      }
      foreach (CelestialBody child in celestial.orbitingBodies) {
        if (expanded_[celestial] || AnyDescendantPinned(child)) {
          RenderSubtree(child, depth + 1);
        }
      }
    }
  }

  private void RenderSubtreeToggleGrid(CelestialBody celestial) {
    int column_width = 2;
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      if (ToggleButton(
              SelectedFrameIs(celestial, FrameType.BODY_CENTRED_NON_ROTATING),
              new UnityEngine.GUIContent(
                  SelectorText(FrameType.BODY_CENTRED_NON_ROTATING, celestial),
                  SelectorTooltip(FrameType.BODY_CENTRED_NON_ROTATING, celestial)),
              GUILayoutWidth(column_width))) {
        EffectChange(() => {
          target_frame_selected = false;
          selected_celestial = celestial;
          frame_type = FrameType.BODY_CENTRED_NON_ROTATING;
        });
      }
      if (ToggleButton(
              SelectedFrameIs(celestial, FrameType.BODY_SURFACE),
              new UnityEngine.GUIContent(
                  SelectorText(FrameType.BODY_SURFACE, celestial),
                  SelectorTooltip(FrameType.BODY_SURFACE, celestial)),
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
              new UnityEngine.GUIContent(
                  SelectorText(FrameType.BODY_CENTRED_PARENT_DIRECTION, celestial),
                  SelectorTooltip(FrameType.BODY_CENTRED_PARENT_DIRECTION, celestial)),
              GUILayoutWidth(column_width))) {
        EffectChange(() => {
          target_frame_selected = false;
          selected_celestial = celestial;
          frame_type = FrameType.BODY_CENTRED_PARENT_DIRECTION;
        });
      }
    }
    if (!celestial.is_leaf(target)) {
      if ((expanded_[celestial] || target_pinned_) &&
          target?.orbit.referenceBody == celestial) {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Button("", UnityEngine.GUI.skin.label, GUILayoutWidth(column_width));
          UnityEngine.GUILayout.Button("", UnityEngine.GUI.skin.label, GUILayoutWidth(column_width));
          if (ToggleButton(
                  target_frame_selected,
                  new UnityEngine.GUIContent(
                    TargetFrameSelectorText(target),
                    TargetFrameSelectorTooltip(target)),
                  GUILayoutWidth(column_width))) {
            EffectChange(() => {
              target_frame_selected = true;
            });
          }
        }
      }
      foreach (CelestialBody child in celestial.orbitingBodies) {
        if (expanded_[celestial] || AnyDescendantPinned(child)) {
          RenderSubtreeToggleGrid(child);
        }
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
    if (is_freshly_constructed_) {
      pinned_[selected_celestial] = true;
      if (!selected_celestial.is_leaf(target)) {
        expanded_[selected_celestial] = true;
      }
    }
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
  private bool ToggleButton(bool value, UnityEngine.GUIContent content,
                            params UnityEngine.GUILayoutOption[] options) {
    return UnityEngine.GUILayout.Button(
        content,
        value ? Style.LitToggleButton() : Style.DarkToggleButton(),
        options) ? !value : value;
  }

  private readonly Callback on_change_;
  private readonly string name_;
  private readonly Dictionary<CelestialBody, bool> expanded_;
  private readonly Dictionary<CelestialBody, bool> pinned_;
  private bool target_pinned_ = true;
  private bool is_freshly_constructed_;
  private ReferenceFrameSelector.FrameType last_orbital_type_ =
      ReferenceFrameSelector.FrameType.BODY_CENTRED_NON_ROTATING;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
