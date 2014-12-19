using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

[KSPAddon(startup : KSPAddon.Startup.MainMenu, once : false)]
public partial class PluginAdapter : UnityEngine.MonoBehaviour {
  // This constant can be at most 32766, since Vectrosity imposes a maximum of
  // 65534 vertices, where there are 2 vertices per point on discrete lines.  We
  // want this to be even since we have two points per line segment.
  // NOTE(egg): Things are fairly slow with the maximum number of points.  We
  // have to do with fewer.  10000 is mostly ok, even fewer would be better.
  // TODO(egg): At the moment we store points in the history every
  // 10 n seconds, where n is maximal such that 10 n seconds is less than the
  // length of a |FixedUpdate|. This means we sometimes have very large gaps.
  // We should store *all* points of the history, then decimate for rendering.
  // This means splines are not needed, since 10 s is small enough to give the
  // illusion of continuity on the scales we are dealing with, and cubics are
  // rendered by Vectrosity as line segments, so that a cubic rendered as
  // 10 segments counts 20 towards |kLinePoints| (and probably takes as long to
  // render as 10 segments from the actual data, with extra overhead for
  // the evaluation of the cubic).
  private const int kLinePoints = 10000;

  private const int kGUIQueueSpot = 3;

  private UnityEngine.Rect main_window_position_;
  private IntPtr plugin_ = IntPtr.Zero;
  // TODO(egg): rendering only one trajectory at the moment.
  private VectorLine rendered_trajectory_;
  private IntPtr rendering_frame_ = IntPtr.Zero;
  private int first_selected_celestial_ = 0;
  private int second_selected_celestial_ = 0;

  private bool time_is_advancing_;

  private bool should_initialize_on_next_fixed_update_ = false;

  private static bool an_instance_is_loaded;

  PluginAdapter() {
    // We create this directory here so we do not need to worry about cross-
    // platform problems in C++.
    System.IO.Directory.CreateDirectory("glog/Principia");
    Log.InitGoogleLogging();
  }

  ~PluginAdapter() {
    Cleanup();
  }

  private bool PluginRunning() {
    return plugin_ != IntPtr.Zero;
  }

  private delegate void BodyProcessor(CelestialBody body);
  private delegate void VesselProcessor(Vessel vessel);

  // Applies |process_body| to all bodies but the sun in the tree of celestial
  // bodies, in topological order.
  private void ApplyToBodyTree(BodyProcessor process_body) {
    // Tree traversal (DFS, not that it matters).
    Stack<CelestialBody> stack = new Stack<CelestialBody>();
    foreach (CelestialBody child in Planetarium.fetch.Sun.orbitingBodies) {
      stack.Push(child);
    }
    CelestialBody body;
    while (stack.Count > 0) {
      body = stack.Pop();
      process_body(body);
      foreach (CelestialBody child in body.orbitingBodies) {
        stack.Push(child);
      }
    }
  }

  // Applies |process_vessel| to all vessels in space.
  private void ApplyToVesselsOnRailsOrInInertialPhysicsBubbleInSpace(
      VesselProcessor process_vessel) {
    var vessels = from vessel in FlightGlobals.Vessels
                  where is_on_rails_in_space(vessel) || 
                        is_in_inertial_physics_bubble_in_space(vessel)
                  select vessel;
    foreach (Vessel vessel in vessels) {
      process_vessel(vessel);
    }
  }

  private void ApplyToVesselsInPhysicsBubble(VesselProcessor process_vessel) {
    var vessels = from vessel in FlightGlobals.Vessels
                  where is_in_physics_bubble(vessel)
                  select vessel;
    foreach (Vessel vessel in vessels) {
      process_vessel(vessel);
    }
  }

  private void UpdateBody(CelestialBody body, double universal_time) {
    UpdateCelestialHierarchy(plugin_,
                             body.flightGlobalsIndex,
                             body.orbit.referenceBody.flightGlobalsIndex);
    Vector3d position =
        (Vector3d)CelestialDisplacementFromParent(plugin_,
                                                  body.flightGlobalsIndex);
    Vector3d velocity =
        (Vector3d)CelestialParentRelativeVelocity(plugin_,
                                                  body.flightGlobalsIndex);
    // TODO(egg): Some of this might be be superfluous and redundant.
    Orbit original = body.orbit;
    Orbit copy = new Orbit(original.inclination, original.eccentricity,
                           original.semiMajorAxis, original.LAN,
                           original.argumentOfPeriapsis,
                           original.meanAnomalyAtEpoch, original.epoch,
                           original.referenceBody);
    copy.UpdateFromStateVectors(position,
                                velocity,
                                copy.referenceBody,
                                universal_time);
    body.orbit.inclination = copy.inclination;
    body.orbit.eccentricity = copy.eccentricity;
    body.orbit.semiMajorAxis = copy.semiMajorAxis;
    body.orbit.LAN = copy.LAN;
    body.orbit.argumentOfPeriapsis = copy.argumentOfPeriapsis;
    body.orbit.meanAnomalyAtEpoch = copy.meanAnomalyAtEpoch;
    body.orbit.epoch = copy.epoch;
    body.orbit.referenceBody = copy.referenceBody;
    body.orbit.Init();
    body.orbit.UpdateFromUT(universal_time);
    body.CBUpdate();
    body.orbit.UpdateFromStateVectors(position,
                                      velocity,
                                      copy.referenceBody,
                                      universal_time);
  }

  private void UpdateVessel(Vessel vessel, double universal_time) {
    bool inserted = InsertOrKeepVessel(
        plugin_,
        vessel.id.ToString(),
        vessel.orbit.referenceBody.flightGlobalsIndex);
    if (inserted) {
      SetVesselStateOffset(plugin               : plugin_,
                           vessel_guid          : vessel.id.ToString(),
                           from_parent_position : (XYZ)vessel.orbit.pos,
                           from_parent_velocity : (XYZ)vessel.orbit.vel);
    }
    Vector3d position =
        (Vector3d)VesselDisplacementFromParent(plugin_,
                                               vessel.id.ToString());
    Vector3d velocity =
        (Vector3d)VesselParentRelativeVelocity(plugin_,
                                               vessel.id.ToString());
    // NOTE(egg): Here we work around a KSP bug: |Orbit.pos| for a vessel
    // corresponds to the position one timestep in the future.  This is not
    // the case for celestial bodies.
    vessel.orbit.UpdateFromStateVectors(
        pos     : position + velocity * UnityEngine.Time.deltaTime,
        vel     : velocity,
        refBody : vessel.orbit.referenceBody,
        UT      : universal_time);
  }

  private void AddToPhysicsBubble(Vessel vessel) {
    bool inserted = InsertOrKeepVessel(
        plugin_,
        vessel.id.ToString(),
        vessel.orbit.referenceBody.flightGlobalsIndex);
    if (inserted) {
      // NOTE(egg): these degrees of freedom are off by one Δt, but they
      // should never actually be used.
      // TODO(egg): we shouldn't have to do this.
      SetVesselStateOffset(plugin               : plugin_,
                           vessel_guid          : vessel.id.ToString(),
                           from_parent_position : (XYZ)vessel.orbit.pos,
                           from_parent_velocity : (XYZ)vessel.orbit.vel);
    }
    Log.Info("vessel has " + vessel.parts.Count() + " parts");
    Vector3d gravity =
        FlightGlobals.getGeeForceAtPosition(vessel.findWorldCenterOfMass());
    Vector3d kraken_velocity = Krakensbane.GetFrameVelocity();
    KSPPart[] parts =
        (from part in vessel.parts
         where part.rb != null  // Physicsless parts have no rigid body.
         select new KSPPart {
             world_position = (XYZ)(Vector3d)part.rb.worldCenterOfMass,
             world_velocity = (XYZ)(kraken_velocity + part.rb.velocity),
             mass = (double)part.mass + (double)part.GetResourceMass(),
             gravitational_acceleration_to_be_applied_by_ksp = (XYZ)gravity,
             id = part.flightID}).ToArray();
    AddVesselToNextPhysicsBubble(plugin      : plugin_,
                                 vessel_guid : vessel.id.ToString(),
                                 parts       : parts,
                                 count       : parts.Count());
  }

  private bool is_in_space(Vessel vessel) {
    return vessel.situation == Vessel.Situations.SUB_ORBITAL ||
           vessel.situation == Vessel.Situations.ORBITING ||
           vessel.situation == Vessel.Situations.ESCAPING;
  }

  private bool is_on_rails_in_space(Vessel vessel) {
    return vessel.packed && is_in_space(vessel);
  }

  private bool is_in_physics_bubble(Vessel vessel) {
    return !vessel.packed && vessel.loaded;
  }

  private bool is_in_inertial_physics_bubble_in_space(Vessel vessel) {
    return is_in_physics_bubble(vessel) &&
           !Planetarium.FrameIsRotating() &&
           is_in_space(vessel);
  }

  private bool has_inertial_physics_bubble_in_space() {
    Vessel active_vessel = FlightGlobals.ActiveVessel;
    return active_vessel != null &&
           is_in_inertial_physics_bubble_in_space(active_vessel);
  }

  #region Unity Lifecycle
  // See the Unity manual on execution order for more information on |Start()|,
  // |OnDestroy()| and |FixedUpdate()|.
  // http://docs.unity3d.com/Manual/ExecutionOrder.html

  // Awake is called once.
  private void Awake() {
    Log.Info("principia.ksp_plugin_adapter.PluginAdapter.Awake()");
    if (an_instance_is_loaded) {
      Log.Info("an instance was loaded");
      UnityEngine.Object.Destroy(gameObject);
    } else {
      UnityEngine.Object.DontDestroyOnLoad(gameObject);
      an_instance_is_loaded = true;
    }
    GameEvents.onGameStateLoad.Add(InitializeOnGameStateLoad);
    main_window_position_ = new UnityEngine.Rect(
        left   : UnityEngine.Screen.width / 2.0f,
        top    : UnityEngine.Screen.height / 2.0f,
        width  : 10,
        height : 10);
    ApplyToBodyTree(body => body.inverseRotThresholdAltitude =
                                body.maxAtmosphereAltitude);
  }

  private void OnGUI() {
    UnityEngine.GUI.skin = HighLogic.Skin;
    main_window_position_ = UnityEngine.GUILayout.Window(
        id         : 1,
        screenRect : main_window_position_,
        func       : DrawMainWindow,
        text       : "Traces of Various Descriptions",
        options    : UnityEngine.GUILayout.MinWidth(500));
  }

  private void FixedUpdate() {
    if (should_initialize_on_next_fixed_update_) {
      ResetPlugin();
      should_initialize_on_next_fixed_update_ = false;
    }
    if (PluginRunning()) {
      double universal_time = Planetarium.GetUniversalTime();
      double plugin_time = current_time(plugin_);
      if (plugin_time > universal_time) {
        Log.Fatal("Closed Timelike Curve");
      } else if (plugin_time == universal_time) {
        time_is_advancing_ = false;
        return;
      }
      time_is_advancing_ = true;
      if (has_inertial_physics_bubble_in_space()) {
        ApplyToVesselsInPhysicsBubble(AddToPhysicsBubble);
      }
      AdvanceTime(plugin_, universal_time, Planetarium.InverseRotAngle);
      ApplyToBodyTree(body => UpdateBody(body, universal_time));
      ApplyToVesselsOnRailsOrInInertialPhysicsBubbleInSpace(vessel => UpdateVessel(vessel, universal_time));
      Vessel active_vessel = FlightGlobals.ActiveVessel;
      if (has_inertial_physics_bubble_in_space()) {
        Vector3d displacement_offset =
            (Vector3d)BubbleDisplacementCorrection(
                          plugin_,
                          (XYZ)Planetarium.fetch.Sun.position);
        Vector3d velocity_offset =
            (Vector3d)BubbleVelocityCorrection(
                          plugin_,
                          active_vessel.orbit.referenceBody.flightGlobalsIndex);
        Krakensbane krakensbane =
            (Krakensbane)FindObjectOfType(typeof(Krakensbane));
        krakensbane.setOffset(displacement_offset);
        krakensbane.FrameVel += velocity_offset;
      }
      if (MapView.MapIsEnabled &&
          active_vessel != null &&
          (is_on_rails_in_space(active_vessel) ||
           is_in_inertial_physics_bubble_in_space(active_vessel))) {
        if (active_vessel.orbitDriver.Renderer.drawMode !=
                OrbitRenderer.DrawMode.OFF ||
            active_vessel.orbitDriver.Renderer.drawIcons !=
                OrbitRenderer.DrawIcons.OBJ ||
            active_vessel.patchedConicRenderer != null) {
          Log.Info("Removing orbit rendering for the active vessel");
          active_vessel.orbitDriver.Renderer.drawMode =
              OrbitRenderer.DrawMode.OFF;
          active_vessel.orbitDriver.Renderer.drawIcons =
              OrbitRenderer.DrawIcons.OBJ;
          active_vessel.DetachPatchedConicsSolver();
          active_vessel.patchedConicRenderer = null;
        }
        IntPtr trajectory_iterator = IntPtr.Zero;
        try {
          trajectory_iterator = RenderedVesselTrajectory(
              plugin_,
              active_vessel.id.ToString(),
              rendering_frame_,
              (XYZ)Planetarium.fetch.Sun.position);

          LineSegment segment;
          int index_in_line_points = kLinePoints -
              NumberOfSegments(trajectory_iterator) * 2;
          while (index_in_line_points < 0) {
            FetchAndIncrement(trajectory_iterator);
            index_in_line_points += 2;
          }
          while (!AtEnd(trajectory_iterator)) {
            segment = FetchAndIncrement(trajectory_iterator);
            // TODO(egg): should we do the |LocalToScaledSpace| conversion in
            // native code?
            // TODO(egg): could we directly assign to
            // |rendered_trajectory_.points3| from C++ using unsafe code and
            // something like the following?
            // |fixed (UnityEngine.Vector3* pts = rendered_trajectory_.points3)|
            rendered_trajectory_.points3[index_in_line_points++] =
                ScaledSpace.LocalToScaledSpace((Vector3d)segment.begin);
            rendered_trajectory_.points3[index_in_line_points++] =
                ScaledSpace.LocalToScaledSpace((Vector3d)segment.end);
          }
        } finally {
          DeleteLineAndIterator(ref trajectory_iterator);
        }
        if (MapView.Draw3DLines) {
          Vector.DrawLine3D(rendered_trajectory_);
        } else {
          Vector.DrawLine(rendered_trajectory_);
        }
      } else {
        ResetRenderedTrajectory();
      }
    }
  }

  #endregion

  private void ResetRenderedTrajectory() {
    DestroyRenderedTrajectory();
    rendered_trajectory_ = new VectorLine(
        lineName     : "rendered_trajectory_",
        linePoints   : new UnityEngine.Vector3[kLinePoints],
        lineMaterial : MapView.OrbitLinesMaterial,
        color        : XKCDColors.AcidGreen,
        width        : 5,
        lineType     : LineType.Discrete);
    rendered_trajectory_.vectorObject.transform.parent =
        ScaledSpace.Instance.transform;
    rendered_trajectory_.vectorObject.renderer.castShadows = false;
    rendered_trajectory_.vectorObject.renderer.receiveShadows = false;
    rendered_trajectory_.layer = 31;
  }

  private void DestroyRenderedTrajectory() {
    if (rendered_trajectory_ != null) {
      Vector.DestroyLine(ref rendered_trajectory_);
    }
  }

  private void Cleanup() {
    DeletePlugin(ref plugin_);
    DeleteRenderingFrame(ref rendering_frame_);
    DestroyRenderedTrajectory();
  }

  private void InitializeOnGameStateLoad(ConfigNode node) {
    should_initialize_on_next_fixed_update_ = true;
  }

  private void DrawMainWindow(int window_id) {
    UnityEngine.GUIStyle style = new UnityEngine.GUIStyle(
        UnityEngine.GUI.skin.button);
    style.normal.textColor = style.focused.textColor = UnityEngine.Color.white;
    style.hover.textColor = style.active.textColor = UnityEngine.Color.yellow;
    style.onNormal.textColor  = UnityEngine.Color.green;
    style.onFocused.textColor = UnityEngine.Color.green;
    style.onHover.textColor   = UnityEngine.Color.green;
    style.onActive.textColor  = UnityEngine.Color.green;
    style.padding             = new UnityEngine.RectOffset(8, 8, 8, 8);

    UnityEngine.GUILayout.BeginVertical();
    String plugin_state;
    if (PluginRunning()) {
      if (UnityEngine.GUILayout.Button(text : "Force Stop")) {
        Cleanup();
      }
    } else {
      if (UnityEngine.GUILayout.Button(text : "Force Start")) {
        ResetPlugin();
      }
    }
    if (!PluginRunning()) {
      plugin_state = "not started";
    } else if (!time_is_advancing_) {
      plugin_state = "holding";
    } else if (!has_inertial_physics_bubble_in_space()) {
      plugin_state = "running";
    } else {
      plugin_state = "managing physics bubble";
    }
    UnityEngine.GUILayout.TextArea(text : "Plugin is " + plugin_state);
    bool barycentric_rotating =
        first_selected_celestial_ != second_selected_celestial_;
    UnityEngine.GUILayout.TextArea(
        text : "Frame is " + (barycentric_rotating ? " " : "not ") + "rotating");
    String reference_frame_description;
    if (barycentric_rotating) {
      reference_frame_description =
          "The reference frame fixing the barycentre of " +
          FlightGlobals.Bodies[first_selected_celestial_].theName + " and " +
          FlightGlobals.Bodies[second_selected_celestial_].theName + ", " +
          "the line through them, and the plane in which they move about the " +
          "barycentre.";
    } else {
      reference_frame_description =
          "The nonrotating reference frame fixing the centre of " +
          FlightGlobals.Bodies[first_selected_celestial_].theName + ".";
    }
    UnityEngine.GUILayout.TextArea(text: reference_frame_description);
    UnityEngine.GUILayout.Label(text : "Reference frame selection:");
    foreach (CelestialBody celestial in FlightGlobals.Bodies) {
      bool changed_reference_frame = false;
      UnityEngine.GUILayout.BeginHorizontal();
      if (UnityEngine.GUILayout.Toggle(
              value : first_selected_celestial_ == celestial.flightGlobalsIndex,
              text  : "") &&
          first_selected_celestial_ != celestial.flightGlobalsIndex) {
        first_selected_celestial_ = celestial.flightGlobalsIndex;
        changed_reference_frame = true;
      }
      if (UnityEngine.GUILayout.Toggle(
              value : second_selected_celestial_ ==
                          celestial.flightGlobalsIndex,
              text  : celestial.name) &&
          second_selected_celestial_ != celestial.flightGlobalsIndex) {
        second_selected_celestial_ = celestial.flightGlobalsIndex;
        changed_reference_frame = true;
      }
      UnityEngine.GUILayout.EndHorizontal();
      if (changed_reference_frame && PluginRunning()) {
        UpdateRenderingFrame();
      }
    }
    UnityEngine.GUILayout.EndVertical();
    UnityEngine.GUI.DragWindow(
        position : new UnityEngine.Rect(left : 0f, top : 0f, width : 10000f,
                                        height : 20f));
  }

  private void UpdateRenderingFrame() {
    DeleteRenderingFrame(ref rendering_frame_);
    if (first_selected_celestial_ == second_selected_celestial_) {
      rendering_frame_ = NewBodyCentredNonRotatingFrame(
                             plugin_,
                             first_selected_celestial_);
    } else {
      rendering_frame_ = NewBarycentricRotatingFrame(
                             plugin_,
                             first_selected_celestial_,
                             second_selected_celestial_);
    }
  }

  private void ResetPlugin() {
    Cleanup();
    ResetRenderedTrajectory();
    plugin_ = NewPlugin(Planetarium.GetUniversalTime(),
                        Planetarium.fetch.Sun.flightGlobalsIndex,
                        Planetarium.fetch.Sun.gravParameter,
                        Planetarium.InverseRotAngle);
    BodyProcessor insert_body = body => {
      Log.Info("Inserting " + body.name + "...");
      InsertCelestial(plugin_,
                      body.flightGlobalsIndex,
                      body.gravParameter,
                      body.orbit.referenceBody.flightGlobalsIndex,
                      (XYZ)body.orbit.pos,
                      (XYZ)body.orbit.vel);
    };
    ApplyToBodyTree(insert_body);
    EndInitialization(plugin_);
    UpdateRenderingFrame();
    VesselProcessor insert_vessel = vessel => {
      Log.Info("Inserting " + vessel.name + "...");
      bool inserted =
          InsertOrKeepVessel(plugin_,
                             vessel.id.ToString(),
                             vessel.orbit.referenceBody.flightGlobalsIndex);
      if (!inserted) {
        Log.Fatal("Plugin initialization: vessel not inserted");
      } else {
        SetVesselStateOffset(plugin_,
                             vessel.id.ToString(),
                             (XYZ)vessel.orbit.pos,
                             (XYZ)vessel.orbit.vel);
      }
    };
    ApplyToVesselsOnRailsOrInInertialPhysicsBubbleInSpace(insert_vessel);
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
