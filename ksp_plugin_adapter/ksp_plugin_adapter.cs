using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

[KSPAddon(startup : KSPAddon.Startup.Flight, once : false)]
public partial class PluginAdapter : UnityEngine.MonoBehaviour {
  // This constant can be at most 32766, since Vectrosity imposes a maximum of
  // 65534 vertices, where there are 2 vertices per point on discrete lines.  We
  // want this to be even since we have two points per line segment.
  // NOTE(egg): Things are fairly slow with the maximum number of points.  We
  // have to do with fewer.  10000 is mostly ok, even fewer would be better.
  // TODO(egg): At the moment we store points in the history  every
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

  private UnityEngine.Rect window_position_;
  private IntPtr plugin_ = IntPtr.Zero;
  // TODO(egg): rendering only one trajectory at the moment.
  private VectorLine rendered_trajectory_;
  private IntPtr rendering_frame_ = IntPtr.Zero;
  private int first_selected_celestial_ = 0;
  private int second_selected_celestial_ = 0;

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
  private void ApplyToVesselsInSpace(VesselProcessor process_vessel) {
    foreach (Vessel vessel in FlightGlobals.Vessels) {
      if (vessel.situation == Vessel.Situations.SUB_ORBITAL ||
          vessel.situation == Vessel.Situations.ORBITING ||
          vessel.situation == Vessel.Situations.ESCAPING) {
        process_vessel(vessel);
      }
    }
  }

  #region Unity Lifecycle
  // See the Unity manual on execution order for more information on |Start()|,
  // |OnDestroy()| and |FixedUpdate()|.
  // http://docs.unity3d.com/Manual/ExecutionOrder.html
  private void Start() {
    Log.Info("principia.ksp_plugin_adapter.PluginAdapter.Start()");
    RenderingManager.AddToPostDrawQueue(queueSpot    : 3,
                                        drawFunction : new Callback(DrawGUI));
    window_position_ = new UnityEngine.Rect(
        left   : UnityEngine.Screen.width / 2.0f,
        top    : UnityEngine.Screen.height / 2.0f,
        width  : 10,
        height : 10);
  }

  private void OnDestroy() {
    Log.Info("principia.ksp_plugin_adapter.PluginAdapter.OnDestroy()");
    RenderingManager.RemoveFromPostDrawQueue(
        queueSpot    : 3,
        drawFunction : new Callback(DrawGUI));
    Cleanup();
  }

  private void FixedUpdate() {
    if (PluginRunning()) {
      double universal_time = Planetarium.GetUniversalTime();
      AdvanceTime(plugin_, universal_time, Planetarium.InverseRotAngle);
      BodyProcessor update_body = body => {
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
      };
      ApplyToBodyTree(update_body);
      VesselProcessor update_vessel = vessel => {
        bool inserted = InsertOrKeepVessel(
            plugin_,
            vessel.id.ToString(),
            vessel.orbit.referenceBody.flightGlobalsIndex);
        if (inserted) {
          SetVesselStateOffset(plugin_,
                               vessel.id.ToString(),
                               (XYZ)vessel.orbit.pos,
                               (XYZ)vessel.orbit.vel);
        }
        Vector3d position =
            (Vector3d)VesselDisplacementFromParent(plugin_,
                                                   vessel.id.ToString());
        Vector3d velocity =
            (Vector3d)VesselParentRelativeVelocity(plugin_,
                                                   vessel.id.ToString());
        vessel.orbit.UpdateFromStateVectors(pos: position, vel: velocity,
                                            refBody: vessel.orbit.referenceBody,
                                            UT: universal_time);
        // Work around a KSP bug: |Orbit.pos| for a vessel in an elliptic orbit
        // corresponds to the position one timestep in the future.
        // TODO(egg): do it.
      };
      ApplyToVesselsInSpace(update_vessel);
      Vessel active_vessel = FlightGlobals.ActiveVessel;
      if (MapView.MapIsEnabled && 
              (active_vessel.situation == Vessel.Situations.SUB_ORBITAL ||
               active_vessel.situation == Vessel.Situations.ORBITING ||
               active_vessel.situation == Vessel.Situations.ESCAPING)) {
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
      }
      LogALot();
    }
  }

  #endregion

  private void Cleanup() {
    DeletePlugin(ref plugin_);
    DeleteRenderingFrame(ref rendering_frame_);
    if (rendered_trajectory_ != null) {
      Vector.DestroyLine(ref rendered_trajectory_);
    }
  }

  private void DrawGUI() {
    UnityEngine.GUI.skin = HighLogic.Skin;
    window_position_ = UnityEngine.GUILayout.Window(
        id         : 1,
        screenRect : window_position_,
        func       : DrawMainWindow,
        text       : "Traces of Various Descriptions",
        options    : UnityEngine.GUILayout.MinWidth(500));
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
    IntPtr hello_ptr = SayHello();
    UnityEngine.GUILayout.TextArea(text : Marshal.PtrToStringAnsi(hello_ptr));
    if (UnityEngine.GUILayout.Button(PluginRunning() ? "Stop plugin"
                                                     : "Start plugin")) {
      if (PluginRunning()) {
        Cleanup();
      } else {
        InitializePlugin();
      }
    }
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
              value : second_selected_celestial_ == celestial.flightGlobalsIndex,
              text  : celestial.name) &&
          second_selected_celestial_ != celestial.flightGlobalsIndex) {
        second_selected_celestial_ = celestial.flightGlobalsIndex;
        changed_reference_frame = true;
      }
      UnityEngine.GUILayout.EndHorizontal();
      if (changed_reference_frame && PluginRunning()) {
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
    }
    if (PluginRunning()) {
      Vessel active_vessel = FlightGlobals.ActiveVessel;
      UnityEngine.GUILayout.TextArea(
          "+ Kraken : " +
              (((Vector3d)active_vessel.rb_velocity) +
                   Krakensbane.GetFrameVelocity()));
      UnityEngine.GUILayout.TextArea(
          "+ Kraken + getRFrmVel: " +
              (((Vector3d)active_vessel.rb_velocity) +
                   Krakensbane.GetFrameVelocity() + 
                   active_vessel.orbit.referenceBody.getRFrmVel(
                       active_vessel.CoM)));
      UnityEngine.GUILayout.TextArea(
          "Principia \"world\", rotating : " +
          (Vector3d)VesselWorldVelocity(
              plugin_,
              active_vessel.id.ToString(),
              new XYZ{x = 0, y = 0, z = 0},
              active_vessel.orbit.referenceBody.rotationPeriod));
      UnityEngine.GUILayout.TextArea(
          "Principia \"world\", no rotation : " +
          (Vector3d)VesselWorldVelocity(
              plugin_,
              active_vessel.id.ToString(),
              new XYZ{x = 0, y = 0, z = 0},
              double.PositiveInfinity));
      UnityEngine.GUILayout.TextArea(
          "Principia \"world\", expected : " +
          (Vector3d)VesselWorldVelocity(
              plugin_,
              active_vessel.id.ToString(),
              new XYZ{x = 0, y = 0, z = 0},
              Planetarium.FrameIsRotating()
                  ? active_vessel.orbit.referenceBody.rotationPeriod
                  : double.PositiveInfinity));
      UnityEngine.GUILayout.TextArea(
          "GetVel : " +
              (Vector3d)active_vessel.orbit.GetVel());
      UnityEngine.GUILayout.TextArea(
          "Root part @ CoM world velocity + Kraken: " +
              (Vector3d)
                  (active_vessel.rootPart.rb.GetPointVelocity(
                       (Vector3d)active_vessel.findWorldCenterOfMass()) +
                   Krakensbane.GetFrameVelocity()));
      UnityEngine.GUILayout.TextArea(
          "CoM : " +
          ((Vector3d)active_vessel.CoM));
      UnityEngine.GUILayout.TextArea(
          "found CoM world position : " +
          ((Vector3d)active_vessel.findWorldCenterOfMass()));
      UnityEngine.GUILayout.TextArea(
          "Principia world : " +
          (Vector3d)VesselWorldPosition(
              plugin_,
              active_vessel.id.ToString(),
              (XYZ)active_vessel.orbit.referenceBody.position));
    }
    UnityEngine.GUILayout.EndVertical();
    UnityEngine.GUI.DragWindow(
        position : new UnityEngine.Rect(left : 0f, top : 0f, width : 10000f,
                                        height : 20f));
  }

  private void LogALot() {
    Vessel active_vessel = FlightGlobals.ActiveVessel;
    Log.Info("UT : " + Planetarium.GetUniversalTime());
    Log.Info(
        "Principia world position : " +
        (Vector3d)VesselWorldPosition(
            plugin_,
            active_vessel.id.ToString(),
            (XYZ)active_vessel.orbit.referenceBody.position));
    Log.Info(
        "Principia world velocity (rotating) : " +
        (Vector3d)VesselWorldVelocity(
            plugin_,
            active_vessel.id.ToString(),
            new XYZ{x = 0, y = 0, z = 0},
            active_vessel.orbit.referenceBody.rotationPeriod));
    Log.Info(
        "Principia world velocity (no rotation) : " +
        (Vector3d)VesselWorldVelocity(
            plugin_,
            active_vessel.id.ToString(),
            new XYZ{x = 0, y = 0, z = 0},
            double.PositiveInfinity));
    Log.Info("reference body position : " +
             active_vessel.orbit.referenceBody.position);
    Log.Info("reference body GetVel : " +
             active_vessel.orbit.referenceBody.orbit.GetVel());
    Log.Info("active vessel found world CoM (32) : " +
             (Vector3d)active_vessel.findWorldCenterOfMass());
    Log.Info("Root part at found world CoM world velocity (32) : " +
             (Vector3d)active_vessel.rootPart.rb.GetPointVelocity(
                 (Vector3d)active_vessel.findWorldCenterOfMass()));
    Log.Info("active vessel orbit.pos : " + active_vessel.orbit.pos);
    Log.Info("active vessel orbit.vel : " + active_vessel.orbit.vel);
    Log.Info("active vessel GetVel : " + active_vessel.orbit.GetVel());
    Log.Info("Principia orbit.pos : " + 
             (Vector3d)VesselDisplacementFromParent(
                 plugin_,
                 active_vessel.id.ToString()));
    Log.Info("Principia orbit.vel : " + 
             (Vector3d)VesselParentRelativeVelocity(
                 plugin_,
                 active_vessel.id.ToString()));
    Log.Info("Kraken : " + Krakensbane.GetFrameVelocity());
    Log.Info("Root rb velocity (32) : " +
             (Vector3d)active_vessel.rootPart.rb.velocity);
    Log.Info("Root part partTransform.position (32) : " +
             (Vector3d)active_vessel.rootPart.partTransform.position);
    Log.Info("Measured orbit.pos : " + active_vessel.orbit.pos);
    Log.Info("Measured orbit.vel : " + active_vessel.orbit.vel);
    Log.Info("reference body reference body position : " +
             active_vessel.orbit.referenceBody.referenceBody.position);
    Log.Info("reference body Principia orbit.pos : " +
             (Vector3d)CelestialDisplacementFromParent(
                 plugin_,
                 active_vessel.orbit.referenceBody.flightGlobalsIndex));
    Log.Info("reference body Principia orbit.vel : " +
             (Vector3d)CelestialParentRelativeVelocity(
                 plugin_,
                 active_vessel.orbit.referenceBody.flightGlobalsIndex));
    Log.Info("reference body measured orbit.pos : " +
             active_vessel.orbit.referenceBody.orbit.pos);
    Log.Info("reference body measured orbit.vel : " +
             active_vessel.orbit.referenceBody.orbit.vel);
    Log.Info("Duna position : " +
             FlightGlobals.Bodies[6].position);
    Log.Info("Ike Principia orbit.pos : " +
             (Vector3d)CelestialDisplacementFromParent(
                 plugin_,
                 FlightGlobals.Bodies[7].flightGlobalsIndex));
    Log.Info("Ike Principia orbit.vel : " +
             (Vector3d)CelestialParentRelativeVelocity(
                 plugin_,
                 FlightGlobals.Bodies[7].flightGlobalsIndex));
    Log.Info("Ike measured orbit.pos : " +
             FlightGlobals.Bodies[7].orbit.pos);
    Log.Info("Ike measured orbit.vel : " +
             FlightGlobals.Bodies[7].orbit.vel);
  }

  private void InitializePlugin() {
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
    first_selected_celestial_ = 0;
    second_selected_celestial_ = 0;
    rendering_frame_ =
        NewBodyCentredNonRotatingFrame(plugin_, first_selected_celestial_);
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
    ApplyToVesselsInSpace(insert_vessel);
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
