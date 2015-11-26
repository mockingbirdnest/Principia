using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

[KSPScenario(createOptions: ScenarioCreationOptions.AddToAllGames,
             tgtScenes: new GameScenes[]{GameScenes.SPACECENTER,
                                         GameScenes.FLIGHT,
                                         GameScenes.TRACKSTATION})]
public partial class PrincipiaPluginAdapter : ScenarioModule {

  private const String kPrincipiaKey = "serialized_plugin";
  private const String kPrincipiaInitialState = "principia_initial_state";
  private const String kPrincipiaGravityModel = "principia_gravity_model";
  private const double kΔt = 10;

  // The number of points in a |VectorLine| can be at most 32766, since
  // Vectrosity imposes a maximum of 65534 vertices, where there are 2 vertices
  // per point on discrete lines.
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
  private const int kMaxVectorLinePoints = 32766;

  private ApplicationLauncherButton toolbar_button_;
  private bool hide_all_gui_ = false;

  // "Persistant" is a KSP typo.
  [KSPField(isPersistant = true)]
  private bool show_main_window_ = true;
  [KSPField(isPersistant = true)]
  private int main_window_x_ = UnityEngine.Screen.width / 2;
  [KSPField(isPersistant = true)]
  private int main_window_y_ = UnityEngine.Screen.height / 3;
  private UnityEngine.Rect main_window_rectangle_;

  private IntPtr plugin_ = IntPtr.Zero;
  // TODO(egg): rendering only one trajectory at the moment.
  private VectorLine rendered_prediction_;
  private VectorLine rendered_trajectory_;
  private IntPtr rendering_frame_ = IntPtr.Zero;

  [KSPField(isPersistant = true)]
  private int first_selected_celestial_ = 0;
  [KSPField(isPersistant = true)]
  private int second_selected_celestial_ = 0;

  [KSPField(isPersistant = true)]
  private bool display_patched_conics_ = false;
  [KSPField(isPersistant = true)]
  private bool fix_navball_in_plotting_frame_ = true;
  [KSPField(isPersistant = true)]
  private bool force_2d_trajectories_ = false;

  private readonly double[] prediction_length_tolerances_ =
      {1E-3, 1E-2, 1E0, 1E1, 1E2, 1E3, 1E4};
  [KSPField(isPersistant = true)]
  private int prediction_length_tolerance_index_ = 1;
  private readonly double[] prediction_lengths_ =
      {1 << 10, 1 << 11, 1 << 12, 1 << 13, 1 << 14, 1 << 15, 1 << 16, 1 << 17,
       1 << 18, 1 << 19, 1 << 20, 1 << 21, 1 << 22, 1 << 23, 1 << 24, 1 << 25,
       1 << 26, 1 << 27, 1 << 28};
  [KSPField(isPersistant = true)]
  private int prediction_length_index_ = 0;
  private readonly double[] history_lengths_ =
      {1 << 10, 1 << 11, 1 << 12, 1 << 13, 1 << 14, 1 << 15, 1 << 16, 1 << 17,
       1 << 18, 1 << 19, 1 << 20, 1 << 21, 1 << 22, 1 << 23, 1 << 24, 1 << 25,
       1 << 26, 1 << 27, 1 << 28, 1 << 29, double.PositiveInfinity};
  [KSPField(isPersistant = true)]
  private int history_length_index_ = 10;

  [KSPField(isPersistant = true)]
  private bool show_reference_frame_selection_ = true;
  [KSPField(isPersistant = true)]
  private bool show_prediction_settings_ = true;
  [KSPField(isPersistant = true)]
  private bool show_ksp_features_ = false;
  [KSPField(isPersistant = true)]
  private bool show_logging_settings_ = false;
  [KSPField(isPersistant = true)]
  private bool show_reset_button_ = false;

  [KSPField(isPersistant = true)]
  private int verbose_logging_ = 0;
  [KSPField(isPersistant = true)]
  private int suppressed_logging_ = 0;
  [KSPField(isPersistant = true)]
  private int stderr_logging_ = 2;
  [KSPField(isPersistant = true)]
  private int buffered_logging_ = 0;
#if CRASH_BUTTON
  [KSPField(isPersistant = true)]
  private bool show_crash_options_ = false;
#endif

  private bool time_is_advancing_;

  private DateTime plugin_construction_;

  private MapRenderer map_renderer_;

  private enum PluginSource {
   SAVED_STATE,
   ORBITAL_ELEMENTS,
   CARTESIAN_CONFIG,
  }
  private PluginSource plugin_source_;

  private Krakensbane krakensbane_;
  private NavBall navball_;
  private UnityEngine.Texture compass_navball_texture_;
  private UnityEngine.Texture inertial_navball_texture_;
  private UnityEngine.Texture barycentric_navball_texture_;
  private bool navball_changed_ = true;

  // The RSAS is the component of the stock KSP autopilot that deals with
  // orienting the vessel towards a specific direction (e.g. prograde).
  // It is, as usual for KSP, an ineffable acronym; it is however likely derived
  // from the name of the SAS, the component of the autopilot that deals with
  // stabilizing the vessel's attitude without fixing it to any particular
  // target.  Note that SAS has several known meanings, none of which are
  // noteworthy.  The interested reader can refer to
  // http://wiki.kerbalspaceprogram.com/wiki/SAS.
  private bool override_rsas_target_ = false;
  private Vector3d rsas_target_;
  private bool reset_rsas_target_ = false;

  PrincipiaPluginAdapter() {
    // We create this directory here so we do not need to worry about cross-
    // platform problems in C++.
    System.IO.Directory.CreateDirectory("glog/Principia");
    Log.InitGoogleLogging();
  }

  ~PrincipiaPluginAdapter() {
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
    QP from_parent = CelestialFromParent(plugin_, body.flightGlobalsIndex);
    // TODO(egg): Some of this might be be superfluous and redundant.
    Orbit original = body.orbit;
    Orbit copy = new Orbit(original.inclination, original.eccentricity,
                           original.semiMajorAxis, original.LAN,
                           original.argumentOfPeriapsis,
                           original.meanAnomalyAtEpoch, original.epoch,
                           original.referenceBody);
    copy.UpdateFromStateVectors((Vector3d)from_parent.q,
                                (Vector3d)from_parent.p,
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
    body.orbit.UpdateFromStateVectors((Vector3d)from_parent.q,
                                      (Vector3d)from_parent.p,
                                      copy.referenceBody,
                                      universal_time);
  }

  private void UpdateVessel(Vessel vessel, double universal_time) {
    bool inserted = InsertOrKeepVessel(
        plugin_,
        vessel.id.ToString(),
        vessel.orbit.referenceBody.flightGlobalsIndex);
    if (inserted) {
      SetVesselStateOffset(plugin      : plugin_,
                           vessel_guid : vessel.id.ToString(),
                           from_parent : new QP{q = (XYZ)vessel.orbit.pos,
                                                p = (XYZ)vessel.orbit.vel});
    }
    QP from_parent = VesselFromParent(plugin_, vessel.id.ToString());
    // NOTE(egg): Here we work around a KSP bug: |Orbit.pos| for a vessel
    // corresponds to the position one timestep in the future.  This is not
    // the case for celestial bodies.
    vessel.orbit.UpdateFromStateVectors(
        pos     : (Vector3d)from_parent.q +
                  (Vector3d)from_parent.p * UnityEngine.Time.deltaTime,
        vel     : (Vector3d)from_parent.p,
        refBody : vessel.orbit.referenceBody,
        UT      : universal_time);
  }

  private void AddToPhysicsBubble(Vessel vessel) {
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
    if (parts.Count() > 0) {
      bool inserted = InsertOrKeepVessel(
          plugin_,
          vessel.id.ToString(),
          vessel.orbit.referenceBody.flightGlobalsIndex);
      if (inserted) {
        // NOTE(egg): this is only used when a (plugin-managed) physics bubble
        // appears with a new vessel (e.g. when exiting the atmosphere).
        // TODO(egg): these degrees of freedom are off by one Δt and we don't
        // compensate for the pos/vel synchronization bug.
        SetVesselStateOffset(plugin      : plugin_,
                             vessel_guid : vessel.id.ToString(),
                             from_parent : new QP{q = (XYZ)vessel.orbit.pos,
                                                  p = (XYZ)vessel.orbit.vel});
      }
      AddVesselToNextPhysicsBubble(plugin      : plugin_,
                                   vessel_guid : vessel.id.ToString(),
                                   parts       : parts,
                                   count       : parts.Count());
    }
  }

  private bool is_in_space(Vessel vessel) {
    return vessel.state != Vessel.State.DEAD &&
           (vessel.situation == Vessel.Situations.SUB_ORBITAL ||
            vessel.situation == Vessel.Situations.ORBITING ||
            vessel.situation == Vessel.Situations.ESCAPING);
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

  private bool has_active_vessel_in_space() {
    Vessel active_vessel = FlightGlobals.ActiveVessel;
    return active_vessel != null &&
           (is_on_rails_in_space(active_vessel) ||
           is_in_inertial_physics_bubble_in_space(active_vessel));
  }

  private bool draw_active_vessel_trajectory() {
    return MapView.MapIsEnabled &&
           has_active_vessel_in_space();
  }

  private void OverrideRSASTarget(FlightCtrlState state) {
    if (override_rsas_target_ && FlightGlobals.ActiveVessel.Autopilot.Enabled) {
      FlightGlobals.ActiveVessel.Autopilot.RSAS.SetTargetOrientation(
          rsas_target_,
          reset_rsas_target_);
    }
    reset_rsas_target_ = false;
  }

  // Returns false and nulls |texture| if the file does not exist.
  private bool LoadTextureIfExists(out UnityEngine.Texture texture,
                                   String path) {
    string full_path =
        KSPUtil.ApplicationRootPath + Path.DirectorySeparatorChar +
        "GameData" + Path.DirectorySeparatorChar +
        "Principia" + Path.DirectorySeparatorChar +
        "assets" + Path.DirectorySeparatorChar +
        path;
    if (File.Exists(full_path)) {
      var texture2d = new UnityEngine.Texture2D(2, 2);
      bool success = texture2d.LoadImage(
          File.ReadAllBytes(full_path));
      if (!success) {
        Log.Fatal("Failed to load texture " + full_path);
      }
      texture = texture2d;
      return true;
    } else {
      texture = null;
      return false;
    }
  }

  private void LoadTextureOrDie(out UnityEngine.Texture texture, String path) {
    bool success = LoadTextureIfExists(out texture, path);
    if (!success) {
      Log.Fatal("Missing texture " + path);
    }
  }

  #region ScenarioModule lifecycle
  // These functions override virtual ones from |ScenarioModule|, but it seems
  // that they're actually called by reflection, so that bad things happen
  // if you don't have, e.g., a function called |OnAwake()| that calls
  // |base.OnAwake()|.  It doesn't matter whether the functions are public or
  // private, overriding or hiding though.

  public override void OnAwake() {
    base.OnAwake();
    // While we're here, we might as well log.
    Log.Info("principia.ksp_plugin_adapter.PrincipiaPluginAdapter.OnAwake()");

    LoadTextureIfExists(out compass_navball_texture_, "navball_compass.png");
    LoadTextureOrDie(out inertial_navball_texture_, "navball_inertial.png");
    LoadTextureOrDie(out barycentric_navball_texture_,
                     "navball_barycentric.png");

    GameEvents.onShowUI.Add(ShowGUI);
    GameEvents.onHideUI.Add(HideGUI);
  }

  public override void OnSave(ConfigNode node) {
    base.OnSave(node);
    if (PluginRunning()) {
      IntPtr serialization = IntPtr.Zero;
      IntPtr serializer = IntPtr.Zero;
      for (;;) {
        try {
          serialization = SerializePlugin(plugin_, ref serializer);
          if (serialization == IntPtr.Zero) {
            break;
          }
          node.AddValue(kPrincipiaKey, Marshal.PtrToStringAnsi(serialization));
        } finally {
          DeletePluginSerialization(ref serialization);
        }
      }
    }
  }

  public override void OnLoad(ConfigNode node) {
    base.OnLoad(node);
    if (node.HasValue(kPrincipiaKey)) {
      Cleanup();
      SetRotatingFrameThresholds();
      Log.SetBufferedLogging(buffered_logging_);
      Log.SetSuppressedLogging(suppressed_logging_);
      Log.SetStderrLogging(stderr_logging_);
      Log.SetVerboseLogging(verbose_logging_);

      IntPtr deserializer = IntPtr.Zero;
      String[] serializations = node.GetValues(kPrincipiaKey);
      Log.Info("Serialization has " + serializations.Length + " chunks");
      foreach (String serialization in serializations) {
        Log.Info("serialization is " + serialization.Length +
                 " characters long");
        DeserializePlugin(serialization,
                          serialization.Length,
                          ref deserializer,
                          ref plugin_);
      }
      DeserializePlugin("", 0, ref deserializer, ref plugin_);

      UpdateRenderingFrame();
      plugin_construction_ = DateTime.Now;
      plugin_source_ = PluginSource.SAVED_STATE;
    } else {
      Log.Warning("No principia state found, creating one");
      ResetPlugin();
    }
  }

  #endregion

  #region Unity Lifecycle
  // See the Unity manual on execution order for more information.
  // http://docs.unity3d.com/Manual/ExecutionOrder.html

  private void OnGUI() {
    if (ApplicationLauncher.Ready && toolbar_button_ == null) {
      UnityEngine.Texture toolbar_button_texture;
      LoadTextureOrDie(out toolbar_button_texture, "toolbar_button.png");
      toolbar_button_ =
          ApplicationLauncher.Instance.AddModApplication(
              onTrue          : ShowMainWindow,
              onFalse         : HideMainWindow,
              onHover         : null,
              onHoverOut      : null,
              onEnable        : null,
              onDisable       : null,
              visibleInScenes : ApplicationLauncher.AppScenes.ALWAYS,
              texture         : toolbar_button_texture);
    }
    // Make sure the state of the toolbar button remains consistent with the
    // state of the window.
    if (show_main_window_) {
      toolbar_button_.SetTrue(makeCall : false);
    } else {
      toolbar_button_.SetFalse(makeCall : false);
    }

    if (hide_all_gui_) {
      return;
    } else if (show_main_window_) {
      UnityEngine.GUI.skin = HighLogic.Skin;
      main_window_rectangle_.xMin = main_window_x_;
      main_window_rectangle_.yMin = main_window_y_;
      main_window_rectangle_ = UnityEngine.GUILayout.Window(
          id         : 1,
          screenRect : main_window_rectangle_,
          func       : DrawMainWindow,
          text       : "Principia",
          options    : UnityEngine.GUILayout.MinWidth(500));
      main_window_x_ = (int)main_window_rectangle_.xMin;
      main_window_y_ = (int)main_window_rectangle_.yMin;
    }
  }

  private void Update() {
    if (MapView.MapIsEnabled && map_renderer_ == null) {
      map_renderer_ = MapView.MapCamera.gameObject.AddComponent<MapRenderer>();
      map_renderer_.render_on_pre_cull = RenderTrajectories;
    }
    override_rsas_target_ = false;
    Vessel active_vessel = FlightGlobals.ActiveVessel;
    if (active_vessel != null &&
        !FlightGlobals.ActiveVessel.isEVA) {
      if (navball_ == null) {
        navball_ = (NavBall)FindObjectOfType(typeof(NavBall));
      }
      if (compass_navball_texture_ == null) {
        compass_navball_texture_ =
            navball_.navBall.renderer.material.mainTexture;
      }

      if (navball_changed_) {
        // Texture the ball.
        navball_changed_ = false;
        if (!fix_navball_in_plotting_frame_ || !PluginRunning()) {
          navball_.navBall.renderer.material.mainTexture =
              compass_navball_texture_;
        } else if (first_selected_celestial_ == second_selected_celestial_) {
          navball_.navBall.renderer.material.mainTexture =
              inertial_navball_texture_;
        } else {
          navball_.navBall.renderer.material.mainTexture =
              barycentric_navball_texture_;
        }
      }

      if (PluginRunning() && fix_navball_in_plotting_frame_) {
        // Orient the ball.
        navball_.navBall.rotation =
            (UnityEngine.QuaternionD)navball_.attitudeGymbal *  // sic.
                (UnityEngine.QuaternionD)NavballOrientation(
                    plugin_,
                    rendering_frame_,
                    (XYZ)Planetarium.fetch.Sun.position,
                    (XYZ)(Vector3d)active_vessel.ReferenceTransform.position);
        // TODO(egg): the navball should be independent from the frame of the
        // Frenet trihedron (seeing your body-centric velocity with a compass
        // navball like in stock makes sense, so does seeing your velocity in
        // any reference frame with the fixed stars navball), although the
        // Frenet trihedron should be in the same frame as the map view
        // trajectory.  Right now when in space the navball is always linked to
        // the frame of the Frenet trihedron and the trajectory.
        if (has_active_vessel_in_space() &&
            has_vessel(plugin_, active_vessel.id.ToString())) {
          // Orient the Frenet trihedron.
          Vector3d prograde =
              (Vector3d)VesselTangent(plugin_,
                                      active_vessel.id.ToString(),
                                      rendering_frame_);
          Vector3d radial =
              (Vector3d)VesselNormal(plugin_,
                                     active_vessel.id.ToString(),
                                     rendering_frame_);
          // Yes, the astrodynamicist's normal is the mathematician's binormal.
          // Don't ask.
          Vector3d normal =
              (Vector3d)VesselBinormal(plugin_,
                                       active_vessel.id.ToString(),
                                       rendering_frame_);

          navball_.progradeVector.transform.localPosition =
              (UnityEngine.QuaternionD)navball_.attitudeGymbal *
                  prograde * 0.05;
          navball_.radialInVector.transform.localPosition =
              (UnityEngine.QuaternionD)navball_.attitudeGymbal *
                  radial * 0.05;
          navball_.normalVector.transform.localPosition =
              (UnityEngine.QuaternionD)navball_.attitudeGymbal *
                  normal * 0.05;
          navball_.retrogradeVector.transform.localPosition =
              -navball_.progradeVector.transform.localPosition;
          navball_.radialOutVector.transform.localPosition =
              -navball_.radialInVector.transform.localPosition;
          navball_.antiNormalVector.transform.localPosition =
              -navball_.normalVector.transform.localPosition;
          // Make the autopilot target our Frenet trihedron.
          // TODO(egg): just the tangent for now.
          if (active_vessel.OnAutopilotUpdate.GetInvocationList()[0] !=
              (Delegate)(FlightInputCallback)OverrideRSASTarget) {
            Log.Info("Prepending RSAS override");
            active_vessel.OnAutopilotUpdate =
                (FlightInputCallback)Delegate.Combine(
                    new FlightInputCallback(OverrideRSASTarget),
                    active_vessel.OnAutopilotUpdate);
          }
          if (active_vessel.Autopilot.Enabled) {
            override_rsas_target_ = true;
            switch (active_vessel.Autopilot.Mode) {
              case VesselAutopilot.AutopilotMode.Prograde:
                rsas_target_ = prograde;
                break;
              case VesselAutopilot.AutopilotMode.Retrograde:
                rsas_target_ = -prograde;
                break;
              case VesselAutopilot.AutopilotMode.RadialIn:
                rsas_target_ = radial;
                break;
              case VesselAutopilot.AutopilotMode.RadialOut:
                rsas_target_ = -radial;
                break;
              case VesselAutopilot.AutopilotMode.Normal:
                rsas_target_ = normal;
                break;
              case VesselAutopilot.AutopilotMode.Antinormal:
                rsas_target_ = -normal;
                break;
              default:
                override_rsas_target_ = false;
                break;
            }
          }
        }
      }
    }
  }

  private void FixedUpdate() {
    if (PluginRunning()) {
      double universal_time = Planetarium.GetUniversalTime();
      double plugin_time = current_time(plugin_);
      if (plugin_time > universal_time) {
        // TODO(Egg): Make this resistant to bad floating points up to 2ULPs,
        // and make it fatal again.
        Log.Error("Closed Timelike Curve: " +
                  plugin_time + " > " + universal_time +
                  " plugin-universal=" + (plugin_time - universal_time));
        time_is_advancing_ = false;
        return;
      } else if (plugin_time == universal_time) {
        time_is_advancing_ = false;
        return;
      }
      time_is_advancing_ = true;
      if (has_inertial_physics_bubble_in_space()) {
        ApplyToVesselsInPhysicsBubble(AddToPhysicsBubble);
      }
      Vessel active_vessel = FlightGlobals.ActiveVessel;
      bool ready_to_draw_active_vessel_trajectory =
          draw_active_vessel_trajectory() &&
          has_vessel(plugin_, active_vessel.id.ToString());
      if (ready_to_draw_active_vessel_trajectory) {
        set_prediction_length_tolerance(
            plugin_,
            prediction_length_tolerances_[prediction_length_tolerance_index_]);
        // TODO(egg): make the speed tolerance independent.
        set_prediction_speed_tolerance(
            plugin_,
            prediction_length_tolerances_[prediction_length_tolerance_index_]);
        set_prediction_length(plugin_,
                              prediction_lengths_[prediction_length_index_]);
      }
      AdvanceTime(plugin_, universal_time, Planetarium.InverseRotAngle);
      if (ready_to_draw_active_vessel_trajectory) {
        UpdatePrediction(plugin_, active_vessel.id.ToString());
      }
      ForgetAllHistoriesBefore(
          plugin_,
          universal_time - history_lengths_[history_length_index_]);
      ApplyToBodyTree(body => UpdateBody(body, universal_time));
      ApplyToVesselsOnRailsOrInInertialPhysicsBubbleInSpace(
          vessel => UpdateVessel(vessel, universal_time));
      if (!PhysicsBubbleIsEmpty(plugin_)) {
        Vector3d displacement_offset =
            (Vector3d)BubbleDisplacementCorrection(
                          plugin_,
                          (XYZ)Planetarium.fetch.Sun.position);
        Vector3d velocity_offset =
            (Vector3d)BubbleVelocityCorrection(
                          plugin_,
                          active_vessel.orbit.referenceBody.flightGlobalsIndex);
        if (krakensbane_ == null) {
          krakensbane_ = (Krakensbane)FindObjectOfType(typeof(Krakensbane));
        }
        krakensbane_.setOffset(displacement_offset);
        krakensbane_.FrameVel += velocity_offset;
      }
    }
  }

  private void OnDisable() {
    Log.Info("principia.ksp_plugin_adapter.PrincipiaPluginAdapter.OnDisable()");
    if (toolbar_button_ != null) {
      ApplicationLauncher.Instance.RemoveModApplication(toolbar_button_);
    }
    Cleanup();
  }

  #endregion

  private void RenderTrajectories() {
    if (!PluginRunning()) {
      return;
    }
    Vessel active_vessel = FlightGlobals.ActiveVessel;
    bool ready_to_draw_active_vessel_trajectory =
        draw_active_vessel_trajectory() &&
        has_vessel(plugin_, active_vessel.id.ToString()); 
    if (ready_to_draw_active_vessel_trajectory) {
      active_vessel.patchedConicRenderer.relativityMode =
          PatchRendering.RelativityMode.RELATIVE;
      if (active_vessel.orbitDriver.Renderer.drawMode !=
              OrbitRenderer.DrawMode.OFF ||
          active_vessel.orbitDriver.Renderer.drawIcons !=
              OrbitRenderer.DrawIcons.OBJ) {
        Log.Info("Removing orbit rendering for the active vessel");
        active_vessel.orbitDriver.Renderer.drawMode =
            OrbitRenderer.DrawMode.OFF;
        active_vessel.orbitDriver.Renderer.drawIcons =
            OrbitRenderer.DrawIcons.OBJ;
      }
      if (display_patched_conics_) {
        foreach (PatchRendering patch_rendering in
                  active_vessel.patchedConicRenderer.patchRenders) {
          patch_rendering.visible = true;
        }
        // For reasons that are unlikely to become clear again at this time,
        // I think the first element of |flightPlanRenders| should not be
        // visible.
        for (int i = 1;
              i < active_vessel.patchedConicRenderer.flightPlanRenders.Count;
              ++i) {
          active_vessel.patchedConicRenderer.flightPlanRenders[i].visible =
              true;
        }
      } else {
        foreach (PatchRendering patch_rendering in
                  active_vessel.patchedConicRenderer.patchRenders) {
          patch_rendering.visible = false;
        }
        foreach (PatchRendering patch_rendering in
                  active_vessel.patchedConicRenderer.flightPlanRenders) {
          patch_rendering.visible = false;
        }
      }
      if (rendered_trajectory_ == null || rendered_prediction_ == null) {
        ResetRenderedTrajectory();
      }
      IntPtr trajectory_iterator = IntPtr.Zero;
      trajectory_iterator = RenderedVesselTrajectory(
                                plugin_,
                                active_vessel.id.ToString(),
                                rendering_frame_,
                                (XYZ)Planetarium.fetch.Sun.position);
      RenderAndDeleteTrajectory(ref trajectory_iterator,
                                rendered_trajectory_);
      if (HasPrediction(plugin_, active_vessel.id.ToString())) {
        trajectory_iterator = RenderedPrediction(
                                  plugin_,
                                  active_vessel.id.ToString(),
                                  rendering_frame_,
                                  (XYZ)Planetarium.fetch.Sun.position);
        RenderAndDeleteTrajectory(ref trajectory_iterator,
                                  rendered_prediction_);
      }
    } else {
      DestroyRenderedTrajectory();
    }
  }

  private void RenderAndDeleteTrajectory(ref IntPtr trajectory_iterator,
                                         VectorLine vector_line) {
    int new_min_draw_index = 0;
    try {
      LineSegment segment;
      int index_in_line_points = vector_line.points3.Length -
          2 * NumberOfSegments(trajectory_iterator);
      // If the |VectorLine| is too big, make sure we're not keeping garbage.
      for (int i = vector_line.minDrawIndex; i < index_in_line_points; ++i) {
        vector_line.points3[i] = UnityEngine.Vector3.zero;
      }
      while (index_in_line_points < 0) {
        FetchAndIncrement(trajectory_iterator);
        index_in_line_points += 2;
      }
      new_min_draw_index = index_in_line_points;
      vector_line.minDrawIndex = Math.Min(vector_line.minDrawIndex,
                                          new_min_draw_index);
      while (!AtEnd(trajectory_iterator)) {
        segment = FetchAndIncrement(trajectory_iterator);
        // TODO(egg): should we do the |LocalToScaledSpace| conversion in
        // native code?
        // TODO(egg): could we directly assign to
        // |vector_line.points3| from C++ using unsafe code and
        // something like the following?
        // |fixed (UnityEngine.Vector3* pts = vector_line.points3)|
        vector_line.points3[index_in_line_points++] =
            ScaledSpace.LocalToScaledSpace((Vector3d)segment.begin);
        vector_line.points3[index_in_line_points++] =
            ScaledSpace.LocalToScaledSpace((Vector3d)segment.end);
      }
    } finally {
      DeleteLineAndIterator(ref trajectory_iterator);
    }
    if (MapView.Draw3DLines && !force_2d_trajectories_) {
      Vector.DrawLine3D(vector_line);
    } else {
      Vector.DrawLine(vector_line);
    }
    vector_line.minDrawIndex = new_min_draw_index;
  }

  private void ResetRenderedTrajectory() {
    DestroyRenderedTrajectory();
    rendered_trajectory_ = new VectorLine(
        lineName     : "rendered_trajectory_",
        linePoints   : new UnityEngine.Vector3[kMaxVectorLinePoints],
        lineMaterial : MapView.OrbitLinesMaterial,
        color        : XKCDColors.AcidGreen,
        width        : 5,
        lineType     : LineType.Discrete);
    rendered_trajectory_.vectorObject.transform.parent =
        ScaledSpace.Instance.transform;
    rendered_trajectory_.vectorObject.renderer.castShadows = false;
    rendered_trajectory_.vectorObject.renderer.receiveShadows = false;
    rendered_trajectory_.layer = 31;
    rendered_prediction_ = new VectorLine(
        lineName     : "rendered_prediction_",
        linePoints   : new UnityEngine.Vector3[kMaxVectorLinePoints],
        lineMaterial : MapView.OrbitLinesMaterial,
        color        : XKCDColors.Fuchsia,
        width        : 5,
        lineType     : LineType.Discrete);
    rendered_prediction_.vectorObject.transform.parent =
        ScaledSpace.Instance.transform;
    rendered_prediction_.vectorObject.renderer.castShadows = false;
    rendered_prediction_.vectorObject.renderer.receiveShadows = false;
    rendered_prediction_.layer = 31;
  }

  private void DestroyRenderedTrajectory() {
    if (rendered_trajectory_ != null) {
      Vector.DestroyLine(ref rendered_trajectory_);
    }
    if (rendered_prediction_ != null) {
      Vector.DestroyLine(ref rendered_prediction_);
    }
  }

  private void Cleanup() {
    UnityEngine.Object.Destroy(map_renderer_);
    map_renderer_ = null;
    DeletePlugin(ref plugin_);
    DeleteRenderingFrame(ref rendering_frame_);
    DestroyRenderedTrajectory();
    navball_changed_ = true;
  }

  private void ShowGUI() {
    hide_all_gui_ = false;
  }

  private void HideGUI() {
    hide_all_gui_ = true;
  }

  private void ShowMainWindow() {
    show_main_window_ = true;
  }

  private void HideMainWindow() {
    show_main_window_ = false;
  }

  private void DrawMainWindow(int window_id) {
    UnityEngine.GUILayout.BeginVertical();
    String plugin_state;
    if (!PluginRunning()) {
      plugin_state = "not started";
    } else if (!time_is_advancing_) {
      plugin_state = "holding";
    } else if (PhysicsBubbleIsEmpty(plugin_)) {
      plugin_state = "running";
    } else {
      plugin_state = "managing physics bubble";
    }
    UnityEngine.GUILayout.TextArea(text : "Plugin is " + plugin_state);
    String last_reset_information;
    if (!PluginRunning()) {
      last_reset_information = "";
    } else {
      String plugin_source = "";
      switch (plugin_source_) {
        case (PluginSource.SAVED_STATE):
           plugin_source = "a saved state";
          break;
        case (PluginSource.ORBITAL_ELEMENTS):
           plugin_source = "KSP orbital elements";
          break;
        case (PluginSource.CARTESIAN_CONFIG):
           plugin_source = "a cartesian configuration file";
          break;
      }
      last_reset_information =
          "Plugin was constructed at " +
          plugin_construction_.ToUniversalTime().ToString("O") + " from " +
          plugin_source;
    }
    UnityEngine.GUILayout.TextArea(last_reset_information);
    bool changed_history_length = false;
    Selector(history_lengths_,
             ref history_length_index_,
             "Max history length",
             ref changed_history_length,
             "{0:0.00e00} s");
    force_2d_trajectories_ =
        UnityEngine.GUILayout.Toggle(force_2d_trajectories_,
                                     "Force 2D trajectories");
    ToggleableSection(name   : "Reference Frame Selection",
                      show   : ref show_reference_frame_selection_,
                      render : ReferenceFrameSelection);
    ToggleableSection(name   : "Prediction Settings",
                      show   : ref show_prediction_settings_,
                      render : PredictionSettings);
    ToggleableSection(name   : "KSP features",
                      show   : ref show_ksp_features_,
                      render : KSPFeatures);
    ToggleableSection(name   : "Logging Settings",
                      show   : ref show_logging_settings_,
                      render : LoggingSettings);
    ToggleableSection(name   : "Reset Principia",
                      show   : ref show_reset_button_,
                      render : ResetButton);
#if CRASH_BUTTON
    ToggleableSection(name   : "CRASH",
                      show   : ref show_crash_options_,
                      render : CrashOptions);
#endif
    UnityEngine.GUILayout.EndVertical();
    UnityEngine.GUI.DragWindow(
        position : new UnityEngine.Rect(left : 0f, top : 0f, width : 10000f,
                                        height : 10000f));
  }

  delegate void GUIRenderer();

  private void ToggleableSection(String name,
                                 ref bool show,
                                 GUIRenderer render) {
    String toggle = show ? "↑ " + name + " ↑"
                         : "↓ " + name + " ↓";
    if (UnityEngine.GUILayout.Button(toggle)) {
      show = !show;
      if (!show) {
        ShrinkMainWindow();
      }
    }
    if (show) {
      render();
    }
  }

#if CRASH_BUTTON
  private void CrashOptions() {
    if (UnityEngine.GUILayout.Button(text : "CRASH ON MAP VIEW")) {
      first_selected_celestial_ = second_selected_celestial_;
      DeleteRenderingFrame(ref rendering_frame_);
      rendering_frame_ = NewBarycentricRotatingRenderingFrame(
                             plugin_,
                             first_selected_celestial_,
                             second_selected_celestial_);
    }
    if (UnityEngine.GUILayout.Button(text : "CRASH NOW")) {
      Log.Fatal("You asked for it!");
    }
  }
#endif

  private void ReferenceFrameSelection() {
    bool was_fixing_navball_in_plotting_frame =
        fix_navball_in_plotting_frame_;
    fix_navball_in_plotting_frame_ = 
        UnityEngine.GUILayout.Toggle(
            value : fix_navball_in_plotting_frame_,
            text  : "Fix navball in plotting frame");
    if (PluginRunning() &&
        was_fixing_navball_in_plotting_frame !=
        fix_navball_in_plotting_frame_) {
      navball_changed_ = true;
      reset_rsas_target_ = true;
    }
    bool barycentric_rotating =
        first_selected_celestial_ != second_selected_celestial_;
    String reference_frame_description =
        "The trajectory of the active vessel is plotted in ";
    if (barycentric_rotating) {
      reference_frame_description +=
          "the reference frame fixing the barycentre of " +
          FlightGlobals.Bodies[first_selected_celestial_].theName + " and " +
          FlightGlobals.Bodies[second_selected_celestial_].theName + ", " +
          "the line through them, and the plane in which they move about the " +
          "barycentre.";
    } else {
      reference_frame_description +=
          "the nonrotating reference frame fixing the centre of " +
          FlightGlobals.Bodies[first_selected_celestial_].theName + ".";
    }
    UnityEngine.GUILayout.TextArea(
        text    : reference_frame_description,
        options : UnityEngine.GUILayout.Height(100));
    foreach (CelestialBody celestial in FlightGlobals.Bodies) {
      bool changed_rendering = false;
      UnityEngine.GUILayout.BeginHorizontal();
      if (UnityEngine.GUILayout.Toggle(
              value : first_selected_celestial_ == celestial.flightGlobalsIndex,
              text  : "") &&
          first_selected_celestial_ != celestial.flightGlobalsIndex) {
        first_selected_celestial_ = celestial.flightGlobalsIndex;
        changed_rendering = true;
      }
      if (UnityEngine.GUILayout.Toggle(
              value : second_selected_celestial_ ==
                          celestial.flightGlobalsIndex,
              text  : celestial.name) &&
          second_selected_celestial_ != celestial.flightGlobalsIndex) {
        second_selected_celestial_ = celestial.flightGlobalsIndex;
        changed_rendering = true;
      }
      UnityEngine.GUILayout.EndHorizontal();
      if (changed_rendering && PluginRunning()) {
        UpdateRenderingFrame();
      }
    }
  }

  private void Selector(
      double[] array,
      ref int index,
      String label,
      ref bool changed,
      String format) {
    UnityEngine.GUILayout.BeginHorizontal();
    UnityEngine.GUILayout.Label(text    : label + ":",
                                options : UnityEngine.GUILayout.Width(150));
    if (UnityEngine.GUILayout.Button(
            text    : index == 0 ? "min" : "-",
            options : UnityEngine.GUILayout.Width(50)) &&
        index != 0) {
      --index;
      changed = true;
    }
    UnityEngine.TextAnchor old_alignment =
        UnityEngine.GUI.skin.textArea.alignment;
    UnityEngine.GUI.skin.textArea.alignment =
        UnityEngine.TextAnchor.MiddleRight;
    // Unity/Mono is screwing with the current culture, let's get unambiguous
    // conventions from a copy of the invariant culture.
    CultureInfo culture = new CultureInfo("");
    culture.NumberFormat.NumberGroupSeparator = "'";
    culture.NumberFormat.PositiveInfinitySymbol = "+∞";
    UnityEngine.GUILayout.TextArea(
        text    : String.Format(culture, format, array[index]),
        options : UnityEngine.GUILayout.Width(75));
    UnityEngine.GUI.skin.textArea.alignment = old_alignment;
    if (UnityEngine.GUILayout.Button(
            text    : index == array.Length - 1 ? "max" : "+",
            options : UnityEngine.GUILayout.Width(50)) &&
        index != array.Length - 1) {
      ++index;
      changed = true;
    }
    UnityEngine.GUILayout.EndHorizontal();
  }

  // NOTE(egg): Dummy UI elements for testing purposes, rendered in an
  // irrelevant part of the UI.
  FlightPlanner test_flight_planner_;
  ReferenceFrameSelector test_reference_frame_selector_;

  private void PredictionSettings() {
    bool changed_settings = false;
    Selector(prediction_length_tolerances_,
             ref prediction_length_tolerance_index_,
             "Tolerance",
             ref changed_settings,
             "{0:0.0e0} m");
    Selector(prediction_lengths_,
             ref prediction_length_index_,
             "Length",
             ref changed_settings,
             "{0:0.00e0} s");
    if (test_flight_planner_ == null) {
      test_flight_planner_ = new FlightPlanner();
    }
    test_flight_planner_.Render();
    if (test_reference_frame_selector_ == null) {
      test_reference_frame_selector_ = new ReferenceFrameSelector();
    }
    test_reference_frame_selector_.Render();
  }

  private void KSPFeatures() {
    display_patched_conics_ =
        UnityEngine.GUILayout.Toggle(value : display_patched_conics_,
                                     text  : "Display patched conics");
    Sun.Instance.sunFlare.enabled =
        UnityEngine.GUILayout.Toggle(value : Sun.Instance.sunFlare.enabled,
                                     text  : "Enable Sun lens flare");
  }

  private void LoggingSettings() {
    UnityEngine.GUILayout.BeginHorizontal();
    UnityEngine.GUILayout.Label(text : "Verbose level:");
    if (UnityEngine.GUILayout.Button(
            text    : "←",
            options : UnityEngine.GUILayout.Width(50))) {
      Log.SetVerboseLogging(Math.Max(verbose_logging_ - 1, 0));
      verbose_logging_ = Log.GetVerboseLogging();
    }
    UnityEngine.GUILayout.TextArea(
        text    : Log.GetVerboseLogging().ToString(),
        options : UnityEngine.GUILayout.Width(50));
    if (UnityEngine.GUILayout.Button(
            text    : "→",
            options : UnityEngine.GUILayout.Width(50))) {
      Log.SetVerboseLogging(Math.Min(verbose_logging_ + 1, 4));
      verbose_logging_ = Log.GetVerboseLogging();
    }
    UnityEngine.GUILayout.EndHorizontal();
    int column_width = 75;
    UnityEngine.GUILayout.BeginHorizontal();
    UnityEngine.GUILayout.Space(column_width);
    UnityEngine.GUILayout.Label(
        text    : "Log",
        options : UnityEngine.GUILayout.Width(column_width));
    UnityEngine.GUILayout.Label(
        text    : "stderr",
        options : UnityEngine.GUILayout.Width(column_width));
    UnityEngine.GUILayout.Label(
        text    : "Flush",
        options : UnityEngine.GUILayout.Width(column_width));
    UnityEngine.GUILayout.EndHorizontal();
    UnityEngine.GUILayout.BeginHorizontal();
    UnityEngine.GUILayout.Space(column_width);
    if (UnityEngine.GUILayout.Button(
            text    : "↑",
            options : UnityEngine.GUILayout.Width(column_width))) {
      Log.SetSuppressedLogging(Math.Max(suppressed_logging_ - 1, 0));
      suppressed_logging_ = Log.GetSuppressedLogging();
    }
    if (UnityEngine.GUILayout.Button(
            text    : "↑",
            options : UnityEngine.GUILayout.Width(column_width))) {
      Log.SetStderrLogging(Math.Max(stderr_logging_ - 1, 0));
      stderr_logging_ = Log.GetStderrLogging();
    }
    if (UnityEngine.GUILayout.Button(
            text    : "↑",
            options : UnityEngine.GUILayout.Width(column_width))) {
      Log.SetBufferedLogging(Math.Max(buffered_logging_ - 1, -1));
      buffered_logging_ = Log.GetBufferedLogging();
    }
    UnityEngine.GUILayout.EndHorizontal();
    for (int severity = 0; severity <= 3; ++severity) {
      UnityEngine.GUILayout.BeginHorizontal();
      UnityEngine.GUILayout.Label(
          text    : Log.kSeverityNames[severity],
          options : UnityEngine.GUILayout.Width(column_width));
      UnityEngine.GUILayout.Toggle(
          value   : severity >= Log.GetSuppressedLogging(),
          text    : "",
          options : UnityEngine.GUILayout.Width(column_width));
      UnityEngine.GUILayout.Toggle(
          value   : severity >= Log.GetStderrLogging(),
          text    : "",
          options : UnityEngine.GUILayout.Width(column_width));
      UnityEngine.GUILayout.Toggle(
          value   : severity > Log.GetBufferedLogging(),
          text    : "",
          options : UnityEngine.GUILayout.Width(column_width));
      UnityEngine.GUILayout.EndHorizontal();
    }
    UnityEngine.GUILayout.BeginHorizontal();
    UnityEngine.GUILayout.Space(column_width);
    if (UnityEngine.GUILayout.Button(
            text    : "↓",
            options : UnityEngine.GUILayout.Width(column_width))) {
      Log.SetSuppressedLogging(Math.Min(suppressed_logging_ + 1, 3));
      suppressed_logging_ = Log.GetSuppressedLogging();
    }
    if (UnityEngine.GUILayout.Button(
            text    : "↓",
            options : UnityEngine.GUILayout.Width(column_width))) {
      Log.SetStderrLogging(Math.Min(stderr_logging_ + 1, 3));
      stderr_logging_ = Log.GetStderrLogging();
    }
    if (UnityEngine.GUILayout.Button(
            text    : "↓",
            options : UnityEngine.GUILayout.Width(column_width))) {
      Log.SetBufferedLogging(Math.Min(buffered_logging_ + 1, 3));
      buffered_logging_ = Log.GetBufferedLogging();
    }
    UnityEngine.GUILayout.EndHorizontal();
  }

  private void ResetButton() {
    if (PluginRunning()) {
      if (UnityEngine.GUILayout.Button(text : "Force Stop")) {
        Cleanup();
      }
    } else {
      if (UnityEngine.GUILayout.Button(text : "Force Start")) {
        ResetPlugin();
      }
    }
  }

  private void ShrinkMainWindow() {
    main_window_rectangle_.height = 0.0f;
    main_window_rectangle_.width = 0.0f;
  }

  private void UpdateRenderingFrame() {
    if (fix_navball_in_plotting_frame_) {
      navball_changed_ = true;
      reset_rsas_target_ = true;
    }
    DeleteRenderingFrame(ref rendering_frame_);
    if (first_selected_celestial_ == second_selected_celestial_) {
      rendering_frame_ = NewBodyCentredNonRotatingRenderingFrame(
                             plugin_,
                             first_selected_celestial_);
    } else {
      rendering_frame_ = NewBarycentricRotatingRenderingFrame(
                             plugin_,
                             first_selected_celestial_,
                             second_selected_celestial_);
    }
  }

  private void ResetPlugin() {
    Cleanup();
    SetRotatingFrameThresholds();
    ResetRenderedTrajectory();
    plugin_construction_ = DateTime.Now;
    if (GameDatabase.Instance.GetConfigs(kPrincipiaInitialState).Length > 0) {
      plugin_source_ = PluginSource.CARTESIAN_CONFIG;
      if (GameDatabase.Instance.GetConfigs(
              kPrincipiaGravityModel).Length == 0) {
        Log.Fatal("missing gravity models");
      }
      if (GameDatabase.Instance.GetConfigs(kPrincipiaInitialState).Length > 1 ||
          GameDatabase.Instance.GetConfigs(
              kPrincipiaGravityModel).Length > 1) {
        Log.Fatal("too many configs");
      }
      try {
        ConfigNode initial_states =
            GameDatabase.Instance.GetConfigs(kPrincipiaInitialState)[0].config;
        ConfigNode gravity_models =
            GameDatabase.Instance.GetConfigs(kPrincipiaGravityModel)[0].config;
        plugin_ = NewPlugin(double.Parse(initial_states.GetValue("epoch")),
                            Planetarium.InverseRotAngle);
        var name_to_initial_state = new Dictionary<String, ConfigNode>();
        var name_to_gravity_model = new Dictionary<String, ConfigNode>();
        foreach (ConfigNode node in initial_states.GetNodes("body")) {
          name_to_initial_state.Add(node.GetValue("name"), node);
        }
        foreach (ConfigNode node in gravity_models.GetNodes("body")) {
          name_to_gravity_model.Add(node.GetValue("name"), node);
        }
        ConfigNode sun_gravity_model =
            name_to_gravity_model[Planetarium.fetch.Sun.name];
        ConfigNode sun_initial_state =
            name_to_initial_state[Planetarium.fetch.Sun.name];
        DirectlyInsertCelestial(
            plugin: plugin_,
            celestial_index: Planetarium.fetch.Sun.flightGlobalsIndex,
            parent_index: IntPtr.Zero,
            gravitational_parameter:
                sun_gravity_model.GetValue("gravitational_parameter"),
            axis_right_ascension:
                sun_gravity_model.HasValue("axis_right_ascension") ?
                sun_gravity_model.GetValue("axis_right_ascension") : null,
            axis_declination:
                sun_gravity_model.HasValue("axis_declination") ?
                sun_gravity_model.GetValue("axis_declination") : null,
            j2: sun_gravity_model.HasValue("j2") ?
                sun_gravity_model.GetValue("j2") : null,
            reference_radius:
                sun_gravity_model.HasValue("reference_radius") ?
                sun_gravity_model.GetValue("reference_radius") : null,
            x: sun_initial_state.GetValue("x"),
            y: sun_initial_state.GetValue("y"),
            z: sun_initial_state.GetValue("z"),
            vx: sun_initial_state.GetValue("vx"),
            vy: sun_initial_state.GetValue("vy"),
            vz: sun_initial_state.GetValue("vz"));
        BodyProcessor insert_body = body => {
          Log.Info("Inserting " + body.name + "...");
          ConfigNode gravity_model = name_to_gravity_model[body.name];
          ConfigNode initial_state = name_to_initial_state[body.name];
          int parent_index = body.orbit.referenceBody.flightGlobalsIndex;
          DirectlyInsertCelestial(
              plugin: plugin_,
              celestial_index: body.flightGlobalsIndex,
              parent_index: ref parent_index,
              gravitational_parameter:
                  gravity_model.GetValue("gravitational_parameter"),
              axis_right_ascension:
                  gravity_model.HasValue("axis_right_ascension") ?
                  gravity_model.GetValue("axis_right_ascension") : null,
              axis_declination:
                  gravity_model.HasValue("axis_declination") ?
                  gravity_model.GetValue("axis_declination") : null,
              j2: gravity_model.HasValue("j2") ?
                  gravity_model.GetValue("j2") : null,
              reference_radius:
                  gravity_model.HasValue("reference_radius") ?
                  gravity_model.GetValue("reference_radius") : null,
              x: initial_state.GetValue("x"),
              y: initial_state.GetValue("y"),
              z: initial_state.GetValue("z"),
              vx: initial_state.GetValue("vx"),
              vy: initial_state.GetValue("vy"),
              vz: initial_state.GetValue("vz"));
        };
        ApplyToBodyTree(insert_body);
        EndInitialization(plugin_);
        AdvanceTime(plugin_,
                    Planetarium.GetUniversalTime(),
                    Planetarium.InverseRotAngle);
      } catch (Exception e) {
        Log.Fatal("Exception while reading initial state: " + e.ToString());
      }
    } else {
      plugin_source_ = PluginSource.ORBITAL_ELEMENTS;
      plugin_ = NewPlugin(Planetarium.GetUniversalTime(),
                          Planetarium.InverseRotAngle);
      InsertSun(plugin_,
                Planetarium.fetch.Sun.flightGlobalsIndex,
                Planetarium.fetch.Sun.gravParameter);
      BodyProcessor insert_body = body => {
        Log.Info("Inserting " + body.name + "...");
        InsertCelestial(plugin_,
                        body.flightGlobalsIndex,
                        body.gravParameter,
                        body.orbit.referenceBody.flightGlobalsIndex,
                        new QP{q = (XYZ)body.orbit.pos,
                               p = (XYZ)body.orbit.vel});
      };
      ApplyToBodyTree(insert_body);
      EndInitialization(plugin_);
    }
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
                             new QP{q = (XYZ)vessel.orbit.pos,
                                    p = (XYZ)vessel.orbit.vel});
      }
    };
    ApplyToVesselsOnRailsOrInInertialPhysicsBubbleInSpace(insert_vessel);
  }

  private void SetRotatingFrameThresholds() {
    ApplyToBodyTree(body => body.inverseRotThresholdAltitude =
                                (float)Math.Max(body.timeWarpAltitudeLimits[1],
                                                body.atmosphereDepth));
  }

}

}  // namespace ksp_plugin_adapter
}  // namespace principia
