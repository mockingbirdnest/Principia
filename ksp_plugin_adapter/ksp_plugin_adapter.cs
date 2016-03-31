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
public partial class PrincipiaPluginAdapter
    : ScenarioModule,
      WindowRenderer.ManagerInterface {

  private const String kPrincipiaKey = "serialized_plugin";
  private const String kPrincipiaInitialState = "principia_initial_state";
  private const String kPrincipiaGravityModel = "principia_gravity_model";
  private const double kΔt = 10;

  private KSP.UI.Screens.ApplicationLauncherButton toolbar_button_;
  private bool hide_all_gui_ = false;

  // "Persistant" is a KSP typo.
  [KSPField(isPersistant = true)]
  private bool show_main_window_ = true;
  [KSPField(isPersistant = true)]
  private int main_window_x_ = UnityEngine.Screen.width / 2;
  [KSPField(isPersistant = true)]
  private int main_window_y_ = UnityEngine.Screen.height / 3;
  private UnityEngine.Rect main_window_rectangle_;

  internal Controlled<ReferenceFrameSelector> plotting_frame_selector_;
  private Controlled<FlightPlanner> flight_planner_;
  private MapNodePool map_node_pool_;

  private IntPtr plugin_ = IntPtr.Zero;

  private bool display_patched_conics_ = false;
  [KSPField(isPersistant = true)]
  private bool fix_navball_in_plotting_frame_ = true;

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

  [KSPField(isPersistant = true)]
  private bool must_record_journal_ = false;
  [KSPField(isPersistant = true)]
  private bool must_record_verbose_journal_ = false;
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
  private KSP.UI.Screens.Flight.NavBall navball_;
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

  private static Dictionary<CelestialBody, Orbit> unmodified_orbits_;

  private String bad_installation_popup_;

  public event Action render_windows;

  PrincipiaPluginAdapter() {
    // We create this directory here so we do not need to worry about cross-
    // platform problems in C++.
    System.IO.Directory.CreateDirectory("glog/Principia");
    string load_error = Loader.LoadPrincipiaDllAndInitGoogleLogging();
    if (load_error != null) {
      bad_installation_popup_ =
          "The Principia DLL failed to load.\n" + load_error;
      UnityEngine.Debug.LogError(bad_installation_popup_);
    }
    map_node_pool_ = new MapNodePool();
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
    plugin_.UpdateCelestialHierarchy(
        body.flightGlobalsIndex,
        body.orbit.referenceBody.flightGlobalsIndex);
    QP from_parent = plugin_.CelestialFromParent(body.flightGlobalsIndex);
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
    bool inserted = plugin_.InsertOrKeepVessel(
        vessel.id.ToString(),
        vessel.orbit.referenceBody.flightGlobalsIndex);
    if (inserted) {
      plugin_.SetVesselStateOffset(
          vessel_guid : vessel.id.ToString(),
          from_parent : new QP{q = (XYZ)vessel.orbit.pos,
                               p = (XYZ)vessel.orbit.vel});
    }
    QP from_parent = plugin_.VesselFromParent(vessel.id.ToString());
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
             mass_in_tonnes =
                 (double)part.mass + (double)part.GetResourceMass(),
             gravitational_acceleration_to_be_applied_by_ksp = (XYZ)gravity,
             id = part.flightID}).ToArray();
    if (parts.Count() > 0) {
      bool inserted = plugin_.InsertOrKeepVessel(
          vessel.id.ToString(),
          vessel.orbit.referenceBody.flightGlobalsIndex);
      if (inserted) {
        // NOTE(egg): this is only used when a (plugin-managed) physics bubble
        // appears with a new vessel (e.g. when exiting the atmosphere).
        // TODO(egg): these degrees of freedom are off by one Δt and we don't
        // compensate for the pos/vel synchronization bug.
        plugin_.SetVesselStateOffset(
            vessel_guid : vessel.id.ToString(),
            from_parent : new QP{q = (XYZ)vessel.orbit.pos,
                                 p = (XYZ)vessel.orbit.vel});
      }
      plugin_.AddVesselToNextPhysicsBubble(vessel_guid : vessel.id.ToString(),
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
    if (bad_installation_popup_ != null) {
      return;
    }
    // While we're here, we might as well log.
    Log.Info("principia.ksp_plugin_adapter.PrincipiaPluginAdapter.OnAwake()");

    LoadTextureIfExists(out compass_navball_texture_, "navball_compass.png");
    LoadTextureOrDie(out inertial_navball_texture_, "navball_inertial.png");
    LoadTextureOrDie(out barycentric_navball_texture_,
                     "navball_barycentric.png");

    if (unmodified_orbits_ == null) {
      unmodified_orbits_ = new Dictionary<CelestialBody, Orbit>();
      foreach (CelestialBody celestial in
               FlightGlobals.Bodies.Where(c => c.orbit != null)) {
        unmodified_orbits_.Add(
            celestial,
            new Orbit(inc  : celestial.orbit.inclination,
                      e    : celestial.orbit.eccentricity,
                      sma  : celestial.orbit.semiMajorAxis,
                      lan  : celestial.orbit.LAN,
                      w    : celestial.orbit.argumentOfPeriapsis,
                      mEp  : celestial.orbit.meanAnomalyAtEpoch,
                      t    : celestial.orbit.epoch,
                      body : celestial.orbit.referenceBody));
      }
      Fix631();
    }

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
          serialization = plugin_.SerializePlugin(ref serializer);
          if (serialization == IntPtr.Zero) {
            break;
          }
          node.AddValue(kPrincipiaKey, Marshal.PtrToStringAnsi(serialization));
        } finally {
          Interface.DeletePluginSerialization(ref serialization);
        }
      }
    }
  }

  public override void OnLoad(ConfigNode node) {
    base.OnLoad(node);
    if (must_record_journal_) {
      Log.ActivateRecorder(true, must_record_verbose_journal_);
    }
    if (node.HasValue(kPrincipiaKey)) {
      Cleanup();
      SetRotatingFrameThresholds();
      RemoveBuggyTidalLocking();
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
        Interface.DeserializePlugin(serialization,
                                    serialization.Length,
                                    ref deserializer,
                                    ref plugin_);
      }
      Interface.DeserializePlugin("", 0, ref deserializer, ref plugin_);

      plotting_frame_selector_.reset(
          new ReferenceFrameSelector(this, 
                                     plugin_,
                                     UpdateRenderingFrame,
                                     "Plotting frame"));
      flight_planner_.reset(new FlightPlanner(this, plugin_));

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
    if (bad_installation_popup_ != null) {
      UnityEngine.Debug.LogError("Spawning: " + bad_installation_popup_);
      // No-one seems to understand what |anchorMin| and |anchorMax| do at this
      // time.
      PopupDialog.SpawnPopupDialog(
          anchorMin           : default(UnityEngine.Vector2),
          anchorMax           : default(UnityEngine.Vector2),
          title               : "Principia",
          message             : bad_installation_popup_,
          buttonMessage       : "OK",
          persistAcrossScenes : true,
          skin                : null,
          isModal             : true);
      bad_installation_popup_ = null;
      return;
    }
    if (KSP.UI.Screens.ApplicationLauncher.Ready && toolbar_button_ == null) {
      UnityEngine.Texture toolbar_button_texture;
      LoadTextureOrDie(out toolbar_button_texture, "toolbar_button.png");
      toolbar_button_ =
          KSP.UI.Screens.ApplicationLauncher.Instance.AddModApplication(
              onTrue          : ShowMainWindow,
              onFalse         : HideMainWindow,
              onHover         : null,
              onHoverOut      : null,
              onEnable        : null,
              onDisable       : null,
              visibleInScenes : KSP.UI.Screens.ApplicationLauncher.AppScenes.
                                    ALWAYS,
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
      UnityEngine.GUI.skin = null;
      main_window_rectangle_.xMin = main_window_x_;
      main_window_rectangle_.yMin = main_window_y_;
      main_window_rectangle_ = UnityEngine.GUILayout.Window(
          id         : this.GetHashCode(),
          screenRect : main_window_rectangle_,
          func       : DrawMainWindow,
          text       : "Principia",
          options    : UnityEngine.GUILayout.MinWidth(500));
      main_window_x_ = (int)main_window_rectangle_.xMin;
      main_window_y_ = (int)main_window_rectangle_.yMin;

      render_windows();
    }
  }

  private void Update() {
    if (MapView.MapIsEnabled && map_renderer_ == null) {
      map_renderer_ =
          PlanetariumCamera.Camera.gameObject.AddComponent<MapRenderer>();
      map_renderer_.post_render = RenderTrajectories;
    }
    override_rsas_target_ = false;
    Vessel active_vessel = FlightGlobals.ActiveVessel;
    if (active_vessel != null &&
        !FlightGlobals.ActiveVessel.isEVA) {
      if (navball_ == null) {
        navball_ = (KSP.UI.Screens.Flight.NavBall)FindObjectOfType(
                       typeof(KSP.UI.Screens.Flight.NavBall));
      }
      var navball_material =
          navball_.navBall.GetComponent<UnityEngine.Renderer>().material;

      if (compass_navball_texture_ == null) {
        compass_navball_texture_ = navball_material.GetTexture("_MainTexture");
      }

      Action<UnityEngine.Texture> set_navball_texture = (texture) =>
          navball_material.SetTexture("_MainTexture", texture);

      if (navball_changed_) {
        // Texture the ball.
        navball_changed_ = false;
        // TODO(egg): switch over all frame types and have more navball textures
        // when more frames are available.
        if (!fix_navball_in_plotting_frame_ || !PluginRunning()) {
          set_navball_texture(compass_navball_texture_);
        } else if (plotting_frame_selector_.get().frame_type ==
                   ReferenceFrameSelector.FrameType.BODY_CENTRED_NON_ROTATING) {
          set_navball_texture(inertial_navball_texture_);
        } else {
          set_navball_texture(barycentric_navball_texture_);
        }
      }

      if (PluginRunning() && fix_navball_in_plotting_frame_) {
        // Orient the ball.
        navball_.navBall.rotation =
            (UnityEngine.QuaternionD)navball_.attitudeGymbal *  // sic.
                (UnityEngine.QuaternionD)plugin_.NavballOrientation(
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
          plugin_.HasVessel(active_vessel.id.ToString())) {
          // Orient the Frenet trihedron.
          Vector3d prograde =
              (Vector3d)plugin_.VesselTangent(active_vessel.id.ToString());
          Vector3d radial =
              (Vector3d)plugin_.VesselNormal(active_vessel.id.ToString());
          // Yes, the astrodynamicist's normal is the mathematician's binormal.
          // Don't ask.
          Vector3d normal =
              (Vector3d)plugin_.VesselBinormal(active_vessel.id.ToString());

          double ball_radius =
              navball_.progradeVector.transform.localPosition.magnitude;

          // TODO(egg): deal with the transparency and disappearance of the
          // markers too, not just their position.
          navball_.progradeVector.transform.localPosition =
              (UnityEngine.QuaternionD)navball_.attitudeGymbal *
                  prograde * ball_radius;
          navball_.radialInVector.transform.localPosition =
              (UnityEngine.QuaternionD)navball_.attitudeGymbal *
                  radial * ball_radius;
          navball_.normalVector.transform.localPosition =
              (UnityEngine.QuaternionD)navball_.attitudeGymbal *
                  normal * ball_radius;
          navball_.retrogradeVector.transform.localPosition =
              -navball_.progradeVector.transform.localPosition;
          navball_.radialOutVector.transform.localPosition =
              -navball_.radialInVector.transform.localPosition;
          navball_.antiNormalVector.transform.localPosition =
              -navball_.normalVector.transform.localPosition;
          // Make the autopilot target our Frenet trihedron.
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
      double plugin_time = plugin_.CurrentTime();
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
          plugin_.HasVessel(active_vessel.id.ToString());
      if (ready_to_draw_active_vessel_trajectory) {
        plugin_.SetPredictionLengthTolerance(
            prediction_length_tolerances_[prediction_length_tolerance_index_]);
        // TODO(egg): make the speed tolerance independent.
        plugin_.SetPredictionSpeedTolerance(
            prediction_length_tolerances_[prediction_length_tolerance_index_]);
        plugin_.SetPredictionLength(
            prediction_lengths_[prediction_length_index_]);
      }
      plugin_.AdvanceTime(universal_time, Planetarium.InverseRotAngle);
      if (ready_to_draw_active_vessel_trajectory) {
        plugin_.UpdatePrediction(active_vessel.id.ToString());
      }
      plugin_.ForgetAllHistoriesBefore(
          universal_time - history_lengths_[history_length_index_]);
      ApplyToBodyTree(body => UpdateBody(body, universal_time));
      ApplyToVesselsOnRailsOrInInertialPhysicsBubbleInSpace(
          vessel => UpdateVessel(vessel, universal_time));
      if (!plugin_.PhysicsBubbleIsEmpty()) {
        Vector3d displacement_offset =
            (Vector3d)plugin_.BubbleDisplacementCorrection(
                          (XYZ)Planetarium.fetch.Sun.position);
        Vector3d velocity_offset =
            (Vector3d)plugin_.BubbleVelocityCorrection(
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
    if (bad_installation_popup_ != null) {
      return;
    }
    Log.Info("principia.ksp_plugin_adapter.PrincipiaPluginAdapter.OnDisable()");
    if (toolbar_button_ != null) {
      KSP.UI.Screens.ApplicationLauncher.Instance.RemoveModApplication(
          toolbar_button_);
    }
    Cleanup();
  }

  #endregion

  private void RemoveStockTrajectoriesIfNeeded(Vessel vessel) {
    vessel.patchedConicRenderer.relativityMode =
        PatchRendering.RelativityMode.RELATIVE;
    if (display_patched_conics_) {
      if (!vessel.patchedConicRenderer.enabled) {
        vessel.patchedConicRenderer.enabled = true;
        vessel.orbitTargeter.enabled = true;
      }
    } else {
      if (vessel.orbitDriver.Renderer.drawMode != OrbitRenderer.DrawMode.OFF ||
          vessel.orbitDriver.Renderer.drawIcons !=
              OrbitRenderer.DrawIcons.OBJ) {
        Log.Info("Removing patched conic rendering for the active vessel");
        vessel.orbitDriver.Renderer.drawMode = OrbitRenderer.DrawMode.OFF;
        vessel.orbitDriver.Renderer.drawIcons = OrbitRenderer.DrawIcons.OBJ;
      }
      vessel.patchedConicRenderer.enabled = false;
      vessel.orbitTargeter.enabled = false;
      foreach (PatchRendering patch_rendering in
               vessel.patchedConicRenderer.patchRenders) {
        patch_rendering.DestroyUINodes();
        patch_rendering.DestroyVector();
      }
      foreach (PatchRendering patch_rendering in
               vessel.patchedConicRenderer.flightPlanRenders) {
        patch_rendering.DestroyUINodes();
        patch_rendering.DestroyVector();
      }
      foreach (ManeuverNode node in vessel.patchedConicSolver.maneuverNodes) {
        node.DetachGizmo();
        if (node.scaledSpaceTarget) {
          MapView.MapCamera.RemoveTarget(node.scaledSpaceTarget);
          node.scaledSpaceTarget.Terminate();
        }
      }
    }
  }

  private void RenderTrajectories() {
    if (!PluginRunning()) {
      return;
    }
    Vessel active_vessel = FlightGlobals.ActiveVessel;
    if (active_vessel == null) {
      return;
    }
    string active_vessel_guid = active_vessel.id.ToString();
    bool ready_to_draw_active_vessel_trajectory =
        draw_active_vessel_trajectory() &&
        plugin_.HasVessel(active_vessel_guid); 
    if (ready_to_draw_active_vessel_trajectory) {
      RemoveStockTrajectoriesIfNeeded(active_vessel);

      XYZ sun_world_position = (XYZ)Planetarium.fetch.Sun.position;

      GLLines.Draw(() => {
        GLLines.RenderAndDeleteTrajectory(
            plugin_.RenderedVesselTrajectory(active_vessel_guid,
                                             sun_world_position),
            XKCDColors.AcidGreen);
        GLLines.RenderAndDeleteTrajectory(
            plugin_.RenderedPrediction(active_vessel_guid, sun_world_position),
            XKCDColors.Fuchsia);
        if (plugin_.FlightPlanExists(active_vessel_guid)) {
          RenderFlightPlanApsides(active_vessel_guid, sun_world_position);

          int number_of_segments =
              plugin_.FlightPlanNumberOfSegments(active_vessel_guid);
          for (int i = 0; i < number_of_segments; ++i) {
            bool is_burn = i % 2 == 1;
            GLLines.RenderAndDeleteTrajectory(
                plugin_.FlightPlanRenderedSegment(active_vessel_guid,
                                                  sun_world_position,
                                                  i),
                is_burn ? XKCDColors.OrangeRed : XKCDColors.RoyalBlue);
            if (is_burn) {
              Vector3d position_at_ignition =
                  (Vector3d)plugin_.FlightPlanRenderedSegmentEndpoints(
                                active_vessel_guid,
                                sun_world_position,
                                i).begin;
              int manoeuvre_index = i / 2;
              NavigationManoeuvre manoeuvre = plugin_.FlightPlanGetManoeuvre(
                                                  active_vessel_guid,
                                                  manoeuvre_index);
              double scale = (ScaledSpace.ScaledToLocalSpace(
                                  MapView.MapCamera.transform.position) -
                              position_at_ignition).magnitude * 0.015;
              Action<XYZ, UnityEngine.Color> add_vector =
                  (world_direction, colour) => {
                UnityEngine.GL.Color(colour);
                GLLines.AddSegment(
                    position_at_ignition,
                    position_at_ignition + scale * (Vector3d)world_direction,
                    hide_behind_bodies : false);
              };
              add_vector(manoeuvre.tangent, XKCDColors.NeonYellow);
              add_vector(manoeuvre.normal, XKCDColors.AquaBlue);
              add_vector(manoeuvre.binormal, XKCDColors.PurplePink);
            }
          }
        }
      });
      map_node_pool_.Update();
    } else {
      map_node_pool_.Clear();
    }
  }

  private void RenderFlightPlanApsides(String vessel_guid,
                                       XYZ sun_world_position) {
    foreach (CelestialBody celestial in
             plotting_frame_selector_.get().Bodies()) {
      IntPtr apoapsis_iterator;
      IntPtr periapsis_iterator;
      plugin_.FlightPlanRenderedApsides(vessel_guid,
                                        celestial.flightGlobalsIndex,
                                        sun_world_position,
                                        out apoapsis_iterator,
                                        out periapsis_iterator);
      map_node_pool_.RenderAndDeleteApsides(apoapsis_iterator,
                                            celestial,
                                            MapObject.ObjectType.Apoapsis);
      map_node_pool_.RenderAndDeleteApsides(periapsis_iterator,
                                            celestial,
                                            MapObject.ObjectType.Periapsis);
    }
  }

  private void Cleanup() {
    UnityEngine.Object.Destroy(map_renderer_);
    map_node_pool_.Clear();
    map_renderer_ = null;
    Interface.DeletePlugin(ref plugin_);
    plotting_frame_selector_.reset();
    flight_planner_.reset();
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
    } else if (plugin_.PhysicsBubbleIsEmpty()) {
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
    String version;
    String build_date;
    Interface.GetVersion(build_date: out build_date, version: out version);
    UnityEngine.GUILayout.TextArea(version + " built on " + build_date);
    bool changed_history_length = false;
    Selector(history_lengths_,
             ref history_length_index_,
             "Max history length",
             ref changed_history_length,
             "{0:0.00e00} s");
    ReferenceFrameSelection();
    if (PluginRunning()) {
      flight_planner_.get().RenderButton();
    }
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
        position : new UnityEngine.Rect(x      : 0f,
                                        y      : 0f,
                                        width  : 10000f,
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
    if (PluginRunning()) {
      bool was_fixing_navball_in_plotting_frame =
          fix_navball_in_plotting_frame_;
      fix_navball_in_plotting_frame_ = 
          UnityEngine.GUILayout.Toggle(
              value : fix_navball_in_plotting_frame_,
              text  : "Fix navball in plotting frame");
      if (was_fixing_navball_in_plotting_frame !=
          fix_navball_in_plotting_frame_) {
        navball_changed_ = true;
        reset_rsas_target_ = true;
      }
      plotting_frame_selector_.get().RenderButton();
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
    UnityEngine.GUILayout.TextArea(
        text    : String.Format(Culture.culture, format, array[index]),
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
  }

  private void KSPFeatures() {
    display_patched_conics_ = UnityEngine.GUILayout.Toggle(
        value : display_patched_conics_,
        text  : "Display patched conics (do not use for flight planning!)");
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

  private void UpdateRenderingFrame(
      NavigationFrameParameters frame_parameters) {
    IntPtr new_plotting_frame = plugin_.NewNavigationFrame(frame_parameters);
    plugin_.SetPlottingFrame(ref new_plotting_frame);
    if (fix_navball_in_plotting_frame_) {
     navball_changed_ = true;
      reset_rsas_target_ = true;
    }
  }

  private void ResetPlugin() {
    Cleanup();
    SetRotatingFrameThresholds();
    RemoveBuggyTidalLocking();
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
        plugin_ =
            Interface.NewPlugin(double.Parse(initial_states.GetValue("epoch")),
                                Planetarium.InverseRotAngle);
        var name_to_initial_state = new Dictionary<String, ConfigNode>();
        var name_to_gravity_model = new Dictionary<String, ConfigNode>();
        foreach (ConfigNode node in initial_states.GetNodes("body")) {
          name_to_initial_state.Add(node.GetValue("name"), node);
        }
        foreach (ConfigNode node in gravity_models.GetNodes("body")) {
          name_to_gravity_model.Add(node.GetValue("name"), node);
        }
        BodyProcessor insert_body = body => {
          Log.Info("Inserting " + body.name + "...");
          ConfigNode gravity_model = name_to_gravity_model[body.name];
          ConfigNode initial_state = name_to_initial_state[body.name];
          int? parent_index = null;
          if (body.orbit != null) {
            parent_index = body.orbit.referenceBody.flightGlobalsIndex;
          }
          plugin_.InsertCelestialAbsoluteCartesian(
              celestial_index: body.flightGlobalsIndex,
              parent_index: parent_index,
              gravitational_parameter:
                  gravity_model.GetValue("gravitational_parameter"),
              mean_radius:
                  gravity_model.GetValue("mean_radius"),
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
        insert_body(Planetarium.fetch.Sun);
        ApplyToBodyTree(insert_body);
        plugin_.EndInitialization();
        plugin_.AdvanceTime(Planetarium.GetUniversalTime(),
                            Planetarium.InverseRotAngle);
      } catch (Exception e) {
        Log.Fatal("Exception while reading initial state: " + e.ToString());
      }
    } else {
      plugin_source_ = PluginSource.ORBITAL_ELEMENTS;
      // We create the plugin at time 0, rather than
      // |Planetarium.GetUniversalTime()|, in order to get a deterministic
      // initial state.
      plugin_ = Interface.NewPlugin(0,
                                    Planetarium.InverseRotAngle);
      plugin_.InsertSun(Planetarium.fetch.Sun.flightGlobalsIndex,
                        Planetarium.fetch.Sun.gravParameter,
                        Planetarium.fetch.Sun.Radius);
      BodyProcessor insert_body = body => {
        Log.Info("Inserting " + body.name + "...");
        Orbit orbit = unmodified_orbits_[body];
        double mean_motion = 2 * Math.PI / orbit.period;
        plugin_.InsertCelestialJacobiKeplerian(
            celestial_index             : body.flightGlobalsIndex,
            parent_index                : body.referenceBody.flightGlobalsIndex,
            gravitational_parameter     : body.gravParameter + " m^3/s^2",
            mean_radius                 : body.Radius + " m",
            axis_right_ascension        : null,
            axis_declination            : null,
            j2                          : null,
            reference_radius            : null,
            eccentricity                : orbit.eccentricity,
            mean_motion                 : mean_motion + " rad/s",
            inclination                 : orbit.inclination + " deg",
            longitude_of_ascending_node : orbit.LAN + " deg",
            argument_of_periapsis       : orbit.argumentOfPeriapsis + " deg",
            mean_anomaly                : orbit.meanAnomalyAtEpoch -
                                          orbit.epoch * mean_motion + " rad");
      };
      ApplyToBodyTree(insert_body);
      plugin_.EndInitialization();
    }
    plotting_frame_selector_.reset(
        new ReferenceFrameSelector(this,
                                   plugin_,
                                   UpdateRenderingFrame,
                                   "Plotting frame"));
    flight_planner_.reset(new FlightPlanner(this, plugin_));
    VesselProcessor insert_vessel = vessel => {
      Log.Info("Inserting " + vessel.name + "...");
      bool inserted =
          plugin_.InsertOrKeepVessel(
              vessel.id.ToString(),
              vessel.orbit.referenceBody.flightGlobalsIndex);
      if (!inserted) {
        Log.Fatal("Plugin initialization: vessel not inserted");
      } else {
        plugin_.SetVesselStateOffset(vessel.id.ToString(),
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

  private void RemoveBuggyTidalLocking() {
    ApplyToBodyTree(body => body.tidallyLocked = false);
  }

  // Deals with issue #631, unstability of the Jool system's resonance.
  private void Fix631() {
    // Check whether this looks like stock.
    if (FlightGlobals.Bodies.Count < 15) {
      return;
    }
    Func<CelestialBody, double> mean_longitude =
        (celestial) =>
            (celestial.orbit.LAN + celestial.orbit.argumentOfPeriapsis) *
                180 / Math.PI +
            celestial.orbit.meanAnomalyAtEpoch;
    CelestialBody jool = FlightGlobals.Bodies[8];
    CelestialBody laythe = FlightGlobals.Bodies[9];
    CelestialBody vall = FlightGlobals.Bodies[10];
    CelestialBody bop = FlightGlobals.Bodies[11];
    CelestialBody tylo = FlightGlobals.Bodies[12];
    CelestialBody pol = FlightGlobals.Bodies[14];
    bool is_stock = true;
    is_stock &= jool.orbitingBodies.Count == 5;
    is_stock &= jool.name   == "Jool";
    is_stock &= laythe.name == "Laythe";
    is_stock &= vall.name   == "Vall";
    is_stock &= bop.name    == "Bop";
    is_stock &= tylo.name   == "Tylo";
    is_stock &= pol.name    == "Pol";
    is_stock &= laythe.referenceBody == jool;
    is_stock &= vall.referenceBody   == jool;
    is_stock &= tylo.referenceBody   == jool;
    is_stock &= (float)laythe.orbit.semiMajorAxis == 27184000f;
    is_stock &= (float)vall.orbit.semiMajorAxis   == 43152000f;
    is_stock &= (float)tylo.orbit.semiMajorAxis   == 68500000f;
    is_stock &= (float)laythe.orbit.inclination == 0f;
    is_stock &= (float)vall.orbit.inclination   == 0f;
    is_stock &= (float)tylo.orbit.inclination   == 0.025f;
    is_stock &= (float)mean_longitude(laythe) == 3.14f;
    is_stock &= (float)mean_longitude(vall)   == 0.9f;
    is_stock &= (float)mean_longitude(tylo)   == 3.14f;
    if (!is_stock) {
      return;
    }
    const double φ = 1.61803398875;
    // The |unmodified_orbits_| are unmodified in the sense that they are
    // unaffected by the plugin's computation; they are the orbits from which
    // we can reproducibly construct a fresh plugin.  In stock we modify them
    // from their stock values for stability reasons.
    unmodified_orbits_[vall] = new Orbit(
        vall.orbit.inclination,
        vall.orbit.eccentricity,
        laythe.orbit.semiMajorAxis * Math.Pow(4 / φ, 2.0 / 3.0),
        vall.orbit.LAN,
        vall.orbit.argumentOfPeriapsis,
        vall.orbit.meanAnomalyAtEpoch,
        vall.orbit.epoch,
        jool);
    unmodified_orbits_[tylo] = new Orbit(
        tylo.orbit.inclination,
        tylo.orbit.eccentricity,
        laythe.orbit.semiMajorAxis *
            Math.Pow(16 / (φ * φ), 2.0 / 3.0),
        tylo.orbit.LAN,
        tylo.orbit.argumentOfPeriapsis,
        tylo.orbit.meanAnomalyAtEpoch,
        tylo.orbit.epoch,
        jool);
    unmodified_orbits_[bop] =
        new Orbit(180.0 - bop.orbit.inclination,
                  bop.orbit.eccentricity,
                  pol.orbit.semiMajorAxis * Math.Pow(0.7, 2.0 / 3.0),
                  bop.orbit.LAN,
                  bop.orbit.argumentOfPeriapsis,
                  bop.orbit.meanAnomalyAtEpoch,
                  bop.orbit.epoch,
                  jool);
    // Vall and Tylo are tidally locked, so KSP will set their rotation period
    // to their orbital period.  Since we disable tidal locking before starting
    // the plugin (because tidal locking is buggy), we set their orbits here
    // to set their rotation period to their orbital period, so that they still
    // appear tidally locked (note that since we set the rotation to the Jacobi
    // osculating orbital period, there remains some noticeable drift; it may be
    // a good idea to compute the mean orbital period offline instead).  We do
    // not do that for Bop (which is tidally locked in stock) to make it look
    // like a more irregular satellite (in any case, Bop orbits retrograde, and
    // it is not clear whether the game supports retrograde rotation).
    foreach (CelestialBody body in new CelestialBody[]{vall, tylo}) {
      body.orbit.inclination = unmodified_orbits_[body].inclination;
      body.orbit.eccentricity = unmodified_orbits_[body].eccentricity;
      body.orbit.semiMajorAxis = unmodified_orbits_[body].semiMajorAxis;
      body.orbit.LAN = unmodified_orbits_[body].LAN;
      body.orbit.argumentOfPeriapsis =
          unmodified_orbits_[body].argumentOfPeriapsis;
      body.orbit.meanAnomalyAtEpoch =
          unmodified_orbits_[body].meanAnomalyAtEpoch;
      body.orbit.epoch = unmodified_orbits_[body].epoch;
      body.orbit.referenceBody = unmodified_orbits_[body].referenceBody;
      body.orbit.Init();
      body.orbit.UpdateFromUT(Planetarium.GetUniversalTime());
      body.CBUpdate();
    }
  }

}

}  // namespace ksp_plugin_adapter
}  // namespace principia
