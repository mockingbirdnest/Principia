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

  private const String principia_key = "serialized_plugin";
  private const String principia_initial_state_config_name =
      "principia_initial_state";
  private const String principia_gravity_model_config_name =
      "principia_gravity_model";
  private const double Δt = 10;

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
      {1e-3, 1e-2, 1e0, 1e1, 1e2, 1e3, 1e4};
  [KSPField(isPersistant = true)]
  private int prediction_length_tolerance_index_ = 1;
  private readonly double[] prediction_steps_ =
      {1 << 2, 1 << 4, 1 << 6, 1 << 8, 1 << 10, 1 << 12, 1 << 14, 1 << 16,
       1 << 18, 1 << 20, 1 << 22, 1 << 24};
  [KSPField(isPersistant = true)]
  private int prediction_steps_index_ = 4;
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

  // Whether a journal will be recorded when the plugin is next constructed.
  [KSPField(isPersistant = true)]
  private bool must_record_journal_ = false;
  // Whether a journal is currently being recorded.
  private static bool journaling_;
#if CRASH_BUTTON
  [KSPField(isPersistant = true)]
  private bool show_crash_options_ = false;
#endif
  // Timing diagnostics.
  private System.Diagnostics.Stopwatch stopwatch_ =
      new System.Diagnostics.Stopwatch();
  private double last_universal_time_;
  private double slowdown_;

  private bool time_is_advancing_;

  private DateTime plugin_construction_;

  private RenderingActions map_renderer_;
  private RenderingActions galaxy_cube_rotator_;

  private enum PluginSource {
   SAVED_STATE,
   ORBITAL_ELEMENTS,
   CARTESIAN_CONFIG,
  }
  private PluginSource plugin_source_;

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
  
  private Krakensbane krakensbane_;
  private Krakensbane krakensbane {
    get {
     return krakensbane_.GetValueOrDefault(
         krakensbane_ = (Krakensbane)FindObjectOfType(typeof(Krakensbane)));
    }
  }

  // The first apocalyptic error message.
  [KSPField(isPersistant = true)]
  private String revelation_ = "";
  // Whether we have encountered an apocalypse already.
  [KSPField(isPersistant = true)]
  private bool is_post_apocalyptic_ = false;
  [KSPField(isPersistant = true)]
  private int apocalypse_window_x_ = UnityEngine.Screen.width / 2;
  [KSPField(isPersistant = true)]
  private int apocalypse_window_y_ = UnityEngine.Screen.height / 3;
  private UnityEngine.Rect apocalypse_window_rectangle_;

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

  private void ApplyToVesselsInSpace(VesselProcessor process_vessel) {
     foreach (Vessel vessel in FlightGlobals.Vessels.Where(is_in_space)) {
       process_vessel(vessel);
    }
  }

  private void ApplyToVesselsOnRails(VesselProcessor process_vessel) {
    foreach (Vessel vessel in
             FlightGlobals.Vessels.Where(is_on_rails_in_space)) {
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
     if (plugin_.HasVessel(vessel.id.ToString())) {
       QP from_parent = plugin_.VesselFromParent(vessel.id.ToString());
       vessel.orbit.UpdateFromStateVectors(pos : (Vector3d)from_parent.q,
                                           vel : (Vector3d)from_parent.p,
                                           refBody : vessel.orbit.referenceBody,
                                           UT : universal_time);
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

  private bool has_active_vessel_in_space() {
    Vessel active_vessel = FlightGlobals.ActiveVessel;
    return active_vessel != null && is_in_space(active_vessel);
  }

  private bool draw_active_vessel_trajectory() {
    return MapView.MapIsEnabled &&
           has_active_vessel_in_space();
  }

  private void OverrideRSASTarget(FlightCtrlState state) {
    if (override_rsas_target_ && FlightGlobals.ActiveVessel.Autopilot.Enabled) {
      FlightGlobals.ActiveVessel.Autopilot.SAS.SetTargetOrientation(
          rsas_target_,
          reset_rsas_target_);
      FlightGlobals.ActiveVessel.Autopilot.SAS.ConnectFlyByWire(
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
            new Orbit(inc   : celestial.orbit.inclination,
                      e     : celestial.orbit.eccentricity,
                      sma   : celestial.orbit.semiMajorAxis,
                      lan   : celestial.orbit.LAN,
                      argPe : celestial.orbit.argumentOfPeriapsis,
                      mEp   : celestial.orbit.meanAnomalyAtEpoch,
                      t     : celestial.orbit.epoch,
                      body  : celestial.orbit.referenceBody));
      }
    }

    GameEvents.onShowUI.Add(ShowGUI);
    GameEvents.onHideUI.Add(HideGUI);
    // Timing0, -8008 on the script execution order page.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.ObscenelyEarly,
                                 DisableVesselPrecalculate);
    // TimingPre, -101 on the script execution order page.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.Precalc,
                                 SetBodyFramesAndPrecalculateVessels);
    // Timing1, -99 on the script execution order page.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.Early,
                                 UpdateVesselOrbits);
    // Timing3, 7.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.FashionablyLate,
                                 ReportVesselsAndParts);
  }

  public override void OnSave(ConfigNode node) {
    base.OnSave(node);
    if (PluginRunning()) {
      String serialization;
      IntPtr serializer = IntPtr.Zero;
      for (;;) {
        serialization = plugin_.SerializePlugin(ref serializer);
        if (serialization == null) {
          break;
        }
        node.AddValue(principia_key, serialization);
      }
    }
  }

  public override void OnLoad(ConfigNode node) {
    base.OnLoad(node);
    if (must_record_journal_) {
      journaling_ = true;
      Log.ActivateRecorder(true);
    }
    if (node.HasValue(principia_key)) {
      Cleanup();
      SetRotatingFrameThresholds();
      RemoveBuggyTidalLocking();
      Log.SetBufferedLogging(buffered_logging_);
      Log.SetSuppressedLogging(suppressed_logging_);
      Log.SetStderrLogging(stderr_logging_);
      Log.SetVerboseLogging(verbose_logging_);

      IntPtr deserializer = IntPtr.Zero;
      String[] serializations = node.GetValues(principia_key);
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

    if (is_post_apocalyptic_) {
      UnityEngine.GUI.skin = null;
      apocalypse_window_rectangle_.xMin = apocalypse_window_x_;
      apocalypse_window_rectangle_.yMin = apocalypse_window_y_;
      apocalypse_window_rectangle_ = UnityEngine.GUILayout.Window(
          id         : this.GetHashCode(),
          screenRect : apocalypse_window_rectangle_,
          func       : (int id) => {
              UnityEngine.GUILayout.BeginVertical();
              UnityEngine.GUILayout.TextArea(revelation_);
              UnityEngine.GUILayout.EndVertical();
            },
          text       : "Principia",
          options    : UnityEngine.GUILayout.MinWidth(500));
      apocalypse_window_x_ = (int)apocalypse_window_rectangle_.xMin;
      apocalypse_window_y_ = (int)apocalypse_window_rectangle_.yMin;
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

  private void LateUpdate() {
    if (map_renderer_ == null) {
      map_renderer_ =
          PlanetariumCamera.Camera.gameObject.AddComponent<RenderingActions>();
      map_renderer_.post_render = RenderTrajectories;
    }

    if (galaxy_cube_rotator_ == null) {
      galaxy_cube_rotator_ = ScaledCamera.Instance.galaxyCamera.gameObject
                                 .AddComponent<RenderingActions>();
      galaxy_cube_rotator_.pre_cull = RotateGalaxyCube;
    }

    // Orient the celestial bodies.
    if (PluginRunning()) {
      foreach (var body in FlightGlobals.Bodies) {
        body.scaledBody.transform.rotation =
            (UnityEngine.QuaternionD)plugin_.CelestialRotation(
                body.flightGlobalsIndex);
      }
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

      // Orient the ball.
      if (PluginRunning() && fix_navball_in_plotting_frame_) {
        navball_.navBall.rotation =
            (UnityEngine.QuaternionD)navball_.attitudeGymbal *  // sic.
            (UnityEngine.QuaternionD)plugin_.NavballOrientation(
                (XYZ)Planetarium.fetch.Sun.position,
                (XYZ)(Vector3d)active_vessel.ReferenceTransform.position);
      }

      if (PluginRunning() &&
          has_active_vessel_in_space() &&
          plugin_.HasVessel(active_vessel.id.ToString()) &&
          FlightGlobals.speedDisplayMode ==
              FlightGlobals.SpeedDisplayModes.Orbit) {
        KSP.UI.Screens.Flight.SpeedDisplay speed_display =
            KSP.UI.Screens.Flight.SpeedDisplay.Instance;
        if (speed_display?.textTitle != null &&
            speed_display?.textSpeed != null) {
          speed_display.textTitle.text =
              plotting_frame_selector_.get().ShortName();
          speed_display.textSpeed.text =
              ((Vector3d)plugin_.VesselVelocity(active_vessel.id.ToString()))
                  .magnitude.ToString("F1") + "m/s";
        }

        // Orient the Frenet trihedron.
        Vector3d prograde =
            (Vector3d)plugin_.VesselTangent(active_vessel.id.ToString());
        Vector3d radial =
            (Vector3d)plugin_.VesselNormal(active_vessel.id.ToString());
        // Yes, the astrodynamicist's normal is the mathematician's binormal.
        // Don't ask.
        Vector3d normal =
            (Vector3d)plugin_.VesselBinormal(active_vessel.id.ToString());

        SetNavballVector(navball_.progradeVector, prograde);
        SetNavballVector(navball_.radialInVector, radial);
        SetNavballVector(navball_.normalVector, normal);
        SetNavballVector(navball_.retrogradeVector, -prograde);
        SetNavballVector(navball_.radialOutVector, -radial);
        SetNavballVector(navball_.antiNormalVector, -normal);

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

  private void FixedUpdate() {
    if (GameSettings.ORBIT_WARP_DOWN_AT_SOI) {
      Log.Info("Setting GameSettings.ORBIT_WARP_DOWN_AT_SOI to false");
      GameSettings.ORBIT_WARP_DOWN_AT_SOI = false;
    }

    if (PluginRunning()) {
      double universal_time = Planetarium.GetUniversalTime();
      slowdown_ = stopwatch_.ElapsedMilliseconds /
                  ((universal_time - last_universal_time_) * 1000.0);
      last_universal_time_ = universal_time;
      stopwatch_.Reset();
      stopwatch_.Start();

      plugin_.SetMainBody(
          FlightGlobals.currentMainBody.GetValueOrDefault(
              FlightGlobals.GetHomeBody()).flightGlobalsIndex);

      Vessel active_vessel = FlightGlobals.ActiveVessel;
      bool ready_to_draw_active_vessel_trajectory =
          draw_active_vessel_trajectory() &&
          plugin_.HasVessel(active_vessel.id.ToString());

      if (ready_to_draw_active_vessel_trajectory) {
        // TODO(egg): make the speed tolerance independent.  Also max_steps.
        AdaptiveStepParameters adaptive_step_parameters =
            new AdaptiveStepParameters {
              max_steps = (Int64)prediction_steps_[prediction_steps_index_],
              length_integration_tolerance =
                  prediction_length_tolerances_[
                      prediction_length_tolerance_index_],
              speed_integration_tolerance =
                  prediction_length_tolerances_[
                      prediction_length_tolerance_index_]};
        plugin_.VesselSetPredictionAdaptiveStepParameters(
            active_vessel.id.ToString(), adaptive_step_parameters);
        plugin_.SetPredictionLength(double.PositiveInfinity);
      }

      if (ready_to_draw_active_vessel_trajectory) {
        plugin_.UpdatePrediction(active_vessel.id.ToString());
      }
      plugin_.ForgetAllHistoriesBefore(universal_time -
                                       history_lengths_[history_length_index_]);
      // TODO(egg): Set the degrees of freedom of the origin of |World| (by
      // toying with Krakensbane and FloatingOrigin) here.

      // Now we let the game and Unity do their thing. among other things,
      // the FashionablyLate callbacks, including ReportNonConservativeForces,
      // then the FlightIntegrator's FixedUpdate will run, then the Vessel's,
      // and eventually the physics simulation.
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
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.ObscenelyEarly,
                                    DisableVesselPrecalculate);
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.Precalc,
                                    SetBodyFramesAndPrecalculateVessels);
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.Early,
                                    UpdateVesselOrbits);
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.FashionablyLate,
                                    ReportVesselsAndParts);
  }

  #endregion

  private void AdvanceTimeAndNudgeVesselsAfterPhysicsSimulation() {
     double universal_time = Planetarium.GetUniversalTime();

     plugin_.PrepareToReportCollisions();

     // The collisions are reported and stored into |currentCollisions| in
     // OnCollisionEnter|Stay|Exit, which occurred while we yielded.
     // Here, the |currentCollisions| are the collisions that occurred in the
     // physics simulation, which is why we report them before calling
     // |AdvanceTime|.
     foreach (Vessel vessel1 in FlightGlobals.VesselsLoaded) {
       if (plugin_.HasVessel(vessel1.id.ToString()) &&
           plugin_.IsLoaded(vessel1.id.ToString()) &&
           !vessel1.packed) {
         if (vessel1.isEVA && vessel1.evaController.OnALadder) {
           var vessel2 = vessel1.evaController.LadderPart.vessel;
           if (vessel2 != null && plugin_.HasVessel(vessel2.id.ToString()) &&
               !vessel2.packed) {
             plugin_.ReportCollision(vessel1.rootPart.flightID,
                                     vessel1.evaController.LadderPart.flightID);
           }
         }
         foreach (Part part1 in vessel1.parts) {
           foreach (var collider in part1.currentCollisions) {
             var part2 =
                 collider.gameObject.GetComponentUpwards<Part>();
             var vessel2 = part2.vessel;
             if (vessel2 != null &&
                 plugin_.HasVessel(vessel2.id.ToString()) &&
                 plugin_.IsLoaded(vessel1.id.ToString())) {
               plugin_.ReportCollision(part1.flightID, part2.flightID);
             }
           }
         }
       }
     }

     double plugin_time = plugin_.CurrentTime();
     if (plugin_time > universal_time) {
       // TODO(egg): Make this resistant to bad floating points up to 2ULPs,
       // and make it fatal again.
       Log.Error("Closed Timelike Curve: " + plugin_time + " > " +
                 universal_time + " plugin-universal=" +
                 (plugin_time - universal_time));
       time_is_advancing_ = false;
       return;
     } else if (plugin_time == universal_time) {
       time_is_advancing_ = false;
       return;
     }
     time_is_advancing_ = true;

     plugin_.FreeVesselsAndPartsAndCollectPileUps();

     foreach (Vessel vessel in FlightGlobals.VesselsLoaded) {
       if (vessel.packed ||
           !plugin_.HasVessel(vessel.id.ToString()) ||
           !plugin_.IsLoaded(vessel.id.ToString())) {
         continue;
       }
       foreach (Part part in vessel.parts) {
         if (part.rb == null) {
           continue;
         }
         plugin_.SetPartApparentDegreesOfFreedom(
             part.flightID,
             // TODO(egg): use the centre of mass.
             new QP{q = (XYZ)(Vector3d)part.rb.position,
                    p = (XYZ)(Vector3d)part.rb.velocity});
       }
     }

     plugin_.AdvanceTime(universal_time, Planetarium.InverseRotAngle);
     is_post_apocalyptic_ |= plugin_.HasEncounteredApocalypse(out revelation_);

     // We don't want to do too many things here, since all the KSP classes
     // still think they're in the preceding step.  We only nudge the Unity
     // transforms of loaded vessels & their parts.
     if (has_active_vessel_in_space() && !FlightGlobals.ActiveVessel.packed) {
       // TODO(egg): move this to the C++, this is just to check that I
       // understand the issue.
       foreach (Vessel vessel in FlightGlobals.VesselsLoaded) {
         // TODO(egg): if I understand anything, there should probably be a
         // special treatment for loaded packed vessels.  I don't understand
         // anything though.
         if (vessel.packed ||
             !plugin_.HasVessel(vessel.id.ToString()) ||
             !plugin_.IsLoaded(vessel.id.ToString())) {
           continue;
         }
         foreach (Part part in vessel.parts) {
           if (part.rb == null) {
             continue;
           }
           QP part_actual_degrees_of_freedom =
               plugin_.GetPartActualDegreesOfFreedom(part.flightID, FlightGlobals.ActiveVessel.rootPart.flightID);
           // TODO(egg): use the centre of mass.  Here it's a bit tedious, some
           // transform nonsense must probably be done.
           part.rb.position = (Vector3d)part_actual_degrees_of_freedom.q;
           part.rb.velocity = (Vector3d)part_actual_degrees_of_freedom.p;
         }
       }
       QP main_body_dof = plugin_.CelestialWorldDegreesOfFreedom(
           FlightGlobals.ActiveVessel.mainBody.flightGlobalsIndex,
           FlightGlobals.ActiveVessel.rootPart.flightID);
       krakensbane.FrameVel = -(Vector3d)main_body_dof.p;
       Vector3d offset = (Vector3d)main_body_dof.q -
                         FlightGlobals.ActiveVessel.mainBody.position;
       // We cannot use FloatingOrigin.SetOffset to move the world here, because
       // as far as I can tell, that does not move the bubble relative to the
       // rest of the universe.
       foreach (CelestialBody celestial in FlightGlobals.Bodies) {
         celestial.position += offset;
       }
       foreach (
           Vessel vessel in FlightGlobals.Vessels.Where(is_on_rails_in_space)) {
         vessel.SetPosition(vessel.transform.position + offset);
       }
       // NOTE(egg): this is almost certainly incorrect, since we give the
       // bodies their positions at the next instant, wherease KSP still expects
       // them at the previous instant, and will propagate them at the beginning
       // of the next frame...
     }
   }

  private void DisableVesselPrecalculate() {
    foreach (var vessel in FlightGlobals.Vessels) {
      vessel.precalc.enabled = false;
    }
  }

  private void SetBodyFramesAndPrecalculateVessels() {
    AdvanceTimeAndNudgeVesselsAfterPhysicsSimulation();
    SetBodyFrames();
    // Unfortunately there is no way to get scheduled between Planetarium and
    // VesselPrecalculate, so we get scheduled after VesselPrecalculate, set the
    // body frames for our weird tilt, and run VesselPrecalculate manually.
    // Sob.
    // NOTE(egg): we cannot use foreach here, and we must iterate downwards,
    // since vessel.precalc.FixedUpdate may remove its vessel.
    for (int i = FlightGlobals.Vessels.Count - 1; i >= 0; --i) {
      var vessel = FlightGlobals.Vessels[i];
      vessel.precalc.enabled = true;
      vessel.precalc.FixedUpdate();
    }
  }

  private void UpdateVesselOrbits() {
    if (PluginRunning()) {
      ApplyToVesselsOnRails(
          vessel => UpdateVessel(vessel, Planetarium.GetUniversalTime()));
    }
  }

  private void ReportVesselsAndParts() {
    // We fetch the forces from the census of nonconservatives here;
    // part.forces, part.force, and part.torque are cleared by the/
    // FlightIntegrator's FixedUpdate (while we are yielding).
    if (PluginRunning()) {
      if (has_active_vessel_in_space() && FlightGlobals.ActiveVessel.packed) {
        if (PhysicsGlobals.GraviticForceMultiplier != 0) {  // sic.
          Log.Info("Killing stock gravity");
          PhysicsGlobals.GraviticForceMultiplier = 0;
        }
      } else if (PhysicsGlobals.GraviticForceMultiplier == 0) {
        Log.Info("Reinstating stock gravity");
        PhysicsGlobals.GraviticForceMultiplier = 1;
      }
      foreach (Vessel vessel in FlightGlobals.Vessels.Where(is_in_space)) {
        bool inserted;
        plugin_.InsertOrKeepVessel(vessel.id.ToString(),
                                   vessel.vesselName,
                                   vessel.mainBody.flightGlobalsIndex,
                                   !vessel.packed,
                                   out inserted);
        if (!vessel.packed) {
          QP main_body_degrees_of_freedom =
              new QP{q = (XYZ)vessel.mainBody.position,
                     p = (XYZ)(-krakensbane.FrameVel)};
          foreach (Part part in vessel.parts.Where((part) => part.rb != null)) {
            plugin_.InsertOrKeepLoadedPart(
                part.flightID,
                part.partName,
                part.physicsMass,
                vessel.id.ToString(),
                vessel.mainBody.flightGlobalsIndex,
                main_body_degrees_of_freedom,
                // TODO(egg): use the centre of mass.
                new QP{q = (XYZ)(Vector3d)part.rb.position,
                       p = (XYZ)(Vector3d)part.rb.velocity});
            plugin_.IncrementPartIntrinsicForce(part.flightID, (XYZ)part.force);
            foreach (var force in part.forces) {
              plugin_.IncrementPartIntrinsicForce(part.flightID,
                                                  (XYZ)force.force);
            }
          }
        } else if (inserted) {
          foreach (ProtoPartSnapshot part in
                   vessel.protoVessel.protoPartSnapshots) {
            plugin_.InsertUnloadedPart(
                part.flightID,
                part.partName,
                vessel.id.ToString(),
                new QP{q = (XYZ)vessel.orbit.pos, p = (XYZ)vessel.orbit.vel});
          }
        }
      }
    }
  }

  private void SetBodyFrames() {
    if (PluginRunning()) {
      if (FlightGlobals.currentMainBody != null) {
        FlightGlobals.currentMainBody.rotationPeriod =
            plugin_.CelestialRotationPeriod(
                FlightGlobals.currentMainBody.flightGlobalsIndex);
        FlightGlobals.currentMainBody.initialRotation =
            plugin_.CelestialInitialRotationInDegrees(
                FlightGlobals.currentMainBody.flightGlobalsIndex);
      }
      ApplyToBodyTree(body => UpdateBody(body, Planetarium.GetUniversalTime()));
      foreach (var body in FlightGlobals.Bodies) {
        // TODO(egg): I have no idea why this |swizzle| thing makes things work.
        // This probably really means something in terms of frames that should
        // be done in the C++ instead---once I figure out what it is.
        var swizzly_body_world_to_world =
            ((UnityEngine.QuaternionD)plugin_.CelestialRotation(
                 body.flightGlobalsIndex)).swizzle;
        body.BodyFrame = new Planetarium.CelestialFrame{
            X = swizzly_body_world_to_world * new Vector3d{x = 1, y = 0, z = 0},
            Y = swizzly_body_world_to_world * new Vector3d{x = 0, y = 1, z = 0},
            Z = swizzly_body_world_to_world * new Vector3d{x = 0, y = 0, z = 1}};
      }
    }
  }

  private void SetNavballVector(UnityEngine.Transform vector,
                                Vector3d direction) {
     vector.localPosition = (UnityEngine.QuaternionD)navball_.attitudeGymbal *
                            direction * navball_.VectorUnitScale;
     vector.gameObject.SetActive(vector.localPosition.z >
                                 navball_.VectorUnitCutoff);
     vector.GetComponent<UnityEngine.MeshRenderer>().materials[0].SetFloat(
         "_Opacity",
         UnityEngine.Mathf.Clamp01(UnityEngine.Vector3.Dot(
             vector.localPosition.normalized,
             UnityEngine.Vector3.forward)));
  }

  private void RotateGalaxyCube() {
    if (PluginRunning()) {
      var initial_rotation =
          UnityEngine.QuaternionD.Inverse(Planetarium.Rotation) *
          (UnityEngine.QuaternionD)
              GalaxyCubeControl.Instance.transform.rotation;
      GalaxyCubeControl.Instance.transform.rotation =
          (UnityEngine.QuaternionD)plugin_.CelestialSphereRotation() *
          initial_rotation;
    }
  }

  private void RemoveStockTrajectoriesIfNeeded(Vessel vessel) {
    vessel.patchedConicRenderer.relativityMode =
        PatchRendering.RelativityMode.RELATIVE;
    if (display_patched_conics_) {
      if (!vessel.patchedConicRenderer.enabled) {
        vessel.patchedConicRenderer.enabled = true;
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
            XKCDColors.AcidGreen,
            GLLines.Style.FADED);
        RenderPredictionApsides(active_vessel_guid, sun_world_position);
        GLLines.RenderAndDeleteTrajectory(
            plugin_.RenderedPrediction(active_vessel_guid, sun_world_position),
            XKCDColors.Fuchsia,
            GLLines.Style.SOLID);
        if (plugin_.FlightPlanExists(active_vessel_guid)) {
          RenderFlightPlanApsides(active_vessel_guid, sun_world_position);

          int number_of_segments =
              plugin_.FlightPlanNumberOfSegments(active_vessel_guid);
          for (int i = 0; i < number_of_segments; ++i) {
            bool is_burn = i % 2 == 1;
            var rendered_segments = plugin_.FlightPlanRenderedSegment(
                active_vessel_guid, sun_world_position, i);
            Vector3d position_at_start =
                (Vector3d)rendered_segments.IteratorGetXYZ();
            GLLines.RenderAndDeleteTrajectory(
                rendered_segments,
                is_burn ? XKCDColors.OrangeRed : XKCDColors.BabyBlue,
                is_burn ? GLLines.Style.SOLID : GLLines.Style.DASHED);
            if (is_burn) {
              int manoeuvre_index = i / 2;
              NavigationManoeuvre manoeuvre = plugin_.FlightPlanGetManoeuvre(
                                                  active_vessel_guid,
                                                  manoeuvre_index);
              double scale = (ScaledSpace.ScaledToLocalSpace(
                                  MapView.MapCamera.transform.position) -
                              position_at_start).magnitude * 0.015;
              Action<XYZ, UnityEngine.Color> add_vector =
                  (world_direction, colour) => {
                UnityEngine.GL.Color(colour);
                GLLines.AddSegment(
                    position_at_start,
                    position_at_start + scale * (Vector3d)world_direction,
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

  private void RenderPredictionApsides(String vessel_guid,
                                       XYZ sun_world_position) {
    foreach (CelestialBody celestial in
             plotting_frame_selector_.get().FixedBodies()) {
      IntPtr apoapsis_iterator;
      IntPtr periapsis_iterator;
      plugin_.RenderedPredictionApsides(vessel_guid,
                                        celestial.flightGlobalsIndex,
                                        sun_world_position,
                                        out apoapsis_iterator,
                                        out periapsis_iterator);
      map_node_pool_.RenderAndDeleteApsides(apoapsis_iterator,
                                            celestial,
                                            MapObject.ObjectType.Apoapsis,
                                            MapNodePool.NodeSource.PREDICTION);
      map_node_pool_.RenderAndDeleteApsides(periapsis_iterator,
                                            celestial,
                                            MapObject.ObjectType.Periapsis,
                                            MapNodePool.NodeSource.PREDICTION);
    }
  }

  private void RenderFlightPlanApsides(String vessel_guid,
                                       XYZ sun_world_position) {
    foreach (CelestialBody celestial in
             plotting_frame_selector_.get().FixedBodies()) {
      IntPtr apoapsis_iterator;
      IntPtr periapsis_iterator;
      plugin_.FlightPlanRenderedApsides(vessel_guid,
                                        celestial.flightGlobalsIndex,
                                        sun_world_position,
                                        out apoapsis_iterator,
                                        out periapsis_iterator);
      map_node_pool_.RenderAndDeleteApsides(apoapsis_iterator,
                                            celestial,
                                            MapObject.ObjectType.Apoapsis,
                                            MapNodePool.NodeSource.FLIGHT_PLAN);
      map_node_pool_.RenderAndDeleteApsides(periapsis_iterator,
                                            celestial,
                                            MapObject.ObjectType.Periapsis,
                                            MapNodePool.NodeSource.FLIGHT_PLAN);
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
    } else {
      plugin_state = "running";
    }
    UnityEngine.GUILayout.TextArea(text : "Plugin is " + plugin_state);
    UnityEngine.GUILayout.TextArea(
        "Time runs slowed by " +
        slowdown_);
    if (FlightGlobals.ActiveVessel != null) {
      UnityEngine.GUILayout.TextArea(FlightGlobals.ActiveVessel.geeForce +
                                     " g0");
    }
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
    Selector(prediction_steps_,
             ref prediction_steps_index_,
             "Steps",
             ref changed_settings,
             "{0:0.00e0}");
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
          text    : Log.severity_names[severity],
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
    UnityEngine.GUILayout.TextArea("Journalling is " +
                                   (journaling_ ? "ON" : "OFF"));
    must_record_journal_ = UnityEngine.GUILayout.Toggle(
        value   : must_record_journal_,
        text    : "Record journal (starts on load)");
    if (journaling_ && !must_record_journal_) {
      // We can deactivate a recorder at any time, but in order for replaying to
      // work, we should only activate one before creating a plugin.
      journaling_ = false;
      Interface.ActivateRecorder(false);
    }
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
    Dictionary<String, ConfigNode> name_to_gravity_model = null;
    var gravity_model_configs =
        GameDatabase.Instance.GetConfigs(principia_gravity_model_config_name);
    var cartesian_configs =
        GameDatabase.Instance.GetConfigs(principia_initial_state_config_name);
    if (gravity_model_configs.Length == 1) {
      name_to_gravity_model =
          gravity_model_configs[0].config.GetNodes("body").
              ToDictionary(node => node.GetValue("name"));
    } else if (gravity_model_configs.Length > 1) {
      Log.Fatal("too many gravity models (" + gravity_model_configs.Length +
                ")");
    }
    if (cartesian_configs.Length > 0) {
      plugin_source_ = PluginSource.CARTESIAN_CONFIG;
      if (cartesian_configs.Length > 1) {
        Log.Fatal("too many Cartesian configs (" + cartesian_configs.Length +
                  ")");
      }
      if (name_to_gravity_model == null) {
        Log.Fatal("Cartesian config without gravity models");
      }
      try {
        ConfigNode initial_states = GameDatabase.Instance.GetConfigs(
            principia_initial_state_config_name)[0].config;
        plugin_ =
            Interface.NewPlugin(initial_states.GetValue("game_epoch"),
                                initial_states.GetValue("solar_system_epoch"),
                                Planetarium.InverseRotAngle);
        var name_to_initial_state =
            initial_states.GetNodes("body").
                ToDictionary(node => node.GetValue("name"));
        BodyProcessor insert_body = body => {
          Log.Info("Inserting " + body.name + "...");
          ConfigNode gravity_model;
          if (!name_to_gravity_model.TryGetValue(body.name,
                                                 out gravity_model)) {
             Log.Fatal("missing gravity model for " + body.name);
          }
          ConfigNode initial_state;
          if (!name_to_initial_state.TryGetValue(body.name,
                                                 out initial_state)) {
             Log.Fatal("missing Cartesian initial state for " + body.name);
          }
          int? parent_index = body.orbit?.referenceBody.flightGlobalsIndex;
          var body_parameters = new BodyParameters{
              name = body.name,
              gravitational_parameter =
                  gravity_model.GetValue("gravitational_parameter"),
              reference_instant       =
                  double.Parse(gravity_model.GetValue("reference_instant")),
              mean_radius             = gravity_model.GetValue("mean_radius"),
              axis_right_ascension    =
                  gravity_model.GetValue("axis_right_ascension"),
              axis_declination        =
                  gravity_model.GetValue("axis_declination"),
              reference_angle         =
                  gravity_model.GetValue("reference_angle"),
              angular_frequency       =
                  gravity_model.GetValue("angular_frequency"),
              j2                      = gravity_model.GetValue("j2"),
              reference_radius        =
                  gravity_model.GetValue("reference_radius")};
          plugin_.InsertCelestialAbsoluteCartesian(
              celestial_index         : body.flightGlobalsIndex,
              parent_index            : parent_index,
              body_parameters         : body_parameters,
              x                       : initial_state.GetValue("x"),
              y                       : initial_state.GetValue("y"),
              z                       : initial_state.GetValue("z"),
              vx                      : initial_state.GetValue("vx"),
              vy                      : initial_state.GetValue("vy"),
              vz                      : initial_state.GetValue("vz"));
        };
        insert_body(Planetarium.fetch.Sun);
        ApplyToBodyTree(insert_body);
        plugin_.EndInitialization();
      } catch (Exception e) {
        Log.Fatal("Exception while reading initial state: " + e.ToString());
      }
    } else {
      plugin_source_ = PluginSource.ORBITAL_ELEMENTS;
      // We create the plugin at time 0, rather than
      // |Planetarium.GetUniversalTime()|, in order to get a deterministic
      // initial state.
      for(;;) {
        plugin_ = Interface.NewPlugin("0 s", "0 s",
                                      Planetarium.InverseRotAngle);
        BodyProcessor insert_body = body => {
          Log.Info("Inserting " + body.name + "...");
          ConfigNode gravity_model = null;
          if (name_to_gravity_model?.TryGetValue(body.name,
                                                 out gravity_model) == true) {
            Log.Info("using custom gravity model");
          }
          Orbit orbit = unmodified_orbits_.GetValueOrNull(body);
          var body_parameters = new BodyParameters{
              name = body.name,
              gravitational_parameter =
                  (gravity_model?.GetValue("gravitational_parameter")).
                      GetValueOrDefault(body.gravParameter + " m^3/s^2"),
              // J2000, because that's when we start non-config games.  We
              // should really parse real-life dates from strings.
              // The origin of rotation in KSP is the x of Barycentric, rather
              // than the y axis as is the case for Earth, so the right
              // ascension is -90 deg.
              reference_instant    = double.Parse(
                  (gravity_model?.GetValue("reference_instant")).
                      GetValueOrDefault("2451545.0")),
              mean_radius          =
                  (gravity_model?.GetValue("mean_radius")).
                      GetValueOrDefault(body.Radius + " m"),
              axis_right_ascension =
                  (gravity_model?.GetValue("axis_right_ascension")).
                      GetValueOrDefault("-90 deg"),
              axis_declination     =
                  (gravity_model?.GetValue("axis_declination")).
                      GetValueOrDefault("90 deg"),
              reference_angle      =
                  (gravity_model?.GetValue("reference_angle")).
                      GetValueOrDefault(body.initialRotation.ToString() +
                                        " deg"),
              angular_frequency    =
                  (gravity_model?.GetValue("angular_frequency")).
                      GetValueOrDefault(body.angularV.ToString() + " rad/s"),
              j2                   = gravity_model?.GetValue("j2"),
              reference_radius     =
                  gravity_model?.GetValue("reference_radius")};
          plugin_.InsertCelestialJacobiKeplerian(
              celestial_index             : body.flightGlobalsIndex,
              parent_index                :
                  orbit?.referenceBody.flightGlobalsIndex,
              body_parameters             : body_parameters,
              keplerian_elements          : orbit?.Elements());
        };
        insert_body(Planetarium.fetch.Sun);
        ApplyToBodyTree(insert_body);
        plugin_.EndInitialization();
        if (plugin_.IsKspStockSystem()) {
          Interface.DeletePlugin(ref plugin_);
          Fix631();
        } else {
          break;
        }
      }
    }
    plugin_.AdvanceTime(Planetarium.GetUniversalTime(),
                        Planetarium.InverseRotAngle);
    plotting_frame_selector_.reset(
        new ReferenceFrameSelector(this,
                                   plugin_,
                                   UpdateRenderingFrame,
                                   "Plotting frame"));
    flight_planner_.reset(new FlightPlanner(this, plugin_));
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
    CelestialBody jool = FlightGlobals.Bodies[8];
    CelestialBody laythe = FlightGlobals.Bodies[9];
    CelestialBody vall = FlightGlobals.Bodies[10];
    CelestialBody bop = FlightGlobals.Bodies[11];
    CelestialBody tylo = FlightGlobals.Bodies[12];
    CelestialBody pol = FlightGlobals.Bodies[14];
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
      body.tidallyLocked = true;
      body.CBUpdate();
    }
    RemoveBuggyTidalLocking();
  }

}

}  // namespace ksp_plugin_adapter
}  // namespace principia
