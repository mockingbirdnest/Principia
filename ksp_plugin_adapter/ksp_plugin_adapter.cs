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

  private const String next_release_name_ = "Fano";
  private const int next_release_lunation_number_ = 238;
  private DateTimeOffset next_release_date_ =
      new DateTimeOffset(2019, 04, 05, 08, 51, 00, TimeSpan.Zero);

  // From https://forum.kerbalspaceprogram.com/index.php?/topic/84273--/,
  // edited 2017-03-09.  Where the name of the layer is not CamelCase, the
  // actual name is commented.
  private enum UnityLayers {
    TransparentFX = 1,
    IgnoreRaycast = 2,  // Ignore Raycast
    Water = 3,
    UI = 5,
    PartListIcons = 8,  // PartsList_Icons
    Atmosphere = 9,
    ScaledScenery = 10,  // Scaled Scenery
    UIDialog = 11,
    UIVectors = 12,
    UIMask = 13,
    Screens = 14,
    LocalScenery = 15,  // Local Scenery
    Kerbals = 16,
    EVA = 17,
    SkySphere = 18,
    PhysicalObjects = 19,
    InternalSpace = 20,  // Internal Space
    PartTriggers = 21,  // Part Triggers
    KerbalInstructors = 22,
    AeroFXIgnore = 23,
    MapFX = 24,
    UIAdditional = 25,
    WheelCollidersIgnore = 26,
    WheelColliders = 27,  // wheelColliders
    TerrainColliders = 28,
    DragRender = 29,
    SurfaceFX = 30,
    Vectors = 31,
  };

  private const String principia_serialized_plugin_ = "serialized_plugin";
  private const String principia_initial_state_config_name_ =
      "principia_initial_state";
  private const String principia_gravity_model_config_name_ =
      "principia_gravity_model";
  private const String principia_numerics_blueprint_config_name_ =
      "principia_numerics_blueprint";

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

#if SELECTABLE_PLOT_METHOD
  [KSPField(isPersistant = true)]
#endif
  private int чебышёв_plotting_method_ = 2;
  private const int чебышёв_plotting_methods_count = 3;

  internal Controlled<ReferenceFrameSelector> plotting_frame_selector_;
  private Controlled<FlightPlanner> flight_planner_;
  private MapNodePool map_node_pool_;

  private bool selecting_active_vessel_target_ = false;
  private bool selecting_target_celestial_ = false;

  private IntPtr plugin_ = IntPtr.Zero;
  internal IntPtr Plugin() {
    return plugin_;
  }

  private bool display_patched_conics_ = false;

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

  // Whether to compress saves.
  [KSPField(isPersistant = true)]
  private string serialization_compression_ = "";
  [KSPField(isPersistant = true)]
  private string serialization_encoding_ = "hexadecimal";

  // Whether the plotting frame must be set to something convenient at the next
  // opportunity.
  private bool must_set_plotting_frame_ = false;

  // Whether a journal is currently being recorded.
  private static bool journaling_;
#if CRASH_BUTTON
  [KSPField(isPersistant = true)]
  private bool show_crash_options_ = false;
#endif

  private bool time_is_advancing_;

  private DateTime plugin_construction_;

  private RenderingActions map_renderer_;
  private RenderingActions galaxy_cube_rotator_;

  private KSP.UI.Screens.Flight.NavBall navball_;
  private UnityEngine.Texture compass_navball_texture_;
  private UnityEngine.Texture inertial_navball_texture_;
  private UnityEngine.Texture barycentric_navball_texture_;
  private UnityEngine.Texture body_direction_navball_texture_;
  private UnityEngine.Texture surface_navball_texture_;
  private UnityEngine.Texture target_navball_texture_;
  private bool navball_changed_ = true;
  private FlightGlobals.SpeedDisplayModes? previous_display_mode_;
  private ReferenceFrameSelector.FrameType last_non_surface_frame_type_ =
      ReferenceFrameSelector.FrameType.BODY_CENTRED_NON_ROTATING;

  List<IntPtr> vessel_futures_ = new List<IntPtr>();

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

  private Krakensbane krakensbane_;
  private Krakensbane krakensbane {
    get {
     if (krakensbane_ == null) {
       krakensbane_ = (Krakensbane)FindObjectOfType(typeof(Krakensbane));
     }
     return krakensbane_;
    }
  }

  private KSP.UI.Screens.SpaceTracking space_tracking_;
  private KSP.UI.Screens.SpaceTracking space_tracking {
    get {
     if (space_tracking_ == null) {
       if (HighLogic.LoadedScene != GameScenes.TRACKSTATION) {
         return null;
       }
       space_tracking_ = (KSP.UI.Screens.SpaceTracking)FindObjectOfType(
                             typeof(KSP.UI.Screens.SpaceTracking));
     }
     return space_tracking_;
    }
  }

  // TODO(egg): these should be moved to the C++; there it can be made a vector
  // because we can test whether we have the part or not. Set in
  // FashionablyLate, before the FlightIntegrator clears the forces.  Used in
  // WaitForFixedUpdate.
  private Dictionary<uint, Part.ForceHolder[]> part_id_to_intrinsic_forces_ =
      new Dictionary<uint, Part.ForceHolder[]>();
  private Dictionary<uint, Vector3d> part_id_to_intrinsic_force_ =
      new Dictionary<uint, Vector3d>();

  // The degrees of freedom at BetterLateThanNever.  Those are used to insert
  // new parts with the correct initial state.
  private Dictionary<uint, QP> part_id_to_degrees_of_freedom_ =
      new Dictionary<uint, QP>();

  // UI for the apocalypse notification.
  // Whether we have encountered an apocalypse already.
  [KSPField(isPersistant = true)]
  private bool is_post_apocalyptic_ = false;
  [KSPField(isPersistant = true)]
  private Dialog apocalypse_dialog_ = new Dialog();

  // UI for the bad installation notification.
  private bool is_bad_installation_ = false;  // Don't persist.
  [KSPField(isPersistant = true)]
  private Dialog bad_installation_dialog_ = new Dialog();

  public event Action render_windows;

  PrincipiaPluginAdapter() {
    // We create this directory here so we do not need to worry about cross-
    // platform problems in C++.
    System.IO.Directory.CreateDirectory("glog/Principia");
    string load_error = Loader.LoadPrincipiaDllAndInitGoogleLogging();
    if (load_error != null) {
      is_bad_installation_ = true;
      bad_installation_dialog_.Message =
          "The Principia DLL failed to load.\n" + load_error;
    }
#if KSP_VERSION_1_3_1
    if (Versioning.version_major != 1 ||
        Versioning.version_minor != 3 ||
        Versioning.Revision != 1) {
      string expected_version = "1.3.1";
#elif KSP_VERSION_1_6_1
    if (!(Versioning.version_major == 1 &&
          (Versioning.version_minor == 4 &&
           (Versioning.Revision >= 1 && Versioning.Revision <= 5)) ||
          (Versioning.version_minor == 5 && Versioning.Revision == 1) ||
          (Versioning.version_minor == 6 && Versioning.Revision == 1))) {
      string expected_version =
          "1.6.1, 1.5.1, 1.4.5, 1.4.4, 1.4.3, 1.4.2, and 1.4.1";
#endif
      Log.Fatal("Unexpected KSP version " + Versioning.version_major + "." +
                Versioning.version_minor + "." + Versioning.Revision +
                "; this build targets " + expected_version + ".");
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

  private void ApplyToManageableVessels(VesselProcessor process_vessel) {
     foreach (Vessel vessel in FlightGlobals.Vessels.Where(is_manageable)) {
       process_vessel(vessel);
    }
  }

  private void ApplyToVesselsOnRails(VesselProcessor process_vessel) {
    foreach (Vessel vessel in
             FlightGlobals.Vessels.Where(is_manageable_on_rails)) {
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
    body.CBUpdate();
  }

  private void UpdatePredictions() {
    Vessel main_vessel = FlightGlobals.ActiveVessel ??
                  space_tracking?.SelectedVessel;
    bool ready_to_draw_active_vessel_trajectory =
        main_vessel != null &&
        MapView.MapIsEnabled &&
        plugin_.HasVessel(main_vessel.id.ToString());

    if (ready_to_draw_active_vessel_trajectory) {
      // TODO(egg): make the speed tolerance independent.  Also max_steps.
      AdaptiveStepParameters adaptive_step_parameters =
          plugin_.VesselGetPredictionAdaptiveStepParameters(
              main_vessel.id.ToString());
      adaptive_step_parameters =
          new AdaptiveStepParameters {
            integrator_kind = adaptive_step_parameters.integrator_kind,
            max_steps = (Int64)prediction_steps_[prediction_steps_index_],
            length_integration_tolerance =
                prediction_length_tolerances_[
                    prediction_length_tolerance_index_],
            speed_integration_tolerance =
                prediction_length_tolerances_[
                    prediction_length_tolerance_index_]};
      plugin_.VesselSetPredictionAdaptiveStepParameters(
          main_vessel.id.ToString(), adaptive_step_parameters);
      plugin_.UpdatePrediction(main_vessel.id.ToString());
      string target_id =
          FlightGlobals.fetch.VesselTarget?.GetVessel()?.id.ToString();
      if (!plotting_frame_selector_.get().target_override &&
          target_id != null && plugin_.HasVessel(target_id)) {
        plugin_.VesselSetPredictionAdaptiveStepParameters(
            target_id, adaptive_step_parameters);
        plugin_.UpdatePrediction(target_id);
      }
    }
  }

  private void UpdateVessel(Vessel vessel, double universal_time) {
     if (plugin_.HasVessel(vessel.id.ToString())) {
       QP from_parent = plugin_.VesselFromParent(
           vessel.mainBody.flightGlobalsIndex,
           vessel.id.ToString());
       vessel.orbit.UpdateFromStateVectors(pos : (Vector3d)from_parent.q,
                                           vel : (Vector3d)from_parent.p,
                                           refBody : vessel.orbit.referenceBody,
                                           UT : universal_time);
     }
  }

  private bool time_is_advancing(double universal_time) {
    double plugin_time = plugin_.CurrentTime();
    if (plugin_time > universal_time) {
      // TODO(egg): Make this resistant to bad floating points up to 2ULPs,
      // and make it fatal again.
      Log.Error("Closed Timelike Curve: " + plugin_time + " > " +
                universal_time +
                " plugin-universal=" + (plugin_time - universal_time));
      return false;
    } else if (plugin_time == universal_time) {
      return false;
    }
    return true;
  }

  private Part closest_physical_parent(Part part) {
    return part.rb ? part : closest_physical_parent(part.parent);
  }

  // Whether the given kerbal is clambering, i.e., climbing something that's not
  // a ladder.
  private bool is_clambering(KerbalEVA kerbal) {
    return kerbal.fsm.CurrentState == kerbal.st_clamber_acquireP1 ||
           kerbal.fsm.CurrentState == kerbal.st_clamber_acquireP2 ||
           kerbal.fsm.CurrentState == kerbal.st_clamber_acquireP3;
  }

  private bool is_manageable(Vessel vessel) {
    return UnmanageabilityReasons(vessel) == null;
  }

  private string UnmanageabilityReasons(Vessel vessel) {
    List<string> reasons = new List<string>(capacity : 3);
    if (vessel.state == Vessel.State.DEAD) {
      reasons.Add("vessel is dead");
    }
    if (vessel.id == Guid.Empty) {
      reasons.Add("vessel has an empty GUID");
    }
    if (!(vessel.situation == Vessel.Situations.SUB_ORBITAL ||
          vessel.situation == Vessel.Situations.ORBITING ||
          vessel.situation == Vessel.Situations.ESCAPING ||
          vessel.situation == Vessel.Situations.FLYING)) {
      reasons.Add("vessel situation is " + vessel.situation);
    }
    double height;
    double vertical_speed;
    if (!vessel.packed &&
        JustAboveTheGround(vessel, out height, out vertical_speed)) {
      reasons.Add("vessel is " + height +
                  " m above ground with a vertical speed of " + vertical_speed +
                  " m/s");
    }
    if (vessel.isEVA && vessel.evaController?.Ready == false) {
      reasons.Add("vessel is an unready Kerbal");
    }
    if (reasons.Count == 0) {
      return null;
    } else {
      return string.Join(", ", reasons.ToArray());
    }
  }

  private bool is_manageable_on_rails(Vessel vessel) {
    return vessel.packed && is_manageable(vessel);
  }

  private bool has_active_manageable_vessel() {
    Vessel active_vessel = FlightGlobals.ActiveVessel;
    return active_vessel != null && is_manageable(active_vessel);
  }

  private bool JustAboveTheGround(Vessel vessel,
                                  out double height,
                                  out double vertical_speed) {
    height = vessel.altitude - vessel.terrainAltitude;
    vertical_speed = vessel.verticalSpeed;
    double Δt = Planetarium.TimeScale * Planetarium.fetch.fixedDeltaTime;
    return height + vertical_speed * Δt < 0;
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
#if KSP_VERSION_1_6_1
      bool success = UnityEngine.ImageConversion.LoadImage(
          texture2d, File.ReadAllBytes(full_path));
#elif KSP_VERSION_1_3_1
      bool success = texture2d.LoadImage(
          File.ReadAllBytes(full_path));
#endif
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
    if (is_bad_installation_) {
      return;
    }
    // While we're here, we might as well log.
    Log.Info("principia.ksp_plugin_adapter.PrincipiaPluginAdapter.OnAwake()");

    LoadTextureIfExists(out compass_navball_texture_, "navball_compass.png");
    LoadTextureOrDie(out inertial_navball_texture_, "navball_inertial.png");
    LoadTextureOrDie(out barycentric_navball_texture_,
                     "navball_barycentric.png");
    LoadTextureOrDie(out body_direction_navball_texture_,
                     "navball_body_direction.png");
    LoadTextureOrDie(out surface_navball_texture_, "navball_surface.png");
    LoadTextureOrDie(out target_navball_texture_, "navball_target.png");

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
                                 ObscenelyEarly);
    // TimingPre, -101 on the script execution order page.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.Precalc,
                                 Precalc);
    // Timing1, -99 on the script execution order page.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.Early,
                                 Early);
    // Timing2, -1 on the script execution order page.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.Earlyish,
                                 Earlyish);
    // Timing3, 7.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.FashionablyLate,
                                 FashionablyLate);
    // TimingFI, 9.  Note that we cannot call the callback FlightIntegrator as
    // that would collide with the type.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.FlightIntegrator,
                                 JaiFailliAttendre);
    // Timing4, 19.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.Late, Late);
    // Timing5, 8008.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.BetterLateThanNever,
                                 BetterLateThanNever);
  }

  public override void OnSave(ConfigNode node) {
    base.OnSave(node);
    if (PluginRunning()) {
      String serialization;
      IntPtr serializer = IntPtr.Zero;
      for (;;) {
        serialization = plugin_.SerializePlugin(
                            ref serializer,
                            serialization_compression_,
                            serialization_encoding_);
        if (serialization == null) {
          break;
        }
        node.AddValue(principia_serialized_plugin_, serialization);
      }
    }
  }

  public override void OnLoad(ConfigNode node) {
    base.OnLoad(node);
    if (is_bad_installation_) {
      return;
    }
    if (must_record_journal_) {
      journaling_ = true;
      Log.ActivateRecorder(true);
    }
    if (node.HasValue(principia_serialized_plugin_)) {
      Cleanup();
      RemoveBuggyTidalLocking();
      Log.SetBufferedLogging(buffered_logging_);
      Log.SetSuppressedLogging(suppressed_logging_);
      Log.SetStderrLogging(stderr_logging_);
      Log.SetVerboseLogging(verbose_logging_);

      IntPtr deserializer = IntPtr.Zero;
      String[] serializations = node.GetValues(principia_serialized_plugin_);
      Log.Info("Serialization has " + serializations.Length + " chunks");
      foreach (String serialization in serializations) {
        Interface.DeserializePlugin(serialization,
                                    serialization.Length,
                                    ref deserializer,
                                    ref plugin_,
                                    serialization_compression_,
                                    serialization_encoding_);
      }
      Interface.DeserializePlugin("",
                                  0,
                                  ref deserializer,
                                  ref plugin_,
                                  serialization_compression_,
                                  serialization_encoding_);
      if (serialization_compression_ == "") {
        serialization_compression_ = "gipfeli";
      }
      if (serialization_encoding_ == "hexadecimal") {
        serialization_encoding_ = "base64";
      }

      plotting_frame_selector_.reset(
          new ReferenceFrameSelector(this, 
                                     plugin_,
                                     UpdateRenderingFrame,
                                     "Plotting frame"));
      previous_display_mode_ = null;
      must_set_plotting_frame_ = true;
      flight_planner_.reset(new FlightPlanner(this, plugin_));

      plugin_construction_ = DateTime.Now;
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
    if (is_bad_installation_) {
      bad_installation_dialog_.Show();
      return;
    }

    if (is_post_apocalyptic_) {
      apocalypse_dialog_.Show();
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
      toolbar_button_?.SetTrue(makeCall : false);
    } else {
      toolbar_button_?.SetFalse(makeCall : false);
    }

    if (hide_all_gui_) {
      WindowUtilities.ClearLock(this);
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
      WindowUtilities.EnsureOnScreen(ref main_window_rectangle_);
      main_window_x_ = (int)main_window_rectangle_.xMin;
      main_window_y_ = (int)main_window_rectangle_.yMin;
      main_window_rectangle_.InputLock(this);

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

    // Handle clicks on planets.
    if (MapView.MapIsEnabled) {
      HandleMapViewClicks();
    }

    override_rsas_target_ = false;
    Vessel active_vessel = FlightGlobals.ActiveVessel;
    if (active_vessel != null) {
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

      if (!PluginRunning()) {
        return;
      }

      var target_vessel = FlightGlobals.fetch.VesselTarget?.GetVessel();
      if (FlightGlobals.speedDisplayMode ==
              FlightGlobals.SpeedDisplayModes.Target &&
          target_vessel != null &&
          plugin_.HasVessel(target_vessel.id.ToString())) {
        plugin_.SetTargetVessel(target_vessel.id.ToString(),
                                plotting_frame_selector_.get()
                                    .selected_celestial.flightGlobalsIndex);
        if (plotting_frame_selector_.get().target_override != target_vessel) {
          navball_changed_ = true;
          plotting_frame_selector_.get().target_override = target_vessel;
        }
      } else {
        plugin_.ClearTargetVessel();
        if (plotting_frame_selector_.get().target_override != null) {
          navball_changed_ = true;
          plotting_frame_selector_.get().target_override = null;
        }
      }

      // Orient the ball.
      navball_.navBall.rotation =
          (UnityEngine.QuaternionD)navball_.attitudeGymbal *  // sic.
          (UnityEngine.QuaternionD)plugin_.NavballOrientation(
              (XYZ)Planetarium.fetch.Sun.position,
              (XYZ)(Vector3d)active_vessel.ReferenceTransform.position);

      if (previous_display_mode_ != FlightGlobals.speedDisplayMode) {
        navball_changed_ = true;
        previous_display_mode_ = FlightGlobals.speedDisplayMode;
        // The navball speed display mode was changed, change the reference
        // frame accordingly.
        switch (FlightGlobals.speedDisplayMode) {
          case FlightGlobals.SpeedDisplayModes.Surface:
            plotting_frame_selector_.get().SetFrameType(
                ReferenceFrameSelector.FrameType.BODY_SURFACE);
            break;
          case FlightGlobals.SpeedDisplayModes.Orbit:
            plotting_frame_selector_.get().SetFrameType(
                last_non_surface_frame_type_);
            break;
        }
      }

      if (navball_changed_ && previous_display_mode_ != null) {
        // Texture the ball.
        navball_changed_ = false;
        if (plotting_frame_selector_.get().target_override) {
          set_navball_texture(target_navball_texture_);
        } else {
          // If we are targeting an unmanageable vessel, keep the navball in
          // target mode; otherwise, put it in the mode that reflects the
          // plotting frame.
          if (FlightGlobals.speedDisplayMode !=
              FlightGlobals.SpeedDisplayModes.Target) {
            if (plotting_frame_selector_.get().frame_type ==
                    ReferenceFrameSelector.FrameType.BODY_SURFACE) {
              if (FlightGlobals.speedDisplayMode !=
                  FlightGlobals.SpeedDisplayModes.Surface) {
                FlightGlobals.SetSpeedMode(
                    FlightGlobals.SpeedDisplayModes.Surface);
              }
            } else {
              if (FlightGlobals.speedDisplayMode !=
                  FlightGlobals.SpeedDisplayModes.Orbit) {
                FlightGlobals.SetSpeedMode(
                    FlightGlobals.SpeedDisplayModes.Orbit);
              }
            }
          }
          switch (plotting_frame_selector_.get().frame_type) {
            case ReferenceFrameSelector.FrameType.BODY_SURFACE:
              set_navball_texture(surface_navball_texture_);
              break;
            case ReferenceFrameSelector.FrameType.BODY_CENTRED_NON_ROTATING:
              set_navball_texture(inertial_navball_texture_);
              break;
            case ReferenceFrameSelector.FrameType.BARYCENTRIC_ROTATING:
              set_navball_texture(barycentric_navball_texture_);
              break;
            case ReferenceFrameSelector.FrameType.BODY_CENTRED_PARENT_DIRECTION:
              set_navball_texture(body_direction_navball_texture_);
              break;
          }
        }
      }

      // Design for compatibility with FAR: if we are in surface mode in an
      // atmosphere, FAR gives options to display the IAS, EAS, and Mach number,
      // in any of a number of units.  In that case, we must force the reference
      // frame to be that of the surface *of the relevant body*, in order to be
      // consistent with the speed display.
      bool ferram_owns_the_speed_display =
          FlightGlobals.speedDisplayMode ==
              FlightGlobals.SpeedDisplayModes.Surface &&
          active_vessel.atmDensity > 0;
      if (ferram_owns_the_speed_display) {
        plotting_frame_selector_.get().SetToSurfaceFrameOf(
            active_vessel.mainBody);
      }

      if (plotting_frame_selector_.get().target_override == null &&
          FlightGlobals.speedDisplayMode ==
              FlightGlobals.SpeedDisplayModes.Target) {
        KSP.UI.Screens.Flight.SpeedDisplay.Instance.textTitle.text = "Target";
      }

      if (FlightGlobals.speedDisplayMode ==
              FlightGlobals.SpeedDisplayModes.Orbit ||
          FlightGlobals.speedDisplayMode ==
              FlightGlobals.SpeedDisplayModes.Surface ||
          plotting_frame_selector_.get().target_override) {
        bool plugin_has_active_manageable_vessel =
            has_active_manageable_vessel() &&
            plugin_.HasVessel(active_vessel.id.ToString());

        KSP.UI.Screens.Flight.SpeedDisplay speed_display =
            KSP.UI.Screens.Flight.SpeedDisplay.Instance;
        if (speed_display?.textTitle != null &&
            speed_display?.textSpeed != null &&
            !ferram_owns_the_speed_display) {
          speed_display.textTitle.text =
              plotting_frame_selector_.get().ShortName();
          var active_vessel_velocity =
              plugin_has_active_manageable_vessel
                  ? (Vector3d)plugin_.VesselVelocity(
                        active_vessel.id.ToString())
                  : (Vector3d)plugin_.UnmanageableVesselVelocity(
                        new QP{q = (XYZ)active_vessel.orbit.pos,
                               p = (XYZ)active_vessel.orbit.vel},
                        active_vessel.orbit.referenceBody.flightGlobalsIndex);
          speed_display.textSpeed.text =
              active_vessel_velocity.magnitude.ToString("F1") + "m/s";
        }

        if (!plugin_has_active_manageable_vessel) {
          // TODO(egg): orient the Frenet trihedron even in the case where the
          // active vessel is unmanageable, similarly to the speed display
          // above.
          return;
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
        if (!active_vessel.OnPreAutopilotUpdate.GetInvocationList().Contains(
                (FlightInputCallback)OverrideRSASTarget)) {
          Log.Info("Prepending RSAS override");
          active_vessel.OnPreAutopilotUpdate += OverrideRSASTarget;
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
            // NOTE(egg): For reasons that are unlikely to become clear again,
            // the button labeled with the radial in icon sets the autopilot
            // mode to |RadialOut|, and vice-versa.  As a result, we must set
            // the target to the outwards radial (negative normal) vector if the
            // mode is |RadialIn|.  Contrast with the navball vectors above,
            // which do not exhibit this inconsistency (thus where
            // |radialInVector| is set to |radial|).
            case VesselAutopilot.AutopilotMode.RadialIn:
              rsas_target_ = -radial;
              break;
            case VesselAutopilot.AutopilotMode.RadialOut:
              rsas_target_ = radial;
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
    if (is_bad_installation_) {
      return;
    }
    if (GameSettings.ORBIT_WARP_DOWN_AT_SOI) {
      Log.Info("Setting GameSettings.ORBIT_WARP_DOWN_AT_SOI to false");
      GameSettings.ORBIT_WARP_DOWN_AT_SOI = false;
    }
    if (must_set_plotting_frame_ && FlightGlobals.currentMainBody != null) {
      must_set_plotting_frame_ = false;
      plotting_frame_selector_.reset(new ReferenceFrameSelector(
          this, plugin_, UpdateRenderingFrame, "Plotting frame"));
      previous_display_mode_ = null;
    }

    if (PluginRunning()) {
      plugin_.SetMainBody(
          (FlightGlobals.currentMainBody
               ?? FlightGlobals.GetHomeBody()).flightGlobalsIndex);

      plugin_.ForgetAllHistoriesBefore(plugin_.CurrentTime() -
                                       history_lengths_[history_length_index_]);
      // TODO(egg): Set the degrees of freedom of the origin of |World| (by
      // toying with Krakensbane and FloatingOrigin) here.

      // Now we let the game and Unity do their thing.  Among other things,
      // the FashionablyLate callbacks, including ReportNonConservativeForces,
      // then the FlightIntegrator's FixedUpdate will run, then the Vessel's,
      // and eventually the physics simulation.
      StartCoroutine(WaitedForFixedUpdate());
    }
  }

  private void OnDisable() {
    if (is_bad_installation_) {
      return;
    }
    Log.Info("principia.ksp_plugin_adapter.PrincipiaPluginAdapter.OnDisable()");
    if (toolbar_button_ != null) {
      KSP.UI.Screens.ApplicationLauncher.Instance.RemoveModApplication(
          toolbar_button_);
    }
    WindowUtilities.ClearLock(this);
    Cleanup();
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.ObscenelyEarly,
                                    ObscenelyEarly);
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.Precalc,
                                    Precalc);
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.Early,
                                    Early);
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.Earlyish,
                                    Earlyish);
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.FashionablyLate,
                                    FashionablyLate);
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.FlightIntegrator,
                                    JaiFailliAttendre);
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.Late, Late);
    TimingManager.FixedUpdateRemove(
        TimingManager.TimingStage.BetterLateThanNever,
        BetterLateThanNever);
  }

  #endregion

  private System.Collections.IEnumerator WaitedForFixedUpdate() {
    yield return new UnityEngine.WaitForFixedUpdate();

  try {
    // Unity's physics has just finished doing its thing.  If we correct the
    // positions here, nobody will know that they're not the ones obtained by
    // Unity.

    if (!time_is_advancing_) {
      yield break;
    }

    double Δt = Planetarium.TimeScale * Planetarium.fetch.fixedDeltaTime;

    QP main_body_degrees_of_freedom =
        new QP{q = (XYZ)(FlightGlobals.currentMainBody ??
                             FlightGlobals.GetHomeBody()).position,
               p = (XYZ)(-krakensbane.FrameVel)};

    // NOTE(egg): Inserting vessels and parts has to occur in
    // |WaitForFixedUpdate|, since some may be destroyed (by collisions) during
    // the physics step.  See also #1281.
    foreach (Vessel vessel in FlightGlobals.Vessels) {
      string unmanageability_reasons = UnmanageabilityReasons(vessel);
      if (unmanageability_reasons != null) {
        if (plugin_.HasVessel(vessel.id.ToString())) {
          Log.Info("Removing vessel " + vessel.vesselName + "(" + vessel.id +
                   ")" + " because " + unmanageability_reasons);
        }
        continue;
      }

      bool inserted;
      plugin_.InsertOrKeepVessel(vessel.id.ToString(),
                                 vessel.vesselName,
                                 vessel.mainBody.flightGlobalsIndex,
                                 !vessel.packed,
                                 out inserted);
      if (!vessel.packed) {
        foreach (Part part in vessel.parts.Where((part) => part.rb != null)) {
          QP degrees_of_freedom;
          if (part_id_to_degrees_of_freedom_.ContainsKey(part.flightID)) {
            degrees_of_freedom = part_id_to_degrees_of_freedom_[part.flightID];
          } else {
            // Assumptions about the invariants of KSP will invariably fail.
            // This is a graceful fallback.
            Log.Error("Unpacked part " + part.name + " (" +
                      part.flightID.ToString("X") +
                      ") appeared between BetterLateThanNever and " +
                      "WaitForFixedUpdate.  Linearly extrapolating its " +
                      "position at the previous frame.");
            degrees_of_freedom =
                new QP{q = (XYZ)((Vector3d)part.rb.position -
                                  Δt * (Vector3d)part.rb.velocity),
                       p = (XYZ)(Vector3d)part.rb.velocity};
          }
          // In the first few frames after spawning a Kerbal, its physicsMass is
          // 0; we use its rb.mass instead.
          // NOTE(egg): the physics engine does not move the celestials, so it
          // is fine to use |main_body_degrees_of_freedom| here rather than to
          // store it during |FixedUpdate| or one of its timings.
          plugin_.InsertOrKeepLoadedPart(
              part.flightID,
              part.name,
              part.physicsMass == 0 ? part.rb.mass : part.physicsMass,
              vessel.id.ToString(),
              vessel.mainBody.flightGlobalsIndex,
              main_body_degrees_of_freedom,
              degrees_of_freedom,
              Δt);
          if (part_id_to_intrinsic_force_.ContainsKey(part.flightID)) {
            // When a Kerbal is doing an EVA and holding on to a ladder, the
            // ladder imbues them with their weight at the location of the
            // vessel to which the ladder is attached.  This leads to fantastic
            // effects where doing an EVA accelerates the vessel, see #1415.
            // Just say no to stupidity.
            if (!(vessel.isEVA && vessel.evaController.OnALadder)) {
              plugin_.IncrementPartIntrinsicForce(
                  part.flightID,
                  (XYZ)part_id_to_intrinsic_force_[part.flightID]);
            }
          }
          if (part_id_to_intrinsic_forces_.ContainsKey(part.flightID)) {
            foreach (
                var force in part_id_to_intrinsic_forces_[part.flightID]) {
              plugin_.IncrementPartIntrinsicForce(part.flightID,
                                                  (XYZ)force.force);
            }
          }
        }
      } else if (inserted) {
        var parts = vessel.protoVessel.protoPartSnapshots;
        // For reasons that are unclear, the asteroid spawning code sometimes
        // generates the same flightID twice; we regenerate the flightID on
        // any asteroid we find.
        if (vessel.vesselType == VesselType.SpaceObject &&
            parts.Count == 1 &&
            parts[0].partName == "PotatoRoid" &&
            parts[0].flightID == parts[0].missionID) {
          var part = parts[0];
          var old_id = part.flightID;
          part.flightID = ShipConstruction.GetUniqueFlightID(
              HighLogic.CurrentGame.flightState);
          Log.Info("Regenerating the part ID of " + vessel.name + ": " +
                   part.flightID + " (was " + old_id + ")");
        }
        foreach (ProtoPartSnapshot part in parts) {
          plugin_.InsertUnloadedPart(
              part.flightID,
              part.partName,
              vessel.id.ToString(),
              new QP{q = (XYZ)vessel.orbit.pos, p = (XYZ)vessel.orbit.vel});
        }
      }
    }

    plugin_.PrepareToReportCollisions();

    // The collisions are reported and stored into |currentCollisions| in
    // OnCollisionEnter|Stay|Exit, which occurred while we yielded.
    // Here, the |currentCollisions| are the collisions that occurred in the
    // physics simulation, which is why we report them before calling
    // |AdvanceTime|.
    foreach (Vessel vessel1 in
             FlightGlobals.Vessels.Where(v => !v.packed && is_manageable(v))) {
      if (plugin_.HasVessel(vessel1.id.ToString())) {
        if (vessel1.isEVA && (vessel1.evaController.OnALadder ||
                              is_clambering(vessel1.evaController))) {
          var vessel2 = vessel1.evaController.LadderPart?.vessel;
          if (vessel2 != null && !vessel2.packed && is_manageable(vessel2)) {
            plugin_.ReportPartCollision(
                vessel1.rootPart.flightID,
                closest_physical_parent(
                    vessel1.evaController.LadderPart).flightID);
          } else {
            // vessel1 is either clambering (quite likely on the ground or on a
            // grounded vessel), or climbing a ladder on an unmanaged vessel.
            Log.Info("Reporting " +
                     (is_clambering(vessel1.evaController)
                          ? "clambering"
                          : "climbing an unmanaged ladder"));
            plugin_.ReportGroundCollision(vessel1.rootPart.flightID);
          }
        }
        foreach (Part part1 in vessel1.parts) {
          if (part1.Modules.OfType<ModuleWheelBase>()
                           .Where(wheel => wheel.isGrounded)
                           .Any()) {
            Log.Info("Reporting grounded wheel");
            plugin_.ReportGroundCollision(
                closest_physical_parent(part1).flightID);
          }
          foreach (var collider in part1.currentCollisions) {
            if (collider == null) {
              // This happens, albeit quite rarely, see #1447.  When it happens,
              // the null collider remains in |currentCollisions| until the next
              // scene change, so we do not log, otherwise we would spam.
              continue;
            }
            if (collider.gameObject.layer == (int)UnityLayers.LocalScenery) {
              Log.Info("Reporting collision with local scenery " +
                       collider.name);
              plugin_.ReportGroundCollision(
                  closest_physical_parent(part1).flightID);
            }
            var part2 = collider.gameObject.GetComponentUpwards<Part>();
            var vessel2 = part2?.vessel;
            if (vessel2 == vessel1) {
              // All parts in a vessel are in the same pile up, so there is no
              // point in reporting this collision; this also causes issues
              // where disappearing kerbals collide with themselves.
              continue;
            }
            if (part1.State == PartStates.DEAD ||
                part2?.State == PartStates.DEAD) {
              continue;
            }
            if (vessel2 != null) {
              if (is_manageable(vessel2)) {
                plugin_.ReportPartCollision(
                    closest_physical_parent(part1).flightID,
                    closest_physical_parent(part2).flightID);
              } else {
                Log.Info("Reporting collision with the unmanageable vessel " +
                         vessel2.vesselName);
                plugin_.ReportGroundCollision(
                    closest_physical_parent(part1).flightID);
              }
            }
          }
        }
      }
    }

    plugin_.FreeVesselsAndPartsAndCollectPileUps(Δt);

    foreach (Vessel vessel in FlightGlobals.Vessels.Where(v => !v.packed)) {
      if (!plugin_.HasVessel(vessel.id.ToString())) {
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
                   p = (XYZ)(Vector3d)part.rb.velocity},
            main_body_degrees_of_freedom);
      }
    }

    // Advance the lagging vessels and kill those which collided with a
    // celestial.
    {
      DisposableIterator collided_vessels;
      plugin_.CatchUpLaggingVessels(out collided_vessels);
      for (; !collided_vessels.IteratorAtEnd();
            collided_vessels.IteratorIncrement()) {
        Guid vessel_guid = new Guid(collided_vessels.IteratorGetVesselGuid());
        Vessel vessel = FlightGlobals.FindVessel(vessel_guid);
        vessel?.Die();
      }
    }

    UpdatePredictions();

    // We don't want to do too many things here, since all the KSP classes
    // still think they're in the preceding step.  We only nudge the Unity
    // transforms of loaded vessels & their parts.
    if (has_active_manageable_vessel() && !FlightGlobals.ActiveVessel.packed &&
        plugin_.HasVessel(FlightGlobals.ActiveVessel.id.ToString())) {
      Vector3d q_correction_at_root_part = Vector3d.zero;
      Vector3d v_correction_at_root_part = Vector3d.zero;
      foreach (Vessel vessel in FlightGlobals.Vessels.Where(v => !v.packed)) {
        // TODO(egg): if I understand anything, there should probably be a
        // special treatment for loaded packed vessels.  I don't understand
        // anything though.
        if (!plugin_.HasVessel(vessel.id.ToString())) {
          continue;
        }
        foreach (Part part in vessel.parts) {
          if (part.rb == null) {
            continue;
          }
          QP part_actual_degrees_of_freedom =
              plugin_.GetPartActualDegreesOfFreedom(
                  part.flightID,
          new Origin{reference_part_is_at_origin  =
                         FloatingOrigin.fetch.continuous,
                     reference_part_is_unmoving =
                         krakensbane_.FrameVel != Vector3d.zero,
                     main_body_centre_in_world =
                         (XYZ)FlightGlobals.ActiveVessel.mainBody.position,
                     reference_part_id =
                         FlightGlobals.ActiveVessel.rootPart.flightID});
          if (part == FlightGlobals.ActiveVessel.rootPart) {
            q_correction_at_root_part =
                (Vector3d)part_actual_degrees_of_freedom.q - part.rb.position;
            v_correction_at_root_part =
                (Vector3d)part_actual_degrees_of_freedom.p - part.rb.velocity;
          }

          // TODO(egg): use the centre of mass.  Here it's a bit tedious, some
          // transform nonsense must probably be done.
          part.rb.position = (Vector3d)part_actual_degrees_of_freedom.q;
          part.rb.transform.position =
              (Vector3d)part_actual_degrees_of_freedom.q;
          part.rb.velocity = (Vector3d)part_actual_degrees_of_freedom.p;
        }
      }
      foreach (
          physicalObject physical_object in FlightGlobals.physicalObjects.Where(
              o => o != null && o.rb != null)) {
        physical_object.rb.position += q_correction_at_root_part;
        physical_object.rb.transform.position += q_correction_at_root_part;
        physical_object.rb.velocity += v_correction_at_root_part;
      }
      QP main_body_dof = plugin_.CelestialWorldDegreesOfFreedom(
          FlightGlobals.ActiveVessel.mainBody.flightGlobalsIndex,
          new Origin{reference_part_is_at_origin  =
                         FloatingOrigin.fetch.continuous,
                     reference_part_is_unmoving =
                         krakensbane_.FrameVel != Vector3d.zero,
                     main_body_centre_in_world =
                         (XYZ)FlightGlobals.ActiveVessel.mainBody.position,
                     reference_part_id =
                         FlightGlobals.ActiveVessel.rootPart.flightID},
          plugin_.CurrentTime());
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
          Vessel vessel in FlightGlobals.Vessels.Where(is_manageable_on_rails)) {
        vessel.SetPosition(vessel.transform.position + offset);
      }
      // NOTE(egg): this is almost certainly incorrect, since we give the
      // bodies their positions at the next instant, wherease KSP still expects
      // them at the previous instant, and will propagate them at the beginning
      // of the next frame...
    }
  } catch (Exception e) { Log.Fatal(e.ToString()); }
  }

  private void ObscenelyEarly() {
    foreach (var vessel in
             FlightGlobals.Vessels.Where(vessel => vessel.precalc != null)) {
      vessel.precalc.enabled = false;
    }
  }

  private void Precalc() {
    if (FlightGlobals.ActiveVessel?.situation == Vessel.Situations.PRELAUNCH &&
        FlightGlobals.ActiveVessel?.orbitDriver?.lastMode ==
            OrbitDriver.UpdateMode.TRACK_Phys &&
        FlightGlobals.ActiveVessel?.orbitDriver?.updateMode ==
            OrbitDriver.UpdateMode.IDLE) {
      Log.Info("Skipping AdvanceTime and SetBodyFrames while waiting for the " +
               "vessel to be fully ready (see #1421).");
    } else {
      double universal_time = Planetarium.GetUniversalTime();
      time_is_advancing_ = time_is_advancing(universal_time);
      if (time_is_advancing_) {
        plugin_.AdvanceTime(universal_time, Planetarium.InverseRotAngle);
        if (!is_post_apocalyptic_) {
          String revelation = "";
          if (plugin_.HasEncounteredApocalypse(out revelation)) {
            is_post_apocalyptic_ = true;
            apocalypse_dialog_.Message = revelation;
          }
        }
        foreach (var vessel in FlightGlobals.Vessels) {
          if (vessel.packed && plugin_.HasVessel(vessel.id.ToString())) {
            vessel_futures_.Add(
                plugin_.FutureCatchUpVessel(vessel.id.ToString()));
          }
        }
      }
      SetBodyFrames();
    }
    // Unfortunately there is no way to get scheduled between Planetarium and
    // VesselPrecalculate, so we get scheduled after VesselPrecalculate, set the
    // body frames for our weird tilt, and run VesselPrecalculate manually.
    // Sob.
    // NOTE(egg): we cannot use foreach here, and we must iterate downwards,
    // since vessel.precalc.FixedUpdate may remove its vessel.
    for (int i = FlightGlobals.Vessels.Count - 1; i >= 0; --i) {
      var vessel = FlightGlobals.Vessels[i];
      if (vessel.precalc == null) {
        continue;
      }
      vessel.precalc.enabled = true;
      // In stock this is equivalent to |FixedUpdate()|.  With
      // ModularFlightIntegrator's ModularVesselPrecalculate, which gets run
      // in TimingPre like us, this comes with a flag that ensures it only gets
      // run once.
      vessel.precalc.MainPhysics(true);
    }
  }

  private void Early() {
    if (PluginRunning()) {
      // Wait for all the asynchronous integrations to complete and kill the
      // vessels that collided with a celestial.  Note that a given vessel may
      // be returned by several calls to FutureWaitForVesselToCatchUp because
      // what we integrate are really pile-ups.
      var all_collided_vessels = new HashSet<Vessel>();
      foreach (var f in vessel_futures_) {
        var future = f;
        DisposableIterator collided_vessels;
        plugin_.FutureWaitForVesselToCatchUp(ref future,
                                             out collided_vessels);
        for (; !collided_vessels.IteratorAtEnd();
             collided_vessels.IteratorIncrement()) {
          Guid vessel_guid =
              new Guid(collided_vessels.IteratorGetVesselGuid());
          Vessel vessel = FlightGlobals.FindVessel(vessel_guid);
          all_collided_vessels.Add(vessel);
        }
      }
      vessel_futures_.Clear();
      foreach (var vessel in all_collided_vessels) {
        vessel?.Die();
      }
      ApplyToVesselsOnRails(
          vessel => UpdateVessel(vessel, Planetarium.GetUniversalTime()));
    }
  }

  private void Earlyish() {}

  private void FashionablyLate() {
    // We fetch the forces from the census of nonconservatives here;
    // part.forces, part.force, and part.torque are cleared by the
    // FlightIntegrator's FixedUpdate (while we are yielding).
    if (PluginRunning()) {
      if (has_active_manageable_vessel() && FlightGlobals.ActiveVessel.packed) {
        if (PhysicsGlobals.GraviticForceMultiplier != 0) {  // sic.
          Log.Info("Killing stock gravity");
          PhysicsGlobals.GraviticForceMultiplier = 0;
        }
      } else if (PhysicsGlobals.GraviticForceMultiplier == 0) {
        Log.Info("Reinstating stock gravity");
        PhysicsGlobals.GraviticForceMultiplier = 1;
      }
      part_id_to_intrinsic_force_.Clear();
      part_id_to_intrinsic_forces_.Clear();
      foreach (Vessel vessel in
               FlightGlobals.Vessels.Where(v => is_manageable(v) &&
                                                !v.packed)) {
        foreach (Part part in vessel.parts.Where((part) => part.rb != null)) {
          if (part.force != Vector3d.zero) {
            part_id_to_intrinsic_force_.Add(part.flightID, part.force);
          }
          if (part.forces.Count > 0) {
            part_id_to_intrinsic_forces_.Add(part.flightID,
                                             part.forces.ToArray());
          }
        }
      }
    }
  }

  private void JaiFailliAttendre() {
    // We fetch the forces from stock aerodynamics, which does not use
    // |Part.AddForce| etc.
    if (PluginRunning()) {
      foreach (Vessel vessel in
               FlightGlobals.Vessels.Where(v => is_manageable(v) &&
                                                !v.packed)) {
        foreach (Part part in vessel.parts) {
          Part physical_parent = closest_physical_parent(part);
          if (part.bodyLiftLocalVector != UnityEngine.Vector3.zero ||
              part.dragVector != UnityEngine.Vector3.zero) {
            if (part_id_to_intrinsic_forces_.ContainsKey(
                    physical_parent.flightID)) {
              var previous_holder =
                  part_id_to_intrinsic_forces_[physical_parent.flightID];
              part_id_to_intrinsic_forces_[physical_parent.flightID] =
                  new Part.ForceHolder[previous_holder.Length + 2];
              previous_holder.CopyTo(
                  part_id_to_intrinsic_forces_[physical_parent.flightID],
                  0);
            } else {
              part_id_to_intrinsic_forces_.Add(physical_parent.flightID,
                                               new Part.ForceHolder[2]);
            }
            int lift_index = part_id_to_intrinsic_forces_[
                                 physical_parent.flightID].Length - 2;
            int drag_index = lift_index + 1;
            part_id_to_intrinsic_forces_[physical_parent.flightID][lift_index] =
                new Part.ForceHolder {
                  force = part.partTransform.TransformDirection(
                              part.bodyLiftLocalVector),
                  pos = part.partTransform.TransformPoint(
                            part.bodyLiftLocalPosition)};
            part_id_to_intrinsic_forces_[physical_parent.flightID][drag_index] =
                new Part.ForceHolder {
                  force = -part.dragVectorDir * part.dragScalar,
                  pos = (physical_parent != part &&
                         PhysicsGlobals.ApplyDragToNonPhysicsPartsAtParentCoM)
                            ? physical_parent.rb.worldCenterOfMass
                            : part.partTransform.TransformPoint(
                                  part.CoPOffset)};
          }
        }
      }
    }
  }

  private void Late() {
    if (PluginRunning()) {
      if (FlightGlobals.currentMainBody?.inverseRotation == true) {
        plugin_.SetWorldRotationalReferenceFrame(
            FlightGlobals.currentMainBody.flightGlobalsIndex);
      } else {
        plugin_.ClearWorldRotationalReferenceFrame();
      }
    }
  }

  private void BetterLateThanNever() {
    if (PluginRunning()) {
      part_id_to_degrees_of_freedom_.Clear();
      foreach (Vessel vessel in
               FlightGlobals.Vessels.Where(v => is_manageable(v) &&
                                                !v.packed)) {
        foreach (Part part in vessel.parts.Where((part) => part.rb != null)) {
          // TODO(egg): use the centre of mass.
          part_id_to_degrees_of_freedom_.Add(
              part.flightID,
              new QP{q = (XYZ)(Vector3d)part.rb.position,
                     p = (XYZ)(Vector3d)part.rb.velocity});
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
    if (vessel.patchedConicRenderer != null) {
      vessel.patchedConicRenderer.relativityMode =
          PatchRendering.RelativityMode.RELATIVE;
    }

    if (display_patched_conics_ || !is_manageable(vessel)) {
      vessel.orbitDriver.Renderer.drawMode =
          vessel.PatchedConicsAttached
              ? OrbitRenderer.DrawMode.OFF
              : OrbitRenderer.DrawMode.REDRAW_AND_RECALCULATE;
      vessel.orbitDriver.Renderer.drawIcons =
          ((vessel.isActiveVessel ||
            vessel == FlightGlobals.ActiveVessel?.targetObject as Vessel) &&
           !vessel.PatchedConicsAttached)
              ? OrbitRenderer.DrawIcons.ALL
              : OrbitRenderer.DrawIcons.OBJ;

      if (vessel.patchedConicRenderer != null &&
          !vessel.patchedConicRenderer.enabled) {
        vessel.patchedConicRenderer.enabled = true;
      }
    } else {
      vessel.orbitDriver.Renderer.drawMode = OrbitRenderer.DrawMode.OFF;
      vessel.orbitDriver.Renderer.drawIcons = OrbitRenderer.DrawIcons.OBJ;

      // vessel.patchedConicRenderer may be null when in career mode, if
      // patched conics have yet to be unlocked.  In that case there is
      // nothing to remove.  Note that Principia doesn't care about unlocks
      // and always displays all the information it can.
      if (vessel.patchedConicRenderer != null) {
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

      if (vessel.orbitTargeter != null) {
        vessel.orbitTargeter.enabled = false;
      }
    }
  }

  private void OnCelestialNodeClick(KSP.UI.Screens.Mapview.MapNode node,
                                    Mouse.Buttons buttons) {
    if (buttons == Mouse.Buttons.Left) {
      if (selecting_target_celestial_) {
        FlightGlobals.fetch.SetVesselTarget(node.mapObject.celestialBody);
        selecting_target_celestial_ = false;
      } else if (PlanetariumCamera.fetch.target != node.mapObject) {
        PlanetariumCamera.fetch.SetTarget(node.mapObject);
      }
    }
  }

  private void OnVesselNodeClick(KSP.UI.Screens.Mapview.MapNode node,
                                 Mouse.Buttons buttons) {
    if (selecting_active_vessel_target_) {
      FlightGlobals.fetch.SetVesselTarget(node.mapObject.vessel);
      selecting_active_vessel_target_ = false;
    } else if (buttons == Mouse.Buttons.Left &&
               PlanetariumCamera.fetch.target != node.mapObject) {
      PlanetariumCamera.fetch.SetTarget(node.mapObject);
    }
  }

  private void HandleMapViewClicks() {
    if (InputLockManager.IsUnlocked(ControlTypes.MAP_UI) &&
        !UnityEngine.EventSystems.EventSystem.current
             .IsPointerOverGameObject() &&
        Mouse.Left.GetClick() && !ManeuverGizmo.HasMouseFocus &&
        !selecting_active_vessel_target_) {
      var ray = PlanetariumCamera.Camera.ScreenPointToRay(
          UnityEngine.Input.mousePosition);
      foreach (var celestial in FlightGlobals.Bodies) {
        double scaled_distance =
            Vector3d.Cross(ray.direction,
                           ScaledSpace.LocalToScaledSpace(celestial.position) -
                               ray.origin).magnitude;
        if (scaled_distance * ScaledSpace.ScaleFactor < celestial.Radius) {
          if (selecting_target_celestial_) {
            FlightGlobals.fetch.SetVesselTarget(celestial);
            selecting_target_celestial_ = false;
          } else if (PlanetariumCamera.fetch.target != celestial.MapObject) {
            PlanetariumCamera.fetch.SetTarget(celestial.MapObject);
          }
        }
      }
    }
  }

  private void RenderTrajectories() {
    if (!PluginRunning()) {
      return;
    }
    foreach (var celestial in FlightGlobals.Bodies.Where(
                                  c => c.MapObject?.uiNode != null)) {
      celestial.MapObject.uiNode.OnClick -= OnCelestialNodeClick;
      celestial.MapObject.uiNode.OnClick += OnCelestialNodeClick;
    }
    foreach (var vessel in FlightGlobals.Vessels.Where(
                               v => v.mapObject?.uiNode != null)) {
      // There is no way to check if we have already added a callback to an
      // event...
      vessel.mapObject.uiNode.OnClick -= OnVesselNodeClick;
      vessel.mapObject.uiNode.OnClick += OnVesselNodeClick;
      RemoveStockTrajectoriesIfNeeded(vessel);
    }
    Vessel main_vessel = FlightGlobals.ActiveVessel ??
                         space_tracking?.SelectedVessel;
    if (main_vessel == null) {
      return;
    }
    string main_vessel_guid = main_vessel.id.ToString();
    bool ready_to_draw_active_vessel_trajectory =
        MapView.MapIsEnabled &&
        plugin_.HasVessel(main_vessel_guid);
    if (ready_to_draw_active_vessel_trajectory) {
      XYZ sun_world_position = (XYZ)Planetarium.fetch.Sun.position;
      using (DisposablePlanetarium planetarium =
                GLLines.NewPlanetarium(plugin_, sun_world_position)) {
        GLLines.Draw(() => {
          using (DisposableIterator rp2_lines_iterator =
                    planetarium.PlanetariumPlotPsychohistory(
                        plugin_,
                        чебышёв_plotting_method_,
                        main_vessel_guid)) {
            GLLines.PlotRP2Lines(rp2_lines_iterator,
                                 XKCDColors.Lime,
                                 GLLines.Style.FADED);
          }
          RenderPredictionMarkers(main_vessel_guid, sun_world_position);
          using (DisposableIterator rp2_lines_iterator =
                    planetarium.PlanetariumPlotPrediction(
                        plugin_,
                        чебышёв_plotting_method_,
                        main_vessel_guid)) {
            GLLines.PlotRP2Lines(rp2_lines_iterator,
                                 XKCDColors.Fuchsia,
                                 GLLines.Style.SOLID);
          }
          string target_id =
              FlightGlobals.fetch.VesselTarget?.GetVessel()?.id.ToString();
          if (FlightGlobals.ActiveVessel != null &&
              !plotting_frame_selector_.get().target_override &&
              target_id != null && plugin_.HasVessel(target_id)) {
            using (DisposableIterator rp2_lines_iterator =
                      planetarium.PlanetariumPlotPsychohistory(
                          plugin_,
                          чебышёв_plotting_method_,
                          target_id)) {
              GLLines.PlotRP2Lines(rp2_lines_iterator,
                                   XKCDColors.Goldenrod,
                                   GLLines.Style.FADED);
            }
            RenderPredictionMarkers(target_id, sun_world_position);
            using (DisposableIterator rp2_lines_iterator =
                      planetarium.PlanetariumPlotPrediction(
                          plugin_,
                          чебышёв_plotting_method_,
                          target_id)) {
              GLLines.PlotRP2Lines(rp2_lines_iterator,
                                   XKCDColors.LightMauve,
                                   GLLines.Style.SOLID);
            }
          }
          if (plugin_.FlightPlanExists(main_vessel_guid)) {
            RenderFlightPlanMarkers(main_vessel_guid, sun_world_position);

            int number_of_segments =
                plugin_.FlightPlanNumberOfSegments(main_vessel_guid);
            for (int i = 0; i < number_of_segments; ++i) {
              bool is_burn = i % 2 == 1;
              using (DisposableIterator rendered_segments =
                        plugin_.FlightPlanRenderedSegment(main_vessel_guid,
                                                          sun_world_position,
                                                          i)) {
                if (rendered_segments.IteratorAtEnd()) {
                  Log.Info("Skipping segment " + i);
                  continue;
                }
                Vector3d position_at_start =
                    (Vector3d)rendered_segments.
                        IteratorGetDiscreteTrajectoryXYZ();
                using (DisposableIterator rp2_lines_iterator =
                          planetarium.PlanetariumPlotFlightPlanSegment(
                              plugin_,
                              чебышёв_plotting_method_,
                              main_vessel_guid,
                              i)) {
                  GLLines.PlotRP2Lines(
                      rp2_lines_iterator,
                      is_burn ? XKCDColors.Pink : XKCDColors.PeriwinkleBlue,
                      is_burn ? GLLines.Style.SOLID : GLLines.Style.DASHED);
                }
                if (is_burn) {
                  int manoeuvre_index = i / 2;
                  NavigationManoeuvreFrenetTrihedron manoeuvre =
                      plugin_.FlightPlanGetManoeuvreFrenetTrihedron(
                          main_vessel_guid,
                          manoeuvre_index);
                  double scale = (ScaledSpace.ScaledToLocalSpace(
                                      MapView.MapCamera.transform.position) -
                                  position_at_start).magnitude * 0.015;
                  Action<XYZ, UnityEngine.Color> add_vector =
                      (world_direction, colour) => {
                        UnityEngine.GL.Color(colour);
                        GLLines.AddSegment(
                            position_at_start,
                            position_at_start +
                                scale * (Vector3d)world_direction);
                      };
                  add_vector(manoeuvre.tangent, XKCDColors.NeonYellow);
                  add_vector(manoeuvre.normal, XKCDColors.AquaBlue);
                  add_vector(manoeuvre.binormal, XKCDColors.PurplePink);
                }
              }
            }
          }
        });
      }
      map_node_pool_.Update();
    } else {
      map_node_pool_.Clear();
    }
  }

  private void RenderPredictionMarkers(String vessel_guid,
                                       XYZ sun_world_position) {
    if (plotting_frame_selector_.get().target_override) {
      Vessel target = plotting_frame_selector_.get().target_override;
      DisposableIterator ascending_nodes_iterator;
      DisposableIterator descending_nodes_iterator;
      DisposableIterator approaches_iterator;
      plugin_.RenderedPredictionNodes(vessel_guid,
                                      sun_world_position,
                                      out ascending_nodes_iterator,
                                      out descending_nodes_iterator);
      plugin_.RenderedPredictionClosestApproaches(vessel_guid,
                                                  sun_world_position,
                                                  out approaches_iterator);
      map_node_pool_.RenderMarkers(
          ascending_nodes_iterator,
          MapObject.ObjectType.AscendingNode,
          MapNodePool.NodeSource.PREDICTION,
          vessel    : target,
          celestial : plotting_frame_selector_.get().selected_celestial);
      map_node_pool_.RenderMarkers(
          descending_nodes_iterator,
          MapObject.ObjectType.DescendingNode,
          MapNodePool.NodeSource.PREDICTION,
          vessel    : target,
          celestial : plotting_frame_selector_.get().selected_celestial);
      map_node_pool_.RenderMarkers(
          approaches_iterator,
          MapObject.ObjectType.ApproachIntersect,
          MapNodePool.NodeSource.PREDICTION,
          vessel    : target,
          celestial : plotting_frame_selector_.get().selected_celestial);
    } else {
      foreach (CelestialBody celestial in
               plotting_frame_selector_.get().FixedBodies()) {
        DisposableIterator apoapsis_iterator;
        DisposableIterator periapsis_iterator;
        plugin_.RenderedPredictionApsides(vessel_guid,
                                          celestial.flightGlobalsIndex,
                                          sun_world_position,
                                          out apoapsis_iterator,
                                          out periapsis_iterator);
        map_node_pool_.RenderMarkers(
            apoapsis_iterator,
            MapObject.ObjectType.Apoapsis,
            MapNodePool.NodeSource.PREDICTION,
            vessel    : null,
            celestial : celestial);
        map_node_pool_.RenderMarkers(
            periapsis_iterator,
            MapObject.ObjectType.Periapsis,
            MapNodePool.NodeSource.PREDICTION,
            vessel    : null,
            celestial : celestial);
      }
      var frame_type = plotting_frame_selector_.get().frame_type;
      if (frame_type ==
              ReferenceFrameSelector.FrameType.BARYCENTRIC_ROTATING ||
          frame_type == ReferenceFrameSelector.FrameType
                            .BODY_CENTRED_PARENT_DIRECTION) {
        var primary =
            plotting_frame_selector_.get().selected_celestial.referenceBody;
        DisposableIterator ascending_nodes_iterator;
        DisposableIterator descending_nodes_iterator;
        plugin_.RenderedPredictionNodes(vessel_guid,
                                        sun_world_position,
                                        out ascending_nodes_iterator,
                                        out descending_nodes_iterator);
        map_node_pool_.RenderMarkers(
            ascending_nodes_iterator,
            MapObject.ObjectType.AscendingNode,
            MapNodePool.NodeSource.PREDICTION,
            vessel    : null,
            celestial : primary);
        map_node_pool_.RenderMarkers(
            descending_nodes_iterator,
            MapObject.ObjectType.DescendingNode,
            MapNodePool.NodeSource.PREDICTION,
            vessel    : null,
            celestial : primary);
      }
    }
  }

  private void RenderFlightPlanMarkers(String vessel_guid,
                                       XYZ sun_world_position) {
    if (plotting_frame_selector_.get().target_override) {
      Vessel target = plotting_frame_selector_.get().target_override;
      DisposableIterator ascending_nodes_iterator;
      DisposableIterator descending_nodes_iterator;
      DisposableIterator approaches_iterator;
      plugin_.FlightPlanRenderedNodes(vessel_guid,
                                      sun_world_position,
                                      out ascending_nodes_iterator,
                                      out descending_nodes_iterator);
      plugin_.FlightPlanRenderedClosestApproaches(vessel_guid,
                                                  sun_world_position,
                                                  out approaches_iterator);
      map_node_pool_.RenderMarkers(
          ascending_nodes_iterator,
          MapObject.ObjectType.AscendingNode,
          MapNodePool.NodeSource.FLIGHT_PLAN,
          vessel    : target,
          celestial : plotting_frame_selector_.get().selected_celestial);
      map_node_pool_.RenderMarkers(
          descending_nodes_iterator,
          MapObject.ObjectType.DescendingNode,
          MapNodePool.NodeSource.FLIGHT_PLAN,
          vessel    : target,
          celestial : plotting_frame_selector_.get().selected_celestial);
      map_node_pool_.RenderMarkers(
          approaches_iterator,
          MapObject.ObjectType.ApproachIntersect,
          MapNodePool.NodeSource.FLIGHT_PLAN,
          vessel    : target,
          celestial : plotting_frame_selector_.get().selected_celestial);
    } else {
      foreach (CelestialBody celestial in
               plotting_frame_selector_.get().FixedBodies()) {
        DisposableIterator apoapsis_iterator;
        DisposableIterator periapsis_iterator;
        plugin_.FlightPlanRenderedApsides(vessel_guid,
                                          celestial.flightGlobalsIndex,
                                          sun_world_position,
                                          out apoapsis_iterator,
                                          out periapsis_iterator);
        map_node_pool_.RenderMarkers(
            apoapsis_iterator,
            MapObject.ObjectType.Apoapsis,
            MapNodePool.NodeSource.FLIGHT_PLAN,
            vessel    : null,
            celestial : celestial);
        map_node_pool_.RenderMarkers(
            periapsis_iterator,
            MapObject.ObjectType.Periapsis,
            MapNodePool.NodeSource.FLIGHT_PLAN,
            vessel    : null,
            celestial : celestial);
      }
      var frame_type = plotting_frame_selector_.get().frame_type;
      if (frame_type ==
              ReferenceFrameSelector.FrameType.BARYCENTRIC_ROTATING ||
          frame_type == ReferenceFrameSelector.FrameType
                            .BODY_CENTRED_PARENT_DIRECTION) {
        var primary =
            plotting_frame_selector_.get().selected_celestial.referenceBody;
        DisposableIterator ascending_nodes_iterator;
        DisposableIterator descending_nodes_iterator;
        plugin_.FlightPlanRenderedNodes(vessel_guid,
                                        sun_world_position,
                                        out ascending_nodes_iterator,
                                        out descending_nodes_iterator);
        map_node_pool_.RenderMarkers(
            ascending_nodes_iterator,
            MapObject.ObjectType.AscendingNode,
            MapNodePool.NodeSource.PREDICTION,
            vessel    : null,
            celestial : primary);
        map_node_pool_.RenderMarkers(
            descending_nodes_iterator,
            MapObject.ObjectType.DescendingNode,
            MapNodePool.NodeSource.PREDICTION,
            vessel    : null,
            celestial : primary);
      }
    }
  }

  private void Cleanup() {
    UnityEngine.Object.Destroy(map_renderer_);
    map_node_pool_.Clear();
    map_renderer_ = null;
    Interface.DeletePlugin(ref plugin_);
    plotting_frame_selector_.reset();
    previous_display_mode_ = null;
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
    using (new UnityEngine.GUILayout.VerticalScope()) {
      if (!PluginRunning()) {
        UnityEngine.GUILayout.TextArea(text : "Plugin is not started");
      }
      if (DateTimeOffset.Now > next_release_date_) {
        UnityEngine.GUILayout.TextArea(
            "Announcement: the new moon of lunation number " +
            next_release_lunation_number_ +
            " has come; please download the latest Principia release, " +
            next_release_name_ + ".");
      }
      String version;
      String unused_build_date;
      Interface.GetVersion(build_date: out unused_build_date,
                           version: out version);
      UnityEngine.GUILayout.TextArea(version);
      bool changed_history_length = false;
      Selector(history_lengths_,
               ref history_length_index_,
               "Max history length",
               ref changed_history_length,
               "{0:0.00e00} s");
      if (MapView.MapIsEnabled &&
          FlightGlobals.ActiveVessel?.orbitTargeter != null) {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          selecting_active_vessel_target_ = UnityEngine.GUILayout.Toggle(
              selecting_active_vessel_target_, "Select target vessel...");
          if (selecting_active_vessel_target_) {
            selecting_target_celestial_ = false;
          }
          if (FlightGlobals.fetch.VesselTarget?.GetVessel()) {
            UnityEngine.GUILayout.Label(
                "Target: " +
                    FlightGlobals.fetch.VesselTarget.GetVessel().vesselName,
                UnityEngine.GUILayout.ExpandWidth(true));
            if (UnityEngine.GUILayout.Button("Clear",
                                             UnityEngine.GUILayout.Width(50))) {
              selecting_active_vessel_target_ = false;
              FlightGlobals.fetch.SetVesselTarget(null);
            }
            if (UnityEngine.GUILayout.Button("Switch To")) {
              var focus_object =
                  new KSP.UI.Screens.Mapview.MapContextMenuOptions.FocusObject(
                      FlightGlobals.fetch.VesselTarget.GetVessel().orbitDriver);
              focus_object.onOptionSelected();
            }
          }
        }
      } else {
        selecting_active_vessel_target_ = false;
      }
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
#if CRASH_BUTTON
      ToggleableSection(name   : "CRASH",
                        show   : ref show_crash_options_,
                        render : CrashOptions);
#endif
    }
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
      plotting_frame_selector_.get().RenderButton();
    }
  }

  private void Selector(
      double[] array,
      ref int index,
      String label,
      ref bool changed,
      String format) {
    using (new UnityEngine.GUILayout.HorizontalScope()) {
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
    }
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
    if (MapView.MapIsEnabled &&
        FlightGlobals.ActiveVessel?.orbitTargeter != null) {
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        selecting_target_celestial_ = UnityEngine.GUILayout.Toggle(
            selecting_target_celestial_, "Select target celestial...");
        if (selecting_target_celestial_) {
          selecting_active_vessel_target_ = false;
        }
        CelestialBody target_celestial =
            FlightGlobals.fetch.VesselTarget as CelestialBody;
        if (target_celestial) {
          UnityEngine.GUILayout.Label("Target: " + target_celestial.name,
                                      UnityEngine.GUILayout.ExpandWidth(true));
          if (UnityEngine.GUILayout.Button("Clear",
                                           UnityEngine.GUILayout.Width(50))) {
            selecting_target_celestial_ = false;
            FlightGlobals.fetch.SetVesselTarget(null);
          }
        }
      }
    } else {
      selecting_target_celestial_ = false;
    }
  }

  private void LoggingSettings() {
#if SELECTABLE_PLOT_METHOD
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label("Чебышёв plotting method:");
      for (int i = 0; i < чебышёв_plotting_methods_count; ++i) {
        if (UnityEngine.GUILayout.Toggle(чебышёв_plotting_method_ == i,
                                         i.ToString())) {
          чебышёв_plotting_method_ = i;
        }
      }
    }
#endif
    using (new UnityEngine.GUILayout.HorizontalScope()) {
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
    }
    int column_width = 75;
    using (new UnityEngine.GUILayout.HorizontalScope()) {
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
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
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
    }
    for (int severity = 0; severity <= 3; ++severity) {
      using (new UnityEngine.GUILayout.HorizontalScope()) {
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
      }
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
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
    }
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

  private void ShrinkMainWindow() {
    main_window_rectangle_.height = 0.0f;
    main_window_rectangle_.width = 0.0f;
  }

  private void UpdateRenderingFrame(
      NavigationFrameParameters frame_parameters) {
    plugin_.SetPlottingFrame(frame_parameters);
    var frame_type =
        (ReferenceFrameSelector.FrameType)frame_parameters.extension;
    if (frame_type != ReferenceFrameSelector.FrameType.BODY_SURFACE) {
      last_non_surface_frame_type_ = frame_type;
    }
    navball_changed_ = true;
    reset_rsas_target_ = true;
  }

  private static void InitializeIntegrators(
      IntPtr plugin,
      ConfigNode numerics_blueprint) {
    if (numerics_blueprint == null) {
      return;
    }
    var ephemeris_parameters =
        numerics_blueprint.GetAtMostOneNode("ephemeris");
    if (ephemeris_parameters != null) {
      plugin.InitializeEphemerisParameters(
          ConfigNodeParsers.NewConfigurationAccuracyParameters(
              ephemeris_parameters),
          ConfigNodeParsers.NewConfigurationFixedStepParameters(
              ephemeris_parameters));
    }

    var history_parameters =
        numerics_blueprint.GetAtMostOneNode("history");
    if (history_parameters != null) {
      plugin.InitializeHistoryParameters(
          ConfigNodeParsers.NewConfigurationFixedStepParameters(
              history_parameters));
    }

    var psychohistory_parameters =
        numerics_blueprint.GetAtMostOneNode("psychohistory");
    if (psychohistory_parameters != null) {
      plugin.InitializePsychohistoryParameters(
          ConfigNodeParsers.NewConfigurationAdaptiveStepParameters(
              psychohistory_parameters));
    }
  }

  private void ResetPlugin() {
  try {
    Cleanup();
    RemoveBuggyTidalLocking();
    plugin_construction_ = DateTime.Now;
    Dictionary<String, ConfigNode> name_to_gravity_model = null;
    ConfigNode gravity_model = GameDatabase.Instance.GetAtMostOneNode(
        principia_gravity_model_config_name_);
    ConfigNode initial_state = GameDatabase.Instance.GetAtMostOneNode(
        principia_initial_state_config_name_);
    ConfigNode numerics_blueprint = GameDatabase.Instance.GetAtMostOneNode(
        principia_numerics_blueprint_config_name_);
    if (gravity_model != null) {
      name_to_gravity_model = gravity_model.GetNodes("body").ToDictionary(
                                  node => node.GetUniqueValue("name"));
    }
    if (initial_state != null) {
      if (name_to_gravity_model == null) {
        Log.Fatal("Cartesian config without gravity models");
      }
      // Note that |game_epoch| is not in the astronomy proto, as it is
      // KSP-specific: it is the |Instant| corresponding to KSP's
      // UniversalTime 0.  It is passed as an argument to
      // |generate_configuration|.
      plugin_ = Interface.NewPlugin(
          initial_state.GetUniqueValue("game_epoch"),
          initial_state.GetUniqueValue("solar_system_epoch"),
          Planetarium.InverseRotAngle);
      InitializeIntegrators(plugin_, numerics_blueprint);
      var name_to_initial_state = initial_state.GetNodes("body").ToDictionary(
                                      node => node.GetUniqueValue("name"));
      BodyProcessor insert_body = body => {
        Log.Info("Inserting " + body.name + "...");
        ConfigNode body_gravity_model;
        if (!name_to_gravity_model.TryGetValue(body.name,
                                               out body_gravity_model)) {
          Log.Fatal("missing gravity model for " + body.name);
        }
        ConfigNode body_initial_state;
        if (!name_to_initial_state.TryGetValue(body.name,
                                               out body_initial_state)) {
          Log.Fatal("missing Cartesian initial state for " + body.name);
        }
        int? parent_index = body.orbit?.referenceBody.flightGlobalsIndex;
        // GetUniqueValue resp. GetAtMostOneValue corresponding to required
        // resp. optional in principia.serialization.GravityModel.Body.
        var body_parameters =
            ConfigNodeParsers.NewCartesianBodyParameters(body,
                                                         body_gravity_model);
        // GetUniqueValue since these are all required fields in
        // principia.serialization.InitialState.Cartesian.Body.
        plugin_.InsertCelestialAbsoluteCartesian(
            celestial_index         : body.flightGlobalsIndex,
            parent_index            : parent_index,
            body_parameters         : body_parameters,
            x                       : body_initial_state.GetUniqueValue("x"),
            y                       : body_initial_state.GetUniqueValue("y"),
            z                       : body_initial_state.GetUniqueValue("z"),
            vx                      : body_initial_state.GetUniqueValue("vx"),
            vy                      : body_initial_state.GetUniqueValue("vy"),
            vz                      : body_initial_state.GetUniqueValue("vz"));
      };
      insert_body(Planetarium.fetch.Sun);
      ApplyToBodyTree(insert_body);
      plugin_.EndInitialization();
    } else {
      // We create the plugin at J2000 (a.k.a. Instant{}), rather than
      // |Planetarium.GetUniversalTime()|, in order to get a deterministic
      // initial state.
      plugin_ = Interface.NewPlugin("JD2451545", "JD2451545",
                                    Planetarium.InverseRotAngle);
      InitializeIntegrators(plugin_, numerics_blueprint);
      BodyProcessor insert_body = body => {
        Log.Info("Inserting " + body.name + "...");
        ConfigNode body_gravity_model = null;
        if (name_to_gravity_model?.TryGetValue(
                body.name,
                out body_gravity_model) == true) {
          Log.Info("using custom gravity model");
        }
        Orbit orbit = unmodified_orbits_.GetValueOrNull(body);
        var body_parameters =
            ConfigNodeParsers.NewKeplerianBodyParameters(body,
                                                         body_gravity_model);
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
    }
    if (Planetarium.GetUniversalTime() > plugin_.CurrentTime()) {
      // Make sure that the plugin has caught up with the game before existing
      // vessels are added.
      plugin_.AdvanceTime(Planetarium.GetUniversalTime(),
                          Planetarium.InverseRotAngle);
    }
    plotting_frame_selector_.reset(
        new ReferenceFrameSelector(this,
                                   plugin_,
                                   UpdateRenderingFrame,
                                   "Plotting frame"));
    must_set_plotting_frame_ = true;
    flight_planner_.reset(new FlightPlanner(this, plugin_));
  } catch (Exception e) {
    Log.Fatal("Exception while resetting plugin: " + e.ToString());
  }
  }

  private void RemoveBuggyTidalLocking() {
    ApplyToBodyTree(body => body.tidallyLocked = false);
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
