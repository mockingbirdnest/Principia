using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace principia {
namespace ksp_plugin_adapter {

[KSPScenario(createOptions: ScenarioCreationOptions.AddToAllGames,
             tgtScenes: new GameScenes[]{GameScenes.FLIGHT,
                                         GameScenes.MAINMENU,
                                         GameScenes.SPACECENTER,
                                         GameScenes.TRACKSTATION})]
public partial class PrincipiaPluginAdapter
    : ScenarioModule,
      SupervisedWindowRenderer.ISupervisor {

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

  private const string principia_serialized_plugin_ = "serialized_plugin";
  private const string principia_initial_state_config_name_ =
      "principia_initial_state";
  private const string principia_gravity_model_config_name_ =
      "principia_gravity_model";
  private const string principia_numerics_blueprint_config_name_ =
      "principia_numerics_blueprint";
  private const string principia_override_version_check_config_name_ =
      "principia_override_version_check";
  private const string principia_flags_ = "principia_flags";

  private KSP.UI.Screens.ApplicationLauncherButton toolbar_button_;
  // Whether the user has hidden the UI.
  private bool hide_all_gui_ = false;
  // Whether we are in a scene/building where we wish to show our UI.
  private bool in_principia_scene_ = true;
  // Whether we are in the main menu.
  private readonly bool in_main_menu_ =
      HighLogic.LoadedScene == GameScenes.MAINMENU;

  private const int чебышёв_plotting_method_ = 2;

  private IntPtr plugin_ = IntPtr.Zero;
  internal IntPtr Plugin() {
    return plugin_;
  }

  // Whether to compress saves.
  [KSPField(isPersistant = true)]
  private string serialization_compression_ = "";
  [KSPField(isPersistant = true)]
  private string serialization_encoding_ = "hexadecimal";

  // Whether the plotting frame must be set to something convenient at the next
  // opportunity.
  private bool must_set_plotting_frame_ = false;

  private bool time_is_advancing_;

  // Used to detect changes of SOI and skip two frames during which we use the
  // EulerSolver.
  CelestialBody last_main_body_;
  private int main_body_change_countdown_ = 1;

  private PlanetariumCameraAdjuster planetarium_camera_adjuster_;

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

  private readonly List<IntPtr> vessel_futures_ = new List<IntPtr>();

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

  private KSP.UI.Screens.DebugToolbar.DebugScreen debug_screen_;
  private KSP.UI.Screens.DebugToolbar.Screens.Cheats.HackGravity hack_gravity_;
  private KSP.UI.Screens.DebugToolbar.Screens.Cheats.HackGravity hack_gravity {
    get {
      if (hack_gravity_ == null) {
        if (debug_screen_ == null) {
          // Force a debug screen to be instantiated, by showing it.
          KSP.UI.Screens.DebugToolbar.DebugScreenSpawner.ShowDebugScreen();
          debug_screen_ = (KSP.UI.Screens.DebugToolbar.DebugScreen)
              FindObjectOfType(typeof(KSP.UI.Screens.DebugToolbar.DebugScreen));
          // Now hide it.
          debug_screen_.Hide();
        }
        // Since we have the debug screen, we can restrict our search for the
        // gravity-hacking control within it, so as not to slow things down to
        // a grind until it is instantiated.
        hack_gravity_ = debug_screen_.contentTransform.GetComponentInChildren<
            KSP.UI.Screens.DebugToolbar.Screens.Cheats.HackGravity>();
      }
      return hack_gravity_;
    }
  }

  // TODO(egg): these should be moved to the C++; there it can be made a vector
  // because we can test whether we have the part or not. Set in
  // FashionablyLate, before the FlightIntegrator clears the forces.  Used in
  // WaitForFixedUpdate.
  private readonly Dictionary<uint, Vector3d> part_id_to_intrinsic_torque_ =
      new Dictionary<uint, Vector3d>();
  private readonly Dictionary<uint, Vector3d> part_id_to_intrinsic_force_ =
      new Dictionary<uint, Vector3d>();

  // Work around the launch backflip issue encountered while releasing
  // Frobenius.
  // The issue is that the position of a  |ForceHolder| collected in
  // |FashionablyLate| or constructed in |JaiFailliAttendre| is made invalid by
  // changes to world coordinates in |FloatingOrigin| by the time it is used in
  // |WaitForFixedUpdate|.
  // TODO(egg): This should be cleaned up a bit, e.g., by making the plugin and
  // interface take the |World| lever arm directly, since the lever arm is what
  // we use in |principia::ksp_plugin::Part|, and it is what we compute below:
  // this would obviate the need to keep track of the |Part| and to adjust its
  // position.
  class PartCentredForceHolder {
    public static PartCentredForceHolder FromPartForceHolder(
        Part part,
        Part.ForceHolder holder) {
      return new PartCentredForceHolder{
          force_in_world_coordinates_ = holder.force,
          lever_arm_in_world_coordinates_ = holder.pos - part.rb.position};
    }

    public Vector3d force => force_in_world_coordinates_;
    public Vector3d lever_arm => lever_arm_in_world_coordinates_;

    private Vector3d force_in_world_coordinates_;
    private Vector3d lever_arm_in_world_coordinates_;
  }
  private readonly Dictionary<uint, PartCentredForceHolder[]>
      part_id_to_intrinsic_forces_ = new Dictionary<uint,
                                                    PartCentredForceHolder[]>();

  // The degrees of freedom at BetterLateThanNever.  Those are used to insert
  // new parts with the correct initial state.
  private readonly Dictionary<uint, QP> part_id_to_degrees_of_freedom_ =
      new Dictionary<uint, QP>();

  private readonly MapNodePool map_node_pool_;
  private ManeuverNode guidance_node_;

  // UI for the apocalypse notification.
  [KSPField(isPersistant = true)]
  private readonly Dialog apocalypse_dialog_ = new Dialog(persist_state: true);

  // UI for the bad installation notification.
  private readonly bool is_bad_installation_ = false;  // Don't persist.
  [KSPField(isPersistant = true)]
  private readonly Dialog bad_installation_dialog_ =
      new Dialog(persist_state: false);

  // The game windows.
  [KSPField(isPersistant = true)]
  private readonly FlightPlanner flight_planner_;
  [KSPField(isPersistant = true)]
  private readonly OrbitAnalyser orbit_analyser_;
  [KSPField(isPersistant = true)]
  internal ReferenceFrameSelector plotting_frame_selector_;
  [KSPField(isPersistant = true)]
  internal MainWindow main_window_;

  public event Action LockClearing;
  public event Action WindowsDisposal;
  public event Action WindowsRendering;

  PrincipiaPluginAdapter() {
    // We create this directory here so we do not need to worry about cross-
    // platform problems in C++.
    Directory.CreateDirectory("glog/Principia");
    string load_error = Loader.LoadPrincipiaDllAndInitGoogleLogging();
    if (load_error == null) {
      is_bad_installation_ = false;
      bad_installation_dialog_.Hide();
    } else {
      is_bad_installation_ = true;
      bad_installation_dialog_.message =
          "The Principia DLL failed to load.\n" + load_error +
          "\n\nWarning: don't load a Principia save before you have fixed this " +
          "error; it might get damaged.";
      bad_installation_dialog_.Show();
    }
#if KSP_VERSION_1_9_1
    if (!(Versioning.version_major == 1 &&
          (Versioning.version_minor == 8 && Versioning.Revision == 1) ||
          (Versioning.version_minor == 9 && Versioning.Revision == 1))) {
      string expected_version = "1.8.1 and 1.9.1";
#elif KSP_VERSION_1_7_3
    if (!(Versioning.version_major == 1 &&
          (Versioning.version_minor == 5 && Versioning.Revision == 1) ||
          (Versioning.version_minor == 6 && Versioning.Revision == 1) ||
          (Versioning.version_minor == 7 && Versioning.Revision <= 3))) {
      string expected_version = "1.7.3, 1.7.2, 1.7.1, 1.7.0, 1.6.1, and 1.5.1";
#endif
      string message = $@"Unexpected KSP version {Versioning.version_major}.{
          Versioning.version_minor}.{Versioning.Revision}; this build targets {
          expected_version}.";
      if (GameDatabase.Instance.GetAtMostOneNode(
              principia_override_version_check_config_name_) == null) {
        Log.Fatal(message);
      } else {
        Log.Error(message);
      }
    }

    map_node_pool_ = new MapNodePool();
    flight_planner_ = new FlightPlanner(this, PredictedVessel);
    orbit_analyser_ = new OrbitAnalyser(this, PredictedVessel);
    plotting_frame_selector_ = new ReferenceFrameSelector(this,
                                                          UpdateRenderingFrame,
                                                          "Plotting frame");
    main_window_ = new MainWindow(this,
                                  flight_planner_,
                                  orbit_analyser_,
                                  plotting_frame_selector_,
                                  PredictedVessel);
  }

  ~PrincipiaPluginAdapter() {
    // We should not get here without deleting the plugin, but just for safety.
    Interface.DeletePlugin(ref plugin_);
  }

  public bool PluginRunning() {
    return plugin_ != IntPtr.Zero;
  }

  private Vessel PredictedVessel() {
    if (!PluginRunning()) {
      return null;
    }
    Vessel vessel =
        FlightGlobals.ActiveVessel ?? space_tracking?.SelectedVessel;
    string vessel_guid = vessel?.id.ToString();
    if (vessel_guid != null && plugin_.HasVessel(vessel_guid)) {
      return vessel;
    } else {
      return null;
    }
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
    while (stack.Count > 0) {
      CelestialBody body = stack.Pop();
      process_body(body);
      foreach (CelestialBody child in body.orbitingBodies) {
        stack.Push(child);
      }
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
    Orbit copy = new Orbit(original.inclination,
                           original.eccentricity,
                           original.semiMajorAxis,
                           original.LAN,
                           original.argumentOfPeriapsis,
                           original.meanAnomalyAtEpoch,
                           original.epoch,
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
    Vessel main_vessel = PredictedVessel();
    bool ready_to_draw_active_vessel_trajectory =
        main_vessel != null && MapView.MapIsEnabled;

    if (ready_to_draw_active_vessel_trajectory) {
      plugin_.UpdatePrediction(main_vessel.id.ToString());
      string target_id =
          FlightGlobals.fetch.VesselTarget?.GetVessel()?.id.ToString();
      if (!plotting_frame_selector_.target_override &&
          target_id != null && plugin_.HasVessel(target_id)) {
        // TODO(phl): It's not nice that we are overriding the target vessel
        // parameters.
        AdaptiveStepParameters adaptive_step_parameters =
            plugin_.VesselGetPredictionAdaptiveStepParameters(
                main_vessel.id.ToString());
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
      if (vessel.loaded) {
        foreach (Part part in vessel.parts.Where(part => part.rb != null &&
                                                         plugin_.PartIsTruthful(
                                                             part.flightID))) {
          // TODO(egg): What if the plugin doesn't have the part? this seems
          // brittle.
          // NOTE(egg): I am not sure what the origin is here, as we are
          // running before the floating origin and krakensbane.  Do everything
          // with respect to the root part, since the overall linear motion of
          // the vessel is handled with by the orbit anyway.
          // TODO(egg): check that the vessel is moved *after* this.  Shouldn't
          // we be calling vessel.orbitDriver.updateFromParameters() after
          // setting the orbit anyway?
          QPRW part_actual_motion = plugin_.PartGetActualRigidMotion(
              part.flightID,
              new Origin{
                  reference_part_is_at_origin = true,
                  reference_part_is_unmoving = true,
                  main_body_centre_in_world =
                      (XYZ)FlightGlobals.ActiveVessel.mainBody.position,
                  reference_part_id = vessel.rootPart.flightID
              });
          part.rb.position = vessel.rootPart.rb.position +
                             (Vector3d)part_actual_motion.qp.q;
          part.rb.transform.position = vessel.rootPart.rb.position +
                                       (Vector3d)part_actual_motion.qp.q;
          part.rb.rotation = (UnityEngine.QuaternionD)part_actual_motion.r;
          part.rb.transform.rotation =
              (UnityEngine.QuaternionD)part_actual_motion.r;
          part.rb.velocity = vessel.rootPart.rb.velocity +
                             (Vector3d)part_actual_motion.qp.p;
          part.rb.angularVelocity = (Vector3d)part_actual_motion.w;
        }
      }
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

  // Æthelred Kerman.
  private bool is_unready_kerbal(Vessel vessel) {
    return vessel.isEVA && vessel.evaController?.Ready == false;
  }

  private bool is_manageable(Vessel vessel) {
    return UnmanageabilityReasons(vessel) == null;
  }

  public static bool is_lit_solid_booster(ModuleEngines module) {
    return module != null &&
        module.EngineIgnited &&
        module.engineType == EngineType.SolidBooster;
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
    if (!vessel.packed &&
        JustAboveTheGround(vessel,
                           out double height,
                           out double vertical_speed)) {
      reasons.Add("vessel is " + height +
                  " m above ground with a vertical speed of " + vertical_speed +
                  " m/s");
    }
    if (is_unready_kerbal(vessel)) {
      reasons.Add("vessel is an unready Kerbal");
    }
    if (vessel.loaded && hack_gravity?.toggle.isOn == true) {
      reasons.Add("vessel is loaded and gravity is being hacked");
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

  private static bool JustAboveTheGround(Vessel vessel,
                                         out double height,
                                         out double vertical_speed) {
    height = vessel.altitude - vessel.terrainAltitude;
    vertical_speed = vessel.verticalSpeed;
    double Δt = Planetarium.TimeScale * Planetarium.fetch.fixedDeltaTime;
    return height + vertical_speed * Δt < 0;
  }

  private static bool IsNaN(Vector3d v) {
    return double.IsNaN(v.x + v.y + v.z);
  }

  private static bool IsNaN(XYZ xyz) {
    return double.IsNaN(xyz.x + xyz.y + xyz.z);
  }

  // It seems that parts sometimes have NaN position, velocity or angular
  // velocity, presumably because they are being destroyed.  Just skip these
  // unfaithful parts as if they had no rigid body.
  private static bool PartIsFaithful(Part part) {
    return part.rb != null &&
           !IsNaN(part.rb.position) &&
           !IsNaN(part.rb.velocity) &&
           !IsNaN(part.rb.angularVelocity);
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
                                   string path) {
    string full_path =
        KSPUtil.ApplicationRootPath + Path.DirectorySeparatorChar +
        "GameData" + Path.DirectorySeparatorChar +
        "Principia" + Path.DirectorySeparatorChar +
        "assets" + Path.DirectorySeparatorChar +
        path;
    if (File.Exists(full_path)) {
      var texture2d = new UnityEngine.Texture2D(2, 2);
      bool success = UnityEngine.ImageConversion.LoadImage(
          texture2d, File.ReadAllBytes(full_path));
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

  private void LoadTextureOrDie(out UnityEngine.Texture texture, string path) {
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
    if (is_bad_installation_ || in_main_menu_) {
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

    GameEvents.onShowUI.Add(() => { hide_all_gui_ = false; });
    GameEvents.onHideUI.Add(() => { hide_all_gui_ = true; });
    GameEvents.onGUIAdministrationFacilitySpawn.Add(() => {
      in_principia_scene_ = false;
    });
    GameEvents.onGUIAdministrationFacilityDespawn.Add(() => {
      in_principia_scene_ = true;
    });
    GameEvents.onGUIAstronautComplexSpawn.Add(() => {
      in_principia_scene_ = false;
    });
    GameEvents.onGUIAstronautComplexDespawn.Add(() => {
      in_principia_scene_ = true;
    });
    GameEvents.onGUIMissionControlSpawn.Add(() => {
      in_principia_scene_ = false;
    });
    GameEvents.onGUIMissionControlDespawn.Add(() => {
      in_principia_scene_ = true;
    });
    GameEvents.onGUIRnDComplexSpawn.Add(() => {
      in_principia_scene_ = false;
    });
    GameEvents.onGUIRnDComplexDespawn.Add(() => {
      in_principia_scene_ = true;
    });
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
    // Custom timing, 301.
    planetarium_camera_adjuster_ =
        gameObject.AddComponent<PlanetariumCameraAdjuster>();
    planetarium_camera_adjuster_.adapter = this;
    // Timing5, 8008.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.BetterLateThanNever,
                                 BetterLateThanNever);
    TimingManager.LateUpdateAdd(
        TimingManager.TimingStage.BetterLateThanNever,
        BetterLateThanNeverLateUpdate);
  }

  public override void OnSave(ConfigNode node) {
    base.OnSave(node);
    if (PluginRunning()) {
      IntPtr serializer = IntPtr.Zero;
      for (;;) {
        string serialization = plugin_.SerializePlugin(
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
    if (is_bad_installation_ || in_main_menu_) {
      return;
    }
    if (node.HasValue(principia_serialized_plugin_)) {
      Cleanup();
      RemoveBuggyTidalLocking();

      IntPtr deserializer = IntPtr.Zero;
      string[] serializations = node.GetValues(principia_serialized_plugin_);
      Log.Info("Serialization has " + serializations.Length + " chunks");
      foreach (string serialization in serializations) {
        Interface.DeserializePlugin(serialization,
                                    ref deserializer,
                                    ref plugin_,
                                    serialization_compression_,
                                    serialization_encoding_);
      }
      Interface.DeserializePlugin("",
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

      previous_display_mode_ = null;
      must_set_plotting_frame_ = true;
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
    if (is_bad_installation_ || in_main_menu_) {
      bad_installation_dialog_.RenderWindow();
      return;
    }

    apocalypse_dialog_.RenderWindow();

    if (KSP.UI.Screens.ApplicationLauncher.Ready && toolbar_button_ == null) {
      LoadTextureOrDie(out UnityEngine.Texture toolbar_button_texture,
                       "toolbar_button.png");
      toolbar_button_ =
          KSP.UI.Screens.ApplicationLauncher.Instance.AddModApplication(
              onTrue          : () => main_window_.Show(),
              onFalse         : () => main_window_.Hide(),
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
    if (main_window_.Shown()) {
      toolbar_button_?.SetTrue(makeCall : false);
    } else {
      toolbar_button_?.SetFalse(makeCall : false);
    }

    if (hide_all_gui_ || !in_principia_scene_) {
      LockClearing();
    } else if (main_window_.Shown()) {
      WindowsRendering();
    } else {
      LockClearing();
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
      RenderNavball(active_vessel);
      if (!PluginRunning()) {
        return;
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
        plotting_frame_selector_.SetToSurfaceFrameOf(
            active_vessel.mainBody);
      }

      if (plotting_frame_selector_.target_override == null &&
          FlightGlobals.speedDisplayMode ==
              FlightGlobals.SpeedDisplayModes.Target) {
        KSP.UI.Screens.Flight.SpeedDisplay.Instance.textTitle.text = "Target";
      }

      if (FlightGlobals.speedDisplayMode ==
              FlightGlobals.SpeedDisplayModes.Orbit ||
          FlightGlobals.speedDisplayMode ==
              FlightGlobals.SpeedDisplayModes.Surface ||
          plotting_frame_selector_.target_override) {
        bool plugin_has_active_manageable_vessel =
            has_active_manageable_vessel() &&
            plugin_.HasVessel(active_vessel.id.ToString());

        KSP.UI.Screens.Flight.SpeedDisplay speed_display =
            KSP.UI.Screens.Flight.SpeedDisplay.Instance;
        if (speed_display?.textTitle != null &&
            speed_display?.textSpeed != null &&
            !ferram_owns_the_speed_display) {
          speed_display.textTitle.text = plotting_frame_selector_.ShortName();
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
        var prograde =
            (Vector3d)plugin_.VesselTangent(active_vessel.id.ToString());
        var radial =
            (Vector3d)plugin_.VesselNormal(active_vessel.id.ToString());
        // Yes, the astrodynamicist's normal is the mathematician's binormal.
        // Don't ask.
        var normal =
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
    if (is_bad_installation_ || in_main_menu_) {
      return;
    }
    if (GameSettings.ORBIT_WARP_DOWN_AT_SOI) {
      Log.Info("Setting GameSettings.ORBIT_WARP_DOWN_AT_SOI to false");
      GameSettings.ORBIT_WARP_DOWN_AT_SOI = false;
    }
    if (must_set_plotting_frame_) {
      must_set_plotting_frame_ = false;
      plotting_frame_selector_.UpdateMainBody();
      previous_display_mode_ = null;
    }

    if (PluginRunning()) {
      plugin_.SetMainBody(
          (FlightGlobals.currentMainBody
               ?? FlightGlobals.GetHomeBody()).flightGlobalsIndex);

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
    if (is_bad_installation_ || in_main_menu_) {
      return;
    }
    Log.Info("principia.ksp_plugin_adapter.PrincipiaPluginAdapter.OnDisable()");
    if (toolbar_button_ != null) {
      KSP.UI.Screens.ApplicationLauncher.Instance.RemoveModApplication(
          toolbar_button_);
    }
    Cleanup();
    WindowsDisposal();
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
    TimingManager.LateUpdateRemove(
        TimingManager.TimingStage.BetterLateThanNever,
        BetterLateThanNeverLateUpdate);
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

      plugin_.InsertOrKeepVessel(vessel.id.ToString(),
                                 vessel.vesselName,
                                 vessel.mainBody.flightGlobalsIndex,
                                 !vessel.packed,
                                 out bool inserted);
      if (!vessel.packed) {
        foreach (Part part in vessel.parts.Where(PartIsFaithful)) {
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
              (XYZ)(Vector3d)part.rb.inertiaTensor,
              (WXYZ)(UnityEngine.QuaternionD)part.rb.inertiaTensorRotation,
              (from PartModule module in part.Modules
               select module as ModuleEngines).Any(is_lit_solid_booster),
              vessel.id.ToString(),
              vessel.mainBody.flightGlobalsIndex,
              main_body_degrees_of_freedom,
              degrees_of_freedom,
              (WXYZ)(UnityEngine.QuaternionD)part.rb.rotation,
              (XYZ)(Vector3d)part.rb.angularVelocity,
              Δt);
          if (part_id_to_intrinsic_torque_.ContainsKey(part.flightID)) {
            plugin_.PartApplyIntrinsicTorque(
                part.flightID,
                (XYZ)part_id_to_intrinsic_torque_[part.flightID]);
          }
          if (part_id_to_intrinsic_force_.ContainsKey(part.flightID)) {
            // When a Kerbal is doing an EVA and holding on to a ladder, the
            // ladder imbues them with their weight at the location of the
            // vessel to which the ladder is attached.  This leads to fantastic
            // effects where doing an EVA accelerates the vessel, see #1415.
            // Just say no to stupidity.
            if (!(vessel.isEVA && vessel.evaController.OnALadder)) {
              plugin_.PartApplyIntrinsicForce(
                  part.flightID,
                  (XYZ)part_id_to_intrinsic_force_[part.flightID]);
            }
          }
          if (part_id_to_intrinsic_forces_.ContainsKey(part.flightID)) {
            foreach (var part_centred_force in
                     part_id_to_intrinsic_forces_[part.flightID]) {
              plugin_.PartApplyIntrinsicForceAtPosition(
                  part.flightID,
                  (XYZ)part_centred_force.force,
                  (XYZ)part_centred_force.lever_arm);
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
          uint old_id = part.flightID;
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
          if (part1.Modules.OfType<ModuleWheelBase>().Any(
                  wheel => wheel.isGrounded)) {
            Log.Info("Reporting grounded wheel");
            plugin_.ReportGroundCollision(
                closest_physical_parent(part1).flightID);
          }
#if KSP_VERSION_1_9_1
          foreach (var collider in part1.currentCollisions.Keys) {
#elif KSP_VERSION_1_7_3
          foreach (var collider in part1.currentCollisions) {
#endif
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
              if (is_unready_kerbal(vessel2)) {
                // A Kerbal that just started an EVA and is not ready yet.  If
                // we let it collide with |vessel1|, it will cause that vessel
                // to become unmanageable and be removed from the plugin.  Of
                // course, the vessel will come back once the Kerbal is ready,
                // but at that point we'll have to trust the position given by
                // KSP, and in practice that might cause the vessel to jump
                // by hundreds of metres (#2590).  The bottom line is: we
                // better ignore the collision, the Kerbal will soon become
                // ready anyway.
              } else if (is_manageable(vessel2)) {
                plugin_.ReportPartCollision(
                    closest_physical_parent(part1).flightID,
                    closest_physical_parent(part2).flightID);
              } else {
                Log.Info("Reporting collision with the unmanageable vessel " +
                         vessel2.vesselName + " (reason: " +
                         UnmanageabilityReasons(vessel2) + ")");
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
      foreach (Part part in vessel.parts.Where(PartIsFaithful)) {
        if (main_body_change_countdown_ == 0 &&
            last_main_body_ == FlightGlobals.ActiveVessel?.mainBody) {
          plugin_.PartSetApparentRigidMotion(
                part.flightID,
                // TODO(egg): use the centre of mass.
                new QP{q = (XYZ)(Vector3d)part.rb.position,
                       p = (XYZ)(Vector3d)part.rb.velocity},
                (WXYZ)(UnityEngine.QuaternionD)part.rb.rotation,
                (XYZ)(Vector3d)part.rb.angularVelocity);
        }
      }
    }

    // Advance the lagging vessels and kill those which collided with a
    // celestial.
    {
      plugin_.CatchUpLaggingVessels(out DisposableIterator collided_vessels);
      for (; !collided_vessels.IteratorAtEnd();
            collided_vessels.IteratorIncrement()) {
        var vessel_guid = new Guid(collided_vessels.IteratorGetVesselGuid());
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
        foreach (Part part in vessel.parts.Where(PartIsFaithful)) {
          QPRW part_actual_motion =
              plugin_.PartGetActualRigidMotion(
                  part.flightID,
                  new Origin{
                      reference_part_is_at_origin  =
                                 FloatingOrigin.fetch.continuous,
                      reference_part_is_unmoving =
                          krakensbane_.FrameVel != Vector3d.zero,
                      main_body_centre_in_world =
                          (XYZ)FlightGlobals.ActiveVessel.mainBody.position,
                      reference_part_id =
                          FlightGlobals.ActiveVessel.rootPart.flightID});
          if (part == FlightGlobals.ActiveVessel.rootPart) {
            QP part_actual_degrees_of_freedom = part_actual_motion.qp;
            q_correction_at_root_part =
                (Vector3d)part_actual_degrees_of_freedom.q - part.rb.position;
            v_correction_at_root_part =
                (Vector3d)part_actual_degrees_of_freedom.p - part.rb.velocity;
          }

          // TODO(egg): use the centre of mass.  Here it's a bit tedious, some
          // transform nonsense must probably be done.
          // NOTE(egg): we must set the position and rotation of the |Transform|
          // as well as that of the |RigidBody| because we are performing this
          // correction after the physics step.
          // See https://github.com/mockingbirdnest/Principia/pull/1427,
          // https://github.com/mockingbirdnest/Principia/issues/1307#issuecomment-478337241.
          part.rb.position = (Vector3d)part_actual_motion.qp.q;
          part.rb.transform.position = (Vector3d)part_actual_motion.qp.q;
          part.rb.rotation = (UnityEngine.QuaternionD)part_actual_motion.r;
          part.rb.transform.rotation =
              (UnityEngine.QuaternionD)part_actual_motion.r;

          part.rb.velocity = (Vector3d)part_actual_motion.qp.p;
          part.rb.angularVelocity = (Vector3d)part_actual_motion.w;
        }
      }
      foreach (
          physicalObject physical_object in FlightGlobals.physicalObjects.Where(
              o => o != null && o.rb != null)) {
        // TODO(egg): This is no longer sensible.
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
      foreach (Vessel vessel in FlightGlobals.Vessels.Where(
          is_manageable_on_rails)) {
        vessel.SetPosition(vessel.transform.position + offset);
      }
      // NOTE(egg): this is almost certainly incorrect, since we give the
      // bodies their positions at the next instant, whereas KSP still expects
      // them at the previous instant, and will propagate them at the beginning
      // of the next frame...
    }

    if (last_main_body_ != FlightGlobals.ActiveVessel?.mainBody) {
      main_body_change_countdown_ = 1;
      last_main_body_ = FlightGlobals.ActiveVessel?.mainBody;
    } else if (main_body_change_countdown_ > 0) {
      --main_body_change_countdown_;
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
        if (!apocalypse_dialog_.Shown()) {
          if (plugin_.HasEncounteredApocalypse(out string revelation)) {
            apocalypse_dialog_.message = revelation;
            apocalypse_dialog_.Show();
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
        plugin_.FutureWaitForVesselToCatchUp(
            ref future,
            out DisposableIterator collided_vessels);
        for (; !collided_vessels.IteratorAtEnd();
             collided_vessels.IteratorIncrement()) {
          var vessel_guid = new Guid(collided_vessels.IteratorGetVesselGuid());
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
      part_id_to_intrinsic_torque_.Clear();
      part_id_to_intrinsic_force_.Clear();
      part_id_to_intrinsic_forces_.Clear();
      foreach (Vessel vessel in
               FlightGlobals.Vessels.Where(v => is_manageable(v) &&
                                                !v.packed)) {
        foreach (Part part in vessel.parts.Where(PartIsFaithful)) {
          if (part.torque != Vector3d.zero) {
            part_id_to_intrinsic_torque_.Add(part.flightID, part.torque);
          }
          if (part.force != Vector3d.zero) {
            part_id_to_intrinsic_force_.Add(part.flightID, part.force);
          }
          if (part.forces.Count > 0) {
            part_id_to_intrinsic_forces_.Add(
                part.flightID,
                (from force in part.forces
                 select PartCentredForceHolder.FromPartForceHolder(
                     part, force)).ToArray());
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
                  new PartCentredForceHolder[previous_holder.Length + 2];
              previous_holder.CopyTo(
                  part_id_to_intrinsic_forces_[physical_parent.flightID],
                  0);
            } else {
              part_id_to_intrinsic_forces_.Add(physical_parent.flightID,
                                               new PartCentredForceHolder[2]);
            }
            int lift_index = part_id_to_intrinsic_forces_[
                                 physical_parent.flightID].Length - 2;
            int drag_index = lift_index + 1;
            part_id_to_intrinsic_forces_[physical_parent.flightID][lift_index] =
                PartCentredForceHolder.FromPartForceHolder(
                    physical_parent,
                    new Part.ForceHolder {
                      force = part.partTransform.TransformDirection(
                                  part.bodyLiftLocalVector),
                      pos = part.partTransform.TransformPoint(
                                part.bodyLiftLocalPosition)});
            part_id_to_intrinsic_forces_[physical_parent.flightID][drag_index] =
                PartCentredForceHolder.FromPartForceHolder(
                    physical_parent,
                    new Part.ForceHolder {
                      force = -part.dragVectorDir * part.dragScalar,
                      pos = (physical_parent != part && PhysicsGlobals.
                                 ApplyDragToNonPhysicsPartsAtParentCoM)
                                ? physical_parent.rb.worldCenterOfMass
                                : part.partTransform.TransformPoint(
                                    part.CoPOffset)});
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
        foreach (Part part in vessel.parts.Where(PartIsFaithful)) {
          // TODO(egg): use the centre of mass.
          part_id_to_degrees_of_freedom_.Add(
              part.flightID,
              new QP{q = (XYZ)(Vector3d)part.rb.position,
                     p = (XYZ)(Vector3d)part.rb.velocity});
        }
      }
    }
  }

  private void BetterLateThanNeverLateUpdate() {
    // While we draw the trajectories directly (and thus do so after everything
    // else has been rendered), we rely on the game to render its map nodes.
    // Since the screen position is determined in |MapNode.NodeUpdate|, it must
    // be called before rendering occurs, but after the cameras have moved;
    // otherwise, the map nodes will lag behind when the camera is moved.
    // The only timing that satisfies these constraints is BetterLateThanNever
    // in LateUpdate.
    string main_vessel_guid = PredictedVessel()?.id.ToString();
    if (MapView.MapIsEnabled && main_vessel_guid != null) {
      XYZ sun_world_position = (XYZ)Planetarium.fetch.Sun.position;
      RenderPredictionMarkers(main_vessel_guid, sun_world_position);
      string target_id =
          FlightGlobals.fetch.VesselTarget?.GetVessel()?.id.ToString();
      if (FlightGlobals.ActiveVessel != null &&
          !plotting_frame_selector_.target_override && target_id != null &&
          plugin_.HasVessel(target_id)) {
        RenderPredictionMarkers(target_id, sun_world_position);
      }
      if (plugin_.FlightPlanExists(main_vessel_guid)) {
        RenderFlightPlanMarkers(main_vessel_guid, sun_world_position);
      }
    }
    map_node_pool_.Update();
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
            Z = swizzly_body_world_to_world * new Vector3d{x = 0, y = 0, z = 1}
        };
      }
    }
  }

  private void RenderGuidance(Vessel active_vessel) {
    string vessel_guid = active_vessel.id.ToString();
    if (plugin_.HasVessel(vessel_guid) &&
        plugin_.FlightPlanExists(vessel_guid)) {
      // Here the vessel is known to the plugin and has a flight plan.
      // This duplicates a bit of code in FlightPlanner.
      // UpdateVesselAndBurnEditors but it's probably not worth factoring out.
      double current_time = plugin_.CurrentTime();
      // Note that we don't want to look at the anomalous manœuvres as we may
      // not even be able to build a Frenet frame for them.
      int number_of_nomalous_manœuvres =
          plugin_.FlightPlanNumberOfManoeuvres(vessel_guid) -
          plugin_.FlightPlanNumberOfAnomalousManoeuvres(vessel_guid);
      int? first_future_manœuvre_index = null;
      for (int i = 0; i < number_of_nomalous_manœuvres; ++i) {
        NavigationManoeuvre manœuvre =
            plugin_.FlightPlanGetManoeuvre(vessel_guid, i);
        if (current_time < manœuvre.final_time) {
          first_future_manœuvre_index = i;
          break;
        }
      }
      if (first_future_manœuvre_index.HasValue) {
        // Here the flight plan has a manœuvre in the future.
        XYZ guidance = plugin_.FlightPlanGetGuidance(
                           vessel_guid,
                           first_future_manœuvre_index.Value);
        Burn burn = plugin_.FlightPlanGetManoeuvre(
                        vessel_guid,
                        first_future_manœuvre_index.Value).burn;
        if (flight_planner_.show_guidance && !IsNaN(guidance)) {
          // The user wants to show the guidance node, and that node was
          // properly computed by the C++ code.
          PatchedConicSolver solver = active_vessel.patchedConicSolver;
          if (guidance_node_ == null ||
              !solver.maneuverNodes.Contains(guidance_node_)) {
            while (solver.maneuverNodes.Count > 0) {
              solver.maneuverNodes.Last().RemoveSelf();
            }
            guidance_node_ = solver.AddManeuverNode(burn.initial_time);
          } else {
            while (solver.maneuverNodes.Count > 1) {
              if (solver.maneuverNodes.First() == guidance_node_) {
                solver.maneuverNodes.Last().RemoveSelf();
              } else {
                solver.maneuverNodes.First().RemoveSelf();
              }
            }
          }
          var stock_orbit = guidance_node_.patch;
          Vector3d stock_velocity_at_node_time =
              stock_orbit.getOrbitalVelocityAtUT(burn.initial_time).xzy;
          Vector3d stock_displacement_from_parent_at_node_time =
              stock_orbit.getRelativePositionAtUT(burn.initial_time).xzy;
          UnityEngine.Quaternion stock_frenet_frame_to_world =
              UnityEngine.Quaternion.LookRotation(
                  stock_velocity_at_node_time,
                  Vector3d.Cross(
                      stock_velocity_at_node_time,
                      stock_displacement_from_parent_at_node_time));
          guidance_node_.DeltaV =
              ((Vector3d)burn.delta_v).magnitude *
               (Vector3d)(UnityEngine.Quaternion.Inverse(
                              stock_frenet_frame_to_world) *
               (Vector3d)guidance);
          guidance_node_.UT = burn.initial_time;
          solver.UpdateFlightPlan();
          // Return here after setting the guidance node.  All other paths will
          // clear the guidance node.
          return;
        }
      }
    }
    if (guidance_node_ != null) {
      guidance_node_.RemoveSelf();
      guidance_node_ = null;
    }
  }

  private void RenderNavball(Vessel active_vessel) {
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
                              plotting_frame_selector_
                                  .selected_celestial.flightGlobalsIndex);
      if (plotting_frame_selector_.target_override != target_vessel) {
        navball_changed_ = true;
        planetarium_camera_adjuster_.should_transfer_camera_coordinates = true;
        plotting_frame_selector_.target_override = target_vessel;
      }
    } else {
      plugin_.ClearTargetVessel();
      if (plotting_frame_selector_.target_override != null) {
        navball_changed_ = true;
        planetarium_camera_adjuster_.should_transfer_camera_coordinates = true;
        plotting_frame_selector_.target_override = null;
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
          plotting_frame_selector_.SetFrameType(
              ReferenceFrameSelector.FrameType.BODY_SURFACE);
          break;
        case FlightGlobals.SpeedDisplayModes.Orbit:
          plotting_frame_selector_.SetFrameType(
              last_non_surface_frame_type_);
          break;
      }
    }

    if (navball_changed_ && previous_display_mode_ != null) {
      // Texture the ball.
      navball_changed_ = false;
      if (plotting_frame_selector_.target_override) {
        set_navball_texture(target_navball_texture_);
      } else {
        // If we are targeting an unmanageable vessel, keep the navball in
        // target mode; otherwise, put it in the mode that reflects the
        // plotting frame.
        if (FlightGlobals.speedDisplayMode !=
            FlightGlobals.SpeedDisplayModes.Target) {
          if (plotting_frame_selector_.frame_type ==
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
        switch (plotting_frame_selector_.frame_type) {
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

    RenderGuidance(active_vessel);
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

  private void RemoveStockTrajectoriesIfNeeded(CelestialBody celestial) {
    if (celestial.orbitDriver == null) {
      return;
    }
    celestial.orbitDriver.Renderer.drawMode =
        main_window_.display_patched_conics
            ? OrbitRenderer.DrawMode.REDRAW_AND_RECALCULATE
            : OrbitRenderer.DrawMode.OFF;
  }

  private void RemoveStockTrajectoriesIfNeeded(Vessel vessel) {
    if (vessel.patchedConicRenderer != null) {
      vessel.patchedConicRenderer.relativityMode =
          PatchRendering.RelativityMode.RELATIVE;
    }

    if (main_window_.display_patched_conics || !is_manageable(vessel)) {
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
      main_window_.SelectTargetCelestial(node.mapObject);
    }
  }

  private void OnVesselNodeClick(KSP.UI.Screens.Mapview.MapNode node,
                                 Mouse.Buttons buttons) {
    main_window_.SelectActiveVesselTarget(
        node.mapObject,
        set_planetarium_camera : buttons == Mouse.Buttons.Left);
  }

  private void HandleMapViewClicks() {
    if (InputLockManager.IsUnlocked(ControlTypes.MAP_UI) &&
        !UnityEngine.EventSystems.EventSystem.current
             .IsPointerOverGameObject() &&
        Mouse.Left.GetClick() && !ManeuverGizmo.HasMouseFocus &&
        !main_window_.selecting_active_vessel_target) {
      var ray = PlanetariumCamera.Camera.ScreenPointToRay(
          UnityEngine.Input.mousePosition);
      foreach (var celestial in FlightGlobals.Bodies) {
        double scaled_distance =
            Vector3d.Cross(ray.direction,
                           ScaledSpace.LocalToScaledSpace(celestial.position) -
                               ray.origin).magnitude;
        if (scaled_distance * ScaledSpace.ScaleFactor < celestial.Radius) {
          main_window_.SelectTargetCelestial(celestial.MapObject);
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
      RemoveStockTrajectoriesIfNeeded(celestial);
    }
    foreach (var vessel in FlightGlobals.Vessels.Where(
                               v => v.mapObject?.uiNode != null)) {
      // There is no way to check if we have already added a callback to an
      // event...
      vessel.mapObject.uiNode.OnClick -= OnVesselNodeClick;
      vessel.mapObject.uiNode.OnClick += OnVesselNodeClick;
      RemoveStockTrajectoriesIfNeeded(vessel);
    }
    string main_vessel_guid = PredictedVessel()?.id.ToString();
    if (MapView.MapIsEnabled) {
      XYZ sun_world_position = (XYZ)Planetarium.fetch.Sun.position;
      using (DisposablePlanetarium planetarium =
                GLLines.NewPlanetarium(plugin_, sun_world_position)) {
        GLLines.Draw(() => {
          PlotCelestialTrajectories(planetarium, main_vessel_guid);

          // Vessel trajectories.
          if (main_vessel_guid == null) {
            return;
          }
          // Main vessel psychohistory and prediction.
          using (DisposableIterator rp2_lines_iterator =
                    planetarium.PlanetariumPlotPsychohistory(
                        plugin_,
                        чебышёв_plotting_method_,
                        main_vessel_guid,
                        main_window_.history_length)) {
            GLLines.PlotRP2Lines(rp2_lines_iterator,
                                 XKCDColors.Lime,
                                 GLLines.Style.Faded);
          }
          using (DisposableIterator rp2_lines_iterator =
                    planetarium.PlanetariumPlotPrediction(
                        plugin_,
                        чебышёв_plotting_method_,
                        main_vessel_guid)) {
            GLLines.PlotRP2Lines(rp2_lines_iterator,
                                 XKCDColors.Fuchsia,
                                 GLLines.Style.Solid);
          }
          // Target psychohistory and prediction.
          string target_id =
              FlightGlobals.fetch.VesselTarget?.GetVessel()?.id.ToString();
          if (FlightGlobals.ActiveVessel != null &&
              !plotting_frame_selector_.target_override &&
              target_id != null && plugin_.HasVessel(target_id)) {
            using (DisposableIterator rp2_lines_iterator =
                      planetarium.PlanetariumPlotPsychohistory(
                          plugin_,
                          чебышёв_plotting_method_,
                          target_id,
                          main_window_.history_length)) {
              GLLines.PlotRP2Lines(rp2_lines_iterator,
                                   XKCDColors.Goldenrod,
                                   GLLines.Style.Faded);
            }
            using (DisposableIterator rp2_lines_iterator =
                      planetarium.PlanetariumPlotPrediction(
                          plugin_,
                          чебышёв_plotting_method_,
                          target_id)) {
              GLLines.PlotRP2Lines(rp2_lines_iterator,
                                   XKCDColors.LightMauve,
                                   GLLines.Style.Solid);
            }
          }
          // Main vessel flight plan.
          if (plugin_.FlightPlanExists(main_vessel_guid)) {
            int number_of_anomalous_manœuvres =
                plugin_.FlightPlanNumberOfAnomalousManoeuvres(main_vessel_guid);
            int number_of_manœuvres =
                plugin_.FlightPlanNumberOfManoeuvres(main_vessel_guid);
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
                      is_burn ? GLLines.Style.Solid : GLLines.Style.Dashed);
                }
                if (is_burn) {
                  int manœuvre_index = i / 2;
                  if (manœuvre_index <
                      number_of_manœuvres - number_of_anomalous_manœuvres) {
                    NavigationManoeuvreFrenetTrihedron manœuvre =
                        plugin_.FlightPlanGetManoeuvreFrenetTrihedron(
                            main_vessel_guid,
                            manœuvre_index);
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
                    add_vector(manœuvre.tangent, Style.Tangent);
                    add_vector(manœuvre.normal, Style.Normal);
                    add_vector(manœuvre.binormal, Style.Binormal);
                  }
                }
              }
            }
          }
        });
      }
    }
  }

  private void PlotCelestialTrajectories(DisposablePlanetarium planetarium,
                                         string main_vessel_guid) {
    foreach (CelestialBody celestial in FlightGlobals.Bodies) {
      if (plotting_frame_selector_.FixedBodies().Contains(celestial)) {
        continue;
      }
      var colour = celestial.MapObject?.uiNode?.VisualIconData.color ??
          XKCDColors.SunshineYellow;
      if (colour.a != 1) {
        // When zoomed into a planetary system, the trajectory of the
        // planet is hidden in stock (because KSP then draws most things
        // in the reference frame centred on that planet).
        // Here we still want to display the trajectory of the primary,
        // e.g., if we are drawing the trajectories of the Jovian system
        // in the heliocentric frame.
        foreach (CelestialBody child in celestial.orbitingBodies) {
          colour.a = Math.Max(
              child.MapObject?.uiNode?.VisualIconData.color.a ?? 1,
              colour.a);
        }
      }
      if (colour.a == 0) {
        continue;
      }
      using (DisposableIterator rp2_lines_iterator =
                planetarium.PlanetariumPlotCelestialTrajectoryForPsychohistory(
                    plugin_,
                    celestial.flightGlobalsIndex,
                    main_vessel_guid,
                    main_window_.history_length)) {
        GLLines.PlotRP2Lines(rp2_lines_iterator,
                              colour,
                              GLLines.Style.Faded);
      }
      if (main_vessel_guid != null) {
        using (DisposableIterator rp2_lines_iterator =
            planetarium.
                PlanetariumPlotCelestialTrajectoryForPredictionOrFlightPlan(
                    plugin_,
                    celestial.flightGlobalsIndex,
                    main_vessel_guid)) {
          GLLines.PlotRP2Lines(rp2_lines_iterator, colour, GLLines.Style.Solid);
        }
      }
    }
  }

  private void RenderPredictionMarkers(string vessel_guid,
                                       XYZ sun_world_position) {
    if (plotting_frame_selector_.target_override) {
      plugin_.RenderedPredictionNodes(
          vessel_guid,
          sun_world_position,
          MapNodePool.MaxRenderedNodes,
          out DisposableIterator ascending_nodes_iterator,
          out DisposableIterator descending_nodes_iterator);
      plugin_.RenderedPredictionClosestApproaches(
          vessel_guid,
          sun_world_position,
          MapNodePool.MaxRenderedNodes,
          out DisposableIterator approaches_iterator);
      map_node_pool_.RenderMarkers(
          ascending_nodes_iterator,
          MapObject.ObjectType.AscendingNode,
          MapNodePool.NodeSource.Prediction,
          plotting_frame_selector_);
      map_node_pool_.RenderMarkers(
          descending_nodes_iterator,
          MapObject.ObjectType.DescendingNode,
          MapNodePool.NodeSource.Prediction,
          plotting_frame_selector_);
      map_node_pool_.RenderMarkers(
          approaches_iterator,
          MapObject.ObjectType.ApproachIntersect,
          MapNodePool.NodeSource.Prediction,
          plotting_frame_selector_);
    } else {
      foreach (CelestialBody celestial in
               plotting_frame_selector_.FixedBodies()) {
        plugin_.RenderedPredictionApsides(
            vessel_guid,
            celestial.flightGlobalsIndex,
            sun_world_position,
            MapNodePool.MaxRenderedNodes,
            out DisposableIterator apoapsis_iterator,
            out DisposableIterator periapsis_iterator);
        map_node_pool_.RenderMarkers(
            apoapsis_iterator,
            MapObject.ObjectType.Apoapsis,
            MapNodePool.NodeSource.Prediction,
            plotting_frame_selector_);
        map_node_pool_.RenderMarkers(
            periapsis_iterator,
            MapObject.ObjectType.Periapsis,
            MapNodePool.NodeSource.Prediction,
            plotting_frame_selector_);
      }
      plugin_.RenderedPredictionNodes(
          vessel_guid,
          sun_world_position,
          MapNodePool.MaxRenderedNodes,
          out DisposableIterator ascending_nodes_iterator,
          out DisposableIterator descending_nodes_iterator);
      map_node_pool_.RenderMarkers(
          ascending_nodes_iterator,
          MapObject.ObjectType.AscendingNode,
          MapNodePool.NodeSource.Prediction,
          plotting_frame_selector_);
      map_node_pool_.RenderMarkers(
          descending_nodes_iterator,
          MapObject.ObjectType.DescendingNode,
          MapNodePool.NodeSource.Prediction,
          plotting_frame_selector_);
    }
  }

  private void RenderFlightPlanMarkers(string vessel_guid,
                                       XYZ sun_world_position) {
    if (plotting_frame_selector_.target_override) {
      plugin_.FlightPlanRenderedNodes(
          vessel_guid,
          sun_world_position,
          MapNodePool.MaxRenderedNodes,
          out DisposableIterator ascending_nodes_iterator,
          out DisposableIterator descending_nodes_iterator);
      plugin_.FlightPlanRenderedClosestApproaches(
          vessel_guid,
          sun_world_position,
          MapNodePool.MaxRenderedNodes,
          out DisposableIterator approaches_iterator);
      map_node_pool_.RenderMarkers(
          ascending_nodes_iterator,
          MapObject.ObjectType.AscendingNode,
          MapNodePool.NodeSource.FlightPlan,
          plotting_frame_selector_);
      map_node_pool_.RenderMarkers(
          descending_nodes_iterator,
          MapObject.ObjectType.DescendingNode,
          MapNodePool.NodeSource.FlightPlan,
          plotting_frame_selector_);
      map_node_pool_.RenderMarkers(
          approaches_iterator,
          MapObject.ObjectType.ApproachIntersect,
          MapNodePool.NodeSource.FlightPlan,
          plotting_frame_selector_);
    } else {
      foreach (CelestialBody celestial in
               plotting_frame_selector_.FixedBodies()) {
        plugin_.FlightPlanRenderedApsides(
            vessel_guid,
            celestial.flightGlobalsIndex,
            sun_world_position,
            MapNodePool.MaxRenderedNodes,
            out DisposableIterator apoapsis_iterator,
            out DisposableIterator periapsis_iterator);
        map_node_pool_.RenderMarkers(
            apoapsis_iterator,
            MapObject.ObjectType.Apoapsis,
            MapNodePool.NodeSource.FlightPlan,
            plotting_frame_selector_);
        map_node_pool_.RenderMarkers(
            periapsis_iterator,
            MapObject.ObjectType.Periapsis,
            MapNodePool.NodeSource.FlightPlan,
            plotting_frame_selector_);
      }
      plugin_.FlightPlanRenderedNodes(
          vessel_guid,
          sun_world_position,
          MapNodePool.MaxRenderedNodes,
          out DisposableIterator ascending_nodes_iterator,
          out DisposableIterator descending_nodes_iterator);
      map_node_pool_.RenderMarkers(
          ascending_nodes_iterator,
          MapObject.ObjectType.AscendingNode,
          MapNodePool.NodeSource.FlightPlan,
          plotting_frame_selector_);
      map_node_pool_.RenderMarkers(
          descending_nodes_iterator,
          MapObject.ObjectType.DescendingNode,
          MapNodePool.NodeSource.FlightPlan,
          plotting_frame_selector_);
    }
  }

  private void Cleanup() {
    UnityEngine.Object.Destroy(map_renderer_);
    map_renderer_ = null;
    map_node_pool_.Clear();
    LockClearing();
    Interface.DeletePlugin(ref plugin_);
    previous_display_mode_ = null;
    navball_changed_ = true;

    // Load the flags.
    Interface.ClearFlags();
    ConfigNode.ValueList flags =
        GameDatabase.Instance.GetAtMostOneNode(principia_flags_)?.values;
    if (flags != null) {
      foreach (ConfigNode.Value flag in flags) {
        Interface.SetFlag(flag.name, flag.value);
      }
    }
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
    planetarium_camera_adjuster_.should_transfer_camera_coordinates = true;
  }

  private static void InitializeIntegrators(
      IntPtr plugin,
      ConfigNode numerics_blueprint) {
    if (numerics_blueprint == null) {
      return;
    }
    ConfigNode ephemeris_parameters =
        numerics_blueprint.GetAtMostOneNode("ephemeris");
    if (ephemeris_parameters != null) {
      plugin.InitializeEphemerisParameters(
          ConfigNodeParsers.NewConfigurationAccuracyParameters(
              ephemeris_parameters),
          ConfigNodeParsers.NewConfigurationFixedStepParameters(
              ephemeris_parameters));
    }

    ConfigNode history_parameters =
        numerics_blueprint.GetAtMostOneNode("history");
    if (history_parameters != null) {
      plugin.InitializeHistoryParameters(
          ConfigNodeParsers.NewConfigurationFixedStepParameters(
              history_parameters));
    }

    ConfigNode psychohistory_parameters =
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
    Dictionary<string, ConfigNode> name_to_gravity_model = null;
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
        if (!name_to_gravity_model.TryGetValue(
                body.name,
                out ConfigNode body_gravity_model)) {
          Log.Fatal("missing gravity model for " + body.name);
        }
        if (!name_to_initial_state.TryGetValue(
                body.name,
                out ConfigNode body_initial_state)) {
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
            celestial_index : body.flightGlobalsIndex,
            parent_index    : parent_index,
            body_parameters : body_parameters,
            x               : body_initial_state.GetUniqueValue("x"),
            y               : body_initial_state.GetUniqueValue("y"),
            z               : body_initial_state.GetUniqueValue("z"),
            vx              : body_initial_state.GetUniqueValue("vx"),
            vy              : body_initial_state.GetUniqueValue("vy"),
            vz              : body_initial_state.GetUniqueValue("vz"));
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
            celestial_index    : body.flightGlobalsIndex,
            parent_index       : orbit?.referenceBody.flightGlobalsIndex,
            body_parameters    : body_parameters,
            keplerian_elements : orbit?.Elements());
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
    must_set_plotting_frame_ = true;
  } catch (Exception e) {
    Log.Fatal($"Exception while resetting plugin: {e}");
  }
  }

  private void RemoveBuggyTidalLocking() {
    ApplyToBodyTree(body => body.tidallyLocked = false);
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
