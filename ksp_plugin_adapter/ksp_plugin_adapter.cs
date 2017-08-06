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
      
  private const String next_release_name = "Чебышёв";
  private const int next_release_lunation_number = 218;
  private DateTimeOffset next_release_date =
      new DateTimeOffset(2017, 08, 21, 18, 30, 11, TimeSpan.Zero);

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

  [KSPField(isPersistant = true)]
  private bool use_cayley_plotting_ = true;
  [KSPField(isPersistant = true)]
  private bool use_чебышёв_plotting_ = false;
  [KSPField(isPersistant = true)]
  private int чебышёв_plotting_method_ = 1;
  private const int чебышёв_plotting_methods_count = 3;

  internal Controlled<ReferenceFrameSelector> plotting_frame_selector_;
  private Controlled<FlightPlanner> flight_planner_;
  private MapNodePool map_node_pool_;

  private bool selecting_active_vessel_target_ = false;
  private bool selecting_target_celestial_ = false;

  private IntPtr plugin_ = IntPtr.Zero;

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
  private FlightGlobals.SpeedDisplayModes previous_display_mode_;
  private ReferenceFrameSelector.FrameType last_non_surface_frame_type_ =
      ReferenceFrameSelector.FrameType.BODY_CENTRED_NON_ROTATING;

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

  private bool is_manageable(Vessel vessel) {
    return UnmanageabilityReasons(vessel) == null;
  }

  private string UnmanageabilityReasons(Vessel vessel) {
    List<string> reasons = new List<string>(capacity : 3);
    if (vessel.state == Vessel.State.DEAD) {
      reasons.Add("vessel is dead");
    }
    if (!(vessel.situation == Vessel.Situations.SUB_ORBITAL ||
          vessel.situation == Vessel.Situations.ORBITING ||
          vessel.situation == Vessel.Situations.ESCAPING ||
          (vessel.situation == Vessel.Situations.FLYING && vessel.packed &&
           (vessel.altitude > vessel.mainBody.inverseRotThresholdAltitude ||
            vessel.altitude == -vessel.mainBody.Radius)))) {
      reasons.Add("vessel situation is " + vessel.situation +
                  " and vessel is " + (vessel.packed ? "packed" : "unpacked") +
                  " at an altitude of " + vessel.altitude + " m above " +
                  vessel.mainBody.NameWithArticle() + " whose threshold is " +
                  vessel.mainBody.inverseRotThresholdAltitude + " m");
    }
    if (!vessel.packed &&
        vessel.altitude <= vessel.mainBody.inverseRotThresholdAltitude) {
      reasons.Add("vessel is unpacked at an altitude of " + vessel.altitude +
                  " m above " + vessel.mainBody.NameWithArticle() +
                  ", below the threshold of " +
                  vessel.mainBody.inverseRotThresholdAltitude + " m");
    }
    if (!vessel.packed && FlightGlobals.RefFrameIsRotating) {
      reasons.Add("vessel is unpacked and frame is rotating");
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
    // Timing5, 8008.
    TimingManager.FixedUpdateAdd(TimingManager.TimingStage.BetterLateThanNever,
                                 StorePartDegreesOfFreedom);
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
    if (bad_installation_popup_ != null) {
      UnityEngine.Debug.LogError("Spawning: " + bad_installation_popup_);
      // No-one seems to understand what |anchorMin| and |anchorMax| do at this
      // time.
      PopupDialog.SpawnPopupDialog(
          anchorMin           : default(UnityEngine.Vector2),
          anchorMax           : default(UnityEngine.Vector2),
#if KSP_VERSION_1_3
          dialogName          : "Principia error",
#endif
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
          id         : this.GetHashCode() + 1,
          screenRect : apocalypse_window_rectangle_,
          func       : (int id) => {
              using (new VerticalLayout()) {
                UnityEngine.GUILayout.TextArea(revelation_);
              }
              UnityEngine.GUI.DragWindow(
                  position : new UnityEngine.Rect(x      : 0f,
                                                  y      : 0f,
                                                  width  : 10000f,
                                                  height : 10000f));
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

      if (navball_changed_) {
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
    if (GameSettings.ORBIT_WARP_DOWN_AT_SOI) {
      Log.Info("Setting GameSettings.ORBIT_WARP_DOWN_AT_SOI to false");
      GameSettings.ORBIT_WARP_DOWN_AT_SOI = false;
    }
    if (must_set_plotting_frame_ && FlightGlobals.currentMainBody != null) {
      must_set_plotting_frame_ = false;
      plotting_frame_selector_.reset(new ReferenceFrameSelector(
          this, plugin_, UpdateRenderingFrame, "Plotting frame"));
    }

    if (PluginRunning()) {
      double universal_time = Planetarium.GetUniversalTime();

      plugin_.SetMainBody(
          (FlightGlobals.currentMainBody
               ?? FlightGlobals.GetHomeBody()).flightGlobalsIndex);

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
        plugin_.SetPredictionLength(double.PositiveInfinity);
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
      plugin_.ForgetAllHistoriesBefore(universal_time -
                                       history_lengths_[history_length_index_]);
      // TODO(egg): Set the degrees of freedom of the origin of |World| (by
      // toying with Krakensbane and FloatingOrigin) here.

      // Now we let the game and Unity do their thing. among other things,
      // the FashionablyLate callbacks, including ReportNonConservativeForces,
      // then the FlightIntegrator's FixedUpdate will run, then the Vessel's,
      // and eventually the physics simulation.
      StartCoroutine(
          AdvanceAndNudgeVesselsAfterPhysicsSimulation(universal_time));
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
    WindowUtilities.ClearLock(this);
    Cleanup();
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.ObscenelyEarly,
                                    DisableVesselPrecalculate);
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.Precalc,
                                    SetBodyFramesAndPrecalculateVessels);
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.Early,
                                    UpdateVesselOrbits);
    TimingManager.FixedUpdateRemove(TimingManager.TimingStage.FashionablyLate,
                                    ReportVesselsAndParts);
    TimingManager.FixedUpdateRemove(
        TimingManager.TimingStage.BetterLateThanNever,
        StorePartDegreesOfFreedom);
  }

  #endregion

  private System.Collections.IEnumerator
  AdvanceAndNudgeVesselsAfterPhysicsSimulation(double universal_time) {
    yield return new UnityEngine.WaitForFixedUpdate();

  try {
    // Unity's physics has just finished doing its thing.  If we correct the
    // positions here, nobody will know that they're not the ones obtained by
    // Unity.

    if (!time_is_advancing_) {
      yield break;
    }

    double Δt = Planetarium.TimeScale * Planetarium.fetch.fixedDeltaTime;

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
        QP main_body_degrees_of_freedom =
            new QP{q = (XYZ)vessel.mainBody.position,
                   p = (XYZ)(-krakensbane.FrameVel)};
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
    foreach (Vessel vessel1 in FlightGlobals.VesselsLoaded) {
      if (plugin_.HasVessel(vessel1.id.ToString()) && !vessel1.packed) {
        if (vessel1.isEVA && vessel1.evaController.OnALadder) {
          var vessel2 = vessel1.evaController.LadderPart.vessel;
          if (vessel2 != null && plugin_.HasVessel(vessel2.id.ToString()) &&
              !vessel2.packed) {
            plugin_.ReportCollision(
                vessel1.rootPart.flightID,
                closest_physical_parent(
                    vessel1.evaController.LadderPart).flightID);
          }
        }
        foreach (Part part1 in vessel1.parts) {
          foreach (var collider in part1.currentCollisions) {
            if (collider == null) {
              // This happens, albeit quite rarely, see #1447.  When it happens,
              // the null collider remains in |currentCollisions| until the next
              // scene change, so we do not log, otherwise we would spam.
              continue;
            }
            var part2 = collider.gameObject.GetComponentUpwards<Part>();
            var vessel2 = part2?.vessel;
            if (vessel2 != null && plugin_.HasVessel(vessel2.id.ToString())) {
              plugin_.ReportCollision(closest_physical_parent(part1).flightID,
                                      closest_physical_parent(part2).flightID);
            }
          }
        }
      }
    }

    plugin_.FreeVesselsAndPartsAndCollectPileUps(Δt);

    foreach (Vessel vessel in FlightGlobals.VesselsLoaded) {
      if (vessel.packed || !plugin_.HasVessel(vessel.id.ToString())) {
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

    plugin_.CatchUpLaggingVessels();

    // We don't want to do too many things here, since all the KSP classes
    // still think they're in the preceding step.  We only nudge the Unity
    // transforms of loaded vessels & their parts.
    if (has_active_manageable_vessel() && !FlightGlobals.ActiveVessel.packed &&
        plugin_.HasVessel(FlightGlobals.ActiveVessel.id.ToString())) {
      Vector3d q_correction_at_root_part = Vector3d.zero;
      Vector3d v_correction_at_root_part = Vector3d.zero;
      foreach (Vessel vessel in FlightGlobals.VesselsLoaded) {
        // TODO(egg): if I understand anything, there should probably be a
        // special treatment for loaded packed vessels.  I don't understand
        // anything though.
        if (vessel.packed || !plugin_.HasVessel(vessel.id.ToString())) {
          continue;
        }
        foreach (Part part in vessel.parts) {
          if (part.rb == null) {
            continue;
          }
          QP part_actual_degrees_of_freedom =
              plugin_.GetPartActualDegreesOfFreedom(
                  part.flightID, FlightGlobals.ActiveVessel.rootPart.flightID);
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
          FlightGlobals.ActiveVessel.rootPart.flightID,
          universal_time);
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

  private void DisableVesselPrecalculate() {
    foreach (var vessel in
             FlightGlobals.Vessels.Where(vessel => vessel.precalc != null)) {
      vessel.precalc.enabled = false;
    }
  }

  private void SetBodyFramesAndPrecalculateVessels() {
    if (FlightGlobals.ActiveVessel?.situation == Vessel.Situations.PRELAUNCH &&
        FlightGlobals.ActiveVessel?.orbitDriver?.lastMode ==
            OrbitDriver.UpdateMode.TRACK_Phys &&
        FlightGlobals.ActiveVessel?.orbitDriver?.updateMode ==
            OrbitDriver.UpdateMode.IDLE) {
      Log.Info("Skipping AdvanceTime and SetBodyFrames while waiting for the " +
               "vessel to be fully ready (see #1421).");
    } else {
      if (PluginRunning()) {
        double plugin_time = plugin_.CurrentTime();
        double universal_time = Planetarium.GetUniversalTime();
        time_is_advancing_ = time_is_advancing(universal_time);
        if (time_is_advancing_) {
          plugin_.AdvanceTime(universal_time, Planetarium.InverseRotAngle);
          if (!is_post_apocalyptic_) {
            is_post_apocalyptic_ |=
                plugin_.HasEncounteredApocalypse(out revelation_);
          }
          foreach (var vessel in FlightGlobals.Vessels) {
            if (vessel.packed && plugin_.HasVessel(vessel.id.ToString())) {
              plugin_.CatchUpVessel(vessel.id.ToString());
            }
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

  private void UpdateVesselOrbits() {
    if (PluginRunning()) {
      ApplyToVesselsOnRails(
          vessel => UpdateVessel(vessel, Planetarium.GetUniversalTime()));
    }
  }

  private void ReportVesselsAndParts() {
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

  private void StorePartDegreesOfFreedom() {
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
      IntPtr planetarium = GLLines.NewPlanetarium(plugin_, sun_world_position);
      try {
        GLLines.Draw(() => {
          if (use_cayley_plotting_) {
            GLLines.RenderAndDeleteTrajectory(
                plugin_.RenderedVesselTrajectory(main_vessel_guid,
                                                 sun_world_position),
                XKCDColors.AcidGreen,
                GLLines.Style.FADED);
          }
          if (use_чебышёв_plotting_) {
            IntPtr rp2_lines_iterator =
                planetarium.PlanetariumPlotPsychohistory(
                    plugin_,
                    чебышёв_plotting_method_,
                    main_vessel_guid);
            GLLines.PlotAndDeleteRP2Lines(rp2_lines_iterator,
                                          XKCDColors.Banana,
                                          GLLines.Style.FADED);
          }
          RenderPredictionMarkers(main_vessel_guid, sun_world_position);
          if (use_cayley_plotting_) {
            GLLines.RenderAndDeleteTrajectory(
                plugin_.RenderedPrediction(main_vessel_guid,
                                           sun_world_position),
                XKCDColors.Fuchsia,
                GLLines.Style.SOLID);
          }
          if (use_чебышёв_plotting_) {
            IntPtr rp2_lines_iterator =
                planetarium.PlanetariumPlotPrediction(
                    plugin_,
                    чебышёв_plotting_method_,
                    main_vessel_guid);
            GLLines.PlotAndDeleteRP2Lines(rp2_lines_iterator,
                                          XKCDColors.Cerise,
                                          GLLines.Style.SOLID);
          }
          string target_id =
              FlightGlobals.fetch.VesselTarget?.GetVessel()?.id.ToString();
          if (!plotting_frame_selector_.get().target_override &&
              target_id != null && plugin_.HasVessel(target_id)) {
            if (use_cayley_plotting_) {
              GLLines.RenderAndDeleteTrajectory(
                  plugin_.RenderedVesselTrajectory(target_id,
                                                   sun_world_position),
                  XKCDColors.Goldenrod,
                  GLLines.Style.FADED);
            }
            if (use_чебышёв_plotting_) {
              IntPtr rp2_lines_iterator =
                  planetarium.PlanetariumPlotPsychohistory(
                      plugin_,
                      чебышёв_plotting_method_,
                      target_id);
              GLLines.PlotAndDeleteRP2Lines(rp2_lines_iterator,
                                            XKCDColors.Orange,
                                            GLLines.Style.FADED);
            }
            RenderPredictionMarkers(target_id, sun_world_position);
            if (use_cayley_plotting_) {
              GLLines.RenderAndDeleteTrajectory(
                  plugin_.RenderedPrediction(target_id, sun_world_position),
                  XKCDColors.PigPink,
                  GLLines.Style.SOLID);
            }
            if (use_чебышёв_plotting_) {
              IntPtr rp2_lines_iterator =
                  planetarium.PlanetariumPlotPrediction(
                      plugin_,
                      чебышёв_plotting_method_,
                      target_id);
              GLLines.PlotAndDeleteRP2Lines(rp2_lines_iterator,
                                            XKCDColors.Raspberry,
                                            GLLines.Style.SOLID);
            }
          }
          if (plugin_.FlightPlanExists(main_vessel_guid)) {
            RenderFlightPlanMarkers(main_vessel_guid, sun_world_position);

            int number_of_segments =
                plugin_.FlightPlanNumberOfSegments(main_vessel_guid);
            for (int i = 0; i < number_of_segments; ++i) {
              bool is_burn = i % 2 == 1;
              var rendered_segments = plugin_.FlightPlanRenderedSegment(
                  main_vessel_guid, sun_world_position, i);
              if (rendered_segments.IteratorAtEnd()) {
                Log.Info("Skipping segment " + i);
                continue;
              }
              Vector3d position_at_start =
                  (Vector3d)rendered_segments.
                      IteratorGetDiscreteTrajectoryXYZ();
              if (use_cayley_plotting_) {
                GLLines.RenderAndDeleteTrajectory(
                    rendered_segments,
                    is_burn ? XKCDColors.OrangeRed : XKCDColors.BabyBlue,
                    is_burn ? GLLines.Style.SOLID : GLLines.Style.DASHED);
              }
              if (use_чебышёв_plotting_) {
                IntPtr rp2_lines_iterator =
                    planetarium.PlanetariumPlotFlightPlanSegment(
                        plugin_,
                        чебышёв_plotting_method_,
                        main_vessel_guid,
                        i);
                GLLines.PlotAndDeleteRP2Lines(
                    rp2_lines_iterator,
                    is_burn ? XKCDColors.Grapefruit : XKCDColors.Blueberry,
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
                    position_at_start + scale * (Vector3d)world_direction,
                    hide_behind_bodies: false);
                    };
                add_vector(manoeuvre.tangent, XKCDColors.NeonYellow);
                add_vector(manoeuvre.normal, XKCDColors.AquaBlue);
                add_vector(manoeuvre.binormal, XKCDColors.PurplePink);
              }
            }
          }
        });
      } finally {
        Interface.PlanetariumDelete(ref planetarium);
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
      IntPtr ascending_nodes_iterator;
      IntPtr descending_nodes_iterator;
      IntPtr approaches_iterator;
      plugin_.RenderedPredictionNodes(vessel_guid,
                                      sun_world_position,
                                      out ascending_nodes_iterator,
                                      out descending_nodes_iterator);
      plugin_.RenderedPredictionClosestApproaches(vessel_guid,
                                                  sun_world_position,
                                                  out approaches_iterator);
      map_node_pool_.RenderAndDeleteMarkers(
          ascending_nodes_iterator,
          MapObject.ObjectType.AscendingNode,
          MapNodePool.NodeSource.PREDICTION,
          vessel    : target,
          celestial : plotting_frame_selector_.get().selected_celestial);
      map_node_pool_.RenderAndDeleteMarkers(
          descending_nodes_iterator,
          MapObject.ObjectType.DescendingNode,
          MapNodePool.NodeSource.PREDICTION,
          vessel    : target,
          celestial : plotting_frame_selector_.get().selected_celestial);
      map_node_pool_.RenderAndDeleteMarkers(
          approaches_iterator,
          MapObject.ObjectType.ApproachIntersect,
          MapNodePool.NodeSource.PREDICTION,
          vessel    : target,
          celestial : plotting_frame_selector_.get().selected_celestial);
    } else {
      foreach (CelestialBody celestial in
               plotting_frame_selector_.get().FixedBodies()) {
        IntPtr apoapsis_iterator;
        IntPtr periapsis_iterator;
        plugin_.RenderedPredictionApsides(vessel_guid,
                                          celestial.flightGlobalsIndex,
                                          sun_world_position,
                                          out apoapsis_iterator,
                                          out periapsis_iterator);
        map_node_pool_.RenderAndDeleteMarkers(
            apoapsis_iterator,
            MapObject.ObjectType.Apoapsis,
            MapNodePool.NodeSource.PREDICTION,
            vessel    : null,
            celestial : celestial);
        map_node_pool_.RenderAndDeleteMarkers(
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
        IntPtr ascending_nodes_iterator;
        IntPtr descending_nodes_iterator;
        plugin_.RenderedPredictionNodes(vessel_guid,
                                        sun_world_position,
                                        out ascending_nodes_iterator,
                                        out descending_nodes_iterator);
        map_node_pool_.RenderAndDeleteMarkers(
            ascending_nodes_iterator,
            MapObject.ObjectType.AscendingNode,
            MapNodePool.NodeSource.PREDICTION,
            vessel    : null,
            celestial : primary);
        map_node_pool_.RenderAndDeleteMarkers(
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
      IntPtr ascending_nodes_iterator;
      IntPtr descending_nodes_iterator;
      IntPtr approaches_iterator;
      plugin_.FlightPlanRenderedNodes(vessel_guid,
                                      sun_world_position,
                                      out ascending_nodes_iterator,
                                      out descending_nodes_iterator);
      plugin_.FlightPlanRenderedClosestApproaches(vessel_guid,
                                                  sun_world_position,
                                                  out approaches_iterator);
      map_node_pool_.RenderAndDeleteMarkers(
          ascending_nodes_iterator,
          MapObject.ObjectType.AscendingNode,
          MapNodePool.NodeSource.FLIGHT_PLAN,
          vessel    : target,
          celestial : plotting_frame_selector_.get().selected_celestial);
      map_node_pool_.RenderAndDeleteMarkers(
          descending_nodes_iterator,
          MapObject.ObjectType.DescendingNode,
          MapNodePool.NodeSource.FLIGHT_PLAN,
          vessel    : target,
          celestial : plotting_frame_selector_.get().selected_celestial);
      map_node_pool_.RenderAndDeleteMarkers(
          approaches_iterator,
          MapObject.ObjectType.ApproachIntersect,
          MapNodePool.NodeSource.FLIGHT_PLAN,
          vessel    : target,
          celestial : plotting_frame_selector_.get().selected_celestial);
    } else {
      foreach (CelestialBody celestial in
               plotting_frame_selector_.get().FixedBodies()) {
        IntPtr apoapsis_iterator;
        IntPtr periapsis_iterator;
        plugin_.FlightPlanRenderedApsides(vessel_guid,
                                          celestial.flightGlobalsIndex,
                                          sun_world_position,
                                          out apoapsis_iterator,
                                          out periapsis_iterator);
        map_node_pool_.RenderAndDeleteMarkers(
            apoapsis_iterator,
            MapObject.ObjectType.Apoapsis,
            MapNodePool.NodeSource.FLIGHT_PLAN,
            vessel    : null,
            celestial : celestial);
        map_node_pool_.RenderAndDeleteMarkers(
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
        IntPtr ascending_nodes_iterator;
        IntPtr descending_nodes_iterator;
        plugin_.FlightPlanRenderedNodes(vessel_guid,
                                        sun_world_position,
                                        out ascending_nodes_iterator,
                                        out descending_nodes_iterator);
        map_node_pool_.RenderAndDeleteMarkers(
            ascending_nodes_iterator,
            MapObject.ObjectType.AscendingNode,
            MapNodePool.NodeSource.PREDICTION,
            vessel    : null,
            celestial : primary);
        map_node_pool_.RenderAndDeleteMarkers(
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
    using (new VerticalLayout()) {
      if (!PluginRunning()) {
        UnityEngine.GUILayout.TextArea(text : "Plugin is not started");
      }
      if (DateTimeOffset.Now > next_release_date) {
        UnityEngine.GUILayout.TextArea(
            "Announcement: the new moon of lunation number " +
            next_release_lunation_number +
            " has come; please download the latest Principia release, " +
            next_release_name + ".");
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
        using (new HorizontalLayout()) {
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
      ToggleableSection(name   : "Reset Principia",
                        show   : ref show_reset_button_,
                        render : ResetButton);
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
    using (new HorizontalLayout()) {
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
      using (new HorizontalLayout()) {
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
    using (new HorizontalLayout()) {
      use_cayley_plotting_ = UnityEngine.GUILayout.Toggle(
          use_cayley_plotting_, "Cayley plotting");
      use_чебышёв_plotting_ = UnityEngine.GUILayout.Toggle(
          use_чебышёв_plotting_, "Чебышёв plotting");
    }
    using (new HorizontalLayout()) {
      UnityEngine.GUILayout.Label("Чебышёв plotting method:");
      for (int i = 0; i < чебышёв_plotting_methods_count; ++i) {
        if (UnityEngine.GUILayout.Toggle(чебышёв_plotting_method_ == i,
                                         i.ToString())) {
          чебышёв_plotting_method_ = i;
        }
      }
    }
    using (new HorizontalLayout()) {
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
    using (new HorizontalLayout()) {
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
    using (new HorizontalLayout()) {
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
      using (new HorizontalLayout()) {
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
    using (new HorizontalLayout()) {
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
    var frame_type =
        (ReferenceFrameSelector.FrameType)frame_parameters.extension;
    if (frame_type != ReferenceFrameSelector.FrameType.BODY_SURFACE) {
      last_non_surface_frame_type_ = frame_type;
    }
    navball_changed_ = true;
    reset_rsas_target_ = true;
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
                  gravity_model.GetValue("reference_instant"),
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
      // We create the plugin at J2000 (a.k.a. Instant{}), rather than
      // |Planetarium.GetUniversalTime()|, in order to get a deterministic
      // initial state.
      plugin_ = Interface.NewPlugin("JD2451545", "JD2451545",
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
                gravity_model?.GetValue("gravitational_parameter")
                     ?? (body.gravParameter + " m^3/s^2"),
            // The origin of rotation in KSP is the x of Barycentric, rather
            // than the y axis as is the case for Earth, so the right
            // ascension is -90 deg.
            reference_instant    = 
                gravity_model?.GetValue("reference_instant") ?? "JD2451545.0",
            mean_radius          =
                gravity_model?.GetValue("mean_radius") ?? (body.Radius + " m"),
            axis_right_ascension =
                gravity_model?.GetValue("axis_right_ascension") ?? "-90 deg",
            axis_declination     =
                gravity_model?.GetValue("axis_declination") ?? "90 deg",
            reference_angle      =
                gravity_model?.GetValue("reference_angle")
                    ?? (body.initialRotation.ToString() + " deg"),
            angular_frequency    =
                gravity_model?.GetValue("angular_frequency")
                    ??  (body.angularV.ToString() + " rad/s"),
            j2                   = gravity_model?.GetValue("j2"),
            reference_radius     = gravity_model?.GetValue("reference_radius")};
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
  }

  private void SetRotatingFrameThresholds() {
    ApplyToBodyTree(
        body => {
          double timewarp_limit = body.timeWarpAltitudeLimits[1];
          if (timewarp_limit == 0) {
            Log.Error("The timewarp limit for " + body.NameWithArticle() +
                      " vanishes");
            if (body.atmosphereDepth == 0) {
              Log.Error(body.NameWithArticle() +
                        " is airless, setting the manageability" +
                        " threshold to 10 km to allow landings");
              timewarp_limit = 10e3;
            }
          }
          body.inverseRotThresholdAltitude =
                (float)Math.Max(timewarp_limit,
                                body.atmosphereDepth);
          Log.Info("Set manageability threshold for " + body.NameWithArticle() +
                   " to " + body.inverseRotThresholdAltitude +
                   " m; its atmosphere extends to " + body.atmosphereDepth +
                   " m and timewarp is allowed above " +
                   body.timeWarpAltitudeLimits[1] + " m");
        });
  }

  private void RemoveBuggyTidalLocking() {
    ApplyToBodyTree(body => body.tidallyLocked = false);
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
