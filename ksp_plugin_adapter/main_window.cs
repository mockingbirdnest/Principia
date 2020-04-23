using System;

namespace principia {
namespace ksp_plugin_adapter {

internal class MainWindow : VesselSupervisedWindowRenderer {
  // Update this section before each release.
  private const string next_release_name_ = "Fuchs";
  private const int next_release_lunation_number_ = 252;
  private readonly DateTimeOffset next_release_date_ =
      new DateTimeOffset(2020, 05, 22, 17, 39, 00, TimeSpan.Zero);

  public MainWindow(PrincipiaPluginAdapter adapter,
                    FlightPlanner flight_planner,
                    OrbitAnalyser orbit_analyser,
                    ReferenceFrameSelector plotting_frame_selector,
                    PredictedVessel predicted_vessel)
      : base(adapter, predicted_vessel) {
    adapter_ = adapter;
    flight_planner_ = flight_planner;
    orbit_analyser_ = orbit_analyser;
    plotting_frame_selector_ = plotting_frame_selector;
    Show();
  }

  public void SelectTargetCelestial(MapObject map_object) {
    if (selecting_target_celestial_) {
      FlightGlobals.fetch.SetVesselTarget(map_object.celestialBody);
      selecting_target_celestial_ = false;
    } else if (PlanetariumCamera.fetch.target != map_object) {
      PlanetariumCamera.fetch.SetTarget(map_object);
    }
  }

  public void SelectActiveVesselTarget(MapObject map_object,
                                       bool set_planetarium_camera) {
    if (selecting_active_vessel_target) {
      FlightGlobals.fetch.SetVesselTarget(map_object.vessel);
      selecting_active_vessel_target = false;
    } else if (set_planetarium_camera &&
               PlanetariumCamera.fetch.target != map_object) {
      PlanetariumCamera.fetch.SetTarget(map_object);
    }
  }

  public bool selecting_active_vessel_target { get; private set; } = false;

  public bool display_patched_conics { get; private set; } = false;

  public double history_length => history_length_.value;

  public override void Load(ConfigNode node) {
    base.Load(node);

    string show_ksp_features_value =
        node.GetAtMostOneValue("show_ksp_features");
    if (show_ksp_features_value != null) {
      show_ksp_features_ = Convert.ToBoolean(show_ksp_features_value);
    }
    string show_logging_settings_value =
        node.GetAtMostOneValue("show_logging_settings");
    if (show_logging_settings_value != null) {
      show_logging_settings_ = Convert.ToBoolean(show_logging_settings_value);
    }
    string show_prediction_settings_value =
        node.GetAtMostOneValue("show_prediction_settings");
    if (show_prediction_settings_value != null) {
      show_prediction_settings_ =
          Convert.ToBoolean(show_prediction_settings_value);
    }

    string history_length_value =
        node.GetAtMostOneValue("history_length");
    if (history_length_value != null) {
      history_length_.value = Convert.ToDouble(history_length_value);
    }

    string buffered_logging_value =
        node.GetAtMostOneValue("buffered_logging");
    if (buffered_logging_value != null) {
      buffered_logging_ = Convert.ToInt32(buffered_logging_value);
    }
    string stderr_logging_value = node.GetAtMostOneValue("stderr_logging");
    if (stderr_logging_value != null) {
      stderr_logging_ = Convert.ToInt32(stderr_logging_value);
    }
    string suppressed_logging_value =
        node.GetAtMostOneValue("suppressed_logging");
    if (suppressed_logging_value != null) {
      suppressed_logging_ = Convert.ToInt32(suppressed_logging_value);
    }
    string verbose_logging_value = node.GetAtMostOneValue("verbose_logging");
    if (verbose_logging_value != null) {
      verbose_logging_ = Convert.ToInt32(verbose_logging_value);
    }

    string must_record_journal_value =
        node.GetAtMostOneValue("must_record_journal");
    if (must_record_journal_value != null) {
      must_record_journal_ = Convert.ToBoolean(must_record_journal_value);
    }

    Log.SetBufferedLogging(buffered_logging_);
    Log.SetSuppressedLogging(suppressed_logging_);
    Log.SetStderrLogging(stderr_logging_);
    Log.SetVerboseLogging(verbose_logging_);

    if (must_record_journal_) {
      journaling_ = true;
      Log.ActivateRecorder(true);
    }
  }

  public override void Save(ConfigNode node) {
    base.Save(node);

    node.SetValue("show_ksp_features",
                  show_ksp_features_,
                  createIfNotFound : true);
    node.SetValue("show_logging_settings",
                  show_logging_settings_,
                  createIfNotFound : true);
    node.SetValue("show_prediction_settings",
                  show_prediction_settings_,
                  createIfNotFound : true);

    node.SetValue("history_length",
                  history_length,
                  createIfNotFound : true);
    node.SetValue("buffered_logging",
                  buffered_logging_,
                  createIfNotFound : true);
    node.SetValue("stderr_logging",
                  stderr_logging_,
                  createIfNotFound : true);
    node.SetValue("suppressed_logging",
                  suppressed_logging_,
                  createIfNotFound : true);
    node.SetValue("verbose_logging",
                  verbose_logging_,
                  createIfNotFound : true);

    node.SetValue("must_record_journal",
                  must_record_journal_,
                  createIfNotFound : true);
  }

  protected override string Title => "Principia";

  protected override void RenderWindow(int window_id) {
    using (new UnityEngine.GUILayout.VerticalScope()) {
      if (!adapter_.PluginRunning()) {
        UnityEngine.GUILayout.Label(
            text : "Plugin is not started",
            style : Style.Warning(UnityEngine.GUI.skin.label));
      }
      if (DateTimeOffset.Now > next_release_date_) {
        if (Versioning.version_minor <= 4) {
          UnityEngine.GUILayout.TextArea(
              "Announcement: the new moon of lunation number " +
              next_release_lunation_number_ +
              " has come; please update KSP to version 1.6.1 and download " +
              "the latest Principia release, " + next_release_name_ + ". " +
              "Note that RealismOverhaul and RealSolarSystem now support " +
              "KSP 1.6.1.",
              style : Style.Multiline(UnityEngine.GUI.skin.textArea));
        } else {
          UnityEngine.GUILayout.TextArea(
              "Announcement: the new moon of lunation number " +
              next_release_lunation_number_ +
              " has come; please download the latest Principia release, " +
              next_release_name_ + ".",
              style: Style.Multiline(UnityEngine.GUI.skin.textArea));
        }
      }
      Interface.GetVersion(build_date : out string _,
                           version    : out string version);
      UnityEngine.GUILayout.Label(
          version,
          style : Style.Info(UnityEngine.GUI.skin.label));
      history_length_.Render(enabled : true);
      if (MapView.MapIsEnabled &&
          FlightGlobals.ActiveVessel?.orbitTargeter != null) {
        show_selection_ui_ = true;
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          selecting_active_vessel_target = UnityEngine.GUILayout.Toggle(
              selecting_active_vessel_target, "Select target vessel...");
          if (selecting_active_vessel_target) {
            selecting_target_celestial_ = false;
          }
          if (FlightGlobals.fetch.VesselTarget?.GetVessel()) {
            UnityEngine.GUILayout.Label(
                "Target: " +
                    FlightGlobals.fetch.VesselTarget.GetVessel().vesselName,
                UnityEngine.GUILayout.ExpandWidth(true));
            if (UnityEngine.GUILayout.Button("Clear", GUILayoutWidth(2))) {
              selecting_active_vessel_target = false;
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
        // This will remove the "Select" UI so it must shrink.
        if (show_selection_ui_) {
          show_selection_ui_ = false;
          Shrink();
        }
        selecting_active_vessel_target = false;
      }

      // The plugin is destroyed, e.g., when using the "Switch To" button.  Wait
      // until it's back alive to display information that requires to cross the
      // interface.
      if (adapter_.PluginRunning()) {
        plotting_frame_selector_.RenderButton();
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          flight_planner_.RenderButton();
          orbit_analyser_.RenderButton();
        }
        RenderToggleableSection(name   : "Prediction Settings",
                                show   : ref show_prediction_settings_,
                                render : RenderPredictionSettings);
      }
      RenderToggleableSection(name   : "KSP Features",
                              show   : ref show_ksp_features_,
                              render : RenderKSPFeatures);
      RenderToggleableSection(name   : "Logging Settings",
                              show   : ref show_logging_settings_,
                              render : RenderLoggingSettings);
    }
    UnityEngine.GUI.DragWindow();
  }

  private void RenderKSPFeatures() {
    if (show_2519_debugging_ui) {
      correct_orientation = UnityEngine.GUILayout.Toggle(
          correct_orientation,
          "Correct orientation");
      correct_angular_velocity = UnityEngine.GUILayout.Toggle(
          correct_angular_velocity,
          "Correct angular velocity");
      thresholding = UnityEngine.GUILayout.Toggle(
          thresholding,
          "Only correct orientation slower than ω");
      Interface.SetAngularMomentumConservation(
          correct_orientation, correct_angular_velocity, thresholding);
      string trace = null;
      if (FlightGlobals.ActiveVessel != null &&
          plugin.HasVessel(FlightGlobals.ActiveVessel.id.ToString())) {
        trace = plugin.VesselGetPileUpTrace(
            FlightGlobals.ActiveVessel.id.ToString());
      }
      UnityEngine.GUILayout.TextArea(
          trace ?? "No managed active vessel",
          style : Style.Multiline(UnityEngine.GUI.skin.textArea));
    }
    display_patched_conics = UnityEngine.GUILayout.Toggle(
        value : display_patched_conics,
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
          selecting_active_vessel_target = false;
        }
        CelestialBody target_celestial =
            FlightGlobals.fetch.VesselTarget as CelestialBody;
        if (target_celestial != null) {
          UnityEngine.GUILayout.Label("Target: " + target_celestial.name,
                                      UnityEngine.GUILayout.ExpandWidth(true));
          if (UnityEngine.GUILayout.Button("Clear", GUILayoutWidth(2))) {
            selecting_target_celestial_ = false;
            FlightGlobals.fetch.SetVesselTarget(null);
          }
        }
      }
    } else {
      selecting_target_celestial_ = false;
    }
  }

  private void RenderLoggingSettings() {
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label(text : "Verbose level:");
      if (UnityEngine.GUILayout.Button(text    : "←",
                                       options : GUILayoutWidth(2))) {
        Log.SetVerboseLogging(Math.Max(verbose_logging_ - 1, 0));
        verbose_logging_ = Log.GetVerboseLogging();
      }
      UnityEngine.GUILayout.TextArea(
          text    : Log.GetVerboseLogging().ToString(),
          options : GUILayoutWidth(2));
      if (UnityEngine.GUILayout.Button(text    : "→",
                                       options : GUILayoutWidth(2))) {
        Log.SetVerboseLogging(Math.Min(verbose_logging_ + 1, 4));
        verbose_logging_ = Log.GetVerboseLogging();
      }
    }
    float column_width = Width(3);
    var gui_layout_column_width = GUILayoutWidth(3);
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Space(column_width);
      UnityEngine.GUILayout.Label(text    : "Log",
                                  options : gui_layout_column_width);
      UnityEngine.GUILayout.Label(text    : "stderr",
                                  options : gui_layout_column_width);
      UnityEngine.GUILayout.Label(text    : "Flush",
                                  options : gui_layout_column_width);
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Space(column_width);
      if (UnityEngine.GUILayout.Button(text    : "↑",
                                       options : gui_layout_column_width)) {
        Log.SetSuppressedLogging(Math.Max(suppressed_logging_ - 1, 0));
        suppressed_logging_ = Log.GetSuppressedLogging();
      }
      if (UnityEngine.GUILayout.Button(text    : "↑",
                                       options : gui_layout_column_width)) {
        Log.SetStderrLogging(Math.Max(stderr_logging_ - 1, 0));
        stderr_logging_ = Log.GetStderrLogging();
      }
      if (UnityEngine.GUILayout.Button(text    : "↑",
                                       options : gui_layout_column_width)) {
        Log.SetBufferedLogging(Math.Max(buffered_logging_ - 1, -1));
        buffered_logging_ = Log.GetBufferedLogging();
      }
    }
    for (int severity = 0; severity <= 3; ++severity) {
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        UnityEngine.GUILayout.Label(
            text    : Log.severity_names[severity],
            options : gui_layout_column_width);
        UnityEngine.GUILayout.Toggle(
            value   : severity >= Log.GetSuppressedLogging(),
            text    : "",
            options : gui_layout_column_width);
        UnityEngine.GUILayout.Toggle(
            value   : severity >= Log.GetStderrLogging(),
            text    : "",
            options : gui_layout_column_width);
        UnityEngine.GUILayout.Toggle(
            value   : severity > Log.GetBufferedLogging(),
            text    : "",
            options : gui_layout_column_width);
      }
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Space(column_width);
      if (UnityEngine.GUILayout.Button(text    : "↓",
                                       options : gui_layout_column_width)) {
        Log.SetSuppressedLogging(Math.Min(suppressed_logging_ + 1, 3));
        suppressed_logging_ = Log.GetSuppressedLogging();
      }
      if (UnityEngine.GUILayout.Button(text    : "↓",
                                       options : gui_layout_column_width)) {
        Log.SetStderrLogging(Math.Min(stderr_logging_ + 1, 3));
        stderr_logging_ = Log.GetStderrLogging();
      }
      if (UnityEngine.GUILayout.Button(text    : "↓",
                                       options : gui_layout_column_width)) {
        Log.SetBufferedLogging(Math.Min(buffered_logging_ + 1, 3));
        buffered_logging_ = Log.GetBufferedLogging();
      }
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      must_record_journal_ = UnityEngine.GUILayout.Toggle(
          value   : must_record_journal_,
          text    : "Record journal (starts on load)");
      UnityEngine.GUILayout.Label(
          "Journaling is " + (journaling_ ? "ON" : "OFF"),
          style : Style.Info(Style.RightAligned(UnityEngine.GUI.skin.label)));
    }
    if (journaling_ && !must_record_journal_) {
      // We can deactivate a recorder at any time, but in order for replaying to
      // work, we should only activate one before creating a plugin.
      journaling_ = false;
      Interface.ActivateRecorder(false);
    }
  }

  private void RenderPredictionSettings() {
    AdaptiveStepParameters? adaptive_step_parameters = null;
    string vessel_guid = predicted_vessel?.id.ToString();
    if (vessel_guid != null) {
      adaptive_step_parameters =
          plugin.VesselGetPredictionAdaptiveStepParameters(vessel_guid);
      prediction_length_tolerance_index_ = Array.FindIndex(
          prediction_length_tolerances_,
          (double tolerance) =>
              tolerance >=
                  adaptive_step_parameters.Value.length_integration_tolerance);
      if (prediction_length_tolerance_index_ < 0) {
        prediction_length_tolerance_index_ =
            default_prediction_length_tolerance_index_;
      }
      prediction_steps_index_ = Array.FindIndex(
          prediction_steps_,
          (long step) => step >= adaptive_step_parameters.Value.max_steps);
      if (prediction_steps_index_ < 0) {
        prediction_steps_index_ = default_prediction_steps_index_;
      }
    }

    // TODO(egg): make the speed tolerance independent.
    if (RenderSelector(prediction_length_tolerances_,
                       ref prediction_length_tolerance_index_,
                       "Tolerance",
                       "{0:0.0e0} m",
                       enabled: adaptive_step_parameters.HasValue)) {
      AdaptiveStepParameters new_adaptive_step_parameters =
          new AdaptiveStepParameters{
            integrator_kind = adaptive_step_parameters.Value.integrator_kind,
            max_steps = prediction_steps,
            length_integration_tolerance = prediction_length_tolerance,
            speed_integration_tolerance = prediction_length_tolerance};
      plugin.VesselSetPredictionAdaptiveStepParameters(
          vessel_guid, new_adaptive_step_parameters);
    }
    if (RenderSelector(prediction_steps_,
                       ref prediction_steps_index_,
                       "Steps",
                       "{0:0.00e0}",
                       enabled: adaptive_step_parameters.HasValue)) {
      AdaptiveStepParameters new_adaptive_step_parameters =
          new AdaptiveStepParameters{
            integrator_kind = adaptive_step_parameters.Value.integrator_kind,
            max_steps = prediction_steps,
            length_integration_tolerance = prediction_length_tolerance,
            speed_integration_tolerance = prediction_length_tolerance};
      plugin.VesselSetPredictionAdaptiveStepParameters(
          vessel_guid, new_adaptive_step_parameters);
    }
  }

  private bool RenderSelector<T>(T[] array,
                                 ref int index,
                                 string label,
                                 string format,
                                 bool enabled) {
    bool changed = false;
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label(text    : label + ":",
                                  options : GUILayoutWidth(6));
      if (UnityEngine.GUILayout.Button(text    : index == 0 ? "min" : "-",
                                       options : GUILayoutWidth(2)) &&
          enabled &&
          index != 0) {
        --index;
        changed = true;
      }
      UnityEngine.GUILayout.TextArea(
          text    : enabled
                        ? string.Format(Culture.culture, format, array[index])
                        : "",
          style   : Style.RightAligned(UnityEngine.GUI.skin.textArea),
          options : GUILayoutWidth(3));
      if (UnityEngine.GUILayout.Button(
              text    : index == array.Length - 1 ? "max" : "+",
              options : GUILayoutWidth(2)) &&
          enabled &&
          index != array.Length - 1) {
        ++index;
        changed = true;
      }
    }
    return changed;
  }

  private void RenderToggleableSection(string name,
                                       ref bool show,
                                       Action render) {
    string toggle = show ? "↑ " + name + " ↑"
                         : "↓ " + name + " ↓";
    if (UnityEngine.GUILayout.Button(toggle)) {
      show = !show;
      if (!show) {
        Shrink();
      }
    }
    if (show) {
      render();
    }
  }

  private IntPtr plugin => adapter_.Plugin();
  private double prediction_length_tolerance =>
      prediction_length_tolerances_[prediction_length_tolerance_index_];
  private long prediction_steps => prediction_steps_[prediction_steps_index_];

  private DifferentialSlider history_length_ = new DifferentialSlider(
      label            : "Max history length",
      unit             : null,
      log10_lower_rate : 0,
      log10_upper_rate : 7,
      min_value        : 10,
      max_value        : double.PositiveInfinity,
      formatter        : Formatters.FormatMissionDuration,
      parser           : Formatters.TryParseMissionDuration,
      label_width      : 5,
      field_width      : 5) {
      value = 7 * 24 * 60 * 60
  };

  // These flags exist to facilitate investigation of #2519.
  // They must not be serialized: their non-default values can lead to absurd
  // behaviour.
  private static bool correct_orientation = true;
  private static bool correct_angular_velocity = true;
  private static bool thresholding = true;
  private static readonly bool show_2519_debugging_ui = false;

  private static readonly double[] prediction_length_tolerances_ =
      {1e-3, 1e-2, 1e0, 1e1, 1e2, 1e3, 1e4};
  private static readonly long[] prediction_steps_ =
      {1 << 2, 1 << 4, 1 << 6, 1 << 8, 1 << 10, 1 << 12, 1 << 14, 1 << 16,
       1 << 18, 1 << 20, 1 << 22, 1 << 24};

  private const int default_prediction_length_tolerance_index_ = 1;
  private const int default_prediction_steps_index_ = 4;

  private readonly PrincipiaPluginAdapter adapter_;
  private readonly FlightPlanner flight_planner_;
  private readonly OrbitAnalyser orbit_analyser_;
  private readonly ReferenceFrameSelector plotting_frame_selector_;

  private bool selecting_target_celestial_ = false;

  private bool show_ksp_features_ = false;
  private bool show_logging_settings_ = false;
  private bool show_prediction_settings_ = true;

  private bool show_selection_ui_ = false;

  private int prediction_length_tolerance_index_ =
      default_prediction_length_tolerance_index_;
  private int prediction_steps_index_ = default_prediction_steps_index_;

  private int buffered_logging_ = 0;
  private int stderr_logging_ = 2;
  private int suppressed_logging_ = 0;
  private int verbose_logging_ = 0;

  // Whether a journal will be recorded when the plugin is next constructed.
  private bool must_record_journal_ = false;
  // Whether a journal is currently being recorded.
  private static bool journaling_ = false;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
