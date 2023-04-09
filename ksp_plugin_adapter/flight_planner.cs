using System;
using System.Collections.Generic;
using System.Linq;
using KSP.Localization;

namespace principia {
namespace ksp_plugin_adapter {

class FlightPlanner : VesselSupervisedWindowRenderer {
  public FlightPlanner(PrincipiaPluginAdapter adapter,
                       PredictedVessel predicted_vessel) : base(
      adapter,
      predicted_vessel) {
    adapter_ = adapter;
    predicted_vessel_ = predicted_vessel;
    final_time_ = new DifferentialSlider(
        label            : L10N.CacheFormat("#Principia_FlightPlan_PlanLength"),
        unit             : null,
        log10_lower_rate : log10_time_lower_rate,
        log10_upper_rate : log10_time_upper_rate,
        min_value        : 10,
        max_value        : double.PositiveInfinity,
        formatter        : FormatPlanLength,
        parser           : TryParsePlanLength,
        field_width      : 7);
    final_trajectory_analyser_ =
        new PlannedOrbitAnalyser(adapter, predicted_vessel);
  }

  public void RenderButton() {
    RenderButton(L10N.CacheFormat("#Principia_FlightPlan_ToggleButton"),
                 GUILayoutWidth(4));
  }

  public bool show_guidance => show_guidance_;

  public override void Load(ConfigNode node) {
    base.Load(node);

    string show_guidance_value = node.GetAtMostOneValue("show_guidance");
    if (show_guidance_value != null) {
      show_guidance_ = Convert.ToBoolean(show_guidance_value);
    }
  }

  public override void Save(ConfigNode node) {
    base.Save(node);

    node.SetValue("show_guidance", show_guidance_, createIfNotFound : true);
  }

  protected override string Title =>
      L10N.CacheFormat("#Principia_FlightPlan_Title");

  protected override void RenderWindowContents(int window_id) {
    // We must ensure that the GUI elements don't change between Layout and
    // Repaint.  This means that any state change must occur before Layout or
    // after Repaint.  This if statement implements the former.  It updates the
    // vessel and the editors to reflect the current state of the plugin and
    // then proceeds with the UI code.
    if (UnityEngine.Event.current.type == UnityEngine.EventType.Layout) {
      UpdateVesselAndBurnEditors();
    }

    // The UI code proper, executed identically for Layout and Repaint.  We
    // can freely change the state in events like clicks (e.g., in if statements
    // for buttons) as these don't happen between Layout and Repaint.
    string vessel_guid = predicted_vessel?.id.ToString();
    if (vessel_guid == null) {
      return;
    }

    using (new UnityEngine.GUILayout.HorizontalScope()) {
      int flight_plans = plugin.FlightPlanCount(vessel_guid);
      int selected_flight_plan = plugin.FlightPlanSelected(vessel_guid);
      for (int i = 0; i < flight_plans; ++i) {
        var id = new string(L10N.CacheFormat("#Principia_AlphabeticList")[i],
                            1);
        if (UnityEngine.GUILayout.Toggle(i == selected_flight_plan, id,
                                         "Button",
                                         GUILayoutWidth(1)) &&
            i != selected_flight_plan) {
          plugin.FlightPlanSelect(vessel_guid, i);
          final_time_.value_if_different =
              plugin.FlightPlanGetDesiredFinalTime(vessel_guid);
          ClearBurnEditors();
          UpdateVesselAndBurnEditors();
        }
      }
      bool must_create_flight_plan = false;
      if (flight_plans == 0) {
        must_create_flight_plan = UnityEngine.GUILayout.Button(
            L10N.CacheFormat("#Principia_FlightPlan_Create"));
      } else if (flight_plans < max_flight_plans) {
        must_create_flight_plan =
            UnityEngine.GUILayout.Button("+", GUILayoutWidth(1));
      }
      if (must_create_flight_plan) {
        plugin.FlightPlanCreate(vessel_guid,
                                plugin.CurrentTime() + 3600,
                                predicted_vessel.GetTotalMass());
        final_time_.value_if_different =
            plugin.FlightPlanGetDesiredFinalTime(vessel_guid);
        ClearBurnEditors();
        UpdateVesselAndBurnEditors();
        return;
      }
    }

    if (plugin.FlightPlanExists(vessel_guid)) {
      RenderFlightPlan(vessel_guid);
    }
    UnityEngine.GUI.DragWindow();
  }

  private void ClearBurnEditors() {
    if (burn_editors_ != null) {
      foreach (BurnEditor editor in burn_editors_) {
        editor.Close();
      }
      burn_editors_ = null;
      Shrink();
    }
  }

  private void UpdateVesselAndBurnEditors() {
    {
      string vessel_guid = predicted_vessel?.id.ToString();
      if (vessel_guid == null ||
          previous_predicted_vessel_ != predicted_vessel ||
          !plugin.FlightPlanExists(vessel_guid) ||
          plugin.FlightPlanNumberOfManoeuvres(vessel_guid) !=
          burn_editors_?.Count) {
        ClearBurnEditors();
        previous_predicted_vessel_ = predicted_vessel;
      }
    }

    if (burn_editors_ == null) {
      string vessel_guid = predicted_vessel?.id.ToString();
      if (vessel_guid != null &&
          plugin.FlightPlanExists(vessel_guid)) {
        burn_editors_ = new List<BurnEditor>();
        final_time_.value_if_different =
            plugin.FlightPlanGetDesiredFinalTime(vessel_guid);
        for (int i = 0;
             i < plugin.FlightPlanNumberOfManoeuvres(vessel_guid);
             ++i) {
          // Dummy initial time, we call |Reset| immediately afterwards.
          burn_editors_.Add(new BurnEditor(adapter_,
                                           predicted_vessel,
                                           initial_time      : 0,
                                           index             : i,
                                           get_burn_at_index : burn_editors_.
                                               ElementAtOrDefault));
          burn_editors_.Last().Reset(
              plugin.FlightPlanGetManoeuvre(vessel_guid, i));
        }
      }
    }

    if (burn_editors_ != null) {
      string vessel_guid = predicted_vessel?.id.ToString();
      double current_time = plugin.CurrentTime();
      first_future_manœuvre_ = null;
      for (int i = 0; i < burn_editors_.Count; ++i) {
        NavigationManoeuvre manœuvre =
            plugin.FlightPlanGetManoeuvre(vessel_guid, i);
        if (current_time < manœuvre.final_time) {
          first_future_manœuvre_ = i;
          break;
        }
      }
      // Must be computed during layout as it affects the layout of some of the
      // differential sliders.
      number_of_anomalous_manœuvres_ =
          plugin.FlightPlanNumberOfAnomalousManoeuvres(vessel_guid);
    }
  }

  private void RenderFlightPlan(string vessel_guid) {
    using (new UnityEngine.GUILayout.VerticalScope()) {
      if (final_time_.Render(enabled : true)) {
        var status =
            plugin.FlightPlanSetDesiredFinalTime(
                vessel_guid,
                final_time_.value);
        UpdateStatus(status, null);
      }
      // Always refresh the final time from C++ as it may have changed because
      // the last burn changed.
      final_time_.value_if_different =
          plugin.FlightPlanGetDesiredFinalTime(vessel_guid);

      FlightPlanAdaptiveStepParameters parameters =
          plugin.FlightPlanGetAdaptiveStepParameters(vessel_guid);
      length_integration_tolerance_index_ = Math.Max(
          0,
          Array.FindIndex(integration_tolerances_,
                          (double tolerance) => tolerance >=
                                                parameters.
                                                    length_integration_tolerance));
      speed_integration_tolerance_index_ = Math.Max(
          0,
          Array.FindIndex(integration_tolerances_,
                          (double tolerance) => tolerance >=
                                                parameters.
                                                    speed_integration_tolerance));
      max_steps_index_ = Math.Max(0,
                                  Array.FindIndex(
                                      max_steps_,
                                      (long step) =>
                                          step >= parameters.max_steps));

      using (new UnityEngine.GUILayout.HorizontalScope()) {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Label(
              L10N.CacheFormat("#Principia_FlightPlan_MaxSteps"),
              GUILayoutWidth(6));
          if (max_steps_index_ == 0) {
            UnityEngine.GUILayout.Button(
                L10N.CacheFormat("#Principia_DiscreteSelector_Min"));
          } else if (UnityEngine.GUILayout.Button("−")) {
            --max_steps_index_;
            UpdateAdaptiveStepParameters(ref parameters);
            var status =
                plugin.FlightPlanSetAdaptiveStepParameters(
                    vessel_guid,
                    parameters);
            UpdateStatus(status, null);
          }
          UnityEngine.GUILayout.TextArea(
              max_steps_[max_steps_index_].ToString(),
              GUILayoutWidth(3));
          if (max_steps_index_ == max_steps_.Length - 1) {
            UnityEngine.GUILayout.Button(
                L10N.CacheFormat("#Principia_DiscreteSelector_Max"));
          } else if (UnityEngine.GUILayout.Button("+")) {
            ++max_steps_index_;
            UpdateAdaptiveStepParameters(ref parameters);
            var status =
                plugin.FlightPlanSetAdaptiveStepParameters(
                    vessel_guid,
                    parameters);
            UpdateStatus(status, null);
          }
        }
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Label(
              L10N.CacheFormat("#Principia_PredictionSettings_ToleranceLabel"),
              GUILayoutWidth(3));
          // Prior to Ἵππαρχος the tolerances were powers of 2, see #3395.
          if (length_integration_tolerance_index_ == 0 ||
              speed_integration_tolerance_index_ == 0) {
            UnityEngine.GUILayout.Button(
                L10N.CacheFormat("#Principia_DiscreteSelector_Min"));
          } else if (UnityEngine.GUILayout.Button("−")) {
            --length_integration_tolerance_index_;
            --speed_integration_tolerance_index_;
            UpdateAdaptiveStepParameters(ref parameters);
            var status =
                plugin.FlightPlanSetAdaptiveStepParameters(
                    vessel_guid,
                    parameters);
            UpdateStatus(status, null);
          }
          UnityEngine.GUILayout.TextArea(
              length_integration_tolerances_names_[
                  length_integration_tolerance_index_],
              GUILayoutWidth(3));
          if (length_integration_tolerance_index_ ==
              integration_tolerances_.Length - 1 ||
              speed_integration_tolerance_index_ ==
              integration_tolerances_.Length - 1) {
            UnityEngine.GUILayout.Button(
                L10N.CacheFormat("#Principia_DiscreteSelector_Max"));
          } else if (UnityEngine.GUILayout.Button("+")) {
            ++length_integration_tolerance_index_;
            ++speed_integration_tolerance_index_;
            UpdateAdaptiveStepParameters(ref parameters);
            var status =
                plugin.FlightPlanSetAdaptiveStepParameters(
                    vessel_guid,
                    parameters);
            UpdateStatus(status, null);
          }
        }
      }

      double Δv = (from burn_editor in burn_editors_
                   select burn_editor.Δv()).Sum();
      UnityEngine.GUILayout.Label(L10N.CacheFormat(
                                      "#Principia_FlightPlan_TotalΔv",
                                      Δv.ToString("0.000")));

      {
        var style = Style.Warning(Style.Multiline(UnityEngine.GUI.skin.label));
        string message = GetStatusMessage();
        // Size the label explicitly so that it doesn't decrease when the
        // message goes away: that causes annoying flicker.  The enclosing
        // window has a width of 20 units, but not all of that is available,
        // hence 19.
        warning_height_ = Math.Max(warning_height_,
                                   style.CalcHeight(
                                       new UnityEngine.GUIContent(message),
                                       Width(19)));
        UnityEngine.GUILayout.Label(message,
                                    style,
                                    UnityEngine.GUILayout.Height(
                                        warning_height_));
      }

      if (burn_editors_.Count == 0 &&
          UnityEngine.GUILayout.Button(
              L10N.CacheFormat("#Principia_FlightPlan_Delete"))) {
        final_trajectory_analyser_.DisposeWindow();
        final_trajectory_analyser_ =
            new PlannedOrbitAnalyser(adapter_, predicted_vessel_);
        plugin.FlightPlanDelete(vessel_guid);
        ResetStatus();
        Shrink();
        // The state change will happen the next time we go through OnGUI.
      } else {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          if (UnityEngine.GUILayout.Button(
              L10N.CacheFormat("#Principia_FlightPlan_Rebase"))) {
            var status = plugin.FlightPlanRebase(
                vessel_guid,
                predicted_vessel.GetTotalMass());
            UpdateStatus(status, null);
            if (status.ok()) {
              // The final time does not change, but since it is displayed with
              // respect to the beginning of the flight plan, the text must be
              // recomputed.
              final_time_.ResetValue(
                  plugin.FlightPlanGetDesiredFinalTime(vessel_guid));
              return;
            }
          }
          if (plugin.FlightPlanCount(vessel_guid) < max_flight_plans &&
              UnityEngine.GUILayout.Button(
              L10N.CacheFormat("#Principia_FlightPlan_Duplicate"))) {
            plugin.FlightPlanDuplicate(vessel_guid);
          }
        }

        if (burn_editors_.Count > 0) {
          RenderUpcomingEvents();
        }

        // Compute the final times for each manœuvre before displaying them.
        var final_times = new List<double>();
        for (int i = 0; i < burn_editors_.Count - 1; ++i) {
          final_times.Add(plugin.FlightPlanGetManoeuvre(vessel_guid, i + 1).
                              burn.initial_time);
        }
        // Allow extending the flight plan.
        final_times.Add(double.PositiveInfinity);

        for (int i = 0; i < burn_editors_.Count; ++i) {
          Style.HorizontalLine();
          if (RenderCoast(i, out double? orbital_period)) {
            return;
          }
          Style.HorizontalLine();
          BurnEditor burn = burn_editors_[i];
          switch (burn.Render(
              header          :
              L10N.CacheFormat("#Principia_FlightPlan_ManœuvreHeader", i + 1),
              anomalous       : i >=
                                burn_editors_.Count -
                                number_of_anomalous_manœuvres_,
              burn_final_time : final_times[i],
              orbital_period  : orbital_period)) {
            case BurnEditor.Event.Deleted: {
              var status = plugin.FlightPlanRemove(vessel_guid, i);
              UpdateStatus(status, null);
              burn_editors_[i].Close();
              burn_editors_.RemoveAt(i);
              UpdateBurnEditorIndices();
              Shrink();
              return;
            }
            case BurnEditor.Event.Minimized:
            case BurnEditor.Event.Maximized: {
              Shrink();
              return;
            }
            case BurnEditor.Event.Changed: {
              var status =
                  plugin.FlightPlanReplace(vessel_guid, burn.Burn(), i);
              UpdateStatus(status, i);
              burn.Reset(plugin.FlightPlanGetManoeuvre(vessel_guid, i));
              break;
            }
            case BurnEditor.Event.None: {
              break;
            }
          }
        }
        Style.HorizontalLine();
        if (RenderCoast(burn_editors_.Count, orbital_period: out _)) {
          return;
        }
      }
    }
  }

  private void RenderUpcomingEvents() {
    string vessel_guid = predicted_vessel.id.ToString();
    double current_time = plugin.CurrentTime();

    Style.HorizontalLine();
    if (first_future_manœuvre_.HasValue) {
      int first_future_manœuvre = first_future_manœuvre_.Value;
      NavigationManoeuvre manœuvre =
          plugin.FlightPlanGetManoeuvre(vessel_guid, first_future_manœuvre);
      if (manœuvre.burn.initial_time > current_time) {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Label(
              L10N.CacheFormat("#Principia_FlightPlan_UpcomingManœuvre",
                               first_future_manœuvre + 1));
          UnityEngine.GUILayout.Label(
              L10N.CacheFormat("#Principia_FlightPlan_IgnitionCountdown",
                               FormatTimeSpan(
                                   current_time - manœuvre.burn.initial_time)),
              style : Style.RightAligned(UnityEngine.GUI.skin.label));
        }
      } else {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Label(
              L10N.CacheFormat("#Principia_FlightPlan_OngoingManœuvre",
                               first_future_manœuvre + 1));
          UnityEngine.GUILayout.Label(
              L10N.CacheFormat("#Principia_FlightPlan_CutoffCountdown",
                               FormatTimeSpan(
                                   current_time - manœuvre.final_time)),
              style : Style.RightAligned(UnityEngine.GUI.skin.label));
        }
      }
      // In career mode, the patched conic solver may be null.  In that case
      // we do not offer the option of showing the manœuvre on the navball,
      // even though the flight planner is still available to plan it.
      // TODO(egg): We may want to consider setting the burn vector directly
      // rather than going through the solver.
      if (predicted_vessel.patchedConicSolver != null) {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          show_guidance_ = UnityEngine.GUILayout.Toggle(
              show_guidance_,
              L10N.CacheFormat("#Principia_FlightPlan_ShowManœuvreOnNavball"));
          if (UnityEngine.GUILayout.Button(
              L10N.CacheFormat("#Principia_FlightPlan_WarpToManœuvre"))) {
            TimeWarp.fetch.WarpTo(manœuvre.burn.initial_time - 60);
          }
        }
      }
    } else {
      // Reserve some space to avoid the UI changing shape if we have
      // nothing to say.
      UnityEngine.GUILayout.Label(
          L10N.CacheFormat(
              "#Principia_FlightPlan_Warning_AllManœuvresInThePast"),
          Style.Warning(UnityEngine.GUI.skin.label));
      UnityEngine.GUILayout.Space(Width(1));
    }
  }

  private bool RenderCoast(int index, out double? orbital_period) {
    string vessel_guid = predicted_vessel.id.ToString();
    var coast_analysis = plugin.FlightPlanGetCoastAnalysis(
        vessel_guid,
        revolutions_per_cycle   : null,
        days_per_cycle          : null,
        ground_track_revolution : 0,
        index);
    string orbit_description = null;
    orbital_period = coast_analysis.elements?.nodal_period;
    if (coast_analysis.primary_index.HasValue) {
      var primary = FlightGlobals.Bodies[coast_analysis.primary_index.Value];
      int? nodal_revolutions = (int?)(coast_analysis.mission_duration /
                                      coast_analysis.elements?.nodal_period);
      orbit_description = OrbitAnalyser.OrbitDescription(
          primary,
          coast_analysis.mission_duration,
          coast_analysis.elements,
          coast_analysis.recurrence,
          coast_analysis.ground_track_equatorial_crossings,
          coast_analysis.solar_times_of_nodes,
          nodal_revolutions);
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      if (index == burn_editors_.Count) {
        final_trajectory_analyser_.index = index;
        final_trajectory_analyser_.RenderButton();
      } else {
        double start_of_coast = index == 0
                                    ? plugin.FlightPlanGetInitialTime(
                                        vessel_guid)
                                    : burn_editors_[index - 1].final_time;
        string coast_duration =
            (burn_editors_[index].initial_time - start_of_coast).FormatDuration(
                show_seconds: false);
        string coast_description = orbit_description == null
                                       ? L10N.CacheFormat(
                                           "#Principia_FlightPlan_Coast",
                                           coast_duration)
                                       : L10N.CacheFormat(
                                           "#Principia_FlightPlan_CoastInOrbit",
                                           orbit_description,
                                           coast_duration);
        UnityEngine.GUILayout.Label(coast_description);
      }
      if (UnityEngine.GUILayout.Button(
              L10N.CacheFormat("#Principia_FlightPlan_AddManœuvre"),
              GUILayoutWidth(4))) {
        double initial_time;
        if (index == 0) {
          initial_time = plugin.CurrentTime() + 60;
        } else {
          initial_time =
              plugin.FlightPlanGetManoeuvre(vessel_guid,
                                            index - 1).final_time + 60;
        }
        var editor = new BurnEditor(adapter_,
                                    predicted_vessel,
                                    initial_time,
                                    index,
                                    get_burn_at_index : burn_editors_.
                                        ElementAtOrDefault);
        editor.minimized = false;
        Burn candidate_burn = editor.Burn();
        var status = plugin.FlightPlanInsert(
            vessel_guid,
            candidate_burn,
            index);

        // The previous call did not necessarily create a manœuvre.  Check if
        // we need to add an editor.
        int number_of_manœuvres =
            plugin.FlightPlanNumberOfManoeuvres(vessel_guid);
        if (number_of_manœuvres > burn_editors_.Count) {
          editor.Reset(plugin.FlightPlanGetManoeuvre(vessel_guid, index));
          burn_editors_.Insert(index, editor);
          UpdateBurnEditorIndices();
          UpdateStatus(status, index);
          Shrink();
          return true;
        }
        // TODO(phl): The error messaging here will be either confusing or
        // wrong.  The messages should mention the new manœuvre without
        // numbering it, since the numbering has not changed (“the new manœuvre
        // would overlap with manœuvre #1 or manœuvre #2” or something along
        // these lines).
        UpdateStatus(status, index);
      }
    }
    return false;
  }

  internal static string FormatPositiveTimeSpan(double seconds) {
    return new PrincipiaTimeSpan(seconds).FormatPositive(
        with_leading_zeroes: true,
        with_seconds: true);
  }

  internal static string FormatTimeSpan (double seconds) {
    return new PrincipiaTimeSpan(seconds).Format(
        with_leading_zeroes: true,
        with_seconds: true);
  }

  internal string FormatPlanLength(double value) {
    return FormatPositiveTimeSpan(value -
                                  plugin.FlightPlanGetInitialTime(
                                      predicted_vessel.id.ToString()));
  }

  internal bool TryParsePlanLength(string text, out double value) {
    value = 0;
    if (!PrincipiaTimeSpan.TryParse(text,
                                    out PrincipiaTimeSpan ts)) {
      return false;
    }
    value = ts.total_seconds +
            plugin.FlightPlanGetInitialTime(predicted_vessel.id.ToString());
    return true;
  }

  private void ResetStatus() {
    status_ = Status.OK;
    first_error_manœuvre_ = null;
    message_was_displayed_ = false;
  }

  private void UpdateStatus(Status status, int? error_manœuvre) {
    if (message_was_displayed_) {
      ResetStatus();
    }
    if (status_.ok() && !status.ok()) {
      status_ = status;
      first_error_manœuvre_ = error_manœuvre;
    }
  }

  private string GetStatusMessage() {
    string vessel_guid = predicted_vessel?.id.ToString();
    string message = "";
    if (vessel_guid != null && !status_.ok()) {
      int anomalous_manœuvres =
          plugin.FlightPlanNumberOfAnomalousManoeuvres(vessel_guid);
      int manœuvres = plugin.FlightPlanNumberOfManoeuvres(vessel_guid);
      double actual_final_time =
          plugin.FlightPlanGetActualFinalTime(vessel_guid);
      bool timed_out = actual_final_time < final_time_.value;

      string remedy_message =
          L10N.CacheFormat(
              "#Principia_FlightPlan_StatusMessage_ChangeFlightPlan");  // Preceded by "Try".
      string status_message = L10N.CacheFormat(
          "#Principia_FlightPlan_StatusMessage_FailedError",
          status_.error,
          status_.message);
      string time_out_message =
          timed_out ? L10N.CacheFormat(
                          "#Principia_FlightPlan_StatusMessage_TimeOut",
                          FormatPositiveTimeSpan(
                              actual_final_time -
                              plugin.FlightPlanGetInitialTime(vessel_guid)))
                    : "";
      if (status_.is_aborted()) {
        status_message = L10N.CacheFormat(
            "#Principia_FlightPlan_StatusMessage_MaxSteps",
            time_out_message);
        remedy_message =
            L10N.CacheFormat("#Principia_FlightPlan_StatusMessage_MaxSegment");
      } else if (status_.is_failed_precondition()) {
        status_message = L10N.CacheFormat(
            "#Principia_FlightPlan_StatusMessage_Singularity",
            time_out_message);
        remedy_message =
            L10N.CacheFormat(
                "#Principia_FlightPlan_StatusMessage_AvoidingCollision");
      } else if (status_.is_invalid_argument()) {
        status_message = L10N.CacheFormat(
            "#Principia_FlightPlan_StatusMessage_Infinite",
            first_error_manœuvre_.Value + 1);
        remedy_message = L10N.CacheFormat(
            "#Principia_FlightPlan_StatusMessage_Adjust",
            first_error_manœuvre_.Value + 1);
      } else if (status_.is_out_of_range()) {
        if (first_error_manœuvre_.HasValue) {
          status_message = L10N.CacheFormat(
              "#Principia_FlightPlan_StatusMessage_OutRange1",
              first_error_manœuvre_.Value + 1,
              first_error_manœuvre_.Value == 0
                  ? L10N.CacheFormat(
                      "#Principia_FlightPlan_StatusMessage_OutRange2")
                  : L10N.CacheFormat(
                      "#Principia_FlightPlan_StatusMessage_OutRange3",
                      first_error_manœuvre_.Value),
              manœuvres == 0 || first_error_manœuvre_.Value == manœuvres - 1
                  ? L10N.CacheFormat(
                      "#Principia_FlightPlan_StatusMessage_OutRange4")
                  : L10N.CacheFormat(
                      "#Principia_FlightPlan_StatusMessage_OutRange5",
                      first_error_manœuvre_.Value + 2));
          remedy_message =  L10N.CacheFormat(
              "#Principia_FlightPlan_StatusMessage_OutRange6",
              manœuvres == 0 || first_error_manœuvre_.Value == manœuvres - 1
                  ? L10N.CacheFormat(
                      "#Principia_FlightPlan_StatusMessage_OutRange7")
                  : "",
              first_error_manœuvre_.Value + 1);
        } else {
          status_message =
              L10N.CacheFormat("#Principia_FlightPlan_StatusMessage_TooShort");
          remedy_message =
              L10N.CacheFormat("#Principia_FlightPlan_StatusMessage_Increase");
        }
      } else if (status_.is_unavailable()) {
        status_message =
            L10N.CacheFormat("#Principia_FlightPlan_StatusMessage_CantRebase");
        remedy_message =
            L10N.CacheFormat("#Principia_FlightPlan_StatusMessage_WaitFinish");
      }

      if (anomalous_manœuvres > 0) {
        message = L10N.CacheFormat(
            "#Principia_FlightPlan_StatusMessage_Last",
            anomalous_manœuvres,
            status_message ,
            remedy_message,
            (anomalous_manœuvres < manœuvres
                 ? L10N.CacheFormat("#Principia_FlightPlan_StatusMessage_Last2",
                                    manœuvres - anomalous_manœuvres)
                 : ""));
      } else {
        message =
            L10N.CacheFormat("#Principia_FlightPlan_StatusMessage_Result",
                             status_message,
                             remedy_message);
      }
    }
    message_was_displayed_ = true;
    return message;
  }

  private void UpdateBurnEditorIndices() {
    for (int j = 0; j < burn_editors_.Count; ++j) {
      burn_editors_[j].index = j;
    }
  }

  private void UpdateAdaptiveStepParameters(
      ref FlightPlanAdaptiveStepParameters parameters) {
    parameters.length_integration_tolerance =
        integration_tolerances_[length_integration_tolerance_index_];
    parameters.speed_integration_tolerance =
        integration_tolerances_[speed_integration_tolerance_index_];
    parameters.max_steps = max_steps_[max_steps_index_];
  }

  private IntPtr plugin => adapter_.Plugin();

  private static readonly double[] integration_tolerances_ =
      {1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6};
  private static readonly string[] length_integration_tolerances_names_ = {
      "1 µm", "10 µm", "100 µm", "1 mm", "1 cm", "10 cm", "1 m", "10 m",
      "100 m", "1 km", "10 km", "100 km", "1000 km"
  };
  private static readonly long[] max_steps_ = {
      1 << 6, 1 << 8, 1 << 10, 1 << 12, 1 << 14, 1 << 16, 1 << 18, 1 << 20
  };

  private readonly PrincipiaPluginAdapter adapter_;
  private readonly PredictedVessel predicted_vessel_;

  // Because this class is stateful (it holds the burn_editors_) we must detect
  // if the vessel changed.  Hence the caching of the vessel.
  private Vessel previous_predicted_vessel_;

  private List<BurnEditor> burn_editors_;
  private PlannedOrbitAnalyser final_trajectory_analyser_;
  private readonly DifferentialSlider final_time_;
  private int? first_future_manœuvre_;
  private int number_of_anomalous_manœuvres_ = 0;

  private int length_integration_tolerance_index_;
  private int speed_integration_tolerance_index_;
  private int max_steps_index_;
  private bool show_guidance_ = false;
  private float warning_height_ = 1;

  private Status status_ = Status.OK;
  private int? first_error_manœuvre_;  // May exceed the number of manœuvres.
  private bool message_was_displayed_ = false;

  private const double log10_time_lower_rate = 0.0;
  private const double log10_time_upper_rate = 7.0;

  private const int max_flight_plans = 10;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
