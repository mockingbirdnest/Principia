using System;
using System.Collections.Generic;
using System.Linq;

namespace principia {
namespace ksp_plugin_adapter {

class FlightPlanner : VesselSupervisedWindowRenderer {
  public FlightPlanner(PrincipiaPluginAdapter adapter,
                       PredictedVessel predicted_vessel)
      : base(adapter, predicted_vessel) {
    adapter_ = adapter;
    final_time_ = new DifferentialSlider(
                      label            : "Plan length",
                      unit             : null,
                      log10_lower_rate : log10_time_lower_rate_,
                      log10_upper_rate : log10_time_upper_rate_,
                      min_value        : 10,
                      max_value        : double.PositiveInfinity,
                      formatter        : FormatPlanLength,
                      parser           : TryParsePlanLength,
                      field_width      : 7);
  }

  public void RenderButton() {
    RenderButton("Flight plan...");
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

    node.SetValue("show_guidance",
                  show_guidance_,
                  createIfNotFound : true);
  }

  protected override string Title => "Flight plan";

  protected override void RenderWindow(int window_id) {
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

    if (plugin.FlightPlanExists(vessel_guid)) {
      RenderFlightPlan(vessel_guid);
    } else if (UnityEngine.GUILayout.Button("Create flight plan")) {
      plugin.FlightPlanCreate(vessel_guid,
                              plugin.CurrentTime() + 3600,
                              predicted_vessel.GetTotalMass());
      final_time_.value = plugin.FlightPlanGetDesiredFinalTime(vessel_guid);
      Shrink();
    }
    UnityEngine.GUI.DragWindow();
  }

  private void UpdateVesselAndBurnEditors() {
    {
      string vessel_guid = predicted_vessel?.id.ToString();
      if (vessel_guid == null ||
          previous_predicted_vessel_ != predicted_vessel ||
          !plugin.FlightPlanExists(vessel_guid) ||
          plugin.FlightPlanNumberOfManoeuvres(vessel_guid) !=
              burn_editors_?.Count) {
        if (burn_editors_ != null) {
          foreach (BurnEditor editor in burn_editors_) {
            editor.Close();
          }
          burn_editors_ = null;
          Shrink();
        }
        previous_predicted_vessel_ = predicted_vessel;
      }
    }

    if (burn_editors_ == null) {
      string vessel_guid = predicted_vessel?.id.ToString();
      if (vessel_guid != null &&
          plugin.FlightPlanExists(vessel_guid)) {
        burn_editors_ = new List<BurnEditor>();
        final_time_.value = plugin.FlightPlanGetDesiredFinalTime(vessel_guid);
        for (int i = 0;
              i < plugin.FlightPlanNumberOfManoeuvres(vessel_guid);
              ++i) {
          // Dummy initial time, we call |Reset| immediately afterwards.
          burn_editors_.Add(
              new BurnEditor(adapter_,
                             predicted_vessel,
                             initial_time  : 0,
                             index         : burn_editors_.Count,
                             previous_burn : burn_editors_.LastOrDefault()));
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
        var status = plugin.FlightPlanSetDesiredFinalTime(vessel_guid,
                                                          final_time_.value);
        UpdateStatus(status, null);
      }
      // Always refresh the final time from C++ as it may have changed because
      // the last burn changed.
      final_time_.value =
          plugin.FlightPlanGetDesiredFinalTime(vessel_guid);

      FlightPlanAdaptiveStepParameters parameters =
          plugin.FlightPlanGetAdaptiveStepParameters(vessel_guid);

      using (new UnityEngine.GUILayout.HorizontalScope()) {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Label("Max. steps per segment:",
                                      GUILayoutWidth(6));
          const int factor = 4;
          if (parameters.max_steps <= 100) {
            UnityEngine.GUILayout.Button("min");
          } else if (UnityEngine.GUILayout.Button("-")) {
            parameters.max_steps /= factor;
            var status = plugin.FlightPlanSetAdaptiveStepParameters(vessel_guid,
                                                                    parameters);
            UpdateStatus(status, null);
          }
          UnityEngine.GUILayout.TextArea(parameters.max_steps.ToString(),
                                          GUILayoutWidth(3));
          if (parameters.max_steps >= long.MaxValue / factor) {
            UnityEngine.GUILayout.Button("max");
          } else if (UnityEngine.GUILayout.Button("+")) {
            parameters.max_steps *= factor;
            var status = plugin.FlightPlanSetAdaptiveStepParameters(vessel_guid,
                                                                    parameters);
            UpdateStatus(status, null);
          }
        }
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Label("Tolerance:",
                                      GUILayoutWidth(3));
          if (parameters.length_integration_tolerance <= 1e-6) {
            UnityEngine.GUILayout.Button("min");
          } else if (UnityEngine.GUILayout.Button("-")) {
            parameters.length_integration_tolerance /= 2;
            parameters.speed_integration_tolerance /= 2;
            var status = plugin.FlightPlanSetAdaptiveStepParameters(vessel_guid,
                                                                    parameters);
            UpdateStatus(status, null);
          }
          UnityEngine.GUILayout.TextArea(
              parameters.length_integration_tolerance.ToString("0.0e0") + " m",
              GUILayoutWidth(3));
          if (parameters.length_integration_tolerance >= 1e6) {
            UnityEngine.GUILayout.Button("max");
          } else if (UnityEngine.GUILayout.Button("+")) {
            parameters.length_integration_tolerance *= 2;
            parameters.speed_integration_tolerance *= 2;
            var status = plugin.FlightPlanSetAdaptiveStepParameters(vessel_guid,
                                                                    parameters);
            UpdateStatus(status, null);
          }
        }
      }

      double Δv = (from burn_editor in burn_editors_
                   select burn_editor.Δv()).Sum();
      UnityEngine.GUILayout.Label(
          "Total Δv : " + Δv.ToString("0.000") + " m/s");

      {
        var style = Style.Warning(Style.Multiline(UnityEngine.GUI.skin.label));
        string message = GetStatusMessage();
        // Size the label explicitly so that it doesn't decrease when the
        // message goes away: that causes annoying flicker.  The enclosing
        // window has a width of 20 units, but not all of that is available,
        // hence 19.
        warning_height_ = Math.Max(
            warning_height_,
            style.CalcHeight(new UnityEngine.GUIContent(message), Width(19)));
        UnityEngine.GUILayout.Label(
            message, style, UnityEngine.GUILayout.Height(warning_height_));
      }

      if (burn_editors_.Count == 0 &&
          UnityEngine.GUILayout.Button("Delete flight plan")) {
        plugin.FlightPlanDelete(vessel_guid);
        ResetStatus();
        Shrink();
        // The state change will happen the next time we go through OnGUI.
      } else {
        if (burn_editors_.Count > 0) {
          RenderUpcomingEvents();
        }

        // Compute the final times for each manœuvre before displaying them.
        var final_times = new List<double>();
        for (int i = 0; i < burn_editors_.Count - 1; ++i) {
          final_times.Add(
              plugin.FlightPlanGetManoeuvre(vessel_guid, i + 1).
                  burn.initial_time);
        }
        // Allow extending the flight plan.
        final_times.Add(double.PositiveInfinity);

        for (int i = 0; i < burn_editors_.Count; ++i) {
          Style.HorizontalLine();
          BurnEditor burn = burn_editors_[i];
          if (burn.Render(header          : "Manœuvre #" + (i + 1),
                          anomalous       :
                              i >= (burn_editors_.Count -
                                    number_of_anomalous_manœuvres_),
                          burn_final_time : final_times[i])) {
            var status = plugin.FlightPlanReplace(vessel_guid, burn.Burn(), i);
            UpdateStatus(status, i);
            burn.Reset(plugin.FlightPlanGetManoeuvre(vessel_guid, i));
          }
        }

        if (burn_editors_.Count > 0) {
          if (UnityEngine.GUILayout.Button(
                  "Delete last manœuvre",
                  UnityEngine.GUILayout.ExpandWidth(true))) {
            var status = plugin.FlightPlanRemoveLast(vessel_guid);
            UpdateStatus(status, null);
            burn_editors_.Last().Close();
            burn_editors_.RemoveAt(burn_editors_.Count - 1);
            Shrink();
          }
        }
        if (UnityEngine.GUILayout.Button(
                "Add manœuvre",
                UnityEngine.GUILayout.ExpandWidth(true))) {
          double initial_time;
          if (burn_editors_.Count == 0) {
            initial_time = plugin.CurrentTime() + 60;
          } else {
            initial_time = plugin.FlightPlanGetManoeuvre(
                               vessel_guid,
                               burn_editors_.Count - 1).final_time + 60;
          }
          var editor =
              new BurnEditor(adapter_,
                             predicted_vessel,
                             initial_time,
                             index         : burn_editors_.Count,
                             previous_burn : burn_editors_.LastOrDefault());
          Burn candidate_burn = editor.Burn();
          var status = plugin.FlightPlanAppend(vessel_guid, candidate_burn);

          // The previous call did not necessarily create a manœuvre.  Check if
          // we need to add an editor.
          int number_of_manœuvres =
              plugin.FlightPlanNumberOfManoeuvres(vessel_guid);
          if (number_of_manœuvres > burn_editors_.Count) {
            editor.Reset(plugin.FlightPlanGetManoeuvre(
                vessel_guid, number_of_manœuvres - 1));
            burn_editors_.Add(editor);
            UpdateStatus(status, number_of_manœuvres - 1);
          } else {
            UpdateStatus(status, number_of_manœuvres);
          }

          Shrink();
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
          UnityEngine.GUILayout.Label("Upcoming manœuvre #" +
                                      (first_future_manœuvre + 1) + ":");
          UnityEngine.GUILayout.Label(
              "Ignition " +
              FormatTimeSpan(current_time - manœuvre.burn.initial_time),
              style : Style.RightAligned(UnityEngine.GUI.skin.label));
        }
      } else {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Label("Ongoing manœuvre #" +
                                      (first_future_manœuvre + 1) + ":");
          UnityEngine.GUILayout.Label(
              "Cutoff " + FormatTimeSpan(current_time - manœuvre.final_time),
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
          show_guidance_ =
              UnityEngine.GUILayout.Toggle(show_guidance_, "Show on navball");
          if (UnityEngine.GUILayout.Button("Warp to manœuvre")) {
            TimeWarp.fetch.WarpTo(manœuvre.burn.initial_time - 60);
          }
        }
      }
    } else {
      // Reserve some space to avoid the UI changing shape if we have
      // nothing to say.
      UnityEngine.GUILayout.Label("All manœuvres are in the past",
                                  Style.Warning(UnityEngine.GUI.skin.label));
      UnityEngine.GUILayout.Space(Width(1));
    }
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
    return FormatPositiveTimeSpan(
        value - plugin.FlightPlanGetInitialTime(
            predicted_vessel.id.ToString()));
  }

  internal bool TryParsePlanLength(string text, out double value) {
    value = 0;
    if (!PrincipiaTimeSpan.TryParse(text,
                                    with_seconds: true,
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

      string remedy_message = "changing the flight plan";  // Preceded by "Try".
      string status_message = "computation failed";  // Preceded by "The".
      string time_out_message = timed_out
                                    ? " after " +
                                      FormatPositiveTimeSpan(
                                          actual_final_time -
                                          plugin.FlightPlanGetInitialTime(
                                              vessel_guid))
                                    : "";
      if (status_.is_aborted()) {
        status_message = "integrator reached the maximum number of steps" +
                         time_out_message;
        remedy_message = "increasing 'Max. steps per segment' or avoiding " +
                         "collisions with a celestial";
      } else if (status_.is_failed_precondition()) {
        status_message = "integrator encountered a singularity (probably the " +
                         "centre of a celestial)" + time_out_message;
        remedy_message = "avoiding collisions with a celestial";
      } else if (status_.is_invalid_argument()) {
        status_message = "manœuvre #" + (first_error_manœuvre_.Value + 1) +
                         " would result in an infinite or indeterminate " +
                         "velocity";
        remedy_message = "adjusting the duration of manœuvre #" +
                         (first_error_manœuvre_.Value + 1);
      } else if (status_.is_out_of_range()) {
        if (first_error_manœuvre_.HasValue) {
          status_message = "manœuvre #" + (first_error_manœuvre_.Value + 1) +
                           " overlaps with " +
                           ((first_error_manœuvre_.Value == 0)
                                ? "the start of the flight plan"
                                : "manœuvre #" +
                                  first_error_manœuvre_.Value) + " or " +
                           ((manœuvres == 0 ||
                             first_error_manœuvre_.Value == manœuvres - 1)
                                ? "the end of the flight plan"
                                : "manœuvre #" +
                                  (first_error_manœuvre_.Value + 2));
          remedy_message = ((manœuvres == 0 ||
                             first_error_manœuvre_.Value == manœuvres - 1)
                               ? "extending the flight plan or "
                               : "") +
                           "adjusting the initial time or reducing the " +
                           "duration of manœuvre #" +
                           (first_error_manœuvre_.Value + 1);
        } else {
          status_message = "flight plan is too short";
          remedy_message = "increasing the flight plan duration";
        }
      }

      if (anomalous_manœuvres > 0) {
        message = "The last " + anomalous_manœuvres + " manœuvres could " +
                  "not be drawn because the " + status_message + "; try " +
                  remedy_message +
                  (anomalous_manœuvres < manœuvres
                       ?  " or adjusting manœuvre #" +
                          (manœuvres - anomalous_manœuvres)
                       : "") + ".";
      } else {
        message = "The " + status_message + "; try " + remedy_message + ".";
      }
    }
    message_was_displayed_ = true;
    return message;
  }

  private IntPtr plugin => adapter_.Plugin();

  private static bool is_stock_day => GameSettings.KERBIN_TIME &&
                                      KSPUtil.dateTimeFormatter.Day ==
                                      6 * 60 * 60;

  private readonly PrincipiaPluginAdapter adapter_;

  // Because this class is stateful (it holds the burn_editors_) we must detect
  // if the vessel changed.  Hence the caching of the vessel.
  private Vessel previous_predicted_vessel_;

  private List<BurnEditor> burn_editors_;
  private readonly DifferentialSlider final_time_;
  private int? first_future_manœuvre_;
  private int number_of_anomalous_manœuvres_ = 0;

  private bool show_guidance_ = false;
  private float warning_height_ = 1;

  private Status status_ = Status.OK;
  private int? first_error_manœuvre_;  // May exceed the number of manœuvres.
  private bool message_was_displayed_ = false;

  private const double log10_time_lower_rate_ = 0.0;
  private const double log10_time_upper_rate_ = 7.0;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
