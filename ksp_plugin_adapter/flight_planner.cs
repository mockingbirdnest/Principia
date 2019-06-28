using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.RegularExpressions;

namespace principia {
namespace ksp_plugin_adapter {

class FlightPlanner : SupervisedWindowRenderer {
  public FlightPlanner(PrincipiaPluginAdapter adapter) : base(adapter) {
    adapter_ = adapter;
    final_time_ = new DifferentialSlider(
                      label            : "Plan length",
                      unit             : null,
                      log10_lower_rate : log10_time_lower_rate,
                      log10_upper_rate : log10_time_upper_rate,
                      min_value        : 10,
                      max_value        : double.PositiveInfinity,
                      formatter        : FormatPlanLength,
                      parser           : TryParsePlanLength);
  }

  public void RenderButton() {
    if (UnityEngine.GUILayout.Button("Flight plan...")) {
      Toggle();
    }
    // Override the state of the toggle if there is no active vessel.
    string vessel_guid = vessel_?.id.ToString();
    if (vessel_guid == null || !plugin.HasVessel(vessel_guid)) {
      Hide();
      vessel_ = FlightGlobals.ActiveVessel;
    }
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
    string vessel_guid = vessel_?.id.ToString();
    if (vessel_guid == null || !plugin.HasVessel(vessel_guid)) {
      return;
    }

    if (plugin.FlightPlanExists(vessel_guid)) {
      RenderFlightPlan(vessel_guid);
    } else if (UnityEngine.GUILayout.Button("Create flight plan")) {
      plugin.FlightPlanCreate(vessel_guid,
                              plugin.CurrentTime() + 1000,
                              vessel_.GetTotalMass());
      final_time_.value = plugin.FlightPlanGetDesiredFinalTime(vessel_guid);
      Shrink();
    }
    UnityEngine.GUI.DragWindow();
  }

  private void UpdateVesselAndBurnEditors() {
    {
      string vessel_guid = vessel_?.id.ToString();
      if (vessel_guid == null ||
          vessel_ != FlightGlobals.ActiveVessel ||
          !plugin.HasVessel(vessel_guid) ||
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
        vessel_ = FlightGlobals.ActiveVessel;
      }
    }

    if (burn_editors_ == null) {
      string vessel_guid = vessel_?.id.ToString();
      if (vessel_guid != null &&
          plugin.HasVessel(vessel_guid) &&
          plugin.FlightPlanExists(vessel_guid)) {
        burn_editors_ = new List<BurnEditor>();
        final_time_.value = plugin.FlightPlanGetDesiredFinalTime(vessel_guid);
        for (int i = 0;
              i < plugin.FlightPlanNumberOfManoeuvres(vessel_guid);
              ++i) {
          // Dummy initial time, we call |Reset| immediately afterwards.
          burn_editors_.Add(
              new BurnEditor(adapter_,
                             vessel_,
                             initial_time  : 0,
                             index         : burn_editors_.Count,
                             previous_burn : burn_editors_.LastOrDefault()));
          burn_editors_.Last().Reset(
              plugin.FlightPlanGetManoeuvre(vessel_guid, i));
        }
      }
    }

    if (burn_editors_ != null) {
      string vessel_guid = vessel_?.id.ToString();
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
    }
  }

  private void RenderFlightPlan(string vessel_guid) {
    using (new UnityEngine.GUILayout.VerticalScope()) {
      if (final_time_.Render(enabled : true)) {
        var status = plugin.FlightPlanSetDesiredFinalTime(vessel_guid,
                                                          final_time_.value);
        UpdateStatus(status, null);
        final_time_.value =
            plugin.FlightPlanGetDesiredFinalTime(vessel_guid);
      }
      double actual_final_time =
          plugin.FlightPlanGetActualFinalTime(vessel_guid);

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
        final_times.Add(plugin.FlightPlanGetActualFinalTime(vessel_guid));
        int number_of_anomalous_manœuvres =
            plugin.FlightPlanNumberOfAnomalousManoeuvres(vessel_guid);

        for (int i = 0; i < burn_editors_.Count; ++i) {
          Style.HorizontalLine();
          BurnEditor burn = burn_editors_[i];
          if (burn.Render(header     : "Manœuvre #" + (i + 1),
                          anomalous  : i >= (burn_editors_.Count -
                                             number_of_anomalous_manœuvres),
                          final_time : final_times[i])) {
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
                             vessel_,
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
    string vessel_guid = vessel_.id.ToString();
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
              "Ignition " + FormatTimeSpan(TimeSpan.FromSeconds(
                                current_time - manœuvre.burn.initial_time)),
              style : Style.RightAligned(UnityEngine.GUI.skin.label));
        }
      } else {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Label("Ongoing manœuvre #" +
                                      (first_future_manœuvre + 1) + ":");
          UnityEngine.GUILayout.Label(
              "Cutoff " + FormatTimeSpan(TimeSpan.FromSeconds(
                              current_time - manœuvre.final_time)),
              style : Style.RightAligned(UnityEngine.GUI.skin.label));
        }
      }
      // In career mode, the patched conic solver may be null.  In that case
      // we do not offer the option of showing the manœuvre on the navball,
      // even though the flight planner is still available to plan it.
      // TODO(egg): We may want to consider setting the burn vector directly
      // rather than going through the solver.
      if (vessel_.patchedConicSolver != null) {
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

  internal static string FormatPositiveTimeSpan(TimeSpan span) {
    return (GameSettings.KERBIN_TIME
                ? (span.Days * 4 + span.Hours / 6).ToString("0000;0000") +
                      " d6 " + (span.Hours % 6).ToString("0;0") + " h "
                : span.Days.ToString("000;000") + " d " +
                      span.Hours.ToString("00;00") + " h ") +
           span.Minutes.ToString("00;00") + " min " +
           (span.Seconds + span.Milliseconds / 1000m).ToString("00.0;00.0") +
           " s";
  }

  internal static string FormatTimeSpan (TimeSpan span) {
    return span.Ticks.ToString("+;-") + FormatPositiveTimeSpan(span);
  }

  internal static bool TryParseTimeSpan(string str, out TimeSpan value) {
    value = TimeSpan.Zero;
    // Using a technology that is customarily used to parse HTML.
    string pattern = @"^[+]?\s*(\d+)\s*" +
                     (GameSettings.KERBIN_TIME ? "d6" : "d") +
                     @"\s*(\d+)\s*h\s*(\d+)\s*min\s*([0-9.,']+)\s*s$";
    var regex = new Regex(pattern);
    var match = regex.Match(str);
    if (!match.Success) {
      return false;
    }
    string days = match.Groups[1].Value;
    string hours = match.Groups[2].Value;
    string minutes = match.Groups[3].Value;
    string seconds = match.Groups[4].Value;
    if (!int.TryParse(days, out int d) ||
        !int.TryParse(hours, out int h) ||
        !int.TryParse(minutes, out int min) ||
        !double.TryParse(seconds.Replace(',', '.'),
                         NumberStyles.AllowDecimalPoint |
                         NumberStyles.AllowThousands,
                         Culture.culture.NumberFormat,
                         out double s)) {
      return false;
    }
    value = TimeSpan.FromDays((double)d / (GameSettings.KERBIN_TIME ? 4 : 1)) +
            TimeSpan.FromHours(h) +
            TimeSpan.FromMinutes(min) +
            TimeSpan.FromSeconds(s);
    return true;
  }

  internal string FormatPlanLength(double value) {
    return FormatPositiveTimeSpan(TimeSpan.FromSeconds(
               value -
               plugin.FlightPlanGetInitialTime(vessel_.id.ToString())));
  }

  internal bool TryParsePlanLength(string str, out double value) {
    value = 0;
    if (!TryParseTimeSpan(str, out TimeSpan ts)) {
      return false;
    }
    value = ts.TotalSeconds +
            plugin.FlightPlanGetInitialTime(vessel_.id.ToString());
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
    string vessel_guid = vessel_?.id.ToString();
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
      string time_out_message =
          timed_out ? " after " +
                      FormatPositiveTimeSpan(TimeSpan.FromSeconds(
                           actual_final_time -
                           plugin.FlightPlanGetInitialTime(vessel_guid)))
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
                  remedy_message + " or adjusting manœuvre " +
                  (manœuvres - anomalous_manœuvres - 1) + ".";
      } else {
        message = "The " + status_message + "; try " + remedy_message + ".";
      }
    }
    message_was_displayed_ = true;
    return message;
  }

  private IntPtr plugin => adapter_.Plugin();

  private readonly PrincipiaPluginAdapter adapter_;
  private Vessel vessel_;
  private List<BurnEditor> burn_editors_;
  private readonly DifferentialSlider final_time_;
  private int? first_future_manœuvre_;

  private bool show_guidance_ = false;
  private float warning_height_ = 1;

  private Status status_ = Status.OK;
  private int? first_error_manœuvre_;  // May exceed the number of manœuvres.
  private bool message_was_displayed_ = false;
  
  private const double log10_time_lower_rate = 0.0;
  private const double log10_time_upper_rate = 7.0;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
