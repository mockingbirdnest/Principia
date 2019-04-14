using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
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

  protected override String Title => "Flight plan";

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
        for (int i = 0;
              i < plugin.FlightPlanNumberOfManoeuvres(vessel_guid);
              ++i) {
          // Dummy initial time, we call |Reset| immediately afterwards.
          final_time_.value =
              plugin.FlightPlanGetDesiredFinalTime(vessel_guid);
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
      first_future_manoeuvre_ = null;
      for (int i = 0; i < burn_editors_.Count; ++i) {
        NavigationManoeuvre manoeuvre =
            plugin.FlightPlanGetManoeuvre(vessel_guid, i);
        if (current_time < manoeuvre.final_time) {
          first_future_manoeuvre_ = i;
          break;
        }
      }
    }
  }

  private void RenderFlightPlan(string vessel_guid) {
    using (new UnityEngine.GUILayout.VerticalScope()) {
      if (final_time_.Render(enabled : true)) {
        plugin.FlightPlanSetDesiredFinalTime(vessel_guid,
                                              final_time_.value);
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
            plugin.FlightPlanSetAdaptiveStepParameters(vessel_guid,
                                                       parameters);
          }
          UnityEngine.GUILayout.TextArea(parameters.max_steps.ToString(),
                                          GUILayoutWidth(3));
          if (parameters.max_steps >= Int64.MaxValue / factor) {
            UnityEngine.GUILayout.Button("max");
          } else if (UnityEngine.GUILayout.Button("+")) {
            parameters.max_steps *= factor;
            plugin.FlightPlanSetAdaptiveStepParameters(vessel_guid,
                                                       parameters);
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
            plugin.FlightPlanSetAdaptiveStepParameters(vessel_guid,
                                                       parameters);
          }
          UnityEngine.GUILayout.TextArea(
              parameters.length_integration_tolerance.ToString("0.0e0") + " m",
              GUILayoutWidth(3));
          if (parameters.length_integration_tolerance >= 1e6) {
            UnityEngine.GUILayout.Button("max");
          } else if (UnityEngine.GUILayout.Button("+")) {
            parameters.length_integration_tolerance *= 2;
            parameters.speed_integration_tolerance *= 2;
            plugin.FlightPlanSetAdaptiveStepParameters(vessel_guid,
                                                       parameters);
          }
        }
      }

      double Δv = (from burn_editor in burn_editors_
                   select burn_editor.Δv()).Sum();
      UnityEngine.GUILayout.Label(
          "Total Δv : " + Δv.ToString("0.000") + " m/s");

      string message = "";
      if (final_time_.value != actual_final_time) {
        message = "Timed out after " +
                  FormatPositiveTimeSpan(TimeSpan.FromSeconds(
                      actual_final_time -
                      plugin.FlightPlanGetInitialTime(vessel_guid)));
      }
      UnityEngine.GUILayout.Label(
          message, Style.Warning(UnityEngine.GUI.skin.label));

      if (burn_editors_.Count == 0 &&
          UnityEngine.GUILayout.Button("Delete flight plan")) {
        plugin.FlightPlanDelete(vessel_guid);
        Shrink();
        // The state change will happen the next time we go through OnGUI.
      } else {
        if (burn_editors_.Count > 0) {
          RenderUpcomingEvents();
        }
        for (int i = 0; i < burn_editors_.Count - 1; ++i) {
          Style.HorizontalLine();
          burn_editors_[i].Render(header: "Manœuvre #" + (i + 1),
                                  enabled : false);
        }
        if (burn_editors_.Count > 0) {
          Style.HorizontalLine();
          BurnEditor last_burn = burn_editors_.Last();
          if (last_burn.Render(header            : "Editing manœuvre #" +
                                                   (burn_editors_.Count),
                               enabled           : true,
                               actual_final_time : actual_final_time)) {
            plugin.FlightPlanReplaceLast(vessel_guid, last_burn.Burn());
            last_burn.Reset(
                plugin.FlightPlanGetManoeuvre(vessel_guid,
                                              burn_editors_.Count - 1));
          }
          if (UnityEngine.GUILayout.Button(
                  "Delete last manœuvre",
                  UnityEngine.GUILayout.ExpandWidth(true))) {
            plugin.FlightPlanRemoveLast(vessel_guid);
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
          bool inserted = plugin.FlightPlanAppend(vessel_guid,
                                                  candidate_burn);
          if (inserted) {
            editor.Reset(plugin.FlightPlanGetManoeuvre(
                vessel_guid, burn_editors_.Count));
            burn_editors_.Add(editor);
          }
          Shrink();
        }
      }
    }
  }

  private void RenderUpcomingEvents() {
    string vessel_guid = vessel_.id.ToString();
    double current_time = plugin.CurrentTime();
    bool should_clear_guidance = true;

    Style.HorizontalLine();
    if (first_future_manoeuvre_.HasValue) {
      int first_future_manoeuvre = first_future_manoeuvre_.Value;
      NavigationManoeuvre manoeuvre =
          plugin.FlightPlanGetManoeuvre(vessel_guid, first_future_manoeuvre);
      if (manoeuvre.burn.initial_time > current_time) {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Label("Upcoming manœuvre #" +
                                      (first_future_manoeuvre + 1) + ":");
          UnityEngine.GUILayout.Label(
              "Ignition " + FormatTimeSpan(TimeSpan.FromSeconds(
                                current_time - manoeuvre.burn.initial_time)),
              style : Style.RightAligned(UnityEngine.GUI.skin.label));
        }
      } else {
        using (new UnityEngine.GUILayout.HorizontalScope()) {
          UnityEngine.GUILayout.Label("Ongoing manœuvre #" +
                                      (first_future_manoeuvre + 1) + ":");
          UnityEngine.GUILayout.Label(
              "Cutoff " + FormatTimeSpan(TimeSpan.FromSeconds(
                              current_time - manoeuvre.final_time)),
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
            TimeWarp.fetch.WarpTo(manoeuvre.burn.initial_time - 60);
          }
        }
        should_clear_guidance &= ShowGuidance(manoeuvre,
                                              first_future_manoeuvre);
      }
    } else {
      // Reserve some space to avoid the UI changing shape if we have
      // nothing to say.
      UnityEngine.GUILayout.Label("All manœuvres are in the past",
                                  Style.Warning(UnityEngine.GUI.skin.label));
      UnityEngine.GUILayout.Space(Width(1));
    }

    if (should_clear_guidance && guidance_node_ != null) {
      guidance_node_.RemoveSelf();
      guidance_node_ = null;
    }
  }

  // Returns true iff the guidance node should be cleared.
  internal bool ShowGuidance(NavigationManoeuvre manoeuvre,
                             int manoeuvre_index) {
    string vessel_guid = vessel_.id.ToString();
    XYZ guidance = plugin.FlightPlanGetGuidance(vessel_guid, manoeuvre_index);
    if (show_guidance_ &&
        !double.IsNaN(guidance.x + guidance.y + guidance.z)) {
      if (guidance_node_ == null ||
          !vessel_.patchedConicSolver.maneuverNodes.Contains(
              guidance_node_)) {
        while (vessel_.patchedConicSolver.maneuverNodes.Count > 0) {
          vessel_.patchedConicSolver.maneuverNodes.Last().RemoveSelf();
        }
        guidance_node_ = vessel_.patchedConicSolver.AddManeuverNode(
            manoeuvre.burn.initial_time);
      } else if (vessel_.patchedConicSolver.maneuverNodes.Count > 1) {
        while (vessel_.patchedConicSolver.maneuverNodes.Count > 1) {
          if (vessel_.patchedConicSolver.maneuverNodes.First() ==
              guidance_node_) {
            vessel_.patchedConicSolver.maneuverNodes.Last().RemoveSelf();
          } else {
            vessel_.patchedConicSolver.maneuverNodes.First().RemoveSelf();
          }
        }
      }
      var stock_orbit = guidance_node_.patch;
      Vector3d stock_velocity_at_node_time =
          stock_orbit.getOrbitalVelocityAtUT(
                            manoeuvre.burn.initial_time).xzy;
      Vector3d stock_displacement_from_parent_at_node_time =
          stock_orbit.getRelativePositionAtUT(
                            manoeuvre.burn.initial_time).xzy;
      UnityEngine.Quaternion stock_frenet_frame_to_world =
          UnityEngine.Quaternion.LookRotation(
              stock_velocity_at_node_time,
              Vector3d.Cross(
                  stock_velocity_at_node_time,
                  stock_displacement_from_parent_at_node_time));
      guidance_node_.DeltaV =
          ((Vector3d)manoeuvre.burn.delta_v).magnitude *
          (Vector3d)(UnityEngine.Quaternion.Inverse(
                          stock_frenet_frame_to_world) *
                      (Vector3d)guidance);
      guidance_node_.UT = manoeuvre.burn.initial_time;
      vessel_.patchedConicSolver.UpdateFlightPlan();
      return false;
    }
    return true;
  }

  internal static string FormatPositiveTimeSpan (TimeSpan span) {
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
    string pattern = @"^([-+]?)\s*(\d+)\s*"+
                     (GameSettings.KERBIN_TIME ? "d6" : "d") +
                     @"\s*(\d+)\s*h\s*(\d+)\s*min\s*([0-9.,']+)\s*s$";
    Regex regex = new Regex(pattern);
    var match = Regex.Match(str, pattern);
    if (!match.Success) {
      return false;
    }
    string sign = match.Groups[1].Value;
    string days = match.Groups[2].Value;
    string hours = match.Groups[3].Value;
    string minutes = match.Groups[4].Value;
    string seconds = match.Groups[5].Value;
    if (!Int32.TryParse(days, out int d) ||
        !Int32.TryParse(hours, out int h) ||
        !Int32.TryParse(minutes, out int min) ||
        !Double.TryParse(seconds.Replace(',', '.'),
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
    if (sign.Length > 0 && sign[0] == '-') {
      value = value.Negate();
    }
    return true;
  }

  internal string FormatPlanLength(double value) {
    return FormatPositiveTimeSpan(TimeSpan.FromSeconds(
               value -
               plugin.FlightPlanGetInitialTime(vessel_.id.ToString())));
  }

  internal bool TryParsePlanLength(string str, out double value) {
    value = 0;
    TimeSpan ts;
    if (!TryParseTimeSpan(str, out ts)) {
      return false;
    }
    value = ts.TotalSeconds +
            plugin.FlightPlanGetInitialTime(vessel_.id.ToString());
    return true;
  }

  private IntPtr plugin => adapter_.Plugin();

  private readonly PrincipiaPluginAdapter adapter_;
  private Vessel vessel_;
  private List<BurnEditor> burn_editors_;
  private DifferentialSlider final_time_;
  private int? first_future_manoeuvre_;

  private bool show_guidance_ = false;
  private ManeuverNode guidance_node_;
  
  private const double log10_time_lower_rate = 0.0;
  private const double log10_time_upper_rate = 7.0;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
