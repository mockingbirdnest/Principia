using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace principia {
namespace ksp_plugin_adapter {

class FlightPlanner : SupervisedWindowRenderer {
  public FlightPlanner(PrincipiaPluginAdapter adapter) : base(adapter) {
    adapter_ = adapter;
  }

  public void Initialize(IntPtr plugin) {
    plugin_ = plugin;
    final_time_ = new DifferentialSlider(
                label            : "Plan length",
                unit             : null,
                log10_lower_rate : Log10TimeLowerRate,
                log10_upper_rate : Log10TimeUpperRate,
                min_value        : 10,
                max_value        : double.PositiveInfinity,
                formatter        : value =>
                    FormatPositiveTimeSpan(
                        TimeSpan.FromSeconds(
                            value - plugin_.FlightPlanGetInitialTime(
                                        vessel_.id.ToString()))));
  }

  public void RenderButton() {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    if (UnityEngine.GUILayout.Button("Flight plan...")) {
      Toggle();
    }
    UnityEngine.GUI.skin = old_skin;
  }

  protected override String Title { get; } = "Flight plan";

  protected override void RenderWindow(int window_id) {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;

    using (new UnityEngine.GUILayout.VerticalScope()) {

      {
        string vessel_guid = vessel_?.id.ToString();
        if (vessel_guid == null ||
            vessel_ != FlightGlobals.ActiveVessel ||
            !plugin_.HasVessel(vessel_guid) ||
            !plugin_.FlightPlanExists(vessel_guid) ||
            plugin_.FlightPlanNumberOfManoeuvres(vessel_guid) !=
                burn_editors_?.Count) {
          Reset();
        }
      }

      if (vessel_ != null) {
        string vessel_guid = vessel_.id.ToString();
        if (burn_editors_ == null) {
          if (plugin_.HasVessel(vessel_guid)) {
            if (plugin_.FlightPlanExists(vessel_guid)) {
              // TODO(phl): Evil change of state between the two calls to
              // RenderPlanner.
              burn_editors_ = new List<BurnEditor>();
              for (int i = 0;
                   i < plugin_.FlightPlanNumberOfManoeuvres(vessel_guid);
                   ++i) {
                // Dummy initial time, we call |Reset| immediately afterwards.
                final_time_.value =
                    plugin_.FlightPlanGetDesiredFinalTime(vessel_guid);
                burn_editors_.Add(new BurnEditor(adapter_,
                                                 plugin_,
                                                 vessel_,
                                                 initial_time : 0));
                burn_editors_.Last().Reset(
                    plugin_.FlightPlanGetManoeuvre(vessel_guid, i));
              }
            } else {
              if (UnityEngine.GUILayout.Button("Create flight plan")) {
                plugin_.FlightPlanCreate(vessel_guid,
                                         plugin_.CurrentTime() + 1000,
                                         vessel_.GetTotalMass());
                final_time_.value =
                    plugin_.FlightPlanGetDesiredFinalTime(vessel_guid);
                Shrink();
              }
            }
          }
        } else {
          if (final_time_.Render(enabled: true)) {
            plugin_.FlightPlanSetDesiredFinalTime(vessel_guid,
                                                  final_time_.value);
            final_time_.value =
                plugin_.FlightPlanGetDesiredFinalTime(vessel_guid);
          }
          double actual_final_time =
              plugin_.FlightPlanGetActualFinalTime(vessel_guid);
          UnityEngine.GUILayout.TextField(
              (final_time_.value == actual_final_time)
                  ? ""
                  : "Timed out after " +
                        FormatPositiveTimeSpan(TimeSpan.FromSeconds(
                            actual_final_time -
                            plugin_.FlightPlanGetInitialTime(vessel_guid))));

          FlightPlanAdaptiveStepParameters parameters =
              plugin_.FlightPlanGetAdaptiveStepParameters(vessel_guid);

          using (new UnityEngine.GUILayout.HorizontalScope()) {
            using (new UnityEngine.GUILayout.HorizontalScope()) {
              UnityEngine.GUILayout.Label("Max. steps per segment:",
                                          UnityEngine.GUILayout.Width(150));
              const int factor = 4;
              if (parameters.max_steps <= 100) {
                UnityEngine.GUILayout.Button("min");
              } else if (UnityEngine.GUILayout.Button("-")) {
                parameters.max_steps /= factor;
                plugin_.FlightPlanSetAdaptiveStepParameters(vessel_guid,
                                                            parameters);
              }
              UnityEngine.GUILayout.TextArea(parameters.max_steps.ToString(),
                                             UnityEngine.GUILayout.Width(75));
              if (parameters.max_steps >= Int64.MaxValue / factor) {
                UnityEngine.GUILayout.Button("max");
              } else if (UnityEngine.GUILayout.Button("+")) {
                parameters.max_steps *= factor;
                plugin_.FlightPlanSetAdaptiveStepParameters(vessel_guid,
                                                            parameters);
              }
            }
            using (new UnityEngine.GUILayout.HorizontalScope()) {
              UnityEngine.GUILayout.Label("Tolerance:",
                                          UnityEngine.GUILayout.Width(75));
              if (parameters.length_integration_tolerance <= 1e-6) {
                UnityEngine.GUILayout.Button("min");
              } else if (UnityEngine.GUILayout.Button("-")) {
                parameters.length_integration_tolerance /= 2;
                parameters.speed_integration_tolerance /= 2;
                plugin_.FlightPlanSetAdaptiveStepParameters(vessel_guid,
                                                            parameters);
              }
              UnityEngine.GUILayout.TextArea(
                  parameters.length_integration_tolerance.ToString("0.0e0") +
                      " m",
                  UnityEngine.GUILayout.Width(75));
              if (parameters.length_integration_tolerance >= 1e6) {
                UnityEngine.GUILayout.Button("max");
              } else if (UnityEngine.GUILayout.Button("+")) {
                parameters.length_integration_tolerance *= 2;
                parameters.speed_integration_tolerance *= 2;
                plugin_.FlightPlanSetAdaptiveStepParameters(vessel_guid,
                                                            parameters);
              }
            }
          }

          double Δv = (from burn_editor in burn_editors_
                       select burn_editor.Δv()).Sum();
          UnityEngine.GUILayout.Label(
              "Total Δv : " + Δv.ToString("0.000") + " m/s");

          if (burn_editors_.Count == 0 && 
              UnityEngine.GUILayout.Button("Delete flight plan")) {
            plugin_.FlightPlanDelete(vessel_guid);
            Reset();
          } else {
            if (burn_editors_.Count > 0) {
              RenderUpcomingEvents();
            }
            for (int i = 0; i < burn_editors_.Count - 1; ++i) {
              UnityEngine.GUILayout.TextArea("Manœuvre #" + (i + 1) + ":");
              burn_editors_[i].Render(enabled : false);
            }
            if (burn_editors_.Count > 0) {
              BurnEditor last_burn = burn_editors_.Last();
              UnityEngine.GUILayout.TextArea("Editing manœuvre #" +
                                             (burn_editors_.Count) + ":");
              if (last_burn.Render(enabled : true)) {
                plugin_.FlightPlanReplaceLast(vessel_guid, last_burn.Burn());
                last_burn.Reset(
                    plugin_.FlightPlanGetManoeuvre(vessel_guid,
                                                   burn_editors_.Count - 1));
              }
              if (UnityEngine.GUILayout.Button(
                      "Delete last manœuvre",
                      UnityEngine.GUILayout.ExpandWidth(true))) {
                plugin_.FlightPlanRemoveLast(vessel_guid);
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
                initial_time = plugin_.CurrentTime() + 60;
              } else {
                initial_time =
                    plugin_.FlightPlanGetManoeuvre(
                        vessel_guid,
                        burn_editors_.Count - 1).final_time + 60;
              }
              var editor =
                  new BurnEditor(adapter_, plugin_, vessel_, initial_time);
              Burn candidate_burn = editor.Burn();
              bool inserted = plugin_.FlightPlanAppend(vessel_guid,
                                                       candidate_burn);
              if (inserted) {
                editor.Reset(plugin_.FlightPlanGetManoeuvre(
                    vessel_guid, burn_editors_.Count));
                burn_editors_.Add(editor);
              }
              Shrink();
            }
          }
        }
      }
    }

    UnityEngine.GUI.DragWindow(
        position : new UnityEngine.Rect(x      : 0f,
                                        y      : 0f,
                                        width  : 10000f,
                                        height : 10000f));

    UnityEngine.GUI.skin = old_skin;
  }

  private void RenderUpcomingEvents() {
    string vessel_guid = vessel_.id.ToString();
    double current_time = plugin_.CurrentTime();
    bool should_clear_guidance = true;
    for (int i = 0; i < burn_editors_.Count; ++i) {
      NavigationManoeuvre manoeuvre =
          plugin_.FlightPlanGetManoeuvre(vessel_guid, i);
      if (manoeuvre.final_time > current_time) {
        if (manoeuvre.burn.initial_time > current_time) {
          UnityEngine.GUILayout.TextArea("Upcoming manœuvre: #" + (i + 1));
          UnityEngine.GUILayout.Label(
              "Ignition " + FormatTimeSpan(TimeSpan.FromSeconds(
                                current_time - manoeuvre.burn.initial_time)));
        } else {
          UnityEngine.GUILayout.TextArea("Ongoing manœuvre: #" + (i + 1));
          UnityEngine.GUILayout.Label(
              "Cutoff " + FormatTimeSpan(TimeSpan.FromSeconds(
                              current_time - manoeuvre.final_time)));
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
          XYZ guidance = plugin_.FlightPlanGetGuidance(vessel_guid, i);
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
            should_clear_guidance = false;
          }
          break;
        }
      }
    }
    if (should_clear_guidance && guidance_node_ != null) {
      guidance_node_.RemoveSelf();
      guidance_node_ = null;
    }
  }

  private void Reset() {
    if (burn_editors_ != null) {
      foreach (BurnEditor editor in burn_editors_) {
        editor.Close();
      }
      Shrink();
    }
    burn_editors_ = null;
    vessel_ = FlightGlobals.ActiveVessel;
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

  // Not owned.
  private readonly PrincipiaPluginAdapter adapter_;
  private IntPtr plugin_;
  private Vessel vessel_;
  private List<BurnEditor> burn_editors_;

  private DifferentialSlider final_time_;

  private bool show_guidance_ = false;
  private ManeuverNode guidance_node_;
  
  private const double Log10TimeLowerRate = 0.0;
  private const double Log10TimeUpperRate = 7.0;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
