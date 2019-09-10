using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace principia {
namespace ksp_plugin_adapter {

internal class OrbitAnalyser : SupervisedWindowRenderer {
  public OrbitAnalyser(PrincipiaPluginAdapter adapter)
      : base(adapter, UnityEngine.GUILayout.MinWidth(0)) {
    adapter_ = adapter;
  }

  public void RenderButton() {
    if (UnityEngine.GUILayout.Button("Orbit analysis...")) {
      Toggle();
    }
  }

  protected override string Title => "Orbit analysis";
  protected override void RenderWindow(int window_id) {
    using (new UnityEngine.GUILayout.VerticalScope(GUILayoutWidth(8))) {
      Vessel vessel = FlightGlobals.ActiveVessel;
      if (plugin == IntPtr.Zero || vessel == null ||
          !plugin.HasVessel(vessel.id.ToString())) {
        UnityEngine.GUILayout.Label("No active vessel");
        return;
      }
      CelestialBody primary = adapter_.plotting_frame_selector_.selected_celestial;
      OrbitAnalysis analysis = plugin.VesselRefreshAnalysis(
          vessel.id.ToString(),
          primary.flightGlobalsIndex,
          mission_duration_.value);

      mission_duration_.Render(enabled : true);
      UnityEngine.GUILayout.Label(
          $"Analysing orbit of {vessel.name} with respect to {primary.NameWithArticle()}...");
      UnityEngine.GUILayout.HorizontalScrollbar(
          value      : analysis.progress_percentage * 0.5f,
          size       : analysis.progress_percentage,
          leftValue  : 0,
          rightValue : 100);

      Style.HorizontalLine();
      UnityEngine.GUILayout.Label(
          $"Orbit of {vessel.name} with respect to {FlightGlobals.Bodies[analysis.primary_index].NameWithArticle()}");
      Style.HorizontalLine();
      UnityEngine.GUILayout.Label("Orbital elements");
      UnityEngine.GUILayout.Label($"Sidereal period: {analysis.elements.sidereal_period / 60} min");
      UnityEngine.GUILayout.Label($"Nodal period: {analysis.elements.nodal_period / 60} min");
      UnityEngine.GUILayout.Label($"Anomalistic period: {analysis.elements.anomalistic_period / 60} min");
      // TODO(egg): represent the intervals.
      UnityEngine.GUILayout.Label($"Semimajor axis: {analysis.elements.mean_semimajor_axis.min} m");
      UnityEngine.GUILayout.Label($"Eccentricity: {analysis.elements.mean_eccentricity.min}");
      UnityEngine.GUILayout.Label($"Inclination: {analysis.elements.mean_inclination.min}");
      UnityEngine.GUILayout.Label($"Longitude of ascending node: {analysis.elements.mean_longitude_of_ascending_nodes.min}");
      UnityEngine.GUILayout.Label($"Argument of periapsis: {analysis.elements.mean_argument_of_periapsis.min}");

      Style.HorizontalLine();

    }
  }

  static internal bool TryParseMissionDuration(string str, out double value) {
    value = 0;
    if (!FlightPlanner.TryParseTimeSpan(str, out TimeSpan ts)) {
      return false;
    }
    value = ts.TotalSeconds;
    return true;
  }

  private IntPtr plugin => adapter_.Plugin();

  private readonly PrincipiaPluginAdapter adapter_;
  private DifferentialSlider mission_duration_ = new DifferentialSlider(
                      label            : "Plan length",
                      unit             : null,
                      log10_lower_rate : 0,
                      log10_upper_rate : 7,
                      min_value        : 10,
                      max_value        : double.PositiveInfinity,
                      formatter        : value =>
                          FlightPlanner.FormatPositiveTimeSpan(
                              TimeSpan.FromSeconds(value)),
                      parser           : TryParseMissionDuration);
}


}  // namespace ksp_plugin_adapter
}  // namespace principia
