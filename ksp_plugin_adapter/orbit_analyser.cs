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
          value      : 0,
          size       : analysis.progress_percentage,
          leftValue  : 0,
          rightValue : 100);

      Style.HorizontalLine();
      UnityEngine.GUILayout.Label(
          $"Orbit of {vessel.name} with respect to {FlightGlobals.Bodies[analysis.primary_index].NameWithArticle()}");
      Style.HorizontalLine();
      UnityEngine.GUILayout.Label("Orbital elements");
      var elements = analysis.elements;
      UnityEngine.GUILayout.Label(
          "Sidereal period: "+
          FlightPlanner.FormatPositiveTimeSpan(
              TimeSpan.FromSeconds(elements.sidereal_period)));
      UnityEngine.GUILayout.Label(
          "Nodal period: " +
          FlightPlanner.FormatPositiveTimeSpan(
              TimeSpan.FromSeconds(elements.nodal_period)));
      UnityEngine.GUILayout.Label(
          "Anomalistic period: " +
          FlightPlanner.FormatPositiveTimeSpan(
              TimeSpan.FromSeconds(analysis.elements.anomalistic_period)));
      // TODO(egg): represent the intervals.
      UnityEngine.GUILayout.Label(
          $"Semimajor axis: {FormatInterval(elements.mean_semimajor_axis)} m");
      UnityEngine.GUILayout.Label(
          $"Eccentricity: {FormatInterval(elements.mean_eccentricity)}");
      UnityEngine.GUILayout.Label(
          $"Inclination: {FormatInterval(elements.mean_inclination, Math.PI / 180, "°")}");
      UnityEngine.GUILayout.Label(
          $"Longitude of ascending node: {FormatInterval(elements.mean_longitude_of_ascending_nodes, Math.PI / 180, "°")}");
      UnityEngine.GUILayout.Label(
          $"Argument of periapsis: {FormatInterval(elements.mean_argument_of_periapsis, Math.PI / 180, "°")}");

      Style.HorizontalLine();

    }
    UnityEngine.GUI.DragWindow();
  }

  // Displays an interval as midpoint±half-width.
  static private string FormatInterval(Interval interval,
                                       double unit = 1,
                                       string unit_symbol = "") {
    double half_width = (interval.max - interval.min) / 2;
    double midpoint = interval.min + (interval.max - interval.min) / 2;
    int fractional_digits =
        Math.Max(0, 1 - (int)Math.Floor(Math.Log10(half_width / unit)));
    string format = $"N{fractional_digits}";
    string formatted_midpoint =
        (midpoint / unit).ToString(format, Culture.culture) + unit_symbol;
    string formatted_half_width =
        (half_width / unit).ToString(format, Culture.culture) + unit_symbol;
    return $"{formatted_midpoint}±{formatted_half_width}";
  }

  static private bool TryParseMissionDuration(string str, out double value) {
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
      label            : "Mission duration",
      unit             : null,
      log10_lower_rate : 0,
      log10_upper_rate : 7,
      min_value        : 10,
      max_value        : double.PositiveInfinity,
      formatter        : value =>
          FlightPlanner.FormatPositiveTimeSpan(
              TimeSpan.FromSeconds(value)),
      parser           : TryParseMissionDuration) {
      value = 30 * 60
  };
}


}  // namespace ksp_plugin_adapter
}  // namespace principia
