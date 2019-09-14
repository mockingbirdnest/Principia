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
          $@"Orbit of {vessel.name} with respect to {
            FlightGlobals.Bodies[
                analysis.primary_index].NameWithArticle()} over {
            FlightPlanner.FormatPositiveTimeSpan(
                TimeSpan.FromSeconds(analysis.mission_duration))}");
      var elements = analysis.elements;
      UnityEngine.GUILayout.Label(
          $@"({analysis.mission_duration /
               elements.anomalistic_period:N0} anomalistic rev., {
               analysis.mission_duration / elements.nodal_period:N0} nodal rev., {
               analysis.mission_duration / elements.sidereal_period:N0} sidereal rev.)".ToString(Culture.culture));
      Style.HorizontalLine();
      UnityEngine.GUILayout.Label("Orbital elements");
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
          $"Semimajor axis: {FormatLengthInterval(elements.mean_semimajor_axis)}");
      UnityEngine.GUILayout.Label(
          $"Eccentricity: {FormatInterval(elements.mean_eccentricity)}");
      UnityEngine.GUILayout.Label(
          $"Inclination: {FormatAngleInterval(elements.mean_inclination)}");
      UnityEngine.GUILayout.Label(
          $"Longitude of ascending node: {FormatAngleInterval(elements.mean_longitude_of_ascending_nodes) ?? "(precesses)"}");
      UnityEngine.GUILayout.Label(
          $"Argument of periapsis: {FormatAngleInterval(elements.mean_argument_of_periapsis) ?? "(precesses)"}");

      Style.HorizontalLine();
      OrbitRecurrence recurrence = analysis.recurrence;
      UnityEngine.GUILayout.Label(
          $"Ground track recurrence: [{recurrence.nuo}; {recurrence.dto}; {recurrence.cto}]");
      UnityEngine.GUILayout.Label(
          $"({recurrence.number_of_revolutions} revolutions every {recurrence.cto} days)");
      UnityEngine.GUILayout.Label(
          $"Subcycle: {recurrence.subcycle} days");
      UnityEngine.GUILayout.Label($"Equatorial shift: {recurrence.equatorial_shift * (180 / Math.PI):N2}° ({recurrence.equatorial_shift * primary.Radius / 1000:N0} km)".ToString(Culture.culture));
      UnityEngine.GUILayout.Label($"Base interval: {recurrence.base_interval * (180 / Math.PI):N2}° ({recurrence.base_interval * primary.Radius / 1000:N0} km)".ToString(Culture.culture));
      UnityEngine.GUILayout.Label($"Grid interval: {recurrence.grid_interval * (180 / Math.PI):N2}° ({recurrence.grid_interval * primary.Radius / 1000:N0} km)".ToString(Culture.culture));
      UnityEngine.GUILayout.Label("Longitudes of equatorial crossings:");
      double half_width_km(Interval interval) => (interval.max - interval.min) / 2 * primary.Radius / 1000;
      UnityEngine.GUILayout.Label($"Pass 1 (ascending): {FormatAngleInterval(analysis.ground_track.reduced_longitudes_of_equator_crossings_of_ascending_passes)} (±{half_width_km(analysis.ground_track.reduced_longitudes_of_equator_crossings_of_ascending_passes):N0} km)");
      UnityEngine.GUILayout.Label($"Pass 2 (descending): {FormatAngleInterval(analysis.ground_track.reduced_longitudes_of_equator_crossings_of_descending_passes)} (±{half_width_km(analysis.ground_track.reduced_longitudes_of_equator_crossings_of_descending_passes):N0} km)");
    }
    UnityEngine.GUI.DragWindow();
  }

  // Displays an interval as midpoint±half-width.
  static private string FormatInterval(Interval interval) {
    double half_width = (interval.max - interval.min) / 2;
    double midpoint = interval.min + (interval.max - interval.min) / 2;
    int fractional_digits =
        Math.Max(0, 1 - (int)Math.Floor(Math.Log10(half_width)));
    string format = $"N{fractional_digits}";
    string formatted_midpoint = midpoint.ToString(format, Culture.culture);
    string formatted_half_width = half_width.ToString(format, Culture.culture);
    return $"{formatted_midpoint}±{formatted_half_width}";
  }
  
  // Displays an interval of lengths as midpoint±half-width, in km if the
  // half-width is 100 m or more.
  static private string FormatLengthInterval(Interval interval) {
    double half_width = (interval.max - interval.min) / 2;
    double midpoint = interval.min + (interval.max - interval.min) / 2;
    string unit = "m";
    if (half_width >= 100) {
      half_width *= 0.01;
      midpoint *= 0.01;
      unit = "km";
    }
    int fractional_digits =
        Math.Max(0, 1 - (int)Math.Floor(Math.Log10(half_width)));
    string format = $"N{fractional_digits}";
    string formatted_midpoint = midpoint.ToString(format, Culture.culture);
    string formatted_half_width = half_width.ToString(format, Culture.culture);
    return $"{formatted_midpoint}±{formatted_half_width} {unit}";
  }

  // Displays an interval of angles (given in radians) as midpoint°±half-width°,
  // or null if the half-width exceeds 180°.
  static private string FormatAngleInterval(Interval interval) {
    double half_width = (interval.max - interval.min) / 2;
    double midpoint = interval.min + (interval.max - interval.min) / 2;
    if (half_width > Math.PI) {
      return null;
    }
    double degree = Math.PI / 180;
    int fractional_digits =
        Math.Max(0, 1 - (int)Math.Floor(Math.Log10(half_width / degree)));
    string format = $"N{fractional_digits}";
    string formatted_midpoint =
        (midpoint / degree).ToString(format, Culture.culture);
    string formatted_half_width =
        (half_width / degree).ToString(format, Culture.culture);
    return $"{formatted_midpoint}°±{formatted_half_width}°";
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
