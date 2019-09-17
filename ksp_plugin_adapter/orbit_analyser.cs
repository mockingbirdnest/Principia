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
          mission_duration_.value,
          autodetect_recurrence_ ? null : (int?)revolutions_per_cycle_,
          autodetect_recurrence_ ? null : (int?)days_per_cycle_,
          ground_track_revolution_);
      if (autodetect_recurrence_ &&
          analysis.recurrence.number_of_revolutions != 0 &&
          analysis.recurrence.cto != 0) {
        revolutions_per_cycle_ = analysis.recurrence.number_of_revolutions;
        days_per_cycle_ = analysis.recurrence.cto;
      }

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
            FormatDuration(analysis.mission_duration)}");
      var elements = analysis.elements;
      UnityEngine.GUILayout.Label(
          analysis.elements_has_value
              ? $@"({analysis.mission_duration /
                     elements.anomalistic_period:N0} anomalistic rev., {
                     analysis.mission_duration / elements.nodal_period:N0} nodal rev., {
                     analysis.mission_duration / elements.sidereal_period:N0} sidereal rev.)".ToString(Culture.culture)
              : "(mission duration is less than one sidereal revolution)");
      Style.HorizontalLine();
      UnityEngine.GUILayout.Label("Orbital elements");
      string element_string(string s) => analysis.elements_has_value ? s : "—";
      UnityEngine.GUILayout.Label(
          "Sidereal period: " +
          element_string(
              FormatDuration(elements.sidereal_period)));
      UnityEngine.GUILayout.Label(
          "Nodal period: " +
          element_string(
              FormatDuration(elements.nodal_period)));
      UnityEngine.GUILayout.Label(
          "Anomalistic period: " +
          element_string(
              FormatDuration(analysis.elements.anomalistic_period)));
      // TODO(egg): represent the intervals.
      UnityEngine.GUILayout.Label(
          $"Semimajor axis: {element_string(FormatLengthInterval(elements.mean_semimajor_axis))}");
      UnityEngine.GUILayout.Label(
          $"Eccentricity: {element_string(FormatInterval(elements.mean_eccentricity))}");
      UnityEngine.GUILayout.Label(
          $"Inclination: {element_string(FormatAngleInterval(elements.mean_inclination))}");
      UnityEngine.GUILayout.Label(
          $"Longitude of ascending node: {element_string(FormatAngleInterval(elements.mean_longitude_of_ascending_nodes) ?? "(precesses)")}");
      UnityEngine.GUILayout.Label(
          $"Argument of periapsis: {element_string(FormatAngleInterval(elements.mean_argument_of_periapsis) ?? "(precesses)")}");

      Style.HorizontalLine();
      OrbitRecurrence recurrence = analysis.recurrence;
      string recurrence_string<T>(T s) => analysis.recurrence_has_value ? s.ToString() : "—";
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        UnityEngine.GUILayout.Label(
            $"Ground track recurrence: [{recurrence_string(recurrence.nuo)}; {recurrence_string(recurrence.dto)}; {recurrence_string(recurrence.cto)}]");
        autodetect_recurrence_ = UnityEngine.GUILayout.Toggle(autodetect_recurrence_, "Auto-detect");
      }
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        string text = UnityEngine.GUILayout.TextField(
            $"{revolutions_per_cycle_}");
        if (int.TryParse(text, out int revolutions_per_cycle) &&
            revolutions_per_cycle > 0) {
          revolutions_per_cycle_ = revolutions_per_cycle;
        }
        UnityEngine.GUILayout.Label("revolutions =");
        text = UnityEngine.GUILayout.TextField($"{days_per_cycle_}");
        if (int.TryParse(text, out int days_per_cycle) &&
            days_per_cycle != 0) {
          days_per_cycle_ = days_per_cycle;
        }
        UnityEngine.GUILayout.Label("days");

      }
      UnityEngine.GUILayout.Label(
          $"Subcycle: {recurrence_string(recurrence.subcycle)} days");
      UnityEngine.GUILayout.Label($"Equatorial shift: {recurrence.equatorial_shift * (180 / Math.PI):N2}° ({recurrence.equatorial_shift * primary.Radius / 1000:N0} km)".ToString(Culture.culture));
      UnityEngine.GUILayout.Label($"Base interval: {recurrence.base_interval * (180 / Math.PI):N2}° ({recurrence.base_interval * primary.Radius / 1000:N0} km)".ToString(Culture.culture));
      UnityEngine.GUILayout.Label($"Grid interval: {recurrence.grid_interval * (180 / Math.PI):N2}° ({recurrence.grid_interval * primary.Radius / 1000:N0} km)".ToString(Culture.culture));
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        UnityEngine.GUILayout.Label("Longitudes of equatorial crossings of revolution #");
        string text = UnityEngine.GUILayout.TextField($"{ground_track_revolution_}");
        if (int.TryParse(text, out int revolution)) {
          ground_track_revolution_ = revolution;
        }
      }
      double half_width_km(Interval interval) => (interval.max - interval.min) / 2 * primary.Radius / 1000;
      var equatorial_crossings = analysis.ground_track.equatorial_crossings;
      UnityEngine.GUILayout.Label($"Ascending pass: {FormatAngleInterval(equatorial_crossings.longitudes_reduced_to_ascending_pass)} (±{half_width_km(equatorial_crossings.longitudes_reduced_to_ascending_pass):N0} km)");
      UnityEngine.GUILayout.Label($"Descending pass: {FormatAngleInterval(equatorial_crossings.longitudes_reduced_to_descending_pass)} (±{half_width_km(equatorial_crossings.longitudes_reduced_to_descending_pass):N0} km)");
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

  // Formats a duration, omitting leading components if they are 0, and omitting
  // leading 0s on the days;
  private static string FormatDuration(double seconds) {
    var span = TimeSpan.FromSeconds(seconds);
    int days = GameSettings.KERBIN_TIME ? span.Days * 4 + span.Hours / 6
                                        : span.Days;
    int hours = GameSettings.KERBIN_TIME ? span.Hours % 6
                                         : span.Hours;
    var components = new List<string>();
    if (days > 0) {
      components.Add(GameSettings.KERBIN_TIME ? $"{days} d6" : $"{days} d");
    }
    if (components.Count > 0 || hours > 0) {
      components.Add(
          GameSettings.KERBIN_TIME ? $"{hours:0} h" : $"{hours:00} h");
    }
    if (components.Count > 0 || span.Minutes > 0) {
      components.Add($"{span.Minutes:00} min");
    }
    components.Add($"{span.Seconds + span.Milliseconds / 1000m:00.0} s");
    return string.Join(" ", components.ToArray());
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
  private bool autodetect_recurrence_ = true;
  private int revolutions_per_cycle_ = 1;
  private int days_per_cycle_ = 1;
  private int ground_track_revolution_ = 1;
}


}  // namespace ksp_plugin_adapter
}  // namespace principia
