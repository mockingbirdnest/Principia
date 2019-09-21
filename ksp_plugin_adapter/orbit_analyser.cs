using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace principia {
namespace ksp_plugin_adapter {

internal static class Formatters {
  // Displays an interval as midpoint±half-width.
  public static string FormatInterval(this Interval interval) {
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
  public static string FormatLengthInterval(this Interval interval) {
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
  // or the string "(precesses)" if the half-width exceeds 180°.
  public static string FormatAngleInterval(this Interval interval) {
    double half_width = (interval.max - interval.min) / 2;
    double midpoint = interval.min + (interval.max - interval.min) / 2;
    if (half_width > Math.PI) {
      return "(precesses)";
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
  // leading 0s on the days; optionally exclude seconds.
  public static string FormatDuration(this double seconds,
                                      bool show_seconds = true) {
    if (double.IsNaN(seconds)) {
      return seconds.ToString();
    }
    var span = TimeSpan.FromSeconds(seconds);
    int days = GameSettings.KERBIN_TIME ? span.Days * 4 + span.Hours / 6
                                        : span.Days;
    int hours = GameSettings.KERBIN_TIME ? span.Hours % 6
                                         : span.Hours;
    var components = new List<string>();
    const string nbsp = "\xA0";
    if (days > 0) {
      components.Add(GameSettings.KERBIN_TIME ? $"{days}{nbsp}d6"
                                              : $"{days}{nbsp}d");
    }
    if (components.Count > 0 || hours > 0) {
      components.Add(
          GameSettings.KERBIN_TIME ? $"{hours:0}{nbsp}h"
                                   : $"{hours:00}{nbsp}h");
    }
    if (components.Count > 0 || span.Minutes > 0 || !show_seconds) {
      components.Add($"{span.Minutes:00}{nbsp}min");
    }
    if (show_seconds) {
      components.Add($"{span.Seconds + span.Milliseconds / 1000m:00.0}{nbsp}s");
    }
    return string.Join(" ", components.ToArray());
  }

  public static string FormatAngularFrequency(this double radians_per_second) {
    double degree = Math.PI / 180;
    double day = GameSettings.KERBIN_TIME ? 6 * 60 * 60 : 24 * 60 * 60;
    string day_unit = GameSettings.KERBIN_TIME ? "d6" : "d";
    double degrees_per_day = radians_per_second / (degree / day);
    return $"{degrees_per_day:N2}°/{day_unit}";
  }

  // Never omit leading 0s (to make keyboard editing easier) but do not show
  // seconds (they are irrelevant for a selector that shows durations much
  // longer than a revolution).
  public static string FormatMissionDuration(double seconds) {
    var span = TimeSpan.FromSeconds(seconds);
    return (GameSettings.KERBIN_TIME
                ? (span.Days * 4 + span.Hours / 6).ToString("0000") +
                      " d6 " + (span.Hours % 6).ToString("0") + " h "
                : span.Days.ToString("000") + " d " +
                      span.Hours.ToString("00") + " h ") +
           span.Minutes.ToString("00") + " min";
  }

  public static bool TryParseMissionDuration(string str, out double value) {
    value = 0;
    // Using a technology that is customarily used to parse HTML.
    string pattern = @"^[+]?\s*(\d+)\s*" +
                     (GameSettings.KERBIN_TIME ? "d6" : "d") +
                     @"\s*(\d+)\s*h\s*(\d+)\s*min$";
    var regex = new Regex(pattern);
    var match = regex.Match(str);
    if (!match.Success) {
      return false;
    }
    string days = match.Groups[1].Value;
    string hours = match.Groups[2].Value;
    string minutes = match.Groups[3].Value;
    if (!int.TryParse(days, out int d) ||
        !int.TryParse(hours, out int h) ||
        !int.TryParse(minutes, out int min)) {
      return false;
    }
    value = (TimeSpan.FromDays((double)d / (GameSettings.KERBIN_TIME ? 4 : 1)) +
             TimeSpan.FromHours(h) +
             TimeSpan.FromMinutes(min)).TotalSeconds;
    return true;
  }
}

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

      mission_duration_.Render(enabled : true);
      var multiline_style = Style.Multiline(UnityEngine.GUI.skin.label);
      float two_lines = multiline_style.CalcHeight(
          new UnityEngine.GUIContent("1\n2"), Width(1));
      float five_lines = multiline_style.CalcHeight(
          new UnityEngine.GUIContent("1\n2\n3\n4\n5"), Width(1));
      UnityEngine.GUILayout.Label(
          $@"Analysing orbit of {vessel.vesselName} with respect to {
            primary.NameWithArticle()}...",
          multiline_style,
          UnityEngine.GUILayout.Height(two_lines));

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
      UnityEngine.GUILayout.HorizontalScrollbar(
          value      : 0,
          size       : analysis.progress_percentage,
          leftValue  : 0,
          rightValue : 100);

      OrbitalElements? elements = null;
      double mission_duration = analysis.mission_duration;
      if (analysis.elements_has_value) {
        elements = analysis.elements;
      }
      primary = FlightGlobals.Bodies[analysis.primary_index];

      Style.HorizontalLine();
      string duration_in_revolutions;
      if (elements.HasValue) {
        int sidereal_revolutions =
            (int)(mission_duration / elements.Value.sidereal_period);
        int nodal_revolutions =
            (int)(mission_duration / elements.Value.nodal_period);
        int anomalistic_revolutions =
            (int)(mission_duration / elements.Value.anomalistic_period);
        int ground_track_cycles = analysis.recurrence_has_value
          ? nodal_revolutions / analysis.recurrence.number_of_revolutions
          : 0;
        string duration_in_ground_track_cycles = ground_track_cycles > 0
            ? $" ({ground_track_cycles:N0} ground track cycles)"
            : "";
        duration_in_revolutions = $@"{
            sidereal_revolutions:N0} sidereal revolutions{"\n"}{
            nodal_revolutions:N0} nodal revolutions{
            duration_in_ground_track_cycles}{"\n"}{
            anomalistic_revolutions:N0} anomalistic revolutions".ToString(
                Culture.culture);
      } else {
        duration_in_revolutions =
            "mission duration is shorter than one sidereal revolution";
      }
      string analysis_description =
          $@"Orbit of {vessel.vesselName} with respect to {
            primary.NameWithArticle()} over {
            mission_duration.FormatDuration(show_seconds : false)}:{"\n"}{
            duration_in_revolutions}";
      UnityEngine.GUILayout.Label(
          analysis_description,
          multiline_style,
          UnityEngine.GUILayout.Height(five_lines));
      Style.HorizontalLine();
      RenderOrbitalElements(elements);
      Style.HorizontalLine();
      OrbitRecurrence recurrence = analysis.recurrence;
      string recurrence_string<T>(T s) => analysis.recurrence_has_value ? s.ToString() : "—";
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        UnityEngine.GUILayout.Label(
            $"Ground track recurrence: [{recurrence_string(recurrence.nuo)}; {recurrence_string(recurrence.dto)}; {recurrence_string(recurrence.cto)}]");
        autodetect_recurrence_ = UnityEngine.GUILayout.Toggle(autodetect_recurrence_, "Auto-detect", UnityEngine.GUILayout.ExpandWidth(false));
      }
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        UnityEngine.GUILayout.Label("Cycle");
        string text = UnityEngine.GUILayout.TextField(
            $"{revolutions_per_cycle_}",
            Style.RightAligned(UnityEngine.GUI.skin.textField),
            GUILayoutWidth(2));
        if (int.TryParse(text, out int revolutions_per_cycle) &&
            revolutions_per_cycle > 0) {
          revolutions_per_cycle_ = revolutions_per_cycle;
        }
        UnityEngine.GUILayout.Label("revolutions =",
                                    UnityEngine.GUILayout.ExpandWidth(false));
        text = UnityEngine.GUILayout.TextField(
            $"{days_per_cycle_}",
            Style.RightAligned(UnityEngine.GUI.skin.textField),
            GUILayoutWidth(2));
        if (int.TryParse(text, out int days_per_cycle) &&
            days_per_cycle != 0) {
          days_per_cycle_ = days_per_cycle;
        }
        UnityEngine.GUILayout.Label("days", UnityEngine.GUILayout.ExpandWidth(false));
      }
      LabeledField(
          "Subcycle",
          $"{recurrence_string(recurrence.subcycle)} days".ToString(Culture.culture));
      LabeledField("Equatorial shift", $"{recurrence.equatorial_shift * (180 / Math.PI):N2}° ({recurrence.equatorial_shift * primary.Radius / 1000:N0} km)".ToString(Culture.culture));
      LabeledField("Grid interval", $"{recurrence.grid_interval * (180 / Math.PI):N2}° ({recurrence.grid_interval * primary.Radius / 1000:N0} km)".ToString(Culture.culture));
      Style.HorizontalLine();
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        UnityEngine.GUILayout.Label("Longitudes of equatorial crossings of rev. #", UnityEngine.GUILayout.ExpandWidth(false));
        string text = UnityEngine.GUILayout.TextField(
            $"{ground_track_revolution_}",
            GUILayoutWidth(2));
        if (int.TryParse(text, out int revolution)) {
          ground_track_revolution_ = revolution;
        }
      }
      double half_width_km(Interval interval) => (interval.max - interval.min) / 2 * primary.Radius / 1000;
      var equatorial_crossings = analysis.ground_track.equatorial_crossings;
      LabeledField("Ascending pass", $"{equatorial_crossings.longitudes_reduced_to_ascending_pass.FormatAngleInterval()} (±{half_width_km(equatorial_crossings.longitudes_reduced_to_ascending_pass):N0} km)");
      LabeledField("Descending pass", $"{equatorial_crossings.longitudes_reduced_to_descending_pass.FormatAngleInterval()} (±{half_width_km(equatorial_crossings.longitudes_reduced_to_descending_pass):N0} km)");
    }
    UnityEngine.GUI.DragWindow();
  }

  private void RenderOrbitalElements(OrbitalElements? elements) {
      UnityEngine.GUILayout.Label("Orbital elements");
      LabeledField("Sidereal period",
                   elements?.sidereal_period.FormatDuration());
      LabeledField("Nodal period",
                   elements?.nodal_period.FormatDuration());
      LabeledField("Anomalistic period",
                   elements?.anomalistic_period.FormatDuration());
      LabeledField("Semimajor axis",
                   elements?.mean_semimajor_axis.FormatLengthInterval());
      LabeledField("Eccentricity",
                   elements?.mean_eccentricity.FormatInterval());
      LabeledField("Inclination",
                   elements?.mean_inclination.FormatAngleInterval());
      LabeledField(
            "Longitude of ascending node",
            elements?.mean_longitude_of_ascending_nodes.FormatAngleInterval());
      LabeledField(
            "Nodal precession",
            elements?.nodal_precession.FormatAngularFrequency());
      LabeledField(
            "Argument of periapsis",
            elements?.mean_argument_of_periapsis.FormatAngleInterval());
  }

  private void LabeledField(
      string label,
      string value) {
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label(label);
      UnityEngine.GUILayout.Label(value ?? em_dash,
          Style.RightAligned(UnityEngine.GUI.skin.label));
    }
  }

  private const string em_dash = "—";
  private const string nbsp = "\xA0";

  private IntPtr plugin => adapter_.Plugin();

  private readonly PrincipiaPluginAdapter adapter_;
  private DifferentialSlider mission_duration_ = new DifferentialSlider(
      label            : "Duration",
      unit             : null,
      log10_lower_rate : 0,
      log10_upper_rate : 7,
      min_value        : 10,
      max_value        : double.PositiveInfinity,
      formatter        : Formatters.FormatMissionDuration,
      parser           : Formatters.TryParseMissionDuration,
      label_width      : 2,
      field_width      : 5) {
      value = 7 * 24 * 60 * 60
  };
  private bool autodetect_recurrence_ = true;
  private int revolutions_per_cycle_ = 1;
  private int days_per_cycle_ = 1;
  private int ground_track_revolution_ = 1;
}


}  // namespace ksp_plugin_adapter
}  // namespace principia
