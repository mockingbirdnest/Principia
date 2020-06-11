using System;

namespace principia {
namespace ksp_plugin_adapter {

// Miscellaneous formatting utilities; these are extension methods to make null
// handling easier, thanks to null-coalescing operators.
internal static class Formatters {
  // Formats |value| in "N<fractional_digits>" using the Principia culture.
  public static string FormatN(this double value, int fractional_digits) {
    return value.ToString($"N{fractional_digits}", Culture.culture);
  }

  public static string FormatN(this int value, int fractional_digits) {
    return value.ToString($"N{fractional_digits}", Culture.culture);
  }

  // Displays an interval as midpoint±half-width.
  public static string FormatInterval(this Interval interval) {
    double half_width = (interval.max - interval.min) / 2;
    double midpoint = interval.min + half_width;
    int fractional_digits =
        Math.Max(0, 1 - (int)Math.Floor(Math.Log10(half_width)));
    string formatted_midpoint = midpoint.FormatN(fractional_digits);
    string formatted_half_width = half_width.FormatN(fractional_digits);
    return $"{formatted_midpoint}±{formatted_half_width}";
  }

  // Displays an interval of lengths as midpoint±half-width, in km if the
  // half-width is 100 m or more.
  public static string FormatLengthInterval(this Interval interval) {
    double half_width = (interval.max - interval.min) / 2;
    double midpoint = interval.min + half_width;
    string unit = "m";
    if (half_width >= 100) {
      half_width *= 0.001;
      midpoint *= 0.001;
      unit = "km";
    }
    int fractional_digits =
        Math.Max(0, 1 - (int)Math.Floor(Math.Log10(half_width)));
    string formatted_midpoint = midpoint.FormatN(fractional_digits);
    string formatted_half_width = half_width.FormatN(fractional_digits);
    return $"{formatted_midpoint}±{formatted_half_width} {unit}";
  }

  // Displays an interval of angles (given in radians) as midpoint°±half-width°,
  // or the string "(precesses)" if the half-width exceeds 180°.
  public static string FormatAngleInterval(this Interval interval) {
    double half_width = (interval.max - interval.min) / 2;
    double midpoint = interval.min + half_width;
    if (half_width > Math.PI) {
      return "(precesses)";
    }
    const double degree = Math.PI / 180;
    int fractional_digits =
        Math.Max(0, 1 - (int)Math.Floor(Math.Log10(half_width / degree)));
    string formatted_midpoint = (midpoint / degree).FormatN(fractional_digits);
    string formatted_half_width =
        (half_width / degree).FormatN(fractional_digits);
    return $"{formatted_midpoint}°±{formatted_half_width}°";
  }

  // Displays an angle (given in radians) in degrees, with the equivalent in km
  // at the equator given in parentheses.
  public static string FormatEquatorialAngle(this double angle,
                                             CelestialBody primary) {
    const double degree = Math.PI / 180;
    string degrees = (angle / degree).FormatN(1);
    string kilometres = (angle * primary.Radius / 1000).FormatN(1);
    return $"{degrees}° ({kilometres} km)";
  }

  // Similar to |FormatAngleInterval|, but annotated with the equivalent
  // of the half-width as a distance at the equator in parentheses.
  public static string FormatEquatorialAngleInterval(this Interval interval,
                                                     CelestialBody primary) {
    double half_width_angle = (interval.max - interval.min) / 2;
    if (half_width_angle > Math.PI) {
      return "(precesses)";
    }
    double half_width_distance = half_width_angle * primary.Radius;
    string formatted_distance = half_width_distance > 1000
        ? $"{(half_width_distance / 1000).FormatN(1)} km"
        : $"{(half_width_distance).FormatN(0)} m";
    return $"{interval.FormatAngleInterval()} ({formatted_distance})";
  }

  // Formats a duration, omitting leading components if they are 0, and omitting
  // leading 0s on the days; optionally exclude seconds.
  public static string FormatDuration(this double seconds,
                                      bool show_seconds = true) {
    return new PrincipiaTimeSpan(seconds).FormatPositive(
        with_leading_zeroes: false,
        with_seconds: show_seconds);
  }

  // Formats an angular frequency (passed in rad/s), in °/d or °/d6.
  public static string FormatAngularFrequency(this double radians_per_second) {
    const double degree = Math.PI / 180;
    double day = PrincipiaTimeSpan.day_duration;
    string day_unit = PrincipiaTimeSpan.day_symbol;
    double degrees_per_day = radians_per_second / (degree / day);
    return $"{degrees_per_day.FormatN(2)}°/{day_unit}";
  }

  // Never omit leading 0s (to make keyboard editing easier) but do not show
  // seconds (they are irrelevant for a selector that shows durations much
  // longer than a revolution).
  public static string FormatMissionDuration(double seconds) {
    return new PrincipiaTimeSpan(seconds).FormatPositive(
        with_leading_zeroes: true,
        with_seconds: false);
  }

  public static bool TryParseMissionDuration(string str,
                                             out double seconds) {
    seconds = 0;
    if (PrincipiaTimeSpan.TryParse(str,
                                   with_seconds: false,
                                   out PrincipiaTimeSpan ts)) {
      seconds = ts.total_seconds;
      return true;
    } else {
      return false;
    }
  }
}

internal class OrbitAnalyser : VesselSupervisedWindowRenderer {
  public OrbitAnalyser(PrincipiaPluginAdapter adapter,
                       PredictedVessel predicted_vessel)
      : base(adapter, predicted_vessel, UnityEngine.GUILayout.MinWidth(0)) {
    adapter_ = adapter;
  }

  public void RenderButton() {
    RenderButton("Orbit analysis...");
  }

  protected override string Title => "Orbit analysis";

  protected override void RenderWindow(int window_id) {
    string vessel_guid = predicted_vessel?.id.ToString();
    if (vessel_guid == null) {
      return;
    }

    using (new UnityEngine.GUILayout.VerticalScope(GUILayoutWidth(8))) {
      mission_duration_.Render(enabled : true);
      var multiline_style = Style.Multiline(UnityEngine.GUI.skin.label);
      float two_lines = multiline_style.CalcHeight(
          new UnityEngine.GUIContent("1\n2"), Width(1));
      float five_lines = multiline_style.CalcHeight(
          new UnityEngine.GUIContent("1\n2\n3\n4\n5"), Width(1));
      UnityEngine.GUILayout.Label(
          $@"Analysing orbit of {predicted_vessel.vesselName}...",
          multiline_style,
          UnityEngine.GUILayout.Height(two_lines));

      OrbitAnalysis analysis = plugin.VesselRefreshAnalysis(
          predicted_vessel.id.ToString(),
          mission_duration_.value,
          autodetect_recurrence_ ? null : (int?)revolutions_per_cycle_,
          autodetect_recurrence_ ? null : (int?)days_per_cycle_,
          ground_track_revolution_);
      if (autodetect_recurrence_ &&
          analysis.recurrence.HasValue &&
          analysis.recurrence.Value.number_of_revolutions != 0 &&
          analysis.recurrence.Value.cto != 0) {
        revolutions_per_cycle_ =
            analysis.recurrence.Value.number_of_revolutions;
        days_per_cycle_ = analysis.recurrence.Value.cto;
      }
      UnityEngine.GUILayout.HorizontalScrollbar(
          value      : 0,
          size       : (float)analysis.progress_of_next_analysis,
          leftValue  : 0,
          rightValue : 1);

      OrbitalElements? elements = analysis.elements;
      OrbitRecurrence? recurrence = analysis.recurrence;
      OrbitGroundTrack? ground_track = analysis.ground_track;
      double mission_duration = analysis.mission_duration;
      CelestialBody primary =
          analysis.primary_index.HasValue ? FlightGlobals.Bodies[analysis.primary_index.Value]
                                          : null;

      Style.HorizontalLine();
      string duration_in_revolutions;
      if (elements.HasValue) {
        int sidereal_revolutions =
            (int)(mission_duration / elements.Value.sidereal_period);
        int nodal_revolutions =
            (int)(mission_duration / elements.Value.nodal_period);
        int anomalistic_revolutions =
            (int)(mission_duration / elements.Value.anomalistic_period);
        int ground_track_cycles = analysis.recurrence.HasValue
                                      ? nodal_revolutions / analysis.recurrence.
                                            Value.number_of_revolutions
                                      : 0;
        string duration_in_ground_track_cycles = ground_track_cycles > 0
            ? $" ({ground_track_cycles.FormatN(0)} ground track cycles)"
            : "";
        duration_in_revolutions = $@"{
            sidereal_revolutions.FormatN(0)} sidereal revolutions{"\n"}{
            nodal_revolutions.FormatN(0)} nodal revolutions{
            duration_in_ground_track_cycles}{"\n"}{
            anomalistic_revolutions.FormatN(0)} anomalistic revolutions";
      } else {
        duration_in_revolutions =
            "could not determine elements; mission duration may be shorter " +
            "than a revolution, or trajectory may not be gravitationally bound";
        if (primary != null) {
          duration_in_revolutions += $" to {primary.NameWithArticle()}.";
        }
        multiline_style = Style.Warning(multiline_style);
      }
      string analysis_description =
      primary == null 
          ? $@"{predicted_vessel.vesselName} is not gravitationally bound over {
               mission_duration.FormatDuration(show_seconds: false)}"
          : $@"Orbit of {predicted_vessel.vesselName} with respect to {
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
      RenderOrbitRecurrence(recurrence, primary);
      Style.HorizontalLine();
      RenderOrbitGroundTrack(ground_track, primary);
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

  private void RenderOrbitRecurrence(OrbitRecurrence? recurrence,
                                     CelestialBody primary) {
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label(
          $@"Recurrence: [{
            recurrence?.nuo.ToString() ?? em_dash}; {
            recurrence?.dto.ToString("+0;-0") ?? em_dash}; {
            recurrence?.cto.ToString() ?? em_dash}]",
          GUILayoutWidth(8));
      UnityEngine.GUILayout.FlexibleSpace();
      autodetect_recurrence_ = UnityEngine.GUILayout.Toggle(
          autodetect_recurrence_,
          "Auto-detect",
          UnityEngine.GUILayout.ExpandWidth(false));
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label("Cycle");
      string text = UnityEngine.GUILayout.TextField(
          recurrence.HasValue || !autodetect_recurrence_
              ? $"{revolutions_per_cycle_}"
              : em_dash,
          Style.RightAligned(UnityEngine.GUI.skin.textField),
          GUILayoutWidth(2));
      if (!autodetect_recurrence_ &&
          int.TryParse(text, out int revolutions_per_cycle) &&
          revolutions_per_cycle > 0) {
        revolutions_per_cycle_ = revolutions_per_cycle;
      }
      UnityEngine.GUILayout.Label("revolutions =",
                                  UnityEngine.GUILayout.ExpandWidth(false));
      text = UnityEngine.GUILayout.TextField(
          recurrence.HasValue || !autodetect_recurrence_
              ? $"{days_per_cycle_}"
              : em_dash,
          Style.RightAligned(UnityEngine.GUI.skin.textField),
          GUILayoutWidth(2));
      if (!autodetect_recurrence_ &&
          int.TryParse(text, out int days_per_cycle) &&
          days_per_cycle != 0) {
        days_per_cycle_ = days_per_cycle;
      }
      UnityEngine.GUILayout.Label("days",
                                  UnityEngine.GUILayout.ExpandWidth(false));
    }
    LabeledField(
        "Subcycle",
        $@"{recurrence?.subcycle.FormatN(0) ?? em_dash} days");
    LabeledField(
        "Equatorial shift",
        recurrence?.equatorial_shift.FormatEquatorialAngle(primary));
    LabeledField(
        "Grid interval",
        recurrence?.grid_interval.FormatEquatorialAngle(primary));
  }

  private void RenderOrbitGroundTrack(OrbitGroundTrack? ground_track,
                                      CelestialBody primary) {
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label(
          "Longitudes of equatorial crossings of rev. #",
          UnityEngine.GUILayout.ExpandWidth(false));
      string text = UnityEngine.GUILayout.TextField(
          $"{ground_track_revolution_}",
          GUILayoutWidth(2));
      if (int.TryParse(text, out int revolution)) {
        ground_track_revolution_ = revolution;
      }
    }
    var equatorial_crossings = ground_track?.equatorial_crossings;
    LabeledField(
        "Ascending pass",
        equatorial_crossings?.longitudes_reduced_to_ascending_pass.
            FormatEquatorialAngleInterval(primary));
    LabeledField(
        "Descending pass",
        equatorial_crossings?.longitudes_reduced_to_descending_pass.
            FormatEquatorialAngleInterval(primary));
  }

  private void LabeledField(string label, string value) {
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label(label);
      UnityEngine.GUILayout.Label(value ?? em_dash,
          Style.RightAligned(UnityEngine.GUI.skin.label));
    }
  }

  private IntPtr plugin => adapter_.Plugin();

  private const string em_dash = "—";

  private readonly PrincipiaPluginAdapter adapter_;
  private readonly DifferentialSlider mission_duration_ =
      new DifferentialSlider(
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
