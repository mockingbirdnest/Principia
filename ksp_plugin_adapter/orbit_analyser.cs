using System;
using KSP.Localization;

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

  // Formats |value| (given in metres) in kilometres or metres.
  // TODO(egg): think about the threshold for metres.
  public static string FormatAltitude(this double value) {
    return $"{(value / 1000).FormatN(0)} km";
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
  // half-width is 100 m or more.  The midpoint is given with respect to
  // |offset| (for instance, if |interval| is an interval of distances from the
  // primary and |offset| is the radius of the primary, altitudes are shown).
  public static string FormatLengthInterval(this Interval interval,
                                            double offset = 0) {
    double half_width = (interval.max - interval.min) / 2;
    double midpoint = interval.min - offset + half_width;
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
      return Localizer.Format("#Principia_OrbitAnalyser_Precesses");
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
      return Localizer.Format("#Principia_OrbitAnalyser_Precesses");
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

  public static bool TryParseMissionDuration(string str, out double seconds) {
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

internal abstract class OrbitAnalyser : VesselSupervisedWindowRenderer {
  public OrbitAnalyser(PrincipiaPluginAdapter adapter,
                       PredictedVessel predicted_vessel) : base(
      adapter,
      predicted_vessel,
      UnityEngine.GUILayout.MinWidth(0)) {
    adapter_ = adapter;
  }

  protected abstract void RequestAnalysis();
  protected abstract OrbitAnalysis GetAnalysis();
  protected abstract string ButtonText(string orbit_description);
  protected abstract string AnalysingText();

  // Whether |RequestAnalysis()| needs to be called.
  protected abstract bool should_request_analysis { get; }

  public void RenderButton() {
    string vessel_guid = predicted_vessel?.id.ToString();
    var now = DateTime.UtcNow;
    if (vessel_guid == null) {
      orbit_description_ = null;
    } else {
      if (should_request_analysis &&
          (Shown() ||
           !last_background_analysis_time_.HasValue ||
           now - last_background_analysis_time_ > TimeSpan.FromSeconds(2))) {
        last_background_analysis_time_ = now;
        // Keep refreshing the analysis (albeit at a reduced rate) even when the
        // analyser is not shown, so that the analysis button can display an
        // up-to-date one-line summary.
        RequestAnalysis();
      }
      OrbitAnalysis analysis = GetAnalysis();
      CelestialBody primary = analysis.primary_index.HasValue
                                  ? FlightGlobals.Bodies[
                                      analysis.primary_index.Value]
                                  : null;
      orbit_description_ = OrbitDescription(primary,
                                            analysis.elements,
                                            analysis.recurrence,
                                            analysis.ground_track,
                                            (int?)(analysis.mission_duration /
                                                   analysis.elements?.
                                                       nodal_period));
    }
    RenderButton(ButtonText(orbit_description_));
  }

  protected override string Title => orbit_description_ == null
                                         ? Localizer.Format(
                                             "#Principia_OrbitAnalyser_Title")
                                         : orbit_description_[0].ToString().
                                               ToUpper() +
                                           orbit_description_.Substring(1);

  protected override void RenderWindow(int window_id) {
    string vessel_guid = predicted_vessel?.id.ToString();
    if (vessel_guid == null) {
      return;
    }

    using (new UnityEngine.GUILayout.VerticalScope(GUILayoutWidth(12))) {
      if (should_request_analysis) {
        mission_duration_.Render(enabled : true);
      }
      var multiline_style = Style.Multiline(UnityEngine.GUI.skin.label);
      float two_lines = multiline_style.CalcHeight(
          new UnityEngine.GUIContent("1一\n2二"),
          Width(1));
      float five_lines = multiline_style.CalcHeight(
          new UnityEngine.GUIContent("1一\n2二\n3三\n4四\n5五"),
          Width(1));
      UnityEngine.GUILayout.Label(AnalysingText(),
                                  multiline_style,
                                  UnityEngine.GUILayout.Height(two_lines));

      OrbitAnalysis analysis = GetAnalysis();

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
      CelestialBody primary = analysis.primary_index.HasValue
                                  ? FlightGlobals.Bodies[
                                      analysis.primary_index.Value]
                                  : null;

      orbit_description_ = OrbitDescription(primary,
                                            elements,
                                            recurrence,
                                            ground_track,
                                            (int?)(mission_duration /
                                                   elements?.nodal_period));

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
                                      ? nodal_revolutions /
                                        analysis.recurrence.Value.
                                            number_of_revolutions
                                      : 0;
        string duration_in_ground_track_cycles =
            ground_track_cycles > 0
                ?  Localizer.Format(
                    "#Principia_OrbitAnalyser_Duration_GroundTrackCycles",
                    ground_track_cycles.FormatN(0))
                : "";
        duration_in_revolutions = Localizer.Format(
            "#Principia_OrbitAnalyser_Duration_Revolutions",
            sidereal_revolutions.FormatN(0),
            nodal_revolutions.FormatN(0),
            duration_in_ground_track_cycles,
            anomalistic_revolutions.FormatN(0));
      } else {
        if (primary == null) {
          duration_in_revolutions = Localizer.Format(
              "#Principia_OrbitAnalyser_Warning_NoElements1");
        } else {
          duration_in_revolutions = Localizer.Format(
              "#Principia_OrbitAnalyser_Warning_NoElements2",
              primary.NameWithArticle());
        }
        multiline_style = Style.Warning(multiline_style);
      }
      string analysis_description =
          primary == null
              ? Localizer.Format(
                  "#Principia_OrbitAnalyser_Warning_NoPrimary",
                  predicted_vessel.vesselName,
                  mission_duration.FormatDuration(show_seconds: false))
              : Localizer.Format(
                  "#Principia_OrbitAnalyser_AnalysisDescription",
                  predicted_vessel.vesselName,
                  primary.NameWithArticle(),
                  mission_duration.FormatDuration(show_seconds: false),
                  duration_in_revolutions);

      UnityEngine.GUILayout.Label(analysis_description,
                                  multiline_style,
                                  UnityEngine.GUILayout.Height(five_lines));
      Style.HorizontalLine();
      RenderLowestAltitude(elements, primary);
      Style.HorizontalLine();
      UnityEngine.GUILayout.Label(
          Localizer.Format("#Principia_OrbitAnalyser_Elements_MeanElements"));
      RenderOrbitalElements(elements, primary);
      Style.HorizontalLine();
      UnityEngine.GUILayout.Label(
          Localizer.Format("#Principia_OrbitAnalyser_GroundTrack"));
      RenderOrbitRecurrence(recurrence, primary);
      Style.LineSpacing();
      RenderOrbitGroundTrack(ground_track, primary);
    }
    UnityEngine.GUI.DragWindow();
  }

  public static string OrbitDescription(CelestialBody primary,
                                        OrbitalElements? elements,
                                        OrbitRecurrence? recurrence,
                                        OrbitGroundTrack? ground_track,
                                        int? nodal_revolutions) {
    if (!elements.HasValue) {
      return null;
    }
    var primary_string = L10N.CelestialString(
        "#Principia_OrbitAnalyser_OrbitDescriptionPrimary",
        L10N.NameWithoutArticle,
        primary);
    string properties = "";
    bool circular = false;
    bool equatorial = false;
    if (elements.Value.mean_eccentricity.max < 0.01) {
      circular = true;
      properties +=
          Localizer.Format(
              "#Principia_OrbitAnalyser_OrbitDescription_Circular");
    } else if (elements.Value.mean_eccentricity.min > 0.5) {
      circular = true;
      properties +=
          Localizer.Format(
              "#Principia_OrbitAnalyser_OrbitDescription_HighlyElliptical");
    }
    const double degree = Math.PI / 180;
    if (elements.Value.mean_inclination.max < 5 * degree ||
        elements.Value.mean_inclination.min > 175 * degree) {
      equatorial = true;
      properties +=
          Localizer.Format(
              "#Principia_OrbitAnalyser_OrbitDescription_Equatorial");
    } else if (elements.Value.mean_inclination.min > 80 * degree &&
               elements.Value.mean_inclination.max < 100 * degree) {
      properties +=
          Localizer.Format("#Principia_OrbitAnalyser_OrbitDescription_Polar");
    } else if (elements.Value.mean_inclination.min > 90 * degree) {
      properties +=
          Localizer.Format(
              "#Principia_OrbitAnalyser_OrbitDescription_Retrograde");
    }
    if (recurrence.HasValue && ground_track.HasValue) {
      Interval ascending_longitudes = ground_track.Value.equatorial_crossings.
          longitudes_reduced_to_ascending_pass;
      Interval descending_longitudes = ground_track.Value.equatorial_crossings.
          longitudes_reduced_to_descending_pass;
      double drift = Math.Max(
          ascending_longitudes.max - ascending_longitudes.min,
          descending_longitudes.max - descending_longitudes.min);
      double revolutions_per_day =
          (double)recurrence.Value.number_of_revolutions / recurrence.Value.cto;
      double days = nodal_revolutions.Value / revolutions_per_day;
      // We ignore 0 drift as it means that there was only one pass, which is
      // insufficient to assess synchronicity.
      if (drift > 0 && drift / days < 0.1 * degree) {
        if (recurrence.Value.cto == 1) {
          switch (recurrence.Value.nuo) {
            case 1:
              if (circular && equatorial) {
                var stationary_primary_string = L10N.CelestialStringOrNull(
                    "Principia_OrbitAnalyser_OrbitDescription_StationaryPrimary",
                    L10N.NameWithoutArticle,
                    primary);
                if (stationary_primary_string != null) {
                  primary_string = stationary_primary_string;
                  properties = "";
                } else {
                  properties = Localizer.Format(
                      "#Principia_OrbitAnalyser_OrbitDescription_Stationary");
                }
              } else {
                var synchronous_primary_string = L10N.CelestialStringOrNull(
                    "Principia_OrbitAnalyser_OrbitDescription_SynchronousPrimary",
                    L10N.NameWithoutArticle,
                    primary);
                if (synchronous_primary_string != null) {
                  primary_string = synchronous_primary_string;
                } else {
                  properties +=
                      Localizer.Format(
                          "#Principia_OrbitAnalyser_OrbitDescription_Synchronous");
                }
              }
              break;
            case 2:
              properties += Localizer.Format(
                  "#Principia_OrbitAnalyser_OrbitDescription_Semisynchronous");
              break;
            default:
              break;
          }
        }
      }
    }
    return Localizer.Format("#Principia_OrbitAnalyser_OrbitDescription",
                            properties,
                            primary_string);
  }

  private void RenderLowestAltitude(OrbitalElements? elements,
                                    CelestialBody primary) {
    double? lowest_distance = elements?.radial_distance.min;
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_Elements_LowestAltitude"),
        (lowest_distance - primary?.Radius)?.FormatAltitude());
    double? lowest_primary_distance = primary?.ocean == true
                                          ? primary.Radius
                                          : primary?.pqsController?.radiusMin;
    string altitude_warning =
        lowest_distance < lowest_primary_distance
            ?
            Localizer.Format("#Principia_OrbitAnalyser_Warning_Collision")
            : lowest_distance < primary?.pqsController?.radiusMax
                ? Localizer.Format(
                    "#Principia_OrbitAnalyser_Warning_CollisionRisk")
                : lowest_distance < primary?.Radius + primary?.atmosphereDepth
                    ? Localizer.Format(
                        "#Principia_OrbitAnalyser_Warning_Reentry")
                    : "";
    UnityEngine.GUILayout.Label(altitude_warning,
                                Style.Warning(UnityEngine.GUI.skin.label));
  }

  private void RenderOrbitalElements(OrbitalElements? elements,
                                     CelestialBody primary) {
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_Elements_SiderealPeriod"),
        elements?.sidereal_period.FormatDuration());
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_Elements_NodalPeriod"),
        elements?.nodal_period.FormatDuration());
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_Elements_AnomalisticPeriod"),
        elements?.anomalistic_period.FormatDuration());
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_Elements_SemimajorAxis"),
        elements?.mean_semimajor_axis.FormatLengthInterval());
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_Elements_Eccentricity"),
        elements?.mean_eccentricity.FormatInterval());
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_Elements_Inclination"),
        elements?.mean_inclination.FormatAngleInterval());
    LabeledField(
        Localizer.Format(
            "#Principia_OrbitAnalyser_Elements_LongitudeOfAscendingNode"),
        elements?.mean_longitude_of_ascending_nodes.FormatAngleInterval());
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_Elements_NodalPrecession"),
        elements?.nodal_precession.FormatAngularFrequency());
    string periapsis = L10N.CelestialString(
        "#Principia_OrbitAnalyser_Elements_Periapsis",
        L10N.NameWithArticle,
        primary);
    string apoapsis = L10N.CelestialString(
        "#Principia_OrbitAnalyser_Elements_Apoapsis",
        L10N.NameWithArticle,
        primary);
    LabeledField(
        Localizer.Format(
            "#Principia_OrbitAnalyser_Elements_ArgumentOfPeriapsis",
            periapsis),
        elements?.mean_argument_of_periapsis.FormatAngleInterval());
    LabeledField(
        Localizer.Format(
            "#Principia_OrbitAnalyser_Elements_MeanPeriapsisAltitude",
            periapsis),
        elements?.mean_periapsis_distance.FormatLengthInterval(primary.Radius));
    LabeledField(
        Localizer.Format(
            "#Principia_OrbitAnalyser_Elements_MeanApoapsisAltitude",
            apoapsis),
        elements?.mean_apoapsis_distance.FormatLengthInterval(primary.Radius));
  }

  private void RenderOrbitRecurrence(OrbitRecurrence? recurrence,
                                     CelestialBody primary) {
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      string νₒ = recurrence?.nuo.ToString() ?? em_dash;
      string Dᴛₒ = recurrence?.dto.ToString("+0;-0") ?? em_dash;
      string Cᴛₒ = recurrence?.cto.ToString() ?? em_dash;
      UnityEngine.GUILayout.Label(
          Localizer.Format("#Principia_OrbitAnalyser_Recurrence_CapderouTriple",
                           $"[{νₒ}; {Dᴛₒ}; {Cᴛₒ}]"),
          GUILayoutWidth(8));
      UnityEngine.GUILayout.FlexibleSpace();
      autodetect_recurrence_ = UnityEngine.GUILayout.Toggle(
          autodetect_recurrence_,
          Localizer.Format("#Principia_OrbitAnalyser_Recurrence_AutoDetect"),
          UnityEngine.GUILayout.ExpandWidth(false));
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label(
          Localizer.Format("#Principia_OrbitAnalyser_Recurrence_Cycle"));
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
      UnityEngine.GUILayout.Label(
          Localizer.Format(
              "#Principia_OrbitAnalyser_Recurrence_Cycle_RevolutionsEquals"),
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
      UnityEngine.GUILayout.Label(
          Localizer.Format("#Principia_OrbitAnalyser_Recurrence_Cycle_Days"),
          UnityEngine.GUILayout.ExpandWidth(false));
    }
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_Recurrence_Subcycle"),
        Localizer.Format(
            "#Principia_OrbitAnalyser_Recurrence_SubcycleLengthInDays",
            recurrence?.subcycle.FormatN(0) ?? em_dash));
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_Recurrence_EquatorialShift"),
        recurrence?.equatorial_shift.FormatEquatorialAngle(primary));
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_Recurrence_GridInterval"),
        recurrence?.grid_interval.FormatEquatorialAngle(primary));
  }

  private void RenderOrbitGroundTrack(OrbitGroundTrack? ground_track,
                                      CelestialBody primary) {
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label(
          Localizer.Format(
              "#Principia_OrbitAnalyser_GroundTrack_LongitudesOfEquatorialCrossings_Prefix"),
          UnityEngine.GUILayout.ExpandWidth(false));
      string text = UnityEngine.GUILayout.TextField(
          $"{ground_track_revolution_}",
          GUILayoutWidth(2));
      UnityEngine.GUILayout.Label(
          Localizer.Format(
              "#Principia_OrbitAnalyser_GroundTrack_LongitudesOfEquatorialCrossings_Suffix"));
      if (int.TryParse(text, out int revolution)) {
        ground_track_revolution_ = revolution;
      }
    }
    var equatorial_crossings = ground_track?.equatorial_crossings;
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_GroundTrack_AscendingPass"),
        equatorial_crossings?.longitudes_reduced_to_ascending_pass.
            FormatEquatorialAngleInterval(primary));
    LabeledField(
        Localizer.Format("#Principia_OrbitAnalyser_GroundTrack_DescendingPass"),
        equatorial_crossings?.longitudes_reduced_to_descending_pass.
            FormatEquatorialAngleInterval(primary));
  }

  private void LabeledField(string label, string value) {
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label(label);
      UnityEngine.GUILayout.Label(value ?? em_dash,
                                  Style.RightAligned(
                                      UnityEngine.GUI.skin.label));
    }
  }

  protected IntPtr plugin => adapter_.Plugin();

  private const string em_dash = "—";

  private readonly PrincipiaPluginAdapter adapter_;
  private readonly DifferentialSlider mission_duration_ =
      new DifferentialSlider(
          label            :
              Localizer.Format("#Principia_OrbitAnalyser_MissionDuration"),
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

  protected double requested_mission_duration => mission_duration_.value;

  private bool autodetect_recurrence_ = true;
  private int revolutions_per_cycle_ = 1;
  private int days_per_cycle_ = 1;
  private int ground_track_revolution_ = 1;

  protected int? manual_revolutions_per_cycle =>
      autodetect_recurrence_ ? null : (int?)revolutions_per_cycle_;

  protected int? manual_days_per_cycle =>
      autodetect_recurrence_ ? null : (int?)days_per_cycle_;

  protected int ground_track_revolution => ground_track_revolution_;

  private string orbit_description_ = null;
  private DateTime? last_background_analysis_time_ = null;
}

internal class CurrentOrbitAnalyser : OrbitAnalyser {
  public CurrentOrbitAnalyser(PrincipiaPluginAdapter adapter,
                              PredictedVessel predicted_vessel) : base(
      adapter,
      predicted_vessel) {}

  protected override void RequestAnalysis() {
    plugin.VesselRequestAnalysis(predicted_vessel.id.ToString(),
                                 requested_mission_duration);
  }

  protected override OrbitAnalysis GetAnalysis() {
    return plugin.VesselGetAnalysis(predicted_vessel.id.ToString(),
                                    manual_revolutions_per_cycle,
                                    manual_days_per_cycle,
                                    ground_track_revolution);
  }

  protected override string ButtonText(string orbit_description) {
    return orbit_description == null
               ? Localizer.Format(
                   "#Principia_CurrentOrbitAnalyser_ToggleButton")
               : Localizer.Format(
                   "#Principia_CurrentOrbitAnalyser_ToggleButtonWithDescription",
                   orbit_description);
  }

  protected override string AnalysingText() {
    return Localizer.Format("#Principia_CurrentOrbitAnalyser_Analysing",
                            predicted_vessel.vesselName);
  }

  protected override bool should_request_analysis => true;
}

internal class PlannedOrbitAnalyser : OrbitAnalyser {
  public PlannedOrbitAnalyser(PrincipiaPluginAdapter adapter,
                              PredictedVessel predicted_vessel) : base(
      adapter,
      predicted_vessel) {}

  protected override void RequestAnalysis() {
    Log.Fatal("A PlannedOrbitAnalyser cannot request analysis");
  }

  protected override OrbitAnalysis GetAnalysis() {
    if (!plugin.FlightPlanExists(predicted_vessel.id.ToString())) {
      Hide();
      return new OrbitAnalysis();
    }
    return plugin.FlightPlanGetCoastAnalysis(predicted_vessel.id.ToString(),
                                             manual_revolutions_per_cycle,
                                             manual_days_per_cycle,
                                             ground_track_revolution,
                                             index);
  }

  protected override string ButtonText(string orbit_description) {
    return orbit_description == null
               ? Localizer.Format(
                   "#Principia_PlannedOrbitAnalyser_ToggleButton")
               : Localizer.Format(
                   "#Principia_PlannedOrbitAnalyser_ToggleButtonWithDescription",
                   orbit_description);
  }

  protected override string AnalysingText() {
    return Localizer.Format("#Principia_PlannedOrbitAnalyser_Analysing",
                            predicted_vessel.vesselName);
  }

  protected override bool should_request_analysis => false;

  public int index { private get; set; }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
