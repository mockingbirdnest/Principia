using System;
using System.Collections.Generic;

namespace principia {
namespace ksp_plugin_adapter {

// Miscellaneous formatting utilities; these are extension methods to make null
// handling easier, thanks to null-coalescing operators.
internal static class Formatters {
  // Formats `value` in "N<fractional_digits>" using the Principia culture.
  public static string FormatN(this double value, int fractional_digits) {
    return value.ToString($"N{fractional_digits}", Culture.culture);
  }

  public static string FormatN(this int value, int fractional_digits) {
    return value.ToString($"N{fractional_digits}", Culture.culture);
  }

  // Formats `value` (given in metres) in kilometres or metres.
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
  // `offset` (for instance, if `interval` is an interval of distances from the
  // primary and `offset` is the radius of the primary, altitudes are shown).
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
      return L10N.CacheFormat("#Principia_OrbitAnalyser_Precesses");
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

  // Similar to `FormatAngleInterval`, but annotated with the equivalent
  // of the half-width as a distance at the equator in parentheses.
  public static string FormatEquatorialAngleInterval(this Interval interval,
    CelestialBody primary) {
    double half_width_angle = (interval.max - interval.min) / 2;
    if (half_width_angle > Math.PI) {
      return L10N.CacheFormat("#Principia_OrbitAnalyser_Precesses");
    }
    double half_width_distance = half_width_angle * primary.Radius;
    string formatted_distance = half_width_distance > 1000
                                    ? $"{(half_width_distance / 1000).FormatN(1)} km"
                                    : $"{(half_width_distance).FormatN(0)} m";
    return $"{interval.FormatAngleInterval()} ({formatted_distance})";
  }

  // Formats an interval of angles (given in radians) as an interval of times of
  // day, mapping [0, 2π] to one day.
  public static string FormatHourAngleInterval(this Interval interval) {
    double half_width = (interval.max - interval.min) / 2;
    double midpoint = interval.min + half_width;
    if (half_width > Math.PI) {
      return L10N.CacheFormat("#Principia_OrbitAnalyser_Precesses");
    }
    var formatted_midpoint = new PrincipiaTimeSpan(
        KSPUtil.dateTimeFormatter.Day * midpoint / (2 * Math.PI))
        .FormatPositive(with_leading_zeroes: false,
                        with_seconds: false,
                        iau_style: true);
    var formatted_half_width = new PrincipiaTimeSpan(
        KSPUtil.dateTimeFormatter.Day * half_width / (2 * Math.PI))
        .FormatPositive(with_leading_zeroes: false,
                        with_seconds: false,
                        iau_style: true);
    return $"{formatted_midpoint}±{formatted_half_width}";
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
                                   out PrincipiaTimeSpan ts)) {
      seconds = ts.total_seconds;
      return true;
    } else {
      return false;
    }
  }
}

internal abstract class OrbitAnalyser : RequiredVesselSupervisedWindowRenderer {
  protected OrbitAnalyser(PrincipiaPluginAdapter adapter,
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

  // Whether `RequestAnalysis()` needs to be called.
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
                                            analysis.mission_duration,
                                            analysis.elements,
                                            analysis.recurrence,
                                            analysis.ground_track_equatorial_crossings,
                                            analysis.solar_times_of_nodes,
                                            (int?)(analysis.mission_duration /
                                                   analysis.elements?.
                                                       nodal_period));
    }
    RenderButton(ButtonText(orbit_description_));
  }

  protected override string Title => orbit_description_ == null
                                         ? L10N.CacheFormat(
                                             "#Principia_OrbitAnalyser_Title")
                                         : orbit_description_[0].ToString().
                                               ToUpper() +
                                           orbit_description_.Substring(1);

  protected override void RenderWindowContents(int window_id) {
    if (graph_icon_ == null) {
      PrincipiaPluginAdapter.LoadTextureOrDie(out graph_icon_, "graph.png");
    }
    string vessel_guid = predicted_vessel?.id.ToString();
    if (vessel_guid == null) {
      return;
    }

    var multiline_style = Style.Multiline(UnityEngine.GUI.skin.label);
    float two_lines = multiline_style.CalcHeight(
        new UnityEngine.GUIContent("1一\n2二"),
        Width(1));
    float five_lines = multiline_style.CalcHeight(
        new UnityEngine.GUIContent("1一\n2二\n3三\n4四\n5五"),
        Width(1));

    using (new UnityEngine.GUILayout.VerticalScope()) {
      OrbitAnalysis analysis;
      using (new UnityEngine.GUILayout.VerticalScope(GUILayoutWidth(12))) {
        if (should_request_analysis) {
          mission_duration_.Render(enabled : true);
          // If the main window is hidden, make sure that the orbit analyser
          // refreshes (not at a reduced rate).
          last_background_analysis_time_ = null;
          RequestAnalysis();
        }
        UnityEngine.GUILayout.Label(AnalysingText(),
                                    multiline_style,
                                    UnityEngine.GUILayout.Height(two_lines));

        analysis = GetAnalysis();

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
      }

      OrbitalElements elements = analysis.elements;
      DrawAllGraphs(elements);
      OrbitRecurrence? recurrence = analysis.recurrence;
      EquatorialCrossings? equatorial_crossings =
          analysis.ground_track_equatorial_crossings;
      SolarTimesOfNodes? solar_times_of_nodes = analysis.solar_times_of_nodes;
      double mission_duration = analysis.mission_duration;
      CelestialBody primary = analysis.primary_index.HasValue
                                  ? FlightGlobals.Bodies[
                                      analysis.primary_index.Value]
                                  : null;

      orbit_description_ = OrbitDescription(primary,
                                            mission_duration,
                                            elements,
                                            recurrence,
                                            equatorial_crossings,
                                            solar_times_of_nodes,
                                            (int?)(mission_duration /
                                                   elements?.nodal_period));

      Style.HorizontalLine();
      string duration_in_revolutions;
      if (elements != null) {
        int sidereal_revolutions =
            (int)(mission_duration / elements.sidereal_period);
        int nodal_revolutions =
            (int)(mission_duration / elements.nodal_period);
        int anomalistic_revolutions =
            (int)(mission_duration / elements.anomalistic_period);
        int ground_track_cycles = nodal_revolutions /
                                  analysis.recurrence?.number_of_revolutions ??
                                  0;
        string duration_in_ground_track_cycles =
            ground_track_cycles > 0
                ?  L10N.CacheFormat(
                    "#Principia_OrbitAnalyser_Duration_GroundTrackCycles",
                    ground_track_cycles.FormatN(0))
                : "";
        duration_in_revolutions = L10N.CacheFormat(
            "#Principia_OrbitAnalyser_Duration_Revolutions",
            sidereal_revolutions.FormatN(0),
            nodal_revolutions.FormatN(0),
            duration_in_ground_track_cycles,
            anomalistic_revolutions.FormatN(0));
      } else {
        if (primary == null) {
          duration_in_revolutions = L10N.CacheFormat(
              "#Principia_OrbitAnalyser_Warning_NoElements1");
        } else {
          duration_in_revolutions = L10N.CacheFormat(
              "#Principia_OrbitAnalyser_Warning_NoElements2",
              primary.Name());
        }
        multiline_style = Style.Warning(multiline_style);
      }
      string analysis_description =
          primary == null
              ? L10N.CacheFormat(
                  "#Principia_OrbitAnalyser_Warning_NoPrimary",
                  predicted_vessel.vesselName,
                  mission_duration.FormatDuration(show_seconds: false))
              : L10N.CacheFormat(
                  "#Principia_OrbitAnalyser_AnalysisDescription",
                  predicted_vessel.vesselName,
                  primary.Name(),
                  mission_duration.FormatDuration(show_seconds: false),
                  duration_in_revolutions);
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        using (new UnityEngine.GUILayout.VerticalScope()) {
          using (new UnityEngine.GUILayout.HorizontalScope()) {
            using (new UnityEngine.GUILayout.VerticalScope(
                       GUILayoutWidth(12))) {
              UnityEngine.GUILayout.Label(analysis_description,
                                          multiline_style,
                                          UnityEngine.GUILayout.Height(
                                              five_lines));
              Style.HorizontalLine();
              RenderLowestAltitude(elements, primary);
              Style.HorizontalLine();
              using (new UnityEngine.GUILayout.HorizontalScope()) {
                UnityEngine.GUILayout.Label(
                    L10N.CacheFormat(
                        "#Principia_OrbitAnalyser_Elements_MeanElements"));
                if (UnityEngine.GUILayout.Button(graph_icon_,
                                                 GUILayoutWidth(1))) {
                  show_graphs_ = !show_graphs_;
                  ScheduleShrink();
                }
              }
              RenderPeriods(elements);
            }
            if (show_graphs_) {
              Style.VerticallLine();
              using (new UnityEngine.GUILayout.VerticalScope()) {
                eccentricity_vector_graph_.Render();
                UnityEngine.GUILayout.Label(
                    L10N.CacheFormat(
                        "#Principia_OrbitAnalyser_Elements_Graphs_EccentricityVector"),
                    Style.Aligned(UnityEngine.TextAnchor.UpperCenter, UnityEngine.GUI.skin.label),
                    UnityEngine.GUILayout.ExpandWidth(true));
              }
            }
          }
          RenderOrbitalElements(elements, primary);
        }
        if (show_graphs_) {
          RenderЛидовGraph();
        }
      }
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        using (new UnityEngine.GUILayout.VerticalScope(GUILayoutWidth(12))) {
          Style.HorizontalLine();
          UnityEngine.GUILayout.Label(
              L10N.CacheFormat("#Principia_OrbitAnalyser_GroundTrack"));
          RenderOrbitRecurrence(recurrence, primary);
          Style.LineSpacing();
          RenderEquatorialCrossings(equatorial_crossings, primary);
          Style.LineSpacing();
          RenderNodeMeanSolarTimes(solar_times_of_nodes);
        }
        // TODO(egg): Show the ground track?
      }
    }
    UnityEngine.GUI.DragWindow();
  }

  public static string OrbitDescription(CelestialBody primary,
                                        double mission_duration,
                                        OrbitalElements elements,
                                        OrbitRecurrence? recurrence,
                                        EquatorialCrossings? equatorial_crossings,
                                        SolarTimesOfNodes? solar_times_of_nodes,
                                        int? nodal_revolutions) {
    if (elements == null) {
      return null;
    }
    string properties = "";
    bool circular = false;
    bool equatorial = false;
    if (elements.mean_eccentricity.max < 0.01) {
      circular = true;
      properties +=
          L10N.CacheFormat(
              "#Principia_OrbitAnalyser_OrbitDescription_Circular");
    } else if (elements.mean_eccentricity.min > 0.5) {
      circular = true;
      properties +=
          L10N.CacheFormat(
              "#Principia_OrbitAnalyser_OrbitDescription_HighlyElliptical");
    }
    const double degree = Math.PI / 180;
    if (elements.mean_inclination.max < 5 * degree ||
        elements.mean_inclination.min > 175 * degree) {
      equatorial = true;
      properties +=
          L10N.CacheFormat(
              "#Principia_OrbitAnalyser_OrbitDescription_Equatorial");
    } else if (elements.mean_inclination.min > 80 * degree &&
               elements.mean_inclination.max < 100 * degree) {
      properties +=
          L10N.CacheFormat("#Principia_OrbitAnalyser_OrbitDescription_Polar");
    } else if (elements.mean_inclination.min > 90 * degree) {
      properties +=
          L10N.CacheFormat(
              "#Principia_OrbitAnalyser_OrbitDescription_Retrograde");
    }
    if (recurrence.HasValue && equatorial_crossings.HasValue) {
      Interval ascending_longitudes = equatorial_crossings.Value.
          longitudes_reduced_to_ascending_pass;
      Interval descending_longitudes = equatorial_crossings.Value.
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
                var stationary_string = L10N.CelestialStringOrNull(
                    "#Principia_OrbitAnalyser_OrbitDescription_Stationary",
                    new[]{primary});
                if (stationary_string != null) {
                  return stationary_string;
                }
                properties = L10N.CacheFormat(
                    "#Principia_OrbitAnalyser_OrbitDescription_Stationary");
              } else {
                var synchronous_string = L10N.CelestialStringOrNull(
                    "#Principia_OrbitAnalyser_OrbitDescription_Synchronous",
                    new[]{primary},
                    properties);
                if (synchronous_string != null) {
                  return synchronous_string;
                }
                properties += L10N.CacheFormat(
                    "#Principia_OrbitAnalyser_OrbitDescription_Synchronous");
              }
              break;
            case 2:
              properties += L10N.CacheFormat(
                  "#Principia_OrbitAnalyser_OrbitDescription_Semisynchronous");
              break;
            default:
              break;
          }
        }
      }
    }
    if (solar_times_of_nodes.HasValue) {
      Interval ascending_times = solar_times_of_nodes.Value.
          mean_solar_times_of_ascending_nodes;
      Interval descending_times = solar_times_of_nodes.Value.
          mean_solar_times_of_descending_nodes;
      double drift = Math.Max(
          ascending_times.max - ascending_times.min,
          descending_times.max - descending_times.min);
      double revolutions = mission_duration / elements.nodal_period;
      // We ignore 0 drift as it means that there was only one pass, which is
      // insufficient to assess synchronicity.
      // TODO(egg): The criterion should be the drift per year, not the drift
      // per revolution, but we would need to bring the mean sun here.
      if (drift > 0 && drift / revolutions < 0.001 * degree) {
        properties += L10N.CacheFormat(
            "#Principia_OrbitAnalyser_OrbitDescription_SunSynchronous");
      }
    }
    return L10N.CelestialString("#Principia_OrbitAnalyser_OrbitDescription",
                                new[]{primary},
                                properties);
  }

  private void RenderLowestAltitude(OrbitalElements elements,
                                    CelestialBody primary) {
    double? lowest_distance = elements?.radial_distance.min;
    LabeledField(
        L10N.CacheFormat("#Principia_OrbitAnalyser_Elements_LowestAltitude"),
        (lowest_distance - primary?.Radius)?.FormatAltitude());
    string altitude_warning =
        elements?.first_collision_time != null
            ?
            L10N.CacheFormat("#Principia_OrbitAnalyser_Warning_CollisionWithin",
                             (elements.first_collision_time.Value -
                              plugin.CurrentTime()).FormatDuration())
            : elements?.first_collision_risk_time != null
                ? L10N.CacheFormat(
                    "#Principia_OrbitAnalyser_Warning_CollisionRiskAfter",
                    (elements.first_collision_risk_time.Value -
                     plugin.CurrentTime()).FormatDuration())
                : elements?.first_reentry_time != null
                    ? L10N.CacheFormat(
                        "#Principia_OrbitAnalyser_Warning_ReentryIn",
                        (elements.first_reentry_time.Value -
                         plugin.CurrentTime()).FormatDuration())
                    :
                    // Compatibility code.  Old saves don't know about the
                    // atmosphere so they may miss the reentry.
                    lowest_distance < primary?.Radius + primary?.atmosphereDepth
                        ? L10N.CacheFormat(
                            "#Principia_OrbitAnalyser_Warning_Reentry")
                        : "";
    UnityEngine.GUILayout.Label(altitude_warning,
                                Style.Warning(UnityEngine.GUI.skin.label));
  }

  double last_t_min_ = double.PositiveInfinity;

  private void DrawAllGraphs(OrbitalElements elements) {
    if (elements == null) {
      return;
    }
    double t_min = double.PositiveInfinity;
    if (!elements.plottable_elements.IteratorAtEnd()) {
      t_min = elements.plottable_elements.IteratorGetPlottableElements().time;
    }
    if (t_min == last_t_min_ && !must_redraw_graphs_) {
      return;
    }
    must_redraw_graphs_ = false;
    last_t_min_ = t_min;
    DrawElementGraphs(elements);
    DrawEccentricityVectorGraph(elements);
    DrawЛидовGraph(elements);
  }

  private void DrawElementGraphs(OrbitalElements elements) {
    if (a_graph_ == null) {
      a_graph_ = new Graph((int)Width(10), (int)Height(1));
      e_graph_ = new Graph((int)Width(10), (int)Height(1));
      i_graph_ = new Graph((int)Width(10), (int)Height(1));
      Ω_graph_ = new Graph((int)Width(10), (int)Height(1));
      ω_graph_ = new Graph((int)Width(10), (int)Height(1));
      periapsis_graph_ = new Graph((int)Width(10), (int)Height(1));
      apoapsis_graph_ = new Graph((int)Width(10), (int)Height(1));
    }
    Interval t_range = Interval.Empty;
    elements.plottable_elements.IteratorReset();
    t_range.min =
        elements.plottable_elements.IteratorGetPlottableElements().time;
    for (;
         !elements.plottable_elements.IteratorAtEnd();
         elements.plottable_elements.IteratorIncrement()) {
      var elements_at_t =
          elements.plottable_elements.IteratorGetPlottableElements();
      t_range.max = elements_at_t.time;
    }
    foreach (var distance_graph in new[]
                 { a_graph_, periapsis_graph_, apoapsis_graph_ }) {
      distance_graph.PrepareCanvas(t_range,
                                   new Interval{
                                       min =
                                           elements.mean_periapsis_distance.min,
                                       max = elements.mean_apoapsis_distance.max
                                   });
    }
    e_graph_.PrepareCanvas(t_range, elements.mean_eccentricity);
    i_graph_.PrepareCanvas(t_range, elements.mean_inclination);
    Ω_graph_.PrepareCanvas(t_range, elements.mean_longitude_of_ascending_nodes);
    ω_graph_.PrepareCanvas(t_range, elements.mean_argument_of_periapsis);
    for (elements.plottable_elements.IteratorReset();
         !elements.plottable_elements.IteratorAtEnd();
         elements.plottable_elements.IteratorIncrement()) {
      var elements_at_t = elements.plottable_elements.IteratorGetPlottableElements();
      double t = elements_at_t.time;
      a_graph_.PlotPoint(t, elements_at_t.semimajor_axis, XKCDColors.Sunflower);
      e_graph_.PlotPoint(t, elements_at_t.eccentricity, XKCDColors.Cornflower);
      i_graph_.PlotPoint(t, elements_at_t.inclination, XKCDColors.Lavender);
      Ω_graph_.PlotPoint(t,
                         elements_at_t.longitude_of_ascending_node,
                         XKCDColors.LightPink);
      ω_graph_.PlotPoint(t,
                         elements_at_t.argument_of_periapsis,
                         XKCDColors.Cornflower);
      periapsis_graph_.PlotPoint(t,
                                 elements_at_t.periapsis_distance,
                                 XKCDColors.Sunflower);
      apoapsis_graph_.PlotPoint(t,
                                elements_at_t.apoapsis_distance,
                                XKCDColors.Sunflower);
    }
  }


  private void DrawEccentricityVectorGraph(OrbitalElements elements) {
    if (eccentricity_vector_graph_ == null) { 
      eccentricity_vector_graph_ = new Graph((int)Width(10), (int)Height(10));
    }
    Interval e_cos_ω_range = Interval.Empty;
    Interval e_sin_ω_range = Interval.Empty;
    for (elements.plottable_elements.IteratorReset();
         !elements.plottable_elements.IteratorAtEnd();
         elements.plottable_elements.IteratorIncrement()) {
      var elements_at_t =
          elements.plottable_elements.IteratorGetPlottableElements();
      e_cos_ω_range.Include(elements_at_t.
                                eccentricity_cos_argument_of_periapsis);
      e_sin_ω_range.Include(elements_at_t.
                                eccentricity_sin_argument_of_periapsis);
    }
    // Show a square region of the eccentricity vector space.
    if (e_cos_ω_range.measure > e_sin_ω_range.measure) {
      double midpoint = e_sin_ω_range.midpoint;
      e_sin_ω_range.min = midpoint - e_cos_ω_range.measure / 2;
      e_sin_ω_range.max = midpoint + e_cos_ω_range.measure / 2;
    } else {
      double midpoint = e_cos_ω_range.midpoint;
      e_cos_ω_range.min = midpoint - e_sin_ω_range.measure / 2;
      e_cos_ω_range.max = midpoint + e_sin_ω_range.measure / 2;
    }
    eccentricity_vector_graph_.PrepareCanvas(e_cos_ω_range, e_sin_ω_range);
    eccentricity_vector_graph_.PlotHorizontalLine(0, XKCDColors.White);
    eccentricity_vector_graph_.PlotVerticalLine(0, XKCDColors.White);
    for (elements.plottable_elements.IteratorReset();
         !elements.plottable_elements.IteratorAtEnd();
         elements.plottable_elements.IteratorIncrement()) {
      var elements_at_t = elements.plottable_elements.IteratorGetPlottableElements();
      eccentricity_vector_graph_.PlotPoint(
          elements_at_t.eccentricity_cos_argument_of_periapsis,
          elements_at_t.eccentricity_sin_argument_of_periapsis,
          XKCDColors.Cornflower);
    }
  }

  private void DrawЛидовGraph(OrbitalElements elements) {
    if (лидов_graph_ == null) {
      лидов_graph_ = new Graph((int)Width(10), (int)Height(10));
    }
    лидов_graph_.PrepareCanvas(
        new Interval{ min = -3.0 / 5.0, max = 2.0 / 5.0 },
        new Interval{ min = 0, max = 1 });
    лидов_graph_.PlotVerticalLine(0, XKCDColors.White);
    лидов_graph_.PlotFunction(Interface.GraphLidovFrozenLine,
                              new Interval{ min = -3.0 / 5.0, max = 0 },
                              XKCDColors.White);
    if (show_max_e_min_i_lines_) {
      for (int ten_e_max = 1; ten_e_max <= 10; ++ten_e_max) {
        double e_max = ten_e_max / 10.0;
        Interval c2_range =
            Interface.GraphLidovMaximalEccentricityLineC2Range(e_max);
        лидов_graph_.PlotFunction(
            c2 => Interface.GraphLidovMaximalEccentricityLine(e_max, c2),
            c2_range,
            XKCDColors.Cornflower);
      }
      for (int i_min_degrees = 0; i_min_degrees <= 80; i_min_degrees += 10) {
        лидов_graph_.PlotFunction(
            c2 => Interface.GraphLidovMinimalInclinationLine(i_min_degrees, c2),
            Interface.GraphLidovMinimalInclinationLineC2Range(i_min_degrees),
            XKCDColors.Lavender);
      }
      for (int i_min_degrees = 10; i_min_degrees <= 30; i_min_degrees += 10) {
        Interval c2 =
            Interface.GraphLidovMinimalInclinationLineC2Range(i_min_degrees);
        лидов_graph_.AddLabel(c2.max,
                              0,
                              $"{i_min_degrees}°",
                              XKCDColors.Lavender,
                              UnityEngine.TextAnchor.UpperCenter);
      }
      for (int i_min_degrees = 40; i_min_degrees <= 80; i_min_degrees += 10) {
        Interval c2 =
            Interface.GraphLidovMinimalInclinationLineC2Range(i_min_degrees);
        лидов_graph_.AddLabel(c2.min,
                              0,
                              $"{i_min_degrees}°",
                              XKCDColors.Lavender,
                              UnityEngine.TextAnchor.UpperCenter);
      }
    }
    if (show_min_e_max_i_lines_) {
      for (int ten_e_min = 1; ten_e_min <= 9; ++ten_e_min) {
        double e_min = ten_e_min / 10.0;
        Interval c2_range =
            Interface.GraphLidovMinimalEccentricityLeftLineC2Range(e_min);
        лидов_graph_.PlotFunction(
            c2 => Interface.GraphLidovMinimalEccentricityLeftLine(e_min, c2),
            c2_range,
            XKCDColors.Cornflower);
        Interface.GraphLidovMinimalEccentricityRightLineC2AndC1Max(
            e_min,
            out double c2_right,
            out double c1_max);
        лидов_graph_.PlotVerticalLine(c2_right,
                                      XKCDColors.Cornflower,
                                      new Interval{ min = 0, max = c1_max });
      }
      for (int i_max_degrees = 0; i_max_degrees <= 90; i_max_degrees += 10) {
        лидов_graph_.PlotFunction(
            c2 => Interface.GraphLidovMaximalInclinationLine(i_max_degrees, c2),
            Interface.GraphLidovMaximalInclinationLineC2Range(i_max_degrees),
            XKCDColors.Lavender);
      }
      for (int ten_e_min = 4; ten_e_min <= 9; ++ten_e_min) {
        double e_min = ten_e_min / 10.0;
        Interval c2 =
            Interface.GraphLidovMinimalEccentricityLeftLineC2Range(e_min);
        лидов_graph_.AddLabel(c2.min,
                              0,
                              $".{ten_e_min}",
                              XKCDColors.Cornflower,
                              UnityEngine.TextAnchor.UpperCenter);
      }
      for (int ten_e_min = 6; ten_e_min <= 9; ++ten_e_min) {
        double e_min = ten_e_min / 10.0;
        Interface.GraphLidovMinimalEccentricityRightLineC2AndC1Max(
            e_min,
            out double c2,
            out double _);
        лидов_graph_.AddLabel(c2,
                              0,
                              $".{ten_e_min}",
                              XKCDColors.Cornflower,
                              UnityEngine.TextAnchor.UpperCenter);
      }
    }
    // The inclination labels on the frozen curve are the same for both the max
    // and min lines (because the inclination is frozen there); it is easiest to
    // position them based on the maximal inclination lines (because they are
    // then at the minimal c₂ for the line, instead of at either end of the c₂
    // interval depending on the inclination).  Likewise the eccentricity labels
    // on the equatorial curve are the for both max and min e.
    if (show_max_e_min_i_lines_ || show_min_e_max_i_lines_) {
      for (int i_degrees = 10; i_degrees <= 60; i_degrees += 10) {
        Interval c2 =
            Interface.GraphLidovMaximalInclinationLineC2Range(i_degrees);
        double c1 =
            Interface.GraphLidovMaximalInclinationLine(i_degrees, c2.min);
        лидов_graph_.AddLabel(c2.min,
                              c1,
                              $"{i_degrees}°",
                              XKCDColors.Lavender,
                              UnityEngine.TextAnchor.MiddleRight);
      }
      for (int ten_e = 2; ten_e <= 9; ++ten_e) {
        double e = ten_e / 10.0;
        Interface.GraphLidovMinimalEccentricityRightLineC2AndC1Max(
            e,
            out double c2,
            out double c1);
        лидов_graph_.AddLabel(c2,
                              c1,
                              $" .{ten_e}",
                              XKCDColors.Cornflower,
                              UnityEngine.TextAnchor.MiddleLeft);
      }
    } else {
      // Plot the boundaries of the (c₂, c₁) space that are normally covered by
      // the coloured lines.
      foreach (int i_degrees in new[]{ 0, 90 }) {
        лидов_graph_.PlotFunction(
            c2 => Interface.GraphLidovMaximalInclinationLine(i_degrees, c2),
            Interface.GraphLidovMaximalInclinationLineC2Range(i_degrees),
            XKCDColors.White);
      }
    }
    for (elements.plottable_elements.IteratorReset();
         !elements.plottable_elements.IteratorAtEnd();
         elements.plottable_elements.IteratorIncrement()) {
      var elements_at_t = elements.plottable_elements.IteratorGetPlottableElements();
      лидов_graph_.PlotPoint(elements_at_t.lidov_c2,
                             elements_at_t.lidov_c1,
                             XKCDColors.RoseRed);
    }
  }

  private void RenderPeriods(OrbitalElements elements) {
    LabeledField(
        L10N.CacheFormat("#Principia_OrbitAnalyser_Elements_SiderealPeriod"),
        elements?.sidereal_period.FormatDuration());
    LabeledField(
        L10N.CacheFormat("#Principia_OrbitAnalyser_Elements_NodalPeriod"),
        elements?.nodal_period.FormatDuration());
    LabeledField(
        L10N.CacheFormat("#Principia_OrbitAnalyser_Elements_AnomalisticPeriod"),
        elements?.anomalistic_period.FormatDuration());
  }

  private void RenderOrbitalElements(OrbitalElements elements,
                                     CelestialBody primary) {
    // If we `show_graphs_`, elements that have an associated plot have that
    // plot to the right of the element range LabeledField, hence the
    // `HorizontalScope`s.
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      LabeledField(
          L10N.CacheFormat("#Principia_OrbitAnalyser_Elements_SemimajorAxis"),
          elements?.mean_semimajor_axis.FormatLengthInterval());
      if (show_graphs_) {
        Style.VerticalLineSpacing();
        a_graph_.Render();
      }
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      LabeledField(
          L10N.CacheFormat("#Principia_OrbitAnalyser_Elements_Eccentricity"),
          elements?.mean_eccentricity.FormatInterval());
      if (show_graphs_) {
        Style.VerticalLineSpacing();
        e_graph_.Render();
      }
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      LabeledField(
          L10N.CacheFormat("#Principia_OrbitAnalyser_Elements_Inclination"),
          elements?.mean_inclination.FormatAngleInterval());
      if (show_graphs_) {
        Style.VerticalLineSpacing();
        i_graph_.Render();
      }
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      LabeledField(
          L10N.CacheFormat(
              "#Principia_OrbitAnalyser_Elements_LongitudeOfAscendingNode"),
          elements?.mean_longitude_of_ascending_nodes.FormatAngleInterval());
      if (show_graphs_) {
        Style.VerticalLineSpacing();
        Ω_graph_.Render();
      }
    }
    LabeledField(
        L10N.CacheFormat("#Principia_OrbitAnalyser_Elements_NodalPrecession"),
        elements?.nodal_precession.FormatAngularFrequency());
    string periapsis =
        L10N.CelestialString("#Principia_OrbitAnalyser_Elements_Periapsis",
                             new[]{ primary });
    string apoapsis =
        L10N.CelestialString("#Principia_OrbitAnalyser_Elements_Apoapsis",
                             new[]{ primary });

    using (new UnityEngine.GUILayout.HorizontalScope()) {
      LabeledField(
          L10N.CacheFormat(
              "#Principia_OrbitAnalyser_Elements_ArgumentOfPeriapsis",
              periapsis),
          elements?.mean_argument_of_periapsis.FormatAngleInterval());
      if (show_graphs_) {
        Style.VerticalLineSpacing();
        ω_graph_.Render();
      }
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      LabeledField(
          L10N.CacheFormat(
              "#Principia_OrbitAnalyser_Elements_MeanPeriapsisAltitude",
              periapsis),
          elements?.mean_periapsis_distance.FormatLengthInterval(
              primary.Radius));
      if (show_graphs_) {
        Style.VerticalLineSpacing();
        periapsis_graph_.Render();
      }
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      LabeledField(
          L10N.CacheFormat(
              "#Principia_OrbitAnalyser_Elements_MeanApoapsisAltitude",
              apoapsis),
          elements?.mean_apoapsis_distance.FormatLengthInterval(
              primary.Radius));
      if (show_graphs_) {
        Style.VerticalLineSpacing();
        apoapsis_graph_.Render();
      }
    }
  }

  private void RenderЛидовGraph() {
    using (new UnityEngine.GUILayout.VerticalScope()) {
      лидов_graph_.Render();
      UnityEngine.GUILayout.Space(Height(0.5f));
      UnityEngine.GUILayout.Label(
          L10N.CacheFormat(
              "#Principia_OrbitAnalyser_Elements_Graphs_ЛидовParameters"),
              Style.Aligned(UnityEngine.TextAnchor.UpperCenter, UnityEngine.GUI.skin.label),
              UnityEngine.GUILayout.ExpandWidth(true));
      using (new UnityEngine.GUILayout.HorizontalScope()) {
        UnityEngine.GUILayout.FlexibleSpace();
        using (new UnityEngine.GUILayout.VerticalScope()) {
          UnityEngine.GUILayout.Label(
              L10N.CacheFormat(
                  "#Principia_OrbitAnalyser_Elements_Graphs_ЛидовParameters_ShowGrid"));
          if (UnityEngine.GUILayout.Toggle(
                  show_max_e_min_i_lines_,
                  L10N.CacheFormat(
                      "#Principia_OrbitAnalyser_Elements_Graphs_ЛидовParameters_MinEMaxI")) !=
              show_max_e_min_i_lines_) {
            show_max_e_min_i_lines_ = !show_max_e_min_i_lines_;
            if (show_max_e_min_i_lines_) {
              show_min_e_max_i_lines_ = false;
            }
            must_redraw_graphs_ = true;
          }
          if (UnityEngine.GUILayout.Toggle(
                  show_min_e_max_i_lines_,
                  L10N.CacheFormat(
                      "#Principia_OrbitAnalyser_Elements_Graphs_ЛидовParameters_MaxEMinI")) !=
              show_min_e_max_i_lines_) {
            show_min_e_max_i_lines_ = !show_min_e_max_i_lines_;
            if (show_min_e_max_i_lines_) {
              show_max_e_min_i_lines_ = false;
            }
            must_redraw_graphs_ = true;
          }
        }
      }
    }
  }

  private void RenderOrbitRecurrence(OrbitRecurrence? recurrence,
                                     CelestialBody primary) {
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      string νₒ = recurrence?.nuo.ToString() ?? em_dash;
      string Dᴛₒ = recurrence?.dto.ToString("+0;−0") ?? em_dash;
      string Cᴛₒ = recurrence?.cto.ToString() ?? em_dash;
      UnityEngine.GUILayout.Label(
          L10N.CacheFormat("#Principia_OrbitAnalyser_Recurrence_CapderouTriple",
                           $"[{νₒ}; {Dᴛₒ}; {Cᴛₒ}]"),
          GUILayoutWidth(8));
      UnityEngine.GUILayout.FlexibleSpace();
      autodetect_recurrence_ = UnityEngine.GUILayout.Toggle(
          autodetect_recurrence_,
          L10N.CacheFormat("#Principia_OrbitAnalyser_Recurrence_AutoDetect"),
          UnityEngine.GUILayout.ExpandWidth(false));
    }
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label(
          L10N.CacheFormat("#Principia_OrbitAnalyser_Recurrence_Cycle"));
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
          L10N.CacheFormat(
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
          L10N.CacheFormat("#Principia_OrbitAnalyser_Recurrence_Cycle_Days"),
          UnityEngine.GUILayout.ExpandWidth(false));
    }
    LabeledField(
        L10N.CacheFormat("#Principia_OrbitAnalyser_Recurrence_Subcycle"),
        L10N.CacheFormat(
            "#Principia_OrbitAnalyser_Recurrence_SubcycleLengthInDays",
            recurrence?.subcycle.FormatN(0) ?? em_dash));
    LabeledField(
        L10N.CacheFormat("#Principia_OrbitAnalyser_Recurrence_EquatorialShift"),
        recurrence?.equatorial_shift.FormatEquatorialAngle(primary));
    LabeledField(
        L10N.CacheFormat("#Principia_OrbitAnalyser_Recurrence_GridInterval"),
        recurrence?.grid_interval.FormatEquatorialAngle(primary));
  }

  private void RenderEquatorialCrossings(EquatorialCrossings? equatorial_crossings,
                                         CelestialBody primary) {
    using (new UnityEngine.GUILayout.HorizontalScope()) {
      UnityEngine.GUILayout.Label(
          L10N.CacheFormat(
              "#Principia_OrbitAnalyser_GroundTrack_LongitudesOfEquatorialCrossings_Prefix"),
          UnityEngine.GUILayout.ExpandWidth(false));
      string text = UnityEngine.GUILayout.TextField(
          $"{ground_track_revolution_}",
          GUILayoutWidth(2));
      UnityEngine.GUILayout.Label(
          L10N.CacheFormat(
              "#Principia_OrbitAnalyser_GroundTrack_LongitudesOfEquatorialCrossings_Suffix"));
      if (int.TryParse(text, out int revolution)) {
        ground_track_revolution_ = revolution;
      }
    }
    LabeledField(
        L10N.CacheFormat("#Principia_OrbitAnalyser_GroundTrack_AscendingPass"),
        equatorial_crossings?.longitudes_reduced_to_ascending_pass.
            FormatEquatorialAngleInterval(primary));
    LabeledField(
        L10N.CacheFormat("#Principia_OrbitAnalyser_GroundTrack_DescendingPass"),
        equatorial_crossings?.longitudes_reduced_to_descending_pass.
            FormatEquatorialAngleInterval(primary));
  }

  private void RenderNodeMeanSolarTimes(SolarTimesOfNodes? solar_times_of_nodes) {
    LabeledField(
        L10N.CacheFormat(
            "#Principia_OrbitAnalyser_MeanSolarTimeOfAscendingNode"),
        solar_times_of_nodes?.mean_solar_times_of_ascending_nodes.
            FormatHourAngleInterval());
    LabeledField(
        L10N.CacheFormat(
            "#Principia_OrbitAnalyser_MeanSolarTimeOfDescendingNode"),
        solar_times_of_nodes?.mean_solar_times_of_descending_nodes.
            FormatHourAngleInterval());
  }

  private void LabeledField(string label, string value) {
    using (new UnityEngine.GUILayout.HorizontalScope(GUILayoutWidth(12))) {
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
          label               :
          L10N.CacheFormat("#Principia_OrbitAnalyser_MissionDuration"),
          unit                : null,
          log10_lower_rate    : 0,
          log10_upper_rate    : 7,
          min_value           : 10,
          max_value           : double.PositiveInfinity,
          formatter           : Formatters.FormatMissionDuration,
          parser              : Formatters.TryParseMissionDuration,
          label_width         : 3,
          field_width         : 5,
          display_zero_button : false) {
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

  private UnityEngine.Texture graph_icon_;
  private bool show_graphs_;
  private Graph a_graph_;
  private Graph e_graph_;
  private Graph i_graph_;
  private Graph Ω_graph_;
  private Graph ω_graph_;
  private Graph periapsis_graph_;
  private Graph apoapsis_graph_;
  private Graph лидов_graph_;
  private Graph eccentricity_vector_graph_;
  private bool must_redraw_graphs_ = false;
  private bool show_max_e_min_i_lines_ = true;
  private bool show_min_e_max_i_lines_ = false;
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
               ? L10N.CacheFormat(
                   "#Principia_CurrentOrbitAnalyser_ToggleButton")
               : L10N.CacheFormat(
                   "#Principia_CurrentOrbitAnalyser_ToggleButtonWithDescription",
                   orbit_description);
  }

  protected override string AnalysingText() {
    return L10N.CacheFormat("#Principia_CurrentOrbitAnalyser_Analysing",
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
               ? L10N.CacheFormat(
                   "#Principia_PlannedOrbitAnalyser_ToggleButton")
               : L10N.CacheFormat(
                   "#Principia_PlannedOrbitAnalyser_ToggleButtonWithDescription",
                   orbit_description);
  }

  protected override string AnalysingText() {
    return L10N.CacheFormat("#Principia_PlannedOrbitAnalyser_Analysing",
                            predicted_vessel.vesselName);
  }

  protected override bool should_request_analysis => false;

  public int index { private get; set; }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
