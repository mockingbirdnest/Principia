using System.Collections.Generic;

namespace principia {
namespace ksp_plugin_adapter {

internal static class ConfigNodeParsers {

  public static BodyParameters NewCartesianBodyParameters(CelestialBody body,
                                                          ConfigNode node) {
    return new BodyParameters{
        name = body.name,
        gravitational_parameter =
            node.GetUniqueValue("gravitational_parameter"),
        reference_instant       =
            node.GetAtMostOneValue("reference_instant"),
        mean_radius             =
            node.GetAtMostOneValue("mean_radius"),
        axis_right_ascension    =
            node.GetAtMostOneValue("axis_right_ascension"),
        axis_declination        =
            node.GetAtMostOneValue("axis_declination"),
        reference_angle         =
            node.GetAtMostOneValue("reference_angle"),
        angular_frequency       =
            node.GetAtMostOneValue("angular_frequency"),
        reference_radius        =
            node.GetAtMostOneValue("reference_radius"),
        j2                      =
            node.GetAtMostOneValue("j2"),
        geopotential            =
            node.GetBodyGeopotentialElements().ToArray()};
  }

  public static ConfigurationAccuracyParameters
  NewConfigurationAccuracyParameters(ConfigNode node) {
    return new ConfigurationAccuracyParameters{
        fitting_tolerance      =
            node.GetUniqueValue("fitting_tolerance"),
        geopotential_tolerance =
            node.GetUniqueValue("geopotential_tolerance")};
  }

  public static ConfigurationAdaptiveStepParameters
  NewConfigurationAdaptiveStepParameters(ConfigNode node) {
    return new ConfigurationAdaptiveStepParameters{
        adaptive_step_size_integrator =
            node.GetUniqueValue("adaptive_step_size_integrator"),
        length_integration_tolerance  =
            node.GetUniqueValue("length_integration_tolerance"),
        speed_integration_tolerance   =
            node.GetUniqueValue("speed_integration_tolerance")};
  }

  public static ConfigurationFixedStepParameters
  NewConfigurationFixedStepParameters(ConfigNode node) {
    return new ConfigurationFixedStepParameters{
        fixed_step_size_integrator =
            node.GetUniqueValue("fixed_step_size_integrator"),
        integration_step_size      =
            node.GetUniqueValue("integration_step_size")};
  }

  public static BodyParameters NewKeplerianBodyParameters(CelestialBody body,
                                                          ConfigNode node) {
    return new BodyParameters{
        name = body.name,
        gravitational_parameter =
            node?.GetAtMostOneValue("gravitational_parameter")
                ?? (body.gravParameter + " m^3/s^2"),
        // The origin of rotation in KSP is the x of Barycentric, rather
        // than the y axis as is the case for Earth, so the right
        // ascension is -90 deg.
        reference_instant    = 
            node?.GetAtMostOneValue("reference_instant")
                ?? "JD2451545.0",
        mean_radius          =
            node?.GetAtMostOneValue("mean_radius")
                ?? (body.Radius + " m"),
        axis_right_ascension =
            node?.GetAtMostOneValue("axis_right_ascension")
                ?? "-90 deg",
        axis_declination     =
            node?.GetAtMostOneValue("axis_declination")
                ?? "90 deg",
        reference_angle      =
            node?.GetAtMostOneValue("reference_angle")
                ?? (body.initialRotation + " deg"),
        angular_frequency    =
            node?.GetAtMostOneValue("angular_frequency")
                ?? (body.angularV + " rad/s"),
        reference_radius     =
            node?.GetAtMostOneValue("reference_radius"),
        j2                   =
            node?.GetAtMostOneValue("j2"),
        geopotential            =
            node?.GetBodyGeopotentialElements()?.ToArray()};
  }

  private static List<BodyGeopotentialElement> GetBodyGeopotentialElements(
      this ConfigNode node) {
    ConfigNode[] geopotential_rows = node.GetNodes("geopotential_row");
    List<BodyGeopotentialElement> elements =
        new List<BodyGeopotentialElement>();
    foreach (ConfigNode geopotential_row in geopotential_rows) {
      string degree = geopotential_row.GetUniqueValue("degree");
      ConfigNode[] geopotential_columns =
          geopotential_row.GetNodes("geopotential_column");
      foreach (ConfigNode geopotential_column in geopotential_columns) {
        string order = geopotential_column.GetUniqueValue("order");
        string j = geopotential_column.GetAtMostOneValue("j");
        string cos = geopotential_column.GetAtMostOneValue("cos");
        string sin = geopotential_column.GetUniqueValue("sin");
        elements.Add(new BodyGeopotentialElement{degree = degree,
                                                 order = order,
                                                 cos = cos,
                                                 j = j,
                                                 sin = sin});
      }
    }
    return elements;
  }

}

}  // namespace ksp_plugin_adapter
}  // namespace principia
