namespace principia {
namespace ksp_plugin_adapter {

internal static class CompatibilityExtensions {
  public static string DisplayName(this CelestialBody body) {
#if KSP_VERSION_1_3
    return body.displayName;
#elif KSP_VERSION_1_2_2
    return body.theName;
#endif
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
