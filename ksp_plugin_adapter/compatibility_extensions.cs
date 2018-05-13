namespace principia {
namespace ksp_plugin_adapter {

internal static class CompatibilityExtensions {
  public static string NameWithArticle(this CelestialBody body) {
#if KSP_VERSION_1_3_1 || KSP_VERSION_1_4_3
    // We are not writing string templates for this mod, sorry.
    return body.displayName.StartsWith("The ") ? "the " + body.name : body.name;
#elif KSP_VERSION_1_2_2
    return body.theName;
#endif
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
