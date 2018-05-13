namespace principia {
namespace ksp_plugin_adapter {

internal static class CompatibilityExtensions {
  public static string NameWithArticle(this CelestialBody body) {
    // We are not writing string templates for this mod, sorry.
    return body.displayName.StartsWith("The ") ? "the " + body.name : body.name;
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
