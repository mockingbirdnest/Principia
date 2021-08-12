using KSP.Localization;

namespace principia {
namespace ksp_plugin_adapter {

internal static class L10N {
  public static string NameWithoutArticle(this CelestialBody body) {
    // This will need to be adjusted when we add support for other languages
    // with articles.
    return body.displayName.StartsWith("The ") ? body.name : body.displayName;
  }

  public static string NameWithArticle(this CelestialBody body) {
    return body.displayName.StartsWith("The ") ? "the " + body.name
                                               : body.displayName;
  }

  public static string FormatOrNull(string template, params object[] args) {
    // Optional translations include the language name so that they do not fall
    // back to English.
    string qualified_template = $"{template}.{Localizer.CurrentLanguage}";
    if (!Localizer.Tags.ContainsKey(qualified_template)) {
      return null;
    }
    return Localizer.Format(qualified_template, args);
  }
}

}
}
