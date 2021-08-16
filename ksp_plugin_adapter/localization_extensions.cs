using KSP.Localization;
using System;
using System.Linq;
using System.Text.RegularExpressions;

namespace principia {
namespace ksp_plugin_adapter {

internal static class L10N {
  private const string english_us_ = "en-us";
  private const string french_ = "fr-fr";

  public static bool IsCJKV(string text) {
    return Regex.IsMatch(
        text,
        @"^[\p{IsCJKUnifiedIdeographs}\p{IsCJKUnifiedIdeographsExtensionA}]*$");
  }

  public static string ZWSPToHyphenBetweenNonCJK(string text) {
    return Regex.Replace(
        text,
        @"([^\p{IsCJKUnifiedIdeographs}\p{IsCJKUnifiedIdeographsExtensionA}])" +
        "\u200B" +
        @"([^\p{IsCJKUnifiedIdeographs}\p{IsCJKUnifiedIdeographsExtensionA}])",
        "$1-$2");
  }

  public static string Standalone(string name) {
    return Localizer.Format("#Principia_GrammaticalForm_Standalone", name);
  }

  public static string StandaloneName(this CelestialBody celestial) {
    return Standalone(celestial.NameWithoutArticle());
  }

  private static bool StartsWithCapitalizedDefiniteArticle(string s) {
    return (Localizer.CurrentLanguage == english_us_ && s.StartsWith("The ")) ||
           (Localizer.CurrentLanguage == french_ &&
            (s.StartsWith("La ") || s.StartsWith("Le ")));
  }

  private static string LingoonaUnqualified(string s) {
    return s.Split(new[]{'^'}, 2)[0];
  }

  private static string LingoonaQualifiers(string s) {
    if (!s.Contains('^')) {
      return "";
    }
    return s.Split(new[]{'^'}, 2)[1];
  }

  private static string LingoonaQualify(string s, string qualifiers) {
    return qualifiers == "" ? s : $"{s}^{qualifiers}";
  }

  private static string NameWithoutArticle(this CelestialBody body) {
    return StartsWithCapitalizedDefiniteArticle(body.displayName)
        ? LingoonaQualify(
              body.displayName.Split(new[]{' '}, 2)[1],
              LingoonaQualifiers(body.displayName).Replace("d", ""))
        : body.displayName;
  }

  public static string NameWithArticle(this CelestialBody body) {
    return StartsWithCapitalizedDefiniteArticle(body.displayName)
        ? body.displayName.Contains('^')
              ? body.displayName.Split(new[]{' '}, 2)[1] + "d"
              : body.displayName.Split(new[]{' '}, 2)[1] + "^d"
        : body.displayName;
  }

  private static string Initial(this CelestialBody body) {
    return IsCJKV(LingoonaUnqualified(body.displayName))
        ? body.displayName
        : body.NameWithoutArticle()[0].ToString();
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

  private static string CelestialOverride(string template,
                                          string[] names,
                                          CelestialBody[] bodies) {
    return FormatOrNull(
        $"{template}({string.Join(",", from body in bodies select body.name)})",
        names);
  }

  public static string CelestialString(string template,
                                       CelestialBody[] bodies,
                                       params object[] args) {
    string[] names =
      (from body in bodies
       from name in new[]{NameWithArticle(body),
                          NameWithoutArticle(body),
                          Initial(body)}
       select name).Concat(from arg in args select arg.ToString()).ToArray();
    return CelestialOverride(template, names, bodies) ??
        Localizer.Format(template, names);
  }

  public static string CelestialStringOrNull(string template,
                                             CelestialBody[] bodies,
                                             params object[] args) {
    string[] names =
      (from body in bodies
       from name in new[]{NameWithArticle(body),
                          NameWithoutArticle(body),
                          Initial(body)}
       select name).Concat(from arg in args select arg.ToString()).ToArray();
    return CelestialOverride(template, names, bodies) ??
        FormatOrNull(template, names);
  }
}

}
}
