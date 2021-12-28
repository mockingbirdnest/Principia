using KSP.Localization;
using System;
using System.Collections.Generic;
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
    return L10N.CacheFormat("#Principia_GrammaticalForm_Standalone", name);
  }

  public static string StandaloneName(this CelestialBody celestial) {
    return Standalone(celestial.Name());
  }

  private static bool StartsWithCapitalizedDefiniteArticle(string s) {
    return (UILanguage() == english_us_ && s.StartsWith("The ")) ||
           (UILanguage() == french_ &&
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

  // Returns the localized name with the appropriate Lingoona grammatical tags.
  // Note that for bodies that may have an article (the Moon, la Terre), the
  // article is absent, and instead the gender is indicated in lowercase
  // (Moon^n, Terre^f), so that an article may be requested using a placeholder
  // of the form <<A:1>>.
  // Bodies that may not have an article (Saturn, Vénus) have an uppercase
  // gender tag (Saturn^N, Vénus^F).
  public static string Name(this CelestialBody body) {
    if (names_.TryGetValue(body, out string result)) {
      return result;
    }
    string name = LingoonaUnqualified(body.displayName);
    string qualifiers = LingoonaQualifiers(body.displayName);
    if (qualifiers == "") {
      // Stock English has everything as a neuter name (^N).
      // For mods that did not try tagging grammar (which is hardly a problem
      // in English since stock strings do not add articles), tag as neuter name
      // (and switch to neuter noun below if we see an article).
      qualifiers = "N";
    }
    if (StartsWithCapitalizedDefiniteArticle(body.displayName)) {
      name = name.Split(new[]{' '}, 2)[1];
      // Lowercase the gender, allowing for articles.
      qualifiers = char.ToLower(qualifiers[0]) + qualifiers.Substring(1);
    }
    return names_[body] = LingoonaQualify(name, qualifiers);
  }

  private static string Initial(this CelestialBody body) {
    if (initials_.TryGetValue(body, out string result)) {
      return result;
    }
    return initials_[body] = IsCJKV(LingoonaUnqualified(body.displayName))
        ? body.displayName
        : body.Name()[0].ToString();
  }

  private static string UILanguage() {
    return CacheFormat("#Principia_UILanguage");
  }

  private static string FormatOrNull(string template, params object[] args) {
    // Optional translations include the language name so that they do not fall
    // back to English.
    string qualified_template = $"{template}.{UILanguage()}";
    if (!Localizer.Tags.ContainsKey(qualified_template)) {
      return null;
    }
    return L10N.CacheFormat(qualified_template, args);
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
       from name in new[]{body.Name(), body.Initial()}
       select name).Concat(from arg in args select arg.ToString()).ToArray();
    return lru_cache_.Get(template, names,
                          () => CelestialOverride(template, names, bodies) ??
                                Format(template, names));
  }

  public static string CelestialStringOrNull(string template,
                                             CelestialBody[] bodies,
                                             params object[] args) {
    string[] names =
      (from body in bodies
       from name in new[]{body.Name(), body.Initial()}
       select name).Concat(from arg in args select arg.ToString()).ToArray();
    return lru_cache_.Get(template, names,
                          () => CelestialOverride(template, names, bodies) ??
                                FormatOrNull(template, names));
  }

  public static string CacheFormat(string name) {
    return lru_cache_.Get(name, () => Localizer.Format(name));
  }

  public static string CacheFormat(string name, params string[] args) {
    return lru_cache_.Get(name, args, () => Format(name, args));
  }

  private static string Format(string name, params string[] args) {
    // KSP substitutes arguments that are themselves localization keys,
    // e.g., the names of the default vessels.
    var resolved_args = new List<string>();
    foreach (var arg in args) {
      // KSP treats nulls as empty strings here (but not in the overload
      // that takes objects).
      resolved_args.Add(
          arg == null ? "" : (Localizer.Tags.GetValueOrNull(arg) ?? arg));
    }

    // Nested <<>> rules don’t work prior to Lingoona 1.7.0 (2021-07-07).
    // We need them for some of the more advanced grammatical constructs.
    // We manually replace <<X:n>> placeholders and leave the result to
    // Lingoona. This allows for some substitution inside <<>> rules, and is
    // enough for our purposes.  <<X:n>> placeholders are substituted
    // untransformed, so the behaviour is consistent with what Lingoona 1.7.0 or
    // later would do.
    string template = Localizer.GetStringByTag(name);
    for (int n = 1; n <= resolved_args.Count; ++n) {
      template = template.Replace($"<<X:{n}>>", resolved_args[n - 1]);
    }

    Lingoona.Grammar.setLanguage(UILanguage());
    // These escapes are handled by KSP *after* formatting, contrary to
    // the \u ones and the ｢｣ for {}, which are handled when the Localizer
    // loads the tags.
    var result = Lingoona.Grammar.useGrammar(template, resolved_args)
        .Replace(@"\n", "\n").Replace(@"\""", "\"").Replace(@"\t", "\t");
    Lingoona.Grammar.setLanguage(Localizer.CurrentLanguage);
    return result;
  }

  public static string CacheFormat(string name, params object[] args) {
    return CacheFormat(name, (from arg in args select arg.ToString()).ToArray());
  }

  private static LRUCache lru_cache_ = new LRUCache();
  private static Dictionary<CelestialBody, string> names_ =
      new Dictionary<CelestialBody, string>();
  private static Dictionary<CelestialBody, string> initials_ =
      new Dictionary<CelestialBody, string>();
}

}
}
