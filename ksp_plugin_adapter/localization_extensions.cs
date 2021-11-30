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

  private static string FormatOrNull(string template, params object[] args) {
    // Optional translations include the language name so that they do not fall
    // back to English.
    string qualified_template = $"{template}.{Localizer.CurrentLanguage}";
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
    return GetCached(new StringCacheEntry(template, names),
                     () => CelestialOverride(template, names, bodies) ??
                           Localizer.Format(template, names));
  }

  public static string CelestialStringOrNull(string template,
                                             CelestialBody[] bodies,
                                             params object[] args) {
    string[] names =
      (from body in bodies
       from name in new[]{body.Name(), body.Initial()}
       select name).Concat(from arg in args select arg.ToString()).ToArray();
    return GetCached(new StringCacheEntry(template, names),
                     () => CelestialOverride(template, names, bodies) ??
                           FormatOrNull(template, names));
  }

  public static string CacheFormat(string name, params object[] args) {
    return GetCached(new StringCacheEntry(name, args),
                     () => Localizer.Format(name, args));
  }

  private static string GetCached(StringCacheEntry key_only,
                                  Func<string> compute_value) {
    var entry = key_only;
    if (cache_.TryGetValue(entry, out StringCacheEntry actual)) {
      entry.value = actual.value;
      cache_.Remove(actual);
      cache_by_time_.Remove(actual);
    } else {
      entry.value = compute_value();
      if (cache_.Count >= max_cache_length) {
        var evicted = cache_by_time_.First();
        cache_.Remove(evicted);
        cache_by_time_.Remove(evicted);
      }
    }
    cache_.Add(entry);
    cache_by_time_.Add(entry);
    return entry.value;
  }

  public struct StringCacheEntry : IEquatable<StringCacheEntry> {
    public StringCacheEntry(string name, object[] args) {
      string unit_separator = "\x1F";
      key = name + unit_separator +
          string.Join(unit_separator, from arg in args select arg.ToString());
      value = null;
      timestamp = System.Diagnostics.Stopwatch.GetTimestamp();
    }

    public bool Equals(StringCacheEntry other) {
      return key == other.key;
    }

    public override int GetHashCode() {
      return key.GetHashCode();
    }

    public string key { get; set; }
    public string value { get; set; }
    public long timestamp { get; set; }
  }

  private const int max_cache_length = 1024;
  private static HashSet<StringCacheEntry> cache_ =
      new HashSet<StringCacheEntry>();
  private static SortedSet<StringCacheEntry> cache_by_time_ =
      new SortedSet<StringCacheEntry>(Comparer<StringCacheEntry>.Create(
          (x, y) => {
            var time_comparison = x.timestamp.CompareTo(y.timestamp);
            return time_comparison != 0 ? time_comparison
                                        : x.key.CompareTo(y.key);
          }));
  private static Dictionary<CelestialBody, string> names_ =
      new Dictionary<CelestialBody, string>();
  private static Dictionary<CelestialBody, string> initials_ =
      new Dictionary<CelestialBody, string>();
}

}
}
