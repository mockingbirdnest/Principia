using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace principia {
namespace ksp_plugin_adapter {

class LRUCache {
  public string Get(Entry key_only, Func<string> compute_value) {
    var entry = key_only;
    if (cache_.TryGetValue(entry, out Entry actual)) {
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

  public struct Entry : IEquatable<Entry> {
    public Entry(string name, object[] args) {
      string unit_separator = "\x1F";
      key = name + unit_separator +
          string.Join(unit_separator, from arg in args select arg.ToString());
      value = null;
      timestamp = System.Diagnostics.Stopwatch.GetTimestamp();
    }

    public bool Equals(Entry other) {
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
  private HashSet<Entry> cache_ = new HashSet<Entry>();
  private SortedSet<Entry> cache_by_time_ =
      new SortedSet<Entry>(Comparer<Entry>.Create(
          (x, y) => {
            var time_comparison = x.timestamp.CompareTo(y.timestamp);
            return time_comparison != 0 ? time_comparison
                                        : x.key.CompareTo(y.key);
          }));
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
