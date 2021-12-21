using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace principia {
namespace ksp_plugin_adapter {

class LRUCache {
  public string Get(string name, Func<string> compute_value) {
    return Get(name, new string[]{}, compute_value);
  }

  public string Get(string name, string[] args, Func<string> compute_value) {
    string key = MakeKey(name, args);
    Entry entry;
    if (cache_.TryGetValue(key, out entry)) {
      // Note that we must remove the entry from |cache_by_time_| before
      // updating the timestamp, since the timestamp is the key there.
      cache_by_time_.Remove(entry);
      entry.timestamp =
          System.Diagnostics.Stopwatch.GetTimestamp();
    } else {
      entry = new Entry(key, compute_value());
      if (cache_.Count >= max_cache_length) {
        var evicted = cache_by_time_.First();
        cache_.Remove(evicted.key);
        cache_by_time_.Remove(evicted);
      }
      cache_.Add(key, entry);
    }
    cache_by_time_.Add(entry);
    return entry.value;
  }

  public string Get(string name, object[] args, Func<string> compute_value) {
    return Get(name, (from arg in args select arg.ToString()).ToArray(), compute_value);
  }

  public string MakeKey(string name, string[] args) {
    string unit_separator = "\x1F";
    return name + unit_separator + string.Join(unit_separator, args);
  }

  private class Entry {
    public Entry(string key, string value) {
      this.key = key;
      this.value = value;
      timestamp = System.Diagnostics.Stopwatch.GetTimestamp();
    }

    public string key { get; }
    public string value { get; }
    public long timestamp { get; set; }
  }

  private const int max_cache_length = 1024;
  private Dictionary<string, Entry> cache_ = new Dictionary<string, Entry>();
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
