using System.Collections.Generic;

namespace principia {
namespace ksp_plugin_adapter {

internal static class NullExtensions {
  internal static Value GetValueOrNull<Key, Value>(
      this Dictionary<Key, Value> dictionary,
      Key key) where Value : class {
    dictionary.TryGetValue(key, out Value value);
    return value;
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
