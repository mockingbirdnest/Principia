using System.Collections.Generic;

namespace principia {
namespace ksp_plugin_adapter {

internal static class NullExtensions {
  internal static Value GetValueOrNull<Key, Value>(
      this Dictionary<Key, Value> dictionary,
      Key key) where Value : class {
    Value value;
    dictionary.TryGetValue(key, out value);
    return value;
  }

  // Similar to |T T?.GetValueOrDefault(T)| for a struct T; note that there is
  // no useful analogue to |T T?.GetValueOrDefault()| in the case of a class.
  internal static Value GetValueOrDefault<Value>(this Value value,
                                                 Value default_value)
  where Value : class {
    return value == null ? default_value : value;
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
