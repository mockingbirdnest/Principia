using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace principia {
namespace ksp_plugin_adapter {

internal static class ConfigNodeExtensions {
  public static string GetAtMostOneValue(this ConfigNode node, string name) {
    return string.Empty;
  }
}

} // namespace ksp_plugin_adapter
} // namespace principia
