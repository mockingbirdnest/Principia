
using System.Linq;

namespace principia {
namespace ksp_plugin_adapter {

// The Principia config nodes correspond to protocol buffers declared in
// astronomy.proto.  The scalar (string or double) protobuf fields are string
// values in a config node, the message fields are config nodes.
// The |GetUniqueValue|, |GetUniqueNode| extension methods correspond to
// required fields.
// The |GetAtMostOneValue|, |GetAtMostNode| extension methods correspond to
// optional fields.  Note that they are also used when a KSP-side default is
// provided for a required field, e.g. using the stock values as defaults for
// gravity models.
// The KSP |GetValues|, |GetNodes| methods should be used directly for repeated
// fields.
internal static class ConfigNodeExtensions {
  public static ConfigNode? GetAtMostOneNode(this ConfigNode node,
                                             string name) {
    var nodes = node.GetNodes(name);
    if (nodes.Length > 1) {
      Log.Fatal("Duplicate |" + name + "| in node |" + node.name + "|:\n" +
                node.ToString());
    } else if (nodes.Length == 0) {
      return null;
    }
    return nodes[0];
  }

  public static ConfigNode GetUniqueNode(this ConfigNode node, string name) {
    var nodes = node.GetNodes(name);
    if (nodes.Length > 1) {
      Log.Fatal("Duplicate |" + name + "| in node |" + node.name + "|:\n" +
                node.ToString());
    } else if (nodes.Length == 0) {
      Log.Fatal("Missing |" + name + "| in node |" + node.name + "|:\n" +
                node.ToString());
    }
    return nodes[0];
  }

  public static string? GetAtMostOneValue(this ConfigNode node, string name) {
    var values = node.GetValues(name);
    if (values.Length > 1) {
      Log.Fatal("Duplicate |" + name + "| in node |" + node.name + "|:\n" +
                node.ToString());
    } else if (values.Length == 0) {
      return null;
    }
    return values[0];
  }

  public static string GetUniqueValue(this ConfigNode node, string name) {
    var values = node.GetValues(name);
    if (values.Length > 1) {
      Log.Fatal("Duplicate |" + name + "| in node |" + node.name + "|:\n" +
                node.ToString());
    } else if (values.Length == 0) {
      Log.Fatal("Missing |" + name + "| in node |" + node.name + "|:\n" +
                node.ToString());
    }
    return values[0];
  }
}

internal static class GameDatabaseExtensions {
  public static ConfigNode? GetAtMostOneNode(this GameDatabase database,
                                             string name) {
    var configs = database.GetConfigs(name);
    if (configs.Length > 1) {
      Log.Fatal(
          "Duplicate config |" + name + "| (" +
          string.Join(", ",
                      (from config in configs select config.url).ToArray()) +
          ")");
    } else if (configs.Length == 0) {
      return null;
    }
    return configs[0].config;
  }

  public static ConfigNode GetUniqueNode(this GameDatabase database,
                                         string name) {
    var configs = database.GetConfigs(name);
    if (configs.Length > 1) {
      Log.Fatal(
          "Duplicate config |" + name + "| (" +
          string.Join(", ",
                      (from config in configs select config.url).ToArray()) +
          ")");
    } else if (configs.Length == 0) {
      Log.Fatal("Missing config |" + name + "|");
    }
    return configs[0].config;
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
