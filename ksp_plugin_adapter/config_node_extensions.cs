
namespace principia {
namespace ksp_plugin_adapter {

internal static class ConfigNodeExtensions {
  public static ConfigNode GetAtMostOneNode(this ConfigNode node, string name) {
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

  public static string GetAtMostOneValue(this ConfigNode node, string name) {
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

}  // namespace ksp_plugin_adapter
}  // namespace principia
