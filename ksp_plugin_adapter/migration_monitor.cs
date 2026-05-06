using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace principia {
namespace ksp_plugin_adapter {

internal class MigrationMonitor : UnsupervisedWindowRenderer {
  public MigrationMonitor(IntPtr plugin_reader) {
    plugin_reader_ = plugin_reader;
  }

  // TODO(egg): L10N.
  protected override string Title => "Principia save file migration";

  protected override void RenderWindow(int window_id) {
    using (new UnityEngine.GUILayout.VerticalScope(GUILayoutWidth(60))) {
      // TODO(egg): L10N.
      UnityEngine.GUILayout.Label(
          "The Principia save requires reprocessing, see Principia issue #4490.\n" +
          "This may take a while; the game will unpause when done…",
          Style.Warning(Style.Multiline(UnityEngine.GUI.skin.label)));
      UnityEngine.GUILayout.TextArea(plugin_reader_.PluginReaderLogs(),
                                     style: Style.Multiline(
                                         UnityEngine.GUI.skin.textArea));
    }
    UnityEngine.GUI.DragWindow();
  }

  private readonly IntPtr plugin_reader_;
}

} // namespace ksp_plugin_adapter
} // namespace principia
