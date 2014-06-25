using System;
using System.Runtime.InteropServices;

using UnityEngine;

namespace principia {
namespace plugin_adapter {

[KSPAddon(KSPAddon.Startup.Flight, false)]
public class PluginAdapter {
  [StructLayout(LayoutKind.Sequential)]

  private void Start() {
    RenderingManager.AddToPostDrawQueue(3, new Callback(DrawGUI));
    window_position_ = new UnityEngine.Rect(left   : Screen.width / 2.0f,
                                            top    : Screen.height / 2.0f,
                                            width  : 10,
                                            height : 10);
  }

  private void DrawGUI() {
    GUI.skin = HighLogic.Skin;
    window_position_ = GUILayout.Window(
        id         : 1,
        screenRect : window_position_,
        func       : DrawMainWindow,
        text       : "Traces of Various Descriptions",
        options    : GUILayout.MinWidth(500));
  }

  private void DrawMainWindow(int window_id) {
    GUIStyle style            = new GUIStyle(GUI.skin.button);
    style.normal.textColor    = style.focused.textColor = Color.white;
    style.hover.textColor     = style.active.textColor  = Color.yellow;
    style.onNormal.textColor  = Color.green;
    style.onFocused.textColor = Color.green;
    style.onHover.textColor   = Color.green;
    style.onActive.textColor  = Color.green;
    style.padding             = new RectOffset(8, 8, 8, 8);
    GUILayout.BeginVertical();

  }

  [DllImport("test_plugin.dll")]
  private static extern string SayHello();

  private UnityEngine.Rect window_position_;
}

}  // namespace principia
}  // namespace principia
