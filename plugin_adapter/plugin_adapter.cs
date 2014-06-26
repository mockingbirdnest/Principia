using System;
using System.Runtime.InteropServices;

using UnityEngine;

namespace principia {
namespace plugin_adapter {

[KSPAddon(KSPAddon.Startup.EveryScene, false)]
public class PluginAdapter : MonoBehaviour {

  private void Start() {
    RenderingManager.AddToPostDrawQueue(queueSpot    : 3,
                                        drawFunction : new Callback(DrawGUI));
    window_position_ = new UnityEngine.Rect(left   : Screen.width / 2.0f,
                                            top    : Screen.height / 2.0f,
                                            width  : 10,
                                            height : 10);
    Debug.Log("principia log: after Start()");
  }

  private void OnDestroy() {
    RenderingManager.RemoveFromPostDrawQueue(
        queueSpot    : 3,
        drawFunction : new Callback(DrawGUI));
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
    GUILayout.TextArea(text : Say33().ToString());
    IntPtr hello_ptr = SayHello();
    GUILayout.TextArea(text : Marshal.PtrToStringAnsi(hello_ptr));
    GUILayout.EndVertical();

    GUI.DragWindow(
        position : new Rect(left : 0f, top : 0f, width : 10000f, height : 20f));
  }

  [DllImport("test_plugin.dll", CallingConvention = CallingConvention.Cdecl)]
  private static extern int Say33();

  [DllImport("test_plugin.dll", CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr SayHello();

  private UnityEngine.Rect window_position_;
}

}  // namespace principia
}  // namespace principia
