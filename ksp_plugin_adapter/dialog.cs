using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal class Dialog {

  public Dialog(ref int x, ref int y) {
    x = UnityEngine.Screen.width / 2;
    y = UnityEngine.Screen.height / 3;
  }

  public void Show(String message, ref int x, ref int y) {
    UnityEngine.Rect dialog_window_rectangle = UnityEngine.Rect.zero;

    UnityEngine.GUI.skin = null;
    dialog_window_rectangle.xMin = x;
    dialog_window_rectangle.yMin = y;
    dialog_window_rectangle = UnityEngine.GUILayout.Window(
        id: this.GetHashCode() + 1,
        screenRect: dialog_window_rectangle,
        func: (int id) => {
          using (new VerticalLayout())
          {
            UnityEngine.GUILayout.TextArea(message);
          }
          UnityEngine.GUI.DragWindow(
              position: new UnityEngine.Rect(x: 0f,
                                             y: 0f,
                                             width: 10000f,
                                             height: 10000f));
        },
        text: "Principia",
        options: UnityEngine.GUILayout.MinWidth(500));
    WindowUtilities.EnsureOnScreen(ref dialog_window_rectangle);
    x = (int)dialog_window_rectangle.xMin;
    y = (int)dialog_window_rectangle.yMin;
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
