using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal class Dialog : IConfigNode {

  public void Show(String message) {
    UnityEngine.Rect dialog_window_rectangle = UnityEngine.Rect.zero;

    UnityEngine.GUI.skin = null;
    dialog_window_rectangle.xMin = x_;
    dialog_window_rectangle.yMin = y_;
    dialog_window_rectangle = UnityEngine.GUILayout.Window(
        id         : this.GetHashCode(),
        screenRect : dialog_window_rectangle,
        func       : (int id) => {
          using (new VerticalLayout())
          {
            UnityEngine.GUILayout.TextArea(message);
          }
          UnityEngine.GUI.DragWindow(
              position: new UnityEngine.Rect(x      : 0f,
                                             y      : 0f,
                                             width  : 10000f,
                                             height : 10000f));
        },
        text: "Principia",
        options: UnityEngine.GUILayout.MinWidth(500));
    WindowUtilities.EnsureOnScreen(ref dialog_window_rectangle);
    x_ = dialog_window_rectangle.xMin;
    y_ = dialog_window_rectangle.yMin;
  }

  void IConfigNode.Load(ConfigNode node) {
    String x_value = node.GetAtMostOneValue("x");
    if (x_value != null) {
      x_ = System.Convert.ToSingle(x_value);
    }
    String y_value = node.GetAtMostOneValue("y");
    if (y_value != null) {
      y_ = System.Convert.ToSingle(y_value);
    }
  }

  void IConfigNode.Save(ConfigNode node) {
    node.SetValue("x", x_, createIfNotFound: true);
    node.SetValue("y", y_, createIfNotFound: true);
  }

  private Single x_ = UnityEngine.Screen.width / 2;
  private Single y_ = UnityEngine.Screen.height / 3;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
