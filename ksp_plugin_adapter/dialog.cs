using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal class Dialog : WindowRenderer, IConfigNode {
  public Dialog(PrincipiaPluginAdapter adapter) : base(adapter) {}

  public String Message {
    set {
      message_ = value;
      UnityEngine.Debug.LogError(message_);
    }
  }

  public override void RenderWindow() {
    UnityEngine.GUI.skin = null;
    Window(func : (int id) => {
             using (new UnityEngine.GUILayout.VerticalScope()) {
               UnityEngine.GUILayout.TextArea(
                   message_ ?? "SHOW WITHOUT MESSAGE");
             }
             UnityEngine.GUI.DragWindow();
           },
           text : "Principia");
    EnsureOnScreen();
  }

  void IConfigNode.Load(ConfigNode node) {
    String x_value = node.GetAtMostOneValue("x");
    if (x_value != null) {
      rectangle_.x = System.Convert.ToSingle(x_value);
    }
    String y_value = node.GetAtMostOneValue("y");
    if (y_value != null) {
      rectangle_.y = System.Convert.ToSingle(y_value);
    }
    message_ = node.GetAtMostOneValue("message");
  }

  void IConfigNode.Save(ConfigNode node) {
    node.SetValue("x", rectangle_.x, createIfNotFound : true);
    node.SetValue("y", rectangle_.y, createIfNotFound : true);
    if (message_ != null) {
      node.SetValue("message", message_, createIfNotFound : true);
    }
  }

  // The message shown, if any.
  private String message_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
