using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal class Dialog : UnsupervisedWindowRenderer, IConfigNode {
  public String Message {
    set {
      message_ = value;
      UnityEngine.Debug.LogError(message_);
    }
  }

  public void Show() {
    RenderWindow();
  }

  protected override void RenderWindow() {
    UnityEngine.GUI.skin = null;
    Window(func : (int id) => {
             using (new UnityEngine.GUILayout.VerticalScope())
             {
               UnityEngine.GUILayout.TextArea(
                   message_ ?? "SHOW WITHOUT MESSAGE");
             }
             UnityEngine.GUI.DragWindow();
           },
           text : "Principia");
  }

  public new void Load(ConfigNode node) {
    base.Load(node);
    message_ = node.GetAtMostOneValue("message");
  }

  public new void Save(ConfigNode node) {
    base.Save(node);
    if (message_ != null) {
      node.SetValue("message", message_, createIfNotFound : true);
    }
  }

  private String message_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
