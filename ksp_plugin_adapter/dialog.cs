namespace principia {
namespace ksp_plugin_adapter {

internal class Dialog : UnsupervisedWindowRenderer, IConfigNode {
  public string Message {
    set {
      message_ = value;
      UnityEngine.Debug.LogError(message_);
    }
  }

  protected override string Title => "Principia";

  protected override void RenderWindow(int window_id) {
    using (new UnityEngine.GUILayout.VerticalScope())
    {
      UnityEngine.GUILayout.TextArea(message_ ?? "SHOW WITHOUT MESSAGE");
    }
    UnityEngine.GUI.DragWindow();
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

  private string message_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
