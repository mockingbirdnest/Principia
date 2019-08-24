namespace principia {
namespace ksp_plugin_adapter {

internal class Dialog : UnsupervisedWindowRenderer, IConfigNode {
  private Dialog() {}

  public Dialog(bool persist_state) {
    persist_state_ = persist_state;
  }

  public string Message {
    set {
      message_ = value;
      UnityEngine.Debug.LogError(message_);
    }
  }

  protected override string Title => "Principia";

  protected override void RenderWindow(int window_id) {
    using (new UnityEngine.GUILayout.VerticalScope()) {
      UnityEngine.GUILayout.TextArea(message_ ?? "SHOW WITHOUT MESSAGE",
                                     style: Style.Multiline(
                                         UnityEngine.GUI.skin.textArea));
    }
    UnityEngine.GUI.DragWindow();
  }

  public override void Load(ConfigNode node) {
    if (persist_state_) {
      base.Load(node);
      message_ = node.GetAtMostOneValue("message");
    } else {
      // For a dialog whose state is not persisted, we still get the graphic
      // properties from the base class, but we preserve show_ and we don't
      // load the message.
      bool saved_shown = Shown();
      base.Load(node);
      if (saved_shown != Shown()) {
        Toggle();
      }
    }
  }

  public override void Save(ConfigNode node) {
    base.Save(node);
    if (persist_state_ && message_ != null) {
      node.SetValue("message", message_, createIfNotFound: true);
    }
  }

  private readonly bool persist_state_;
  private string message_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
