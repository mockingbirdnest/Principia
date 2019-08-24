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
    using (new UnityEngine.GUILayout.VerticalScope()) {
      UnityEngine.GUILayout.TextArea(message_ ?? "SHOW WITHOUT MESSAGE",
                                     style: Style.Multiline(
                                         UnityEngine.GUI.skin.textArea));
    }
    UnityEngine.GUI.DragWindow();
  }

  public new void Load(ConfigNode node) {
    // Whether a dialog is shown or not should not be a persisted property.
    // It's convenient to save and restore show_ in the base class, but we
    // override this behaviour here.  For the same reason, we don't persist the
    // message.
    bool saved_shown = Shown();
    base.Load(node);
    if (saved_shown != Shown()) {
      Toggle();
    }
  }

  private string message_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
