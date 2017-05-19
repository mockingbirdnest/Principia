using System;

namespace principia {
namespace ksp_plugin_adapter {

internal class HorizontalLayout : IDisposable {
  public HorizontalLayout(params UnityEngine.GUILayoutOption[] options) {
    UnityEngine.GUILayout.BeginHorizontal(options);
  }

  public void Dispose() {
    UnityEngine.GUILayout.EndHorizontal();
  }
}

internal class VerticalLayout : IDisposable {
  public VerticalLayout(params UnityEngine.GUILayoutOption[] options) {
    UnityEngine.GUILayout.BeginVertical(options);
  }

  public void Dispose() {
    UnityEngine.GUILayout.EndVertical();
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
