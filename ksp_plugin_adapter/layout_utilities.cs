using System;

namespace principia {
namespace ksp_plugin_adapter {

internal class HorizontalLayout : IDisposable {
  public HorizontalLayout() {
    UnityEngine.GUILayout.BeginHorizontal();
  }

  public void Dispose() {
    UnityEngine.GUILayout.EndHorizontal();
  }
}

internal class VerticalLayout : IDisposable {
  public VerticalLayout() {
    UnityEngine.GUILayout.BeginVertical();
  }

  public void Dispose() {
    UnityEngine.GUILayout.EndVertical();
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
