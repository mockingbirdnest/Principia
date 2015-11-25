using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace principia {
namespace ksp_plugin_adapter {

class FlightPlanner {

  private List<BurnEditor> burns_;

  public FlightPlanner() {
    burns_ = new List<BurnEditor>();
  }

  public void Render() {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    UnityEngine.GUILayout.BeginVertical();
    foreach (BurnEditor burn_editor in burns_) {
      burn_editor.Render();
    } 
    if (burns_.Count > 0) {
      if (UnityEngine.GUILayout.Button(
              "Delete",
              UnityEngine.GUILayout.ExpandWidth(true))) {
        burns_.RemoveAt(burns_.Count - 1);
      }
    }
    if (UnityEngine.GUILayout.Button("Add",
                                     UnityEngine.GUILayout.ExpandWidth(true))) {
      burns_.Add(new BurnEditor());
    }
    UnityEngine.GUILayout.EndVertical();
    UnityEngine.GUI.skin = old_skin;
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
