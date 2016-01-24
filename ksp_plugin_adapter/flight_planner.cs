using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace principia {
namespace ksp_plugin_adapter {

class FlightPlanner {
  public FlightPlanner(WindowRenderer.ManagerInterface manager,
                       IntPtr plugin,
                       Vessel vessel) {
    manager_ = manager;
    plugin_ = plugin;
    vessel_ = vessel;
    string vessel_guid = vessel.id.ToString();
    if (plugin_.HasVessel(vessel_guid) &&
        plugin_.FlightPlanExists(vessel_guid)) {
      burns_ = new List<BurnEditor>();
      for (int i = 0;
           i < plugin_.FlightPlanNumberOfManoeuvres(vessel_guid);
           ++i) {
        burns_.Add(new BurnEditor(manager_, plugin_, vessel_));
        burns_.Last().Reset(
            plugin_.FlightPlanGetManoeuvre(vessel.id.ToString(), i).burn);
      }
    }
  }

  public void Render() {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    UnityEngine.GUILayout.BeginVertical();
    if (burns_ == null) {
    } else {
      for (int i = 0; i < burns_.Count - 1; ++i) {
        burns_[i].Render(enabled : false);
      }
      burns_.Last().Render(enabled : true);
      if (burns_.Count > 0) {
        if (UnityEngine.GUILayout.Button(
                "Delete",
                UnityEngine.GUILayout.ExpandWidth(true))) {
          plugin_.FlightPlanRemoveLast(vessel_.id.ToString());
          burns_.Last().Close();
          burns_.RemoveAt(burns_.Count - 1);
        }
      }
      if (UnityEngine.GUILayout.Button(
              "Add",
              UnityEngine.GUILayout.ExpandWidth(true))) {
        burns_.Add(new BurnEditor(manager_, plugin_, vessel_));
      }
    }
    UnityEngine.GUILayout.EndVertical();
    UnityEngine.GUI.skin = old_skin;
  }

  private void ComputeEngineProperties() { }

  // Not owned.
  private readonly IntPtr plugin_;
  private readonly WindowRenderer.ManagerInterface manager_;
  private readonly Vessel vessel_;
  private List<BurnEditor> burns_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
