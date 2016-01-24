using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace principia {
namespace ksp_plugin_adapter {

class FlightPlanner {
  public FlightPlanner(WindowRenderer.ManagerInterface manager,
                       IntPtr plugin) {
    manager_ = manager;
    plugin_ = plugin;
  }

  public void Render() {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    UnityEngine.GUILayout.BeginVertical();

    if (vessel_ == null || vessel_ != FlightGlobals.ActiveVessel ||
        !plugin_.HasVessel(vessel_.id.ToString())) {
      Reset();
    }

    if (vessel_ != null) {
      string vessel_guid = vessel_.id.ToString();
      if (burn_editors_ == null) {
        if (plugin_.HasVessel(vessel_guid)) {
          if (plugin_.FlightPlanExists(vessel_guid)) {
            burn_editors_ = new List<BurnEditor>();
            for (int i = 0;
                 i < plugin_.FlightPlanNumberOfManoeuvres(vessel_guid);
                 ++i) {
              // Dummy initial time, we call |Reset| immediately afterwards.
              burn_editors_.Add(
                  new BurnEditor(manager_, plugin_, vessel_, initial_time : 0));
              burn_editors_.Last().Reset(
                  plugin_.FlightPlanGetManoeuvre(vessel_guid, i).burn);
            }
          } else {
            if (UnityEngine.GUILayout.Button("Create flight plan")) {
              plugin_.FlightPlanCreate(vessel_guid,
                                       plugin_.CurrentTime() + 1000,
                                       vessel_.GetTotalMass());
            }
          }
        }
      } else {
        if (UnityEngine.GUILayout.Button("Delete flight plan")) {
          plugin_.FlightPlanDelete(vessel_guid);
        } else {
          for (int i = 0; i < burn_editors_.Count - 1; ++i) {
            burn_editors_[i].Render(enabled : false);
          }
          burn_editors_.Last().Render(enabled : true);
          if (burn_editors_.Count > 0) {
            if (UnityEngine.GUILayout.Button(
                    "Delete",
                    UnityEngine.GUILayout.ExpandWidth(true))) {
              plugin_.FlightPlanRemoveLast(vessel_guid);
              burn_editors_.Last().Close();
              burn_editors_.RemoveAt(burn_editors_.Count - 1);
            }
          }
          if (UnityEngine.GUILayout.Button(
                  "Add",
                  UnityEngine.GUILayout.ExpandWidth(true))) {
            double initial_time =
                (burn_editors_.Count == 0
                     ? plugin_.CurrentTime()
                     : plugin_.FlightPlanGetManoeuvre(
                           vessel_guid,
                           burn_editors_.Count - 1).final_time) + 60;
            burn_editors_.Add(
                new BurnEditor(manager_, plugin_, vessel_, initial_time));
            bool inserted = plugin_.FlightPlanAppend(
                                vessel_guid,
                                burn_editors_.Last().Burn());
            if (!inserted) {
              burn_editors_.RemoveAt(burn_editors_.Count - 1);
            }
          }
        }
      }
    }
    UnityEngine.GUILayout.EndVertical();
    UnityEngine.GUI.skin = old_skin;
  }

  private void Reset() {
    foreach (BurnEditor editor in burn_editors_) {
      editor.Close();
    }
    burn_editors_ = null;
    vessel_ = FlightGlobals.ActiveVessel;
  }

  // Not owned.
  private readonly IntPtr plugin_;
  private readonly WindowRenderer.ManagerInterface manager_;
  private Vessel vessel_;
  private List<BurnEditor> burn_editors_;

  // TODO(egg): make mutable.
  private const double excess_time_ = 60;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
