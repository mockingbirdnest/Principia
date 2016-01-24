using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace principia {
namespace ksp_plugin_adapter {

class FlightPlanner : WindowRenderer {
  public FlightPlanner(WindowRenderer.ManagerInterface manager,
                       IntPtr plugin) : base(manager) {
    manager_ = manager;
    plugin_ = plugin;
    window_rectangle_.x = UnityEngine.Screen.width / 2;
    window_rectangle_.y = UnityEngine.Screen.height / 3;
  }

  protected override void RenderWindow() {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    if (show_planner_) {
      window_rectangle_ = UnityEngine.GUILayout.Window(
                              id         : this.GetHashCode(),
                              screenRect : window_rectangle_,
                              func       : RenderPlanner,
                              text       : "Flight plan",
                              options    : UnityEngine.GUILayout.MinWidth(500));
    }
    UnityEngine.GUI.skin = old_skin;
  }

  public void RenderButton() {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    if (UnityEngine.GUILayout.Button("Flight plan...")) {
      show_planner_ = !show_planner_;
    }
    UnityEngine.GUI.skin = old_skin;
  }

  private void RenderPlanner(int window_id) {
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
              Shrink();
            }
          }
        }
      } else {
        if (UnityEngine.GUILayout.Button("Delete flight plan")) {
          plugin_.FlightPlanDelete(vessel_guid);
          Reset();
        } else {
          for (int i = 0; i < burn_editors_.Count - 1; ++i) {
            burn_editors_[i].Render(enabled : false);
          }
          if (burn_editors_.Count > 0) {
            BurnEditor last_burn = burn_editors_.Last();
            if (last_burn.Render(enabled : true)) {
              plugin_.FlightPlanReplaceLast(vessel_guid, last_burn.Burn());
              last_burn.Reset(
                  plugin_.FlightPlanGetManoeuvre(
                      vessel_guid,
                      burn_editors_.Count - 1).burn);
            }
            if (UnityEngine.GUILayout.Button(
                    "Delete",
                    UnityEngine.GUILayout.ExpandWidth(true))) {
              plugin_.FlightPlanRemoveLast(vessel_guid);
              burn_editors_.Last().Close();
              burn_editors_.RemoveAt(burn_editors_.Count - 1);
              Shrink();
            }
          }
          if (UnityEngine.GUILayout.Button(
                  "Add",
                  UnityEngine.GUILayout.ExpandWidth(true))) {
            double initial_time;
            if (burn_editors_.Count == 0) {
              initial_time = plugin_.CurrentTime() + 60;
            } else {
              initial_time =
                  plugin_.FlightPlanGetManoeuvre(
                      vessel_guid,
                      burn_editors_.Count - 1).final_time + 60;
            }
            burn_editors_.Add(
                new BurnEditor(manager_, plugin_, vessel_, initial_time));
            Burn candidate_burn = burn_editors_.Last().Burn();
            Log.Info("CANDIDATE BURN " + candidate_burn.delta_v.x + " " +
                     candidate_burn.delta_v.y + " " + candidate_burn.delta_v.z +
                     " " + candidate_burn.frame.extension + " " +
                     candidate_burn.frame.centre_index + " " +
                     candidate_burn.frame.primary_index + " " +
                     candidate_burn.frame.secondary_index + " " +
                     candidate_burn.initial_time + " " +
                     candidate_burn.specific_impulse_in_seconds_g0 + " " +
                     candidate_burn.thrust_in_kilonewtons);
            bool inserted = plugin_.FlightPlanAppend(vessel_guid,
                                                     candidate_burn);
            if (!inserted) {
              burn_editors_.RemoveAt(burn_editors_.Count - 1);
            }
            Shrink();
          }
        }
      }
    }
    UnityEngine.GUILayout.EndVertical();

    UnityEngine.GUI.DragWindow(
        position : new UnityEngine.Rect(left : 0f, top : 0f, width : 10000f,
                                        height : 10000f));

    UnityEngine.GUI.skin = old_skin;
  }

  private void Reset() {
    if (burn_editors_ != null) {
      foreach (BurnEditor editor in burn_editors_) {
        editor.Close();
      }
      Shrink();
    }
    burn_editors_ = null;
    vessel_ = FlightGlobals.ActiveVessel;
  }

  private void Shrink() {
    window_rectangle_.height = 0.0f;
    window_rectangle_.width = 0.0f;
  }

  // Not owned.
  private readonly IntPtr plugin_;
  private readonly WindowRenderer.ManagerInterface manager_;
  private Vessel vessel_;
  private List<BurnEditor> burn_editors_;

  private bool show_planner_ = false;
  private UnityEngine.Rect window_rectangle_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
