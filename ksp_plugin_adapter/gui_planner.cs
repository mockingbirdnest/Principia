/*
 * Copyright© (c) 2018 Maarten Maathuis, (aka madman2003).
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
 
using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using UnityEngine;

namespace principia {
namespace ksp_plugin_adapter {

    [KSPAddon(KSPAddon.Startup.Flight, false)]
    public sealed class PlannerGUI: MonoBehaviour
    {
        private const float button_width = 100.0f;
        private const float button_height = 25.0f;
        private const string magic_maneuver_string = "Delete me please";

        //
        // Definitions to support the execution page
        //
        private const string show_on_navball_string = "<color=#ffffffff>Show on navball</color>";
        private const float show_on_navball_string_length = 70f;
        private const string ignition_delta_time_name_string = "<color=#ffffffff>Ignition Δt: </color>";
        private const string cutoff_delta_time_name_string = "<color=#ffffffff>Cutoff Δt: </color>";
        private const float ignition_delta_time_name_string_length = 60f;

        private const string total_delta_velocity_prefix_string = "<color=#ffffffff>Total Δv: </color>";
        private const float total_delta_velocity_prefix_string_length = 60f;
        
        //
        // Definitions to support the planner page
        //
        private const string velocity_tangent_string = "<color=#ffff00ff>Δv tangent</color>";
        private const string velocity_normal_string = "<color=#00ffffff>Δv normal</color>";
        private const string velocity_binormal_string = "<color=#ff00ffff>Δv binormal</color>";
        private const float velocity_name_string_length = 60f;
        private const string velocity_value_string = "{0:E3} m/s";
        private const float velocity_value_string_length = 90f;

        private const string time_name_string = "<color=#ffffffff>t</color>";
        private const float time_name_string_length = 10f;
        private const string delta_time_name_string = "<color=#ffffffff>Δt</color>";
        private const float delta_time_name_string_length = 20f;
        private const float time_value_string_length_two_lines = 70f;
        private const float time_value_string_length_single_line = 130f;

        // Fixed reference frame as opposed to a frenet frame that moves with the object
        private const string inertially_fixed_burn_frame_string = "<color=#ffffffff>Inertially fixed</color>";
        private const float inertially_fixed_burn_frame_string_length = 70f;
        private const string inertial_frame_string = "<color=#ffffffff>Inertial frame</color>";
        private const string frenet_frame_string = "<color=#ffffffff>Frenet frame</color>";
        private const float inertial_or_frenet_frame_string_length = 80f;

        private const string burn_mode_prefix_string = "<color=#ffffffff>Burn mode: </color>";
        private const float burn_mode_prefix_string_length = 50f;
        private const float burn_mode_string_length = 50f;

        private const string burn_delta_velocity_prefix_string = "<color=#ffffffff>Burn Δv: </color>";
        private const float burn_delta_velocity_prefix_string_length = 40f;
        private const string burn_delta_velocity_string = "{0:E3} m/s";
        private const float burn_delta_velocity_string_length = 90f;

        private const string burn_time_prefix_string = "<color=#ffffffff>Burn Δt: </color>";
        private const float burn_time_prefix_string_length = 40f;
        private const string burn_time_string = "{0:E3} s";
        private const float burn_time_string_length = 80f;

        //
        // Defintions to support the settings page
        //
        private const string plan_length_time_name_string = "<color=#ffffffff>Plan length Δt: </color>";
        private const float plan_length_time_name_string_length = 80f;
        private const string max_steps_per_segment_name_string = "<color=#ffffffff>Max steps per segment: </color>";
        private const float max_steps_per_segment_name_string_length = 140f;
        private const string max_steps_per_segment_value_string = "{0}";
        private const float max_steps_per_segment_value_string_length = 30f;
        private const string tolerance_name_string = "<color=#ffffffff>Tolerance: </color>";
        private const float tolerance_name_string_length = 50f;
        private const string tolerance_value_string = "{0:E3} m";
        private const float tolerance_value_string_length = 60f;

        private bool planner_window_visible = false;
        private float x_pos = 0.95f;
        private float y_pos = 0.05f;

        private DialogGUIBase execution_page;
        private DialogGUIBase planning_page;
        private DialogGUIBase settings_page;
        private DialogGUIBase planner_box;
        private MultiOptionDialog multi_page_planner;
        private PopupDialog popup_dialog;

        private KSP.UI.Screens.ApplicationLauncherButton toolbar_button;

        public PlannerGUI()
        {
            InitializePlannerGUI();
        }

        private void Awake()
        {
            GameEvents.onGUIApplicationLauncherReady.Add(InitializeToolbarIcon);
            GameEvents.onGameSceneLoadRequested.Add(TerminateToolbarIcon);
            GameEvents.onVesselChange.Add(OnVesselChange);
        }

        // If we don't remove our event subscriptions we end up getting a growing amount of toolbar icons
        // Might be because the previous objects haven't gone through the garbage collector yet
        public void OnDestroy()
        {
            GameEvents.onGUIApplicationLauncherReady.Remove(InitializeToolbarIcon);
            GameEvents.onGameSceneLoadRequested.Remove(TerminateToolbarIcon);
            GameEvents.onVesselChange.Remove(OnVesselChange);
        }

        //
        // Generic redrawing code
        //
        private void ClearPlannnerBox()
        {
            planner_box.children[0].uiItem.gameObject.DestroyGameObjectImmediate();
            planner_box.children.Clear();
        }

        //
        // Updating to another settings page
        //
        private void OnButtonClick_Execution()
        {
            ClearPlannnerBox();
            planner_box.children.Add(execution_page);
            GUISupport.ForceGUIUpdate(planner_box, execution_page);
        }

        private void OnButtonClick_Planning()
        {
            ClearPlannnerBox();
            planner_box.children.Add(planning_page);
            GUISupport.ForceGUIUpdate(planner_box, planning_page);
        }

        private void OnButtonClick_Settings()
        {
            ClearPlannnerBox();
            planner_box.children.Add(settings_page);
            GUISupport.ForceGUIUpdate(planner_box, settings_page);
        }

        //
        // Execution page support code
        //

        private void OnButtonClick_WarpToManeuver()
        {
            DataServices.WarpToManeuver();
        }

        private string EngineIgnitionCutoffString()
        {
            if (DataServices.IsEngineDeltaTimeIgnition())
                return ignition_delta_time_name_string;
            else
                return cutoff_delta_time_name_string;
        }

        //
        // Planning page support code
        //

        // The following 3 functions may look like code duplication, the reason they exists is because C#
        // does not allow constructing lambda expressions with reference variables
        private DialogGUIBase AddVelocityTangentToManeuver()
        {
            int index = DataServices.GetLastManeuverIndex();
            return new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                new DialogGUILabel(velocity_tangent_string, velocity_name_string_length),
                new DialogGUILabel(() => { return string.Format(velocity_value_string, DataServices.GetManeuverDeltaVelocityTangent(index)); }, velocity_value_string_length),
                new DialogGUIButton("-1000", () => { DataServices.SetManeuverDeltaVelocityTangent(DataServices.GetManeuverDeltaVelocityTangent(index) - 1000.0); }, false),
                new DialogGUIButton("-100", () => { DataServices.SetManeuverDeltaVelocityTangent(DataServices.GetManeuverDeltaVelocityTangent(index) - 100.0); }, false),
                new DialogGUIButton("-10", () => { DataServices.SetManeuverDeltaVelocityTangent(DataServices.GetManeuverDeltaVelocityTangent(index) - 10.0); }, false),
                new DialogGUIButton("-1", () => { DataServices.SetManeuverDeltaVelocityTangent(DataServices.GetManeuverDeltaVelocityTangent(index) - 1.0); }, false),
                new DialogGUIButton("-0.1", () => { DataServices.SetManeuverDeltaVelocityTangent(DataServices.GetManeuverDeltaVelocityTangent(index) - 0.1); }, false),
                new DialogGUIButton("0", () => { DataServices.SetManeuverDeltaVelocityTangent(0.0); }, false),
                new DialogGUIButton("+0.1", () => { DataServices.SetManeuverDeltaVelocityTangent(DataServices.GetManeuverDeltaVelocityTangent(index) + 0.1); }, false),
                new DialogGUIButton("+1", () => { DataServices.SetManeuverDeltaVelocityTangent(DataServices.GetManeuverDeltaVelocityTangent(index) + 1.0); }, false),
                new DialogGUIButton("+10", () => { DataServices.SetManeuverDeltaVelocityTangent(DataServices.GetManeuverDeltaVelocityTangent(index) + 10.0); }, false),
                new DialogGUIButton("+100", () => { DataServices.SetManeuverDeltaVelocityTangent(DataServices.GetManeuverDeltaVelocityTangent(index) + 100.0); }, false),
                new DialogGUIButton("+1000", () => { DataServices.SetManeuverDeltaVelocityTangent(DataServices.GetManeuverDeltaVelocityTangent(index) + 1000.0); }, false));
        }

        private DialogGUIBase AddVelocityNormalToManeuver()
        {
            int index = DataServices.GetLastManeuverIndex();
            return new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                new DialogGUILabel(velocity_normal_string, velocity_name_string_length),
                new DialogGUILabel(() => { return string.Format(velocity_value_string, DataServices.GetManeuverDeltaVelocityNormal(index)); }, velocity_value_string_length),
                new DialogGUIButton("-1000", () => { DataServices.SetManeuverDeltaVelocityNormal(DataServices.GetManeuverDeltaVelocityNormal(index) - 1000.0); }, false),
                new DialogGUIButton("-100", () => { DataServices.SetManeuverDeltaVelocityNormal(DataServices.GetManeuverDeltaVelocityNormal(index) - 100.0); }, false),
                new DialogGUIButton("-10", () => { DataServices.SetManeuverDeltaVelocityNormal(DataServices.GetManeuverDeltaVelocityNormal(index) - 10.0); }, false),
                new DialogGUIButton("-1", () => { DataServices.SetManeuverDeltaVelocityNormal(DataServices.GetManeuverDeltaVelocityNormal(index) - 1.0); }, false),
                new DialogGUIButton("-0.1", () => { DataServices.SetManeuverDeltaVelocityNormal(DataServices.GetManeuverDeltaVelocityNormal(index) - 0.1); }, false),
                new DialogGUIButton("0", () => { DataServices.SetManeuverDeltaVelocityNormal(0.0); }, false),
                new DialogGUIButton("+0.1", () => { DataServices.SetManeuverDeltaVelocityNormal(DataServices.GetManeuverDeltaVelocityNormal(index) + 0.1); }, false),
                new DialogGUIButton("+1", () => { DataServices.SetManeuverDeltaVelocityNormal(DataServices.GetManeuverDeltaVelocityNormal(index) + 1.0); }, false),
                new DialogGUIButton("+10", () => { DataServices.SetManeuverDeltaVelocityNormal(DataServices.GetManeuverDeltaVelocityNormal(index) + 10.0); }, false),
                new DialogGUIButton("+100", () => { DataServices.SetManeuverDeltaVelocityNormal(DataServices.GetManeuverDeltaVelocityNormal(index) + 100.0); }, false),
                new DialogGUIButton("+1000", () => { DataServices.SetManeuverDeltaVelocityNormal(DataServices.GetManeuverDeltaVelocityNormal(index) + 1000.0); }, false));
        }

        private DialogGUIBase AddVelocityBinormalToManeuver()
        {
            int index = DataServices.GetLastManeuverIndex();
            return new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                new DialogGUILabel(velocity_binormal_string, velocity_name_string_length),
                new DialogGUILabel(() => { return string.Format(velocity_value_string, DataServices.GetManeuverDeltaVelocityBinormal(index)); }, velocity_value_string_length),
                new DialogGUIButton("-1000", () => { DataServices.SetManeuverDeltaVelocityBinormal(DataServices.GetManeuverDeltaVelocityBinormal(index) - 1000.0); }, false),
                new DialogGUIButton("-100", () => { DataServices.SetManeuverDeltaVelocityBinormal(DataServices.GetManeuverDeltaVelocityBinormal(index) - 100.0); }, false),
                new DialogGUIButton("-10", () => { DataServices.SetManeuverDeltaVelocityBinormal(DataServices.GetManeuverDeltaVelocityBinormal(index) - 10.0); }, false),
                new DialogGUIButton("-1", () => { DataServices.SetManeuverDeltaVelocityBinormal(DataServices.GetManeuverDeltaVelocityBinormal(index) - 1.0); }, false),
                new DialogGUIButton("-0.1", () => { DataServices.SetManeuverDeltaVelocityBinormal(DataServices.GetManeuverDeltaVelocityBinormal(index) - 0.1); }, false),
                new DialogGUIButton("0", () => { DataServices.SetManeuverDeltaVelocityBinormal(0.0); }, false),
                new DialogGUIButton("+0.1", () => { DataServices.SetManeuverDeltaVelocityBinormal(DataServices.GetManeuverDeltaVelocityBinormal(index) + 0.1); }, false),
                new DialogGUIButton("+1", () => { DataServices.SetManeuverDeltaVelocityBinormal(DataServices.GetManeuverDeltaVelocityBinormal(index) + 1.0); }, false),
                new DialogGUIButton("+10", () => { DataServices.SetManeuverDeltaVelocityBinormal(DataServices.GetManeuverDeltaVelocityBinormal(index) + 10.0); }, false),
                new DialogGUIButton("+100", () => { DataServices.SetManeuverDeltaVelocityBinormal(DataServices.GetManeuverDeltaVelocityBinormal(index) + 100.0); }, false),
                new DialogGUIButton("+1000", () => { DataServices.SetManeuverDeltaVelocityBinormal(DataServices.GetManeuverDeltaVelocityBinormal(index) + 1000.0); }, false));
        }
        
        // TODO: check out how to deal with 6 hour days (Kerbin) vs 24 hour days (Earth)
        private DialogGUIBase AddDeltaTimeToManeuver()
        {
            int index = DataServices.GetLastManeuverIndex();
            return new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                new DialogGUILabel(time_name_string, time_name_string_length),
                new DialogGUILabel(() => { return GUISupport.FormatTimeSpan(TimeSpan.FromSeconds(DataServices.GetManeuverTime(index))); }, time_value_string_length_two_lines),
                // User needs to move forward in time in big steps, but backwards is only
                // for finetuning where a burn needs to be, hence the assymetric number of buttons
                new DialogGUIButton("-1H", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) - 1*3600); }, false),
                new DialogGUIButton("-10M", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) - 10*60); }, false),
                new DialogGUIButton("-1M", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) - 1*60); }, false),
                new DialogGUIButton("-10S", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) - 10); }, false),
                new DialogGUIButton("-1S", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) - 1); }, false),
                new DialogGUIButton("DEF", () => { DataServices.SetManeuverDeltaTime(DataServices.DEFAULT_MANEUVER_DELTA_TIME); }, false),
                new DialogGUIButton("+1S", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) + 1); }, false),
                new DialogGUIButton("+10S", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) + 10); }, false),
                new DialogGUIButton("+1M", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) + 1*60); }, false),
                new DialogGUIButton("+10M", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) + 10*60); }, false),
                new DialogGUIButton("+1H", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) + 1*3600); }, false),
                new DialogGUIButton("+6H", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) + 6*3600); }, false),
                new DialogGUIButton("+1D", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) + 24*3600); }, false),
                new DialogGUIButton("+10D", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) + 10*24*3600); }, false),
                new DialogGUIButton("+100D", () => { DataServices.SetManeuverDeltaTime(DataServices.GetManeuverDeltaTime(index) + 100*24*3600); }, false));
        }

        private float GetBurnMode()
        {
            int index = DataServices.GetLastManeuverIndex();
            switch (DataServices.GetBurnMode(index)) {
                case DataServices.BurnMode.Engine:
                    return 0f;
                case DataServices.BurnMode.RCS:
                    return 1f;
                case DataServices.BurnMode.Instant:
                default:
                    return 2f;
            }
        }

        private void SetBurnMode(float value)
        {
            if (value > 1.5f) { DataServices.SetBurnMode(DataServices.BurnMode.Instant); }
            else if (value > 0.5f) { DataServices.SetBurnMode(DataServices.BurnMode.RCS); }
            else { DataServices.SetBurnMode(DataServices.BurnMode.Engine); }
        }

        private DialogGUIBase CreateNonMutableManeuver(int index)
        {
            return new DialogGUIVerticalLayout(true, true, 0, new RectOffset(), TextAnchor.MiddleCenter,
                new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                    new DialogGUILabel(velocity_tangent_string, velocity_name_string_length),
                    new DialogGUILabel(() => { return string.Format(velocity_value_string, DataServices.GetManeuverDeltaVelocityTangent(index)); }, velocity_value_string_length),
                    new DialogGUILabel(velocity_normal_string, velocity_name_string_length),
                    new DialogGUILabel(() => { return string.Format(velocity_value_string, DataServices.GetManeuverDeltaVelocityNormal(index)); }, velocity_value_string_length),
                    new DialogGUILabel(velocity_binormal_string, velocity_name_string_length),
                    new DialogGUILabel(() => { return string.Format(velocity_value_string, DataServices.GetManeuverDeltaVelocityBinormal(index)); }, velocity_value_string_length),
                    new DialogGUILabel(() => { return burn_delta_velocity_prefix_string + string.Format(burn_delta_velocity_string, DataServices.GetBurnDeltaVelocity(index)); }, burn_delta_velocity_prefix_string_length + burn_delta_velocity_string_length)
                ),
                new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                    new DialogGUILabel(time_name_string, time_name_string_length),
                    new DialogGUILabel(() => { return GUISupport.FormatTimeSpan(TimeSpan.FromSeconds(DataServices.GetManeuverTime(index))); }, time_value_string_length_single_line),
                    new DialogGUILabel(delta_time_name_string, delta_time_name_string_length),
                    // User should see the time relative to now, so he/she can assess (roughly) how long from now the maneuver needs to be executed
                    // Even more detailed information about this, including taking into account burn times only belong in the execution tab of the planner
                    new DialogGUILabel(() => { return GUISupport.FormatTimeSpan(TimeSpan.FromSeconds(DataServices.GetManeuverDeltaTime(index))); }, time_value_string_length_single_line),
                    new DialogGUILabel(DataServices.GetManeuverIntertiallyFixed(index) ? inertial_frame_string : frenet_frame_string, inertial_or_frenet_frame_string_length),
                    new DialogGUILabel(() => { return burn_mode_prefix_string + DataServices.GetBurnModeString(index); }, burn_delta_velocity_prefix_string_length + burn_mode_string_length),
                    new DialogGUILabel(() => { return burn_time_prefix_string + string.Format(burn_time_string, DataServices.GetBurnTime(index)); }, burn_time_prefix_string_length + burn_time_string_length)
                ),
                new DialogGUISpace(5.0f));
        }
        
        private DialogGUIBase CreateMutableManeuver()
        {
            int index = DataServices.GetLastManeuverIndex();
            return new DialogGUIVerticalLayout(true, true, 0, new RectOffset(), TextAnchor.MiddleCenter,
                AddVelocityTangentToManeuver(),
                AddVelocityNormalToManeuver(),
                AddVelocityBinormalToManeuver(),
                AddDeltaTimeToManeuver(),
                new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                    new DialogGUIToggle(DataServices.GetManeuverIntertiallyFixed(index), inertially_fixed_burn_frame_string, (value) => { DataServices.SetManeuverIntertiallyFixed(value); }, inertially_fixed_burn_frame_string_length),
                    new DialogGUISlider(GetBurnMode, 0f, 2f, true, -1, -1, SetBurnMode),
                    new DialogGUILabel(() => { return burn_mode_prefix_string + DataServices.GetBurnModeString(index); }, burn_mode_prefix_string_length + burn_mode_string_length),
                    new DialogGUILabel(() => { return burn_delta_velocity_prefix_string + string.Format(burn_delta_velocity_string, DataServices.GetBurnDeltaVelocity(index)); }, burn_delta_velocity_prefix_string_length + burn_delta_velocity_string_length),
                    new DialogGUILabel(() => { return burn_time_prefix_string + string.Format(burn_time_string, DataServices.GetBurnTime(index)); }, burn_time_prefix_string_length + burn_time_string_length)
                )
            );
        }

        private void OnButtonClick_AddManeuver()
        {
            int burn_index_to_make_non_mutable = DataServices.GetLastManeuverIndex();

            if (DataServices.AddManeuver())
            {
                List<DialogGUIBase> rows = planning_page.children;

                // If needed make the maneuver non-mutable
                int pre_number_of_rows = rows.Count;
                DeleteLastGUIManeuver(planning_page);
                int post_number_of_row = rows.Count;

                if (pre_number_of_rows > post_number_of_row)
                {
                    DialogGUIBase maneuver_non_mutable = CreateNonMutableManeuver(burn_index_to_make_non_mutable);
                    AddGUIManouver(planning_page, maneuver_non_mutable);
                }

                DialogGUIBase maneuver = CreateMutableManeuver();
                AddGUIManouver(planning_page, maneuver);
            }
        }

        private void OnButtonClick_DeleteManeuver()
        {
            DeleteLastGUIManeuver(planning_page);
            DataServices.RemoveLastManeuver();

            // From here on pure-GUI adjustments
            List<DialogGUIBase> rows = planning_page.children;
            // If needed make the new last maneuver mutable
            int pre_number_of_rows = rows.Count;
            DeleteLastGUIManeuver(planning_page);
            int post_number_of_row = rows.Count;

            if (pre_number_of_rows > post_number_of_row)
            {
                DialogGUIBase maneuver = CreateMutableManeuver();
                AddGUIManouver(planning_page, maneuver);
            }
        }

        private void AddGUIManouver(DialogGUIBase parent, DialogGUIBase maneuver)
        {
            List<DialogGUIBase> rows = parent.children;
            maneuver.SetOptionText(magic_maneuver_string);
            rows.Add(maneuver);
            GUISupport.ForceGUIUpdate(planning_page, maneuver);
        }

        private void DeleteLastGUIManeuver(DialogGUIBase parent)
        {
            List<DialogGUIBase> rows = parent.children;
            if (rows.Count > 0) {
                DialogGUIBase child = rows.ElementAt(rows.Count - 1);
                if (child.OptionText == magic_maneuver_string) {
                    rows.RemoveAt(rows.Count - 1);
                    if (child.uiItem && child.uiItem.gameObject)
                    {
                        child.uiItem.gameObject.DestroyGameObjectImmediate(); // ensure memory gets freed
                    }
                }
            }
        }

        //
        // Settings page support code
        //

        private DialogGUIBase AddFlightPlanLength()
        {
            return new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                new DialogGUILabel(plan_length_time_name_string, plan_length_time_name_string_length),
                new DialogGUILabel(() => { return GUISupport.FormatTimeSpan(TimeSpan.FromSeconds(DataServices.GetPlanTimeLength())); }, time_value_string_length_single_line),
                new DialogGUIButton("-1H", () => { DataServices.SetPlanTimeLength(DataServices.GetPlanTimeLength() - 1*3600); }, false),
                new DialogGUIButton("-10M", () => { DataServices.SetPlanTimeLength(DataServices.GetPlanTimeLength() - 10*60); }, false),
                new DialogGUIButton("-1M", () => { DataServices.SetPlanTimeLength(DataServices.GetPlanTimeLength() - 1*60); }, false),
                new DialogGUIButton("-10S", () => { DataServices.SetPlanTimeLength(DataServices.GetPlanTimeLength() - 10); }, false),
                new DialogGUIButton("-1S", () => { DataServices.SetPlanTimeLength(DataServices.GetPlanTimeLength() - 1); }, false),
                new DialogGUIButton("DEF", () => { DataServices.SetPlanTimeLength(DataServices.PLAN_TIME_LENGTH_DEFAULT); }, false),
                new DialogGUIButton("+1S", () => { DataServices.SetPlanTimeLength(DataServices.GetPlanTimeLength() + 1); }, false),
                new DialogGUIButton("+10S", () => { DataServices.SetPlanTimeLength(DataServices.GetPlanTimeLength() + 10); }, false),
                new DialogGUIButton("+1M", () => { DataServices.SetPlanTimeLength(DataServices.GetPlanTimeLength() + 1*60); }, false),
                new DialogGUIButton("+10M", () => { DataServices.SetPlanTimeLength(DataServices.GetPlanTimeLength() + 10*60); }, false),
                new DialogGUIButton("+1H", () => { DataServices.SetPlanTimeLength(DataServices.GetPlanTimeLength() + 1*3600); }, false)
            );
        }

        private DialogGUIBase AddMaxStepsPerSegment()
        {
            return new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                new DialogGUILabel(max_steps_per_segment_name_string, max_steps_per_segment_name_string_length),
                new DialogGUILabel(() => { return string.Format(max_steps_per_segment_value_string, DataServices.GetPlanMaxStepsPerSegment()); }, max_steps_per_segment_value_string_length),
                new DialogGUIButton("-1000", () => { DataServices.SetPlanMaxStepsPerSegment(DataServices.GetPlanMaxStepsPerSegment() - 1000); }, false),
                new DialogGUIButton("-100", () => { DataServices.SetPlanMaxStepsPerSegment(DataServices.GetPlanMaxStepsPerSegment() - 100); }, false),
                new DialogGUIButton("-10", () => { DataServices.SetPlanMaxStepsPerSegment(DataServices.GetPlanMaxStepsPerSegment() - 10); }, false),
                new DialogGUIButton("-1", () => { DataServices.SetPlanMaxStepsPerSegment(DataServices.GetPlanMaxStepsPerSegment() - 1); }, false),
                new DialogGUIButton("DEF", () => { DataServices.SetPlanMaxStepsPerSegment(DataServices.PLAN_MAX_STEPS_PER_SEGMENT_DEFAULT); }, false),
                new DialogGUIButton("+1", () => { DataServices.SetPlanMaxStepsPerSegment(DataServices.GetPlanMaxStepsPerSegment() + 1); }, false),
                new DialogGUIButton("+10", () => { DataServices.SetPlanMaxStepsPerSegment(DataServices.GetPlanMaxStepsPerSegment() + 10); }, false),
                new DialogGUIButton("+100", () => { DataServices.SetPlanMaxStepsPerSegment(DataServices.GetPlanMaxStepsPerSegment() + 100); }, false),
                new DialogGUIButton("+1000", () => { DataServices.SetPlanMaxStepsPerSegment(DataServices.GetPlanMaxStepsPerSegment() + 1000); }, false)
            );
        }

        private DialogGUIBase AddTolerance()
        {
            return new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                new DialogGUILabel(tolerance_name_string, tolerance_name_string_length),
                new DialogGUILabel(() => { return string.Format(tolerance_value_string, DataServices.GetPlanTolerance()); }, tolerance_value_string_length),
                new DialogGUIButton("-10", () => { DataServices.SetPlanTolerance(DataServices.GetPlanTolerance() - 10); }, false),
                new DialogGUIButton("-1", () => { DataServices.SetPlanTolerance(DataServices.GetPlanTolerance() - 1); }, false),
                new DialogGUIButton("-0.1", () => { DataServices.SetPlanTolerance(DataServices.GetPlanTolerance() - 0.1); }, false),
                new DialogGUIButton("-0.01", () => { DataServices.SetPlanTolerance(DataServices.GetPlanTolerance() - 0.01); }, false),
                new DialogGUIButton("DEF", () => { DataServices.SetPlanTolerance(DataServices.PLAN_TOLERANCE_DEFAULT); }, false),
                new DialogGUIButton("+0.01", () => { DataServices.SetPlanTolerance(DataServices.GetPlanTolerance() + 0.01); }, false),
                new DialogGUIButton("+0.1", () => { DataServices.SetPlanTolerance(DataServices.GetPlanTolerance() + 0.1); }, false),
                new DialogGUIButton("+1", () => { DataServices.SetPlanTolerance(DataServices.GetPlanTolerance() + 1); }, false),
                new DialogGUIButton("+10", () => { DataServices.SetPlanTolerance(DataServices.GetPlanTolerance() + 10); }, false)
            );
        }

        //
        // GUI initialization
        //
        private void InitializePlannerGUI()
        {
            execution_page = new DialogGUIVerticalLayout(true, true, 0, new RectOffset(), TextAnchor.UpperCenter,
                new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                    // Note: haven't figured a way to reliably detect if we can actually show on navball, so always add the UI element
                    new DialogGUIToggle(DataServices.GetShowOnNavball(), show_on_navball_string, (value) => { DataServices.SetShowOnNavball(value); }, show_on_navball_string_length),
                    new DialogGUIButton("Warp to Maneuver", OnButtonClick_WarpToManeuver, button_width, button_height, false)),
                new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                    new DialogGUILabel(() => { return EngineIgnitionCutoffString() + GUISupport.FormatTimeSpan(TimeSpan.FromSeconds(DataServices.GetEngineDeltaTime())); }, ignition_delta_time_name_string_length + time_value_string_length_single_line),
                    new DialogGUILabel(() => { return total_delta_velocity_prefix_string + string.Format(burn_delta_velocity_string, DataServices.GetDeltaVelocityOfAllBurns()); }, total_delta_velocity_prefix_string_length + burn_delta_velocity_string_length)
                )
            );

            planning_page = new DialogGUIVerticalLayout(true, true, 0, new RectOffset(), TextAnchor.UpperCenter,
                new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                    new DialogGUIButton("Add Maneuver", OnButtonClick_AddManeuver, button_width, button_height, false),
                    new DialogGUIButton("Delete Maneuver", OnButtonClick_DeleteManeuver, button_width, button_height, false)
                )
            );

            settings_page = new DialogGUIVerticalLayout(true, true, 0, new RectOffset(), TextAnchor.UpperCenter,
                    AddFlightPlanLength(),
                    AddMaxStepsPerSegment(),
                    AddTolerance()
            );

            // Do not use a DialogGUIBox for this, it will not respect automatic resizing
            planner_box = new DialogGUIHorizontalLayout(true, true, 0, new RectOffset(), TextAnchor.UpperCenter, execution_page);

            multi_page_planner = new MultiOptionDialog(
                "PrincipiaPlannerGUI",
                "",
                "Principia Planner",
                HighLogic.UISkin,
                new Rect(x_pos, y_pos, 700.0f, 50.0f), // for reasons beyond me we have to get the width correct ourselves
                new DialogGUIBase[]
                {
                    // buttons to select page
                    new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                        new DialogGUIButton("Execution", OnButtonClick_Execution, button_width, button_height, false),
                        new DialogGUIButton("Planning", OnButtonClick_Planning, button_width, button_height, false),
                        new DialogGUIButton("Settings", OnButtonClick_Settings, button_width, button_height, false)),
                    // box that contains the actual content
                    planner_box
                });
        }

        private void ShowPlannerWindow()
        {
            // Bring GUI in sync with whatever is the current vessel, before making the planner window visible
            // This is done to avoid a race condition with the creation of the GUI game objects that we
            // normally use to force refresh the GUI
            RegenerateGUIManeuvers();

            planner_window_visible = true;
            popup_dialog = PopupDialog.SpawnPopupDialog(new Vector2(x_pos, y_pos), new Vector2(x_pos, y_pos),
                                                        multi_page_planner, true, HighLogic.UISkin, false);
        }

        private void HidePlannerWindow()
        {
            planner_window_visible = false;
            if (popup_dialog) {
                popup_dialog.Dismiss();
            }
        }

        private void OnVesselChange(Vessel value)
        {
            DataServices.UpdateSettingsUponVesselChange();
            RegenerateGUIManeuvers();
        }

        private void RegenerateGUIManeuvers()
        {
            List<DialogGUIBase> rows = planning_page.children;
            int num_rows = rows.Count;

            for (int i = 0; i < num_rows; i++)
            {
                DeleteLastGUIManeuver(planning_page);
            }

            int max_index = DataServices.GetLastManeuverIndex();
            for (int i = 0; i < max_index; i++)
            {
                DialogGUIBase maneuver_non_mutable = CreateNonMutableManeuver(i);
                AddGUIManouver(planning_page, maneuver_non_mutable);
            }
            if (max_index >= 0)
            {
                DialogGUIBase maneuver = CreateMutableManeuver();
                AddGUIManouver(planning_page, maneuver);
            }
        }

        private void InitializeToolbarIcon()
        {
            if (toolbar_button == null) {
                UnityEngine.Texture toolbar_button_texture = GUISupport.LoadTextureOrDie("toolbar_button.png");
                toolbar_button = KSP.UI.Screens.ApplicationLauncher.Instance.AddModApplication(
                      onTrue          : ShowPlannerWindow,
                      onFalse         : HidePlannerWindow,
                      onHover         : null,
                      onHoverOut      : null,
                      onEnable        : null,
                      onDisable       : null,
                      visibleInScenes : KSP.UI.Screens.ApplicationLauncher.AppScenes.FLIGHT,
                      texture         : toolbar_button_texture);
            }
        }

        private void TerminateToolbarIcon(GameScenes scenes) {
            if (toolbar_button != null) {
                KSP.UI.Screens.ApplicationLauncher.Instance.RemoveModApplication(toolbar_button);
                toolbar_button = null;
            }
        }

        public void Start()
        {
            if (planner_window_visible) {
                ShowPlannerWindow();
            }
        }
    }
}  // namespace ksp_plugin_adapter
}  // namespace principia