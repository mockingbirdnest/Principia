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
using System.IO;
using UnityEngine;
using UnityEngine.UI;

namespace principia {
namespace ksp_plugin_adapter {

    [KSPAddon(KSPAddon.Startup.AllGameScenes, false)]
    public sealed class SettingsGUI: MonoBehaviour
    {
        private const float button_width = 50.0f;
        private const float button_height = 25.0f;

        // Plotting frame defines
        private const string plotting_frame_body_name_string = "<color=#ffffffff>Selected celestial body for plotting frame: </color>";
        private const float plotting_frame_body_name_string_length = 200f;
        private const float plotting_frame_body_value_string_length = 60f;

        private const float plotting_frame_string_length = 300f;

        // Logging settings defines
        private const string verbose_level_name_string = "<color=#ffffffff>Verbose level: </color>";
        private const float verbose_level_name_string_length = 60f;
        private const float verbose_level_value_string_length = 60f;

        private const string log_level_name_string = "<color=#ffffffff>Log level: </color>";
        private const string stderr_level_name_string = "<color=#ffffffff>Stderr level: </color>";
        private const string flush_level_name_string = "<color=#ffffffff>Flush level: </color>";
        private const float log_level_name_string_length = 60f;
        private const float log_level_value_string_length = 60f;

        private const string record_journal_at_next_startup_name_string = "<color=#ffffffff>Record journal (enabling requires restart)</color>";
        private const float record_journal_at_next_startup_name_string_length = 350f;

        private const string record_journal_in_progress_name_string = "Journaling is ON";
        private const string record_journal_not_in_progress_name_string = "Journaling is OFF";
        private const float record_journal_in_progress_name_string_length = 100f;

        bool settings_window_visible = false;
        private float x_pos = 0.95f;
        private float y_pos = 0.90f;

        private DialogGUIBase main_settings_page;
        private DialogGUIBase plotting_frame_page;
        private DialogGUIBase logging_settings_page;

        private DialogGUIBase settings_box;
        private MultiOptionDialog multi_page_settings;
        private PopupDialog popup_dialog;

        private KSP.UI.Screens.ApplicationLauncherButton toolbar_button;

        public SettingsGUI()
        {
            InitializeSettingsGUI();
        }

        public void Awake()
        {
            GameEvents.onGUIApplicationLauncherReady.Add(InitializeToolbarIcon);
            GameEvents.onGameSceneLoadRequested.Add(TerminateToolbarIcon);
        }

        // If we don't remove our event subscriptions we end up getting a growing amount of toolbar icons
        // Might be because the previous objects haven't gone through the garbage collector yet
        public void OnDestroy()
        {
            GameEvents.onGUIApplicationLauncherReady.Remove(InitializeToolbarIcon);
            GameEvents.onGameSceneLoadRequested.Remove(TerminateToolbarIcon);
        }

        //
        // Generic redrawing code
        //
        private void ClearSettingsBox()
        {
            DialogGUIBase child = settings_box.children[0];
            settings_box.children.RemoveAt(0);
            child.uiItem.gameObject.DestroyGameObjectImmediate();
        }

        //
        // Updating to another settings page
        //
        private void OnButtonClick_MainSettings()
        {
            ClearSettingsBox();
            settings_box.children.Add(main_settings_page);
            GUISupport.ForceGUIUpdate(settings_box, settings_box.children[0]);
        }

        private void OnButtonClick_PlottingFrameSettings()
        {
            DialogGUIBase backup_reference = plotting_frame_page;

            ClearSettingsBox();
            // Force the complete re-rending of the plotting frame settings
            // otherwise the scroll-bar breaks. And the scroll bar is a hard
            // requirement to make the entire list fit the window.
            plotting_frame_page = AddPlottingFrameSelectionUI();
            settings_box.children.Add(plotting_frame_page);
            GUISupport.ForceGUIUpdate(settings_box, settings_box.children[0]);

            // Delete old GUI elements
            // This is done afterwards, because for reasons I do not know
            // The new GUI must be fully updated before we can do this
            GUISupport.RecursivelyDeleteLayout(backup_reference);
        }
        
        private void OnButtonClick_LoggingSettings()
        {
            ClearSettingsBox();
            settings_box.children.Add(logging_settings_page);
            GUISupport.ForceGUIUpdate(settings_box, settings_box.children[0]);
        }
        
        //
        // Data required for main settings page
        //
        private String GetVersion()
        {
            String version;
            String unused_build_date;
            Interface.GetVersion(build_date: out unused_build_date,
                                    version: out version);
            return version;
        }
        
        private DialogGUIBase AddMainSettingsUI()
        {
            return new DialogGUIVerticalLayout(false, true, 0, new RectOffset(), TextAnchor.UpperCenter,
                // Version
                new DialogGUIHorizontalLayout(TextAnchor.MiddleLeft,
                    new DialogGUILabel("Version: " + GetVersion())
                ),
                // History length
                new DialogGUIHorizontalLayout(TextAnchor.MiddleLeft,
                    new DialogGUILabel(() => { return "Maximum history length: " + string.Format("{0:E2}", DataServices.GetHistoryLength()) + " s"; }),
                    new DialogGUISlider(DataServices.GetHistoryMagnitude, 10f, 30f, true, -1, -1, DataServices.SetHistoryMagnitude)
                ),
                // Prediction settings
                new DialogGUIHorizontalLayout(TextAnchor.MiddleCenter,
                    new DialogGUILabel("<color=#ffffffff>Prediction settings</color>")
                ),
                new DialogGUIHorizontalLayout(TextAnchor.MiddleLeft,
                    new DialogGUILabel(() => { return "Tolerance: " + string.Format("{0:E2}", DataServices.GetPredictionTolerance()) + " m"; }),
                    new DialogGUISlider(DataServices.GetPredictionToleranceMagnitude, -3f, 4f, true, -1, -1, DataServices.SetPredictionToleranceMagnitude)
                ),
                new DialogGUIHorizontalLayout(TextAnchor.MiddleLeft,
                    new DialogGUILabel(() => { return "Steps: " + string.Format("{0:E2}", DataServices.GetPredictionStep()); }),
                    new DialogGUISlider(DataServices.GetPredictionStepMagnitude, 2f, 24f, true, -1, -1, DataServices.SetPredictionStepMagnitude)
                ),
                // KSP settings
                new DialogGUIHorizontalLayout(TextAnchor.MiddleCenter,
                    new DialogGUILabel("<color=#ffffffff>KSP settings</color>")
                ),
                new DialogGUIHorizontalLayout(
                    new DialogGUIToggle(DataServices.GetPatchedConicsEnabled, "Display patched conics (not intended for flight planning)", DataServices.SetPatchedConicsEnabled)
                ),
                new DialogGUIHorizontalLayout(
                    new DialogGUIToggle(DataServices.GetSolarFlareEnabled, "Enable system-star lens flare", DataServices.SetSolarFlareEnabled)
                )
            );
        }

        //
        // Support code for the plotting frame selection
        //
        private float GetFrameType()
        {
            switch (DataServices.GetReferenceFrame())
            {
                case DataServices.FrameType.BODY_SURFACE:
                    return 0f;
                case DataServices.FrameType.BODY_CENTRED_NON_ROTATING:
                    return 1f;
                case DataServices.FrameType.BARYCENTRIC_ROTATING:
                    return 2f;
                case DataServices.FrameType.BODY_CENTRED_PARENT_DIRECTION:
                default:
                    return 3f;
            }
        }

        private void SetFrameType(float value)
        {
            DataServices.FrameType reference_frame;
            if (value > 2.5f) {
                reference_frame = DataServices.FrameType.BODY_CENTRED_PARENT_DIRECTION;
            } else if (value > 1.5f) {
                reference_frame = DataServices.FrameType.BARYCENTRIC_ROTATING;
            } else if (value > 0.5f) {
                reference_frame = DataServices.FrameType.BODY_CENTRED_NON_ROTATING;
            } else {
                reference_frame = DataServices.FrameType.BODY_SURFACE;
            }
            DataServices.SetReferenceFrame(reference_frame);
        }

        private string GetFrameTypeString()
        {
            CelestialBody selected_celestial_body = DataServices.GetSelectedCelestialBody();
            CelestialBody parent = selected_celestial_body.referenceBody;
            switch (DataServices.GetReferenceFrame())
            {
                case DataServices.FrameType.BODY_SURFACE:
                    return string.Format("Reference frame fixing the surface of {0}", selected_celestial_body.name);
                case DataServices.FrameType.BODY_CENTRED_NON_ROTATING:
                    return string.Format("Non-rotating reference frame fixing the center of {0}", selected_celestial_body.name);
                case DataServices.FrameType.BARYCENTRIC_ROTATING:
                    return string.Format("Reference frame fixing the barycenter of {0} and {1}, the plane which they move about the barycenter, and the line between them", selected_celestial_body.name, parent.name);
                case DataServices.FrameType.BODY_CENTRED_PARENT_DIRECTION:
                default:
                    return string.Format("Reference frame fixing the center of {0}, the plane of its orbit around {1}, and the line between them", selected_celestial_body.name, parent.name);
            }
        }

        private void AddOrbitingBodies(CelestialBody body, ref DialogGUILayoutBase gui, string parent_name)
        {
            if (parent_name != null) {
                gui.children.Add(new DialogGUIButton(body.name + String.Format(" (Parent body: {0})", parent_name), () => { DataServices.SetSelectedCelestialBody(body); }, false));
            } else {
                gui.children.Add(new DialogGUIButton(body.name, () => { DataServices.SetSelectedCelestialBody(body); }, false));
            }
            foreach (CelestialBody child_body in body.orbitingBodies)
            {
                AddOrbitingBodies(child_body, ref gui, body.name);
            }
        }

        private DialogGUIBase AddPlottingFrameSelectionUI()
        {
            CelestialBody root_body = null;
            foreach (CelestialBody body in FlightGlobals.Bodies) {
                if (FlightGlobals.GetBodyIndex(body) == 0)
                {
                    root_body = body; // typically the star
                    break;
                }
            }
            if (root_body == null) {
                throw Log.Fatal("No root body of celestials could be found");
            }

            DialogGUILayoutBase celestial_body_list = new DialogGUIVerticalLayout(true, true, 0, new RectOffset(0, 18, 0, 0) /* prevent scroll bar overlap */, TextAnchor.UpperCenter,
                new DialogGUIContentSizer(ContentSizeFitter.FitMode.Unconstrained, ContentSizeFitter.FitMode.PreferredSize, true));
            AddOrbitingBodies(root_body, ref celestial_body_list, null);

            DialogGUILayoutBase gui = new DialogGUIVerticalLayout(true, true, 0, new RectOffset(), TextAnchor.UpperCenter,
                new DialogGUILabel(() => { return plotting_frame_body_name_string + DataServices.GetSelectedCelestialBody().name; }, plotting_frame_body_name_string_length + plotting_frame_body_value_string_length),
                new DialogGUILabel(GetFrameTypeString, plotting_frame_string_length),
                new DialogGUISlider(GetFrameType, 0f, 3f, true, -1, -1, SetFrameType),
                // there will be too many celestial bodies, putting them inside something that can scroll vertically
                new DialogGUIScrollList(new Vector2(450f, 500f), false, true, celestial_body_list)
            );
            return gui;
        }

        //
        // Support code for logging settings
        //
        private float GetVerboseLevel() { return (float)DataServices.GetVerboseLevel(); }
        private void SetVerboseLevel(float value) { DataServices.SetVerboseLevel((int)value); }

        private float GetLogLevel() { return (float)DataServices.GetLogLevel(); }
        private void SetLogLevel(float value) { DataServices.SetLogLevel((int)value); }
        private float GetStderrLevel() { return (float)DataServices.GetStderrLevel(); }
        private void SetStderrLevel(float value) { DataServices.SetStderrLevel((int)value); }
        private float GetFlushLevel() { return (float)DataServices.GetFlushLevel(); }
        private void SetFlushLevel(float value) { DataServices.SetFlushLevel((int)value); }

        private string GetFlushLevelString(int index)
        {
            if (index < 0)
            {
                return "ALL";
            }
            else
            {
                return Log.severity_names[index];
            }
        }

        private DialogGUIBase AddLoggingSettingsUI()
        {
            return new DialogGUIVerticalLayout(true, true, 0, new RectOffset(), TextAnchor.UpperCenter,
                new DialogGUIHorizontalLayout(
                    new DialogGUISlider(GetVerboseLevel, 0f, 4f, true, -1, -1, SetVerboseLevel),
                    new DialogGUILabel(() => { return verbose_level_name_string + DataServices.GetVerboseLevel(); }, verbose_level_name_string_length + verbose_level_value_string_length)),
                new DialogGUIHorizontalLayout(
                    new DialogGUISlider(GetLogLevel, 0f, 3f, true, -1, -1, SetLogLevel),
                    new DialogGUILabel(() => { return log_level_name_string + Log.severity_names[DataServices.GetLogLevel()]; }, log_level_name_string_length + log_level_value_string_length)),
                new DialogGUIHorizontalLayout(
                    new DialogGUISlider(GetStderrLevel, 0f, 3f, true, -1, -1, SetStderrLevel),
                    new DialogGUILabel(() => { return stderr_level_name_string + Log.severity_names[DataServices.GetStderrLevel()]; }, log_level_name_string_length + log_level_value_string_length)),
                new DialogGUIHorizontalLayout(
                    new DialogGUISlider(GetFlushLevel, -1f, 3f, true, -1, -1, SetFlushLevel),
                    new DialogGUILabel(() => { return flush_level_name_string + GetFlushLevelString(DataServices.GetFlushLevel()); }, log_level_name_string_length + log_level_value_string_length)),
                new DialogGUIHorizontalLayout(
                    new DialogGUIToggle(DataServices.GetRecordJournalAtNextStartup(), record_journal_at_next_startup_name_string, (value) => { DataServices.SetRecordJournalAtNextStartup(value); }, record_journal_at_next_startup_name_string_length),
                    new DialogGUILabel(() => { if (DataServices.GetRecordJournalInProgress()) { return record_journal_in_progress_name_string; } else { return record_journal_not_in_progress_name_string; } }, record_journal_in_progress_name_string_length))
            );
        }

        private void InitializeSettingsGUI()
        {
            main_settings_page = AddMainSettingsUI();
            plotting_frame_page = AddPlottingFrameSelectionUI();
            logging_settings_page = AddLoggingSettingsUI();

            // Do not use a DialogGUIBox for this, it will not respect automatic resizing
            settings_box = new DialogGUIVerticalLayout(true, true, 0, new RectOffset(), TextAnchor.UpperCenter, main_settings_page);

            multi_page_settings = new MultiOptionDialog(
                "PrincipiaSettingsGUI",
                "",
                "Principia Settings",
                HighLogic.UISkin,
                new Rect(x_pos, y_pos, 500.0f, 50.0f),
                new DialogGUIBase[]
                {
                    // buttons to select page
                    new DialogGUIHorizontalLayout(true, false, 0, new RectOffset(), TextAnchor.MiddleCenter,
                        new DialogGUIButton("Main", OnButtonClick_MainSettings, button_width, button_height, false),
                        new DialogGUIButton("Plotting frame", OnButtonClick_PlottingFrameSettings, button_width, button_height, false),
                        new DialogGUIButton("Logging", OnButtonClick_LoggingSettings, button_width, button_height, false)),
                    // box that contains actual settings
                    settings_box
                });
        }

        private void ShowSettingsWindow()
        {
            settings_window_visible = true;
            popup_dialog = PopupDialog.SpawnPopupDialog(new Vector2(x_pos, y_pos), new Vector2(x_pos, y_pos),
                                                        multi_page_settings, true, HighLogic.UISkin, false);
        }

        private void HideSettingsWindow()
        {
            settings_window_visible = false;
            if (popup_dialog) {
                popup_dialog.Dismiss();
            }
        }

        private void InitializeToolbarIcon()
        {
            if (toolbar_button == null) {
                UnityEngine.Texture toolbar_button_texture = GUISupport.LoadTextureOrDie("toolbar_button.png");
                toolbar_button = KSP.UI.Screens.ApplicationLauncher.Instance.AddModApplication(
                      onTrue          : ShowSettingsWindow,
                      onFalse         : HideSettingsWindow,
                      onHover         : null,
                      onHoverOut      : null,
                      onEnable        : null,
                      onDisable       : null,
                      visibleInScenes : KSP.UI.Screens.ApplicationLauncher.AppScenes.SPACECENTER | KSP.UI.Screens.ApplicationLauncher.AppScenes.TRACKSTATION | KSP.UI.Screens.ApplicationLauncher.AppScenes.FLIGHT,
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
            if (settings_window_visible) {
                ShowSettingsWindow();
            }
        }
    }
}  // namespace ksp_plugin_adapter
}  // namespace principia