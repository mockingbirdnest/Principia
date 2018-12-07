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

namespace principia {
namespace ksp_plugin_adapter {

    [KSPAddon(KSPAddon.Startup.FlightAndKSC, false)]
    public sealed class SettingsGUI: MonoBehaviour
    {
        private const float button_width = 50.0f;
        private const float button_height = 25.0f;

        bool settings_window_visible = false;

        private DialogGUIVerticalLayout main_settings_page;
        private DialogGUIVerticalLayout plotting_frame_page;
        private DialogGUIVerticalLayout logging_settings_page;

        private DialogGUIBox settings_box;
        private MultiOptionDialog multi_page_settings;
        private PopupDialog popup_dialog;

        private KSP.UI.Screens.ApplicationLauncherButton toolbar_button;

        //
        // Generic redrawing code
        //
        private void ClearSettingsPage()
        {
            settings_box.children[0].uiItem.gameObject.DestroyGameObjectImmediate();
            settings_box.children.Clear();
        }

        private void ForceGUIUpdate()
        {
            Stack<Transform> stack = new Stack<Transform>();
            stack.Push(settings_box.uiItem.gameObject.transform);
            settings_box.children[0].Create(ref stack, HighLogic.UISkin);
        }

        //
        // Updating to another settings page
        //
        private void OnButtonClick_MainSettings()
        {
            ClearSettingsPage();
            settings_box.children.Add(main_settings_page);
            ForceGUIUpdate();
        }

        private void OnButtonClick_PlottingFrameSettings()
        {
            ClearSettingsPage();
            settings_box.children.Add(plotting_frame_page);
            ForceGUIUpdate();
        }
        
        private void OnButtonClick_LoggingSettings()
        {
            ClearSettingsPage();
            settings_box.children.Add(logging_settings_page);
            ForceGUIUpdate();
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

        //
        // History length
        //
        private int history_magnitude = 20; // (1 << index) is the history time in seconds, with the exception of 30, which is +infinity
        private float GetHistoryMagnitude() { return (float)history_magnitude; }
        private void SetHistoryMagnitude(float value) { history_magnitude = (int)value; }
        private double GetHistoryLength()
        {
            if (history_magnitude == 30) {
                return double.PositiveInfinity;
            } else {
                return (1 << history_magnitude);
            }
        }

        //
        // Prediction
        //
        private int prediction_tolerance_magnitude = -2;
        private int prediction_step_magnitude = 8;
        private float GetPredictionToleranceMagnitude() { return (float)prediction_tolerance_magnitude; }
        private void SetPredictionToleranceMagnitude(float value) { prediction_tolerance_magnitude = (int)value; }
        private double GetPredictionTolerance() { return 10^prediction_tolerance_magnitude; }
        private float GetPredictionStepMagnitude() { return (float)prediction_step_magnitude; }
        private void SetPredictionStepMagnitude(float value) { prediction_step_magnitude = (int)value; }
        private double GetPredictionStep() { return (1 << prediction_step_magnitude); }

        //
        // KSP settings
        //
        private bool display_patched_conics = false;
        private bool display_solar_flare = false;
        private bool GetPatchedConics() { return display_patched_conics; }
        private void SetPatchedConics(bool value) { display_patched_conics = value; }
        private bool GetSolarFlare() { return display_solar_flare; }
        private void SetSolarFlare(bool value) { display_solar_flare = value; }

        public SettingsGUI()
        {
            InitializeSettingsGUI();

            GameEvents.onGUIApplicationLauncherReady.Add(onAppLauncherReady);
            GameEvents.onGUIApplicationLauncherDestroyed.Add(onAppLauncherDestroyed);
        }
        
        private void InitializeSettingsGUI()
        {
            main_settings_page = new DialogGUIVerticalLayout(false, true, 0, new RectOffset(), TextAnchor.UpperCenter,
                // Version
                new DialogGUIHorizontalLayout(TextAnchor.MiddleLeft,
                    new DialogGUILabel("Version: " + GetVersion())
                ),
                new DialogGUISpace(20.0f),
                // History length
                new DialogGUIHorizontalLayout(TextAnchor.MiddleLeft,
                    new DialogGUILabel(() => { return "Maximum history length: " + string.Format("{0:E2}", GetHistoryLength()) + " s"; })
                ),
                new DialogGUIHorizontalLayout(
                    new DialogGUISlider(GetHistoryMagnitude, 10f, 30f, true, -1, -1, SetHistoryMagnitude)
                ),
                new DialogGUISpace(20.0f),
                // Prediction settings
                new DialogGUIHorizontalLayout(TextAnchor.MiddleCenter,
                    new DialogGUILabel("Prediction settings")
                ),
                new DialogGUISpace(10.0f),
                new DialogGUIHorizontalLayout(TextAnchor.MiddleLeft,
                    new DialogGUILabel(() => { return "Tolerance: " + string.Format("{0:E2}", GetPredictionTolerance()) + " m"; })
                ),
                new DialogGUIHorizontalLayout(
                    new DialogGUISlider(GetPredictionToleranceMagnitude, -3f, 4f, true, -1, -1, SetPredictionToleranceMagnitude)
                ),
                new DialogGUIHorizontalLayout(TextAnchor.MiddleLeft,
                    new DialogGUILabel(() => { return "Steps: " + string.Format("{0:E2}", GetPredictionStep()); })
                ),
                new DialogGUIHorizontalLayout(
                    new DialogGUISlider(GetPredictionStepMagnitude, 2f, 24f, true, -1, -1, SetPredictionStepMagnitude)
                ),
                new DialogGUISpace(20.0f),
                // KSP settings
                new DialogGUIHorizontalLayout(TextAnchor.MiddleCenter,
                    new DialogGUILabel("KSP settings")
                ),
                new DialogGUISpace(10.0f),
                new DialogGUIHorizontalLayout(
                    new DialogGUIToggle(GetPatchedConics, "Display patched conics (not intended for flight planning)", SetPatchedConics)
                ),
                new DialogGUIHorizontalLayout(
                    new DialogGUIToggle(GetSolarFlare, "Enable system-star lens flare", SetSolarFlare)
                ),
                new DialogGUISpace(100.0f)
            );
            plotting_frame_page = new DialogGUIVerticalLayout(false, true, 0, new RectOffset(), TextAnchor.UpperCenter,
                new DialogGUILabel("testing 456")
            );
            logging_settings_page = new DialogGUIVerticalLayout(false, true, 0, new RectOffset(), TextAnchor.UpperCenter,
                new DialogGUILabel("testing 789")
            );

            settings_box = new DialogGUIBox(null, -1, -1, () => true, main_settings_page);

            multi_page_settings = new MultiOptionDialog(
                "PrincipiaSettingsGUI",
                "",
                "Principia Settings",
                HighLogic.UISkin,
                new Rect(0.5f, 0.5f, 450.0f, 450.0f),
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

        // Returns false and nulls |texture| if the file does not exist.
        private bool LoadTextureIfExists(out UnityEngine.Texture texture,
                                         String path) {
            string full_path =
                KSPUtil.ApplicationRootPath + Path.DirectorySeparatorChar +
                "GameData" + Path.DirectorySeparatorChar +
                "Principia" + Path.DirectorySeparatorChar +
                "assets" + Path.DirectorySeparatorChar +
                path;
            if (File.Exists(full_path)) {
                var texture2d = new UnityEngine.Texture2D(2, 2);
#if KSP_VERSION_1_5_1
                bool success = UnityEngine.ImageConversion.LoadImage(
                    texture2d, File.ReadAllBytes(full_path));
#elif KSP_VERSION_1_3_1
                bool success = texture2d.LoadImage(
                    File.ReadAllBytes(full_path));
#endif
                if (!success) {
                    Log.Fatal("Failed to load texture " + full_path);
                }
                texture = texture2d;
                return true;
            } else {
                texture = null;
                return false;
            }
        }
                
        private UnityEngine.Texture LoadTextureOrDie(String path) {
            UnityEngine.Texture texture;
            bool success = LoadTextureIfExists(out texture, path);
            if (!success) {
                Log.Fatal("Missing texture " + path);
            }
            return texture;
        }        

        private void ShowSettingsWindow()
        {
            settings_window_visible = true;
            popup_dialog = PopupDialog.SpawnPopupDialog(new Vector2(0.5f, 0.5f), new Vector2(0.5f, 0.5f), 
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
                UnityEngine.Texture toolbar_button_texture = LoadTextureOrDie("toolbar_button.png");
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
        private void onAppLauncherReady() {
            InitializeToolbarIcon();
        }
        
        private void onAppLauncherDestroyed() {
            KSP.UI.Screens.ApplicationLauncher.Instance.RemoveModApplication(toolbar_button);
            toolbar_button = null;
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