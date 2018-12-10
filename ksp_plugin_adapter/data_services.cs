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

namespace principia {
namespace ksp_plugin_adapter {

    // For the sake of not having to know the GUI complexity when working with the controllers
    // that actually do stuff with data, we isolate the data into a seperate class
    // It's intended to be as simple as possible, so that at any time both the GUI and control
    // developers can easily read this
    public sealed class DataServices
    {
        //
        // Global settings
        //

        //
        // History length
        //
        private static int history_magnitude = 20; // (1 << index) is the history time in seconds, with the exception of 30, which is +infinity
        public static float GetHistoryMagnitude() { return (float)history_magnitude; }
        public static void SetHistoryMagnitude(float value) { history_magnitude = (int)value; }
        public static double GetHistoryLength()
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
        private static int prediction_tolerance_magnitude = -2;
        private static int prediction_step_magnitude = 8;
        public static float GetPredictionToleranceMagnitude() { return (float)prediction_tolerance_magnitude; }
        public static void SetPredictionToleranceMagnitude(float value) { prediction_tolerance_magnitude = (int)value; }
        public static double GetPredictionTolerance() { return Math.Pow(10, prediction_tolerance_magnitude); }
        public static float GetPredictionStepMagnitude() { return (float)prediction_step_magnitude; }
        public static void SetPredictionStepMagnitude(float value) { prediction_step_magnitude = (int)value; }
        public static int GetPredictionStep() { return (1 << prediction_step_magnitude); }

        //
        // KSP settings
        //
        private static bool display_patched_conics = false;
        public static bool GetPatchedConicsEnabled() { return display_patched_conics; }
        public static void SetPatchedConicsEnabled(bool value) { display_patched_conics = value; }
        // For legacy reasons this setting is directly configured on an KSP object
        public static bool GetSolarFlareEnabled() { return Sun.Instance.sunFlare.enabled; }
        public static void SetSolarFlareEnabled(bool value) { Sun.Instance.sunFlare.enabled = value; }

        //
        // Reference frame settings
        //
        private static CelestialBody selected_celestial_body = FlightGlobals.GetHomeBody();
        // I suspect the hardcoded numbers are to allow the C++ code to understand this as well
        public enum FrameType {
            BARYCENTRIC_ROTATING = 6001,
            BODY_CENTRED_NON_ROTATING = 6000,
            BODY_CENTRED_PARENT_DIRECTION = 6002,
            BODY_SURFACE = 6003
        }
        private static FrameType reference_frame = FrameType.BODY_CENTRED_NON_ROTATING;
        // internal is inherited from NavigationFrameParameters
        internal delegate void NavigationFrameParametersCallback(NavigationFrameParameters frame_parameters);
        private static NavigationFrameParametersCallback on_update_celestial_body_or_reference_frame;

        public static CelestialBody GetSelectedCelestialBody() { return selected_celestial_body; }
        public static void SetSelectedCelestialBody(CelestialBody value)
        {
            CelestialBody prev = selected_celestial_body;
            selected_celestial_body = value;
            if (prev != selected_celestial_body)
            {
                UpdatePluginWithCelestialBodyAndReferenceFrame();
            }
        }
        public static FrameType GetReferenceFrame() { return reference_frame; }
        public static void SetReferenceFrame(FrameType value)
        {
            FrameType prev = reference_frame;
            reference_frame = value;
            // We cannot support frames that require 2 bodies, if we selected the root body
            // (and thus only have 1 body to work with)
            if (selected_celestial_body.is_root() &&
               (reference_frame == FrameType.BARYCENTRIC_ROTATING ||
                reference_frame == FrameType.BODY_CENTRED_PARENT_DIRECTION)) {
                reference_frame = FrameType.BODY_CENTRED_NON_ROTATING;
            }
            if (prev != reference_frame)
            {
                UpdatePluginWithCelestialBodyAndReferenceFrame();
            }
        }

        // internal is inherited from NavigationFrameParameters
        internal static void InitializeSelectedCelestialBodyAndReferenceFrame(NavigationFrameParametersCallback callback)
        {
            on_update_celestial_body_or_reference_frame = callback;
            reference_frame = FrameType.BODY_CENTRED_NON_ROTATING;
            selected_celestial_body = FlightGlobals.currentMainBody ?? FlightGlobals.GetHomeBody();
        }
        
        private static NavigationFrameParameters GenerateNavigationFrameParameters()
        {
            switch (reference_frame) {
                case FrameType.BODY_CENTRED_NON_ROTATING:
                case FrameType.BODY_SURFACE:
                    return new NavigationFrameParameters{
                        extension = (int)reference_frame,
                        centre_index = selected_celestial_body.flightGlobalsIndex};
                case FrameType.BARYCENTRIC_ROTATING:
                    return new NavigationFrameParameters{
                        extension = (int)reference_frame,
                        primary_index = selected_celestial_body.referenceBody.flightGlobalsIndex,
                        secondary_index = selected_celestial_body.flightGlobalsIndex};
                case FrameType.BODY_CENTRED_PARENT_DIRECTION:
                    // We put the primary body as secondary, because the one we want fixed
                    // is the secondary body (which means it has to be the primary in the
                    // terminology of |BodyCentredBodyDirection|).
                    return new NavigationFrameParameters{
                        extension = (int)reference_frame,
                        primary_index = selected_celestial_body.flightGlobalsIndex,
                        secondary_index = selected_celestial_body.referenceBody.flightGlobalsIndex};
                default:
                    throw Log.Fatal("Unexpected reference_frame " + reference_frame.ToString());
            }
        }

        private static void UpdatePluginWithCelestialBodyAndReferenceFrame()
        {
            NavigationFrameParameters frame_parameters = GenerateNavigationFrameParameters();
            on_update_celestial_body_or_reference_frame(frame_parameters);
        }

        //
        // Logging settings
        //
        private static bool record_journal_at_next_startup = false;
        private static bool record_journal_in_progress = false;
        private static int verbose_level = 0;
        private static int supressed_logging_level = 0;
        private static int stderr_logging_level = 2;
        private static int flush_logging_level = 0;

        public static bool GetRecordJournalInProgress() { return record_journal_in_progress; }
        public static bool GetRecordJournalAtNextStartup() { return record_journal_at_next_startup; }
        public static void SetRecordJournalAtNextStartup(bool value) { record_journal_at_next_startup = value; }
        
        public static void InitializeJournaling()
        {
            if (record_journal_at_next_startup) {
                Log.ActivateRecorder(true);
                record_journal_in_progress = true;
            }
        }

        // Implicit requirement: Don't accept a log setting, until the C++ side of principia has accepted it
        public static int GetVerboseLevel() { return verbose_level; }
        public static void SetVerboseLevel(int value) { Log.SetVerboseLogging(value); verbose_level = Log.GetVerboseLogging(); }
        public static int GetLogLevel() { return supressed_logging_level; }
        public static void SetLogLevel(int value) { Log.SetSuppressedLogging(value); supressed_logging_level = Log.GetSuppressedLogging(); }
        public static int GetStderrLevel() { return stderr_logging_level; }
        public static void SetStderrLevel(int value) { Log.SetStderrLogging(value); stderr_logging_level = Log.GetStderrLogging(); }
        public static int GetFlushLevel() { return flush_logging_level; }
        public static void SetFlushLevel(int value) { Log.SetBufferedLogging(value); flush_logging_level = Log.GetBufferedLogging(); }

        public static void InitializeLoggingSettings()
        {
            SetVerboseLevel(verbose_level);
            SetLogLevel(supressed_logging_level);
            SetStderrLevel(stderr_logging_level);
            SetFlushLevel(flush_logging_level);
        }
    }
}  // namespace ksp_plugin_adapter
}  // namespace principia