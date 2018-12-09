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
    public sealed class DataModel
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
        
    }
}  // namespace ksp_plugin_adapter
}  // namespace principia