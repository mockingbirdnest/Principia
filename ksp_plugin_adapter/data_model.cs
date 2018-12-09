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
        
    }
}  // namespace ksp_plugin_adapter
}  // namespace principia