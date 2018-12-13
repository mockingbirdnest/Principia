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
    public sealed class GUISupport {
        // In some rare cases, such as when using a scrolling layout
        // simply re-connecting a previously rendered layout causes
        // the layout to break. For this case we recreate the whole
        // layout.
        public static void RecursivelyDeleteLayout(DialogGUIBase layout)
        {
            RecursivelyDeleteLayoutChildren(layout);

            if (layout.uiItem && layout.uiItem.gameObject)
            {
                layout.uiItem.gameObject.DestroyGameObjectImmediate();
            }
        }

        public static void RecursivelyDeleteLayoutChildren(DialogGUIBase layout)
        {
            if (layout.children.Count < 1) {
                return;
            }
            for (int i = layout.children.Count - 1; i >= 0; i--)
            {
                RecursivelyDeleteLayoutChildren(layout.children[i]);

                DialogGUIBase child = layout.children[i];
                layout.children.RemoveAt(i);
                if (child.uiItem && child.uiItem.gameObject)
                {
                    child.uiItem.gameObject.DestroyGameObjectImmediate();
                }
            }
        }

        public static void ForceGUIUpdate(DialogGUIBase parent, DialogGUIBase child)
        {
            Stack<Transform> stack = new Stack<Transform>(); // some data on hierarchy of the GUI components
            if (parent.uiItem && parent.uiItem.gameObject)
            {
                stack.Push(parent.uiItem.gameObject.transform); // need the reference point of the parent GUI component for position and size
                child.Create(ref stack, HighLogic.UISkin); // required to force the GUI creation
            }
        }

        // Returns false and nulls |texture| if the file does not exist.
        private static bool LoadTextureIfExists(out UnityEngine.Texture texture,
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

        public static UnityEngine.Texture LoadTextureOrDie(String path) {
            UnityEngine.Texture texture;
            bool success = LoadTextureIfExists(out texture, path);
            if (!success) {
                Log.Fatal("Missing texture " + path);
            }
            return texture;
        }

        // TODO: check if this is the correct way of printing time
        private static string FormatPositiveTimeSpan (TimeSpan span) {
            return (GameSettings.KERBIN_TIME
                ? (span.Days * 4 + span.Hours / 6).ToString("0000;0000") +
                      " d6 " + (span.Hours % 6).ToString("0;0") + " h "
                : span.Days.ToString("000;000") + " d " +
                      span.Hours.ToString("00;00") + " h ") +
            span.Minutes.ToString("00;00") + " min " +
            (span.Seconds + span.Milliseconds / 1000m).ToString("00.0;00.0") +
            " s";
        }

        public static string FormatTimeSpan (TimeSpan span) {
            return span.Ticks.ToString("+;-") + FormatPositiveTimeSpan(span);
        }
    }
}  // namespace ksp_plugin_adapter
}  // namespace principia