using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal static class Style {

  public static UnityEngine.Color Tangent { get; } = XKCDColors.NeonYellow;
  public static UnityEngine.Color Normal { get; } = XKCDColors.AquaBlue;
  public static UnityEngine.Color Binormal { get; } = XKCDColors.PurplePink;

  public static UnityEngine.GUIStyle Info(UnityEngine.GUIStyle style) {
    var info_style = new UnityEngine.GUIStyle(style);
    info_style.focused.textColor = ultra_cool_grey_;
    info_style.normal.textColor = ultra_cool_grey_;
    return info_style;
  }

  public static UnityEngine.GUIStyle Warning(UnityEngine.GUIStyle style) {
    var warning_style = new UnityEngine.GUIStyle(style);
    warning_style.focused.textColor = XKCDColors.Orange;
    warning_style.normal.textColor = XKCDColors.Orange;
    return warning_style;
  }

  public static UnityEngine.GUIStyle Error(UnityEngine.GUIStyle style) {
    var error_style = new UnityEngine.GUIStyle(style);
    error_style.focused.textColor = XKCDColors.Red;
    error_style.normal.textColor = XKCDColors.Red;
    return error_style;
  }

  public static UnityEngine.GUIStyle RightAligned(UnityEngine.GUIStyle style) {
    var right_aligned_style = new UnityEngine.GUIStyle(style);
    right_aligned_style.alignment = UnityEngine.TextAnchor.MiddleRight;
    return right_aligned_style;
  }

  public static UnityEngine.GUIStyle Multiline(UnityEngine.GUIStyle style) {
    var multiline_style = new UnityEngine.GUIStyle(style);
    multiline_style.alignment = UnityEngine.TextAnchor.UpperLeft;
    multiline_style.fixedHeight = 0;
    multiline_style.wordWrap = true;
    return multiline_style;
  }

  public static void HorizontalLine() {
    var horizontal_line_style = new UnityEngine.GUIStyle
                                    (UnityEngine.GUI.skin.horizontalSlider);
    horizontal_line_style.fixedHeight /= 5;
    horizontal_line_style.normal.background = ultra_cool_grey_texture;
    UnityEngine.GUILayout.Label("", horizontal_line_style);
  }

  private static UnityEngine.Texture2D ultra_cool_grey_texture {
    get {
      var texture = new UnityEngine.Texture2D(width : 4, height : 4);
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          texture.SetPixel(i, j, ultra_cool_grey_);
        }
      }
      return texture;
    }
  }

  // Close to XKCD "cool grey", but hue-neutral
  private static UnityEngine.Color ultra_cool_grey_ =
      new UnityEngine.Color(0.64f, 0.64f, 0.64f);
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
