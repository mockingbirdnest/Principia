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
    info_style.focused.textColor = XKCDColors.LightGrey;
    info_style.normal.textColor = XKCDColors.LightGrey;
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

  public static void HorizontalLine() {
    var horizontal_line_style = new UnityEngine.GUIStyle
                                    (UnityEngine.GUI.skin.horizontalSlider);
    horizontal_line_style.fixedHeight /= 5;
    UnityEngine.GUILayout.Label("", horizontal_line_style);
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
