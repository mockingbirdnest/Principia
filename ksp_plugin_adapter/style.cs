namespace principia {
namespace ksp_plugin_adapter {

internal static class Style {
  static Style() {
    horizontal_line_style_ =
        new UnityEngine.GUIStyle(UnityEngine.GUI.skin.horizontalSlider);
    horizontal_line_style_.fixedHeight /= 5;
    horizontal_line_style_.normal.background = ultra_cool_grey_texture;
    line_spacing_style_ = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.label);
    line_spacing_style_.fixedHeight /= 5;
  }

  public static UnityEngine.Color Tangent { get; } = XKCDColors.NeonYellow;
  public static UnityEngine.Color Normal { get; } = XKCDColors.AquaBlue;
  public static UnityEngine.Color Binormal { get; } = XKCDColors.PurplePink;

  public static UnityEngine.GUIStyle DarkToggleButton() {
    var style = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.button);
    style.active.textColor = ultra_cool_grey_;
    style.hover.textColor = ultra_cool_grey_;
    style.normal.textColor = ultra_cool_grey_;
    return style;
  }

  public static UnityEngine.GUIStyle LitToggleButton() {
    var style = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.button);
    style.onActive.textColor = XKCDColors.Spearmint;
    style.onHover.textColor = XKCDColors.Spearmint;
    style.onNormal.textColor = XKCDColors.Spearmint;
    style.active = style.onActive;
    style.hover = style.onHover;
    style.normal = style.onNormal;
    return style;
  }

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
    UnityEngine.GUILayout.Label("", horizontal_line_style_);
  }

  public static void LineSpacing() {
    UnityEngine.GUILayout.Label("", line_spacing_style_);
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
  private static readonly UnityEngine.Color ultra_cool_grey_ =
      new UnityEngine.Color(0.64f, 0.64f, 0.64f);

  // It is very important to use *static* styles for these elements.  If we use
  // local variables at the place of use, the stupid Unity ends up burning 9 kiB
  // for each horizontal line we display.  See #3064.
  private static readonly UnityEngine.GUIStyle horizontal_line_style_;
  private static readonly UnityEngine.GUIStyle line_spacing_style_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
