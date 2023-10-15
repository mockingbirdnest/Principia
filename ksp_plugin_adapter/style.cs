namespace principia {
namespace ksp_plugin_adapter {

internal static class Style {
  public static UnityEngine.Color Tangent { get; } = XKCDColors.NeonYellow;
  public static UnityEngine.Color Normal { get; } = XKCDColors.AquaBlue;
  public static UnityEngine.Color Binormal { get; } = XKCDColors.PurplePink;

  public static UnityEngine.Color ManœuvreMarkerBase { get; } = XKCDColors.LightGrey;

  // As used by stock map node markers.
  public static UnityEngine.Color MarkerCaption { get; }
    = new UnityEngine.Color(191f / 255f, 1f, 0f, 1f);
  public static UnityEngine.Color MarkerCaptionPinned { get; }
    = new UnityEngine.Color(191f / 255f, 1f, 0f, 0.6f);

  public static UnityEngine.GUIStyle DarkToggleButton() {
    var style = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.button){
        active = {
            textColor = ultra_cool_grey_
        },
        hover = {
            textColor = ultra_cool_grey_
        },
        normal = {
            textColor = ultra_cool_grey_
        }
    };
    return style;
  }

  public static UnityEngine.GUIStyle LitToggleButton() {
    var style = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.button){
        onActive = {
            textColor = XKCDColors.Spearmint
        },
        onHover = {
            textColor = XKCDColors.Spearmint
        },
        onNormal = {
            textColor = XKCDColors.Spearmint
        }
    };
    style.active = style.onActive;
    style.hover = style.onHover;
    style.normal = style.onNormal;
    return style;
  }

  public static UnityEngine.GUIStyle Info(UnityEngine.GUIStyle style) {
    var info_style = new UnityEngine.GUIStyle(style){
        focused = {
            textColor = ultra_cool_grey_
        },
        normal = {
            textColor = ultra_cool_grey_
        }
    };
    return info_style;
  }

  public static UnityEngine.GUIStyle Warning(UnityEngine.GUIStyle style) {
    var warning_style = new UnityEngine.GUIStyle(style){
        focused = {
            textColor = XKCDColors.Orange
        },
        normal = {
            textColor = XKCDColors.Orange
        }
    };
    return warning_style;
  }

  public static UnityEngine.GUIStyle Error(UnityEngine.GUIStyle style) {
    var error_style = new UnityEngine.GUIStyle(style){
        focused = {
            textColor = XKCDColors.Red
        },
        normal = {
            textColor = XKCDColors.Red
        }
    };
    return error_style;
  }

  public static UnityEngine.GUIStyle Aligned(UnityEngine.TextAnchor alignment,
                                             UnityEngine.GUIStyle style) {
    var aligned_style = new UnityEngine.GUIStyle(style){
        alignment = alignment
    };
    return aligned_style;
  }

  public static UnityEngine.GUIStyle MiddleLeftAligned(
      UnityEngine.GUIStyle style,
      float height) {
    return new UnityEngine.GUIStyle(style){
        alignment = UnityEngine.TextAnchor.MiddleLeft,
        fixedHeight = height
    };
  }

  public static UnityEngine.GUIStyle RightAligned(UnityEngine.GUIStyle style) {
    return Aligned(UnityEngine.TextAnchor.MiddleRight, style);
  }

  public static UnityEngine.GUIStyle Multiline(UnityEngine.GUIStyle style) {
    var multiline_style = new UnityEngine.GUIStyle(style){
        alignment = UnityEngine.TextAnchor.UpperLeft,
        fixedHeight = 0,
        wordWrap = true
    };
    return multiline_style;
  }

  public static void HorizontalLine() {
    if (horizontal_line_style_ == null) {
      horizontal_line_style_ =
          new UnityEngine.GUIStyle(UnityEngine.GUI.skin.horizontalSlider);
      horizontal_line_style_.fixedHeight /= 5;
      horizontal_line_style_.normal.background = ultra_cool_grey_texture;
    }
    UnityEngine.GUILayout.Label("", horizontal_line_style_);
  }

  public static void LineSpacing() {
    if (line_spacing_style_ == null) {
      line_spacing_style_ = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.label);
      line_spacing_style_.fixedHeight /= 5;
    }
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
  private static UnityEngine.GUIStyle horizontal_line_style_;
  private static UnityEngine.GUIStyle line_spacing_style_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
