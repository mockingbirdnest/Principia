namespace principia {
namespace ksp_plugin_adapter {

internal static class Style {
  public static UnityEngine.Color Tangent { get; } = XKCDColors.NeonYellow;
  public static UnityEngine.Color Normal { get; } = XKCDColors.AquaBlue;
  public static UnityEngine.Color Binormal { get; } = XKCDColors.PurplePink;
  
  private static UnityEngine.Texture2D lit_toggle_button_on_active_;
  private static UnityEngine.Texture2D lit_toggle_button_on_hover_;
  private static UnityEngine.Texture2D lit_toggle_button_on_normal_;

  private static UnityEngine.Texture2D MakeGreen(UnityEngine.Texture2D source) {
    var render = UnityEngine.RenderTexture.GetTemporary(
        source.width,
        source.height,
        depthBuffer: 0,
        UnityEngine.RenderTextureFormat.Default,
        UnityEngine.RenderTextureReadWrite.Linear);
    UnityEngine.Graphics.Blit(source, render);
    UnityEngine.RenderTexture old_active = UnityEngine.RenderTexture.active;
    UnityEngine.RenderTexture.active = render;
    UnityEngine.Texture2D result = new UnityEngine.Texture2D(source.width, source.height);
    result.ReadPixels(new UnityEngine.Rect(0, 0, render.width, render.height), 0, 0);
    result.Apply();
    UnityEngine.RenderTexture.active = old_active;
    UnityEngine.RenderTexture.ReleaseTemporary(render);

    UnityEngine.Color32[] pixels = result.GetPixels32();
    for (int i = 0; i < pixels.Length; ++i) {
      pixels[i].r = 0;
      pixels[i].b = 0;
    }
    result.SetPixels32(pixels);
    result.Apply();
    return result;
  }

  public static UnityEngine.GUIStyle LitToggleButton() {
    var style = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.button);
    style.active.textColor = ultra_cool_grey_;
    style.hover.textColor = ultra_cool_grey_;
    style.normal.textColor = ultra_cool_grey_;
    style.onActive.textColor = XKCDColors.Green;
    style.onHover.textColor = XKCDColors.Green;
    style.onNormal.textColor = XKCDColors.Green;/*
    if (lit_toggle_button_on_active_ == null) {
      lit_toggle_button_on_active_ = MakeGreen(style.onActive.background);
      lit_toggle_button_on_hover_ = MakeGreen(style.onHover.background);
      lit_toggle_button_on_normal_ = MakeGreen(style.onNormal.background);
    }
    style.onActive.background = lit_toggle_button_on_active_;
    style.onHover.background = lit_toggle_button_on_hover_;
    style.onNormal.background = lit_toggle_button_on_normal_;*/
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
    var horizontal_line_style =
        new UnityEngine.GUIStyle(UnityEngine.GUI.skin.horizontalSlider);
    horizontal_line_style.fixedHeight /= 5;
    horizontal_line_style.normal.background = ultra_cool_grey_texture;
    UnityEngine.GUILayout.Label("", horizontal_line_style);
  }

  public static void LineSpacing() {
    var horizontal_line_style =
        new UnityEngine.GUIStyle(UnityEngine.GUI.skin.label);
    horizontal_line_style.fixedHeight /= 5;
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
  private static readonly UnityEngine.Color ultra_cool_grey_ =
      new UnityEngine.Color(0.64f, 0.64f, 0.64f);
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
