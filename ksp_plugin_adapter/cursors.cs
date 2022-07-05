using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace principia {
namespace ksp_plugin_adapter {

static class DigitScrollCursor {
  // A hardware cursor will have the size of the system cursor, which is fine for
  // most purposes, but annoying for the digit scroll cursor which should be sized
  // like the digits.
  class TextureSoftwareCursor : Cursors.TextureCursor {
    public override void SetCursor() {
      UnityEngine.Cursor.SetCursor(texture, hotspot,
                                   UnityEngine.CursorMode.ForceSoftware);
    }
  }

  private const string principia_digit_scroll_cursor_name =
      "principia_digit_scroll";

  public static void RequestCursor() {
    digit_scroll_requested = true;
    if (!cursor_is_digit_scroll) {
      cursor_is_digit_scroll = true;
      if (!Cursors.CursorController.Instance.Contains(
          principia_digit_scroll_cursor_name)) {
        // TODO(egg): pick images based on resolution?
        PrincipiaPluginAdapter.LoadTextureOrDie(
            out UnityEngine.Texture digit_scroll_cursor_texture,
            "digit_scroll_cursor.png");
        Cursors.CursorController.Instance.AddCursor(
            principia_digit_scroll_cursor_name,
            new TextureSoftwareCursor{
                texture = (UnityEngine.Texture2D)digit_scroll_cursor_texture,
                hotspot = new UnityEngine.Vector2(4, 15)});
      }
      Cursors.CursorController.Instance.ChangeCursor(
          principia_digit_scroll_cursor_name);
    }
  }

  public static void Update() {
    if (cursor_is_digit_scroll && !digit_scroll_requested) {
      cursor_is_digit_scroll = false;
      Cursors.CursorController.Instance.ForceDefaultCursor();
    }
    digit_scroll_requested = false;
  }
  
  private static bool digit_scroll_requested = false;
  private static bool cursor_is_digit_scroll = false;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
