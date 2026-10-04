using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace principia {
namespace ksp_plugin_adapter {

internal class Graph : ScalingRenderer {
  public unsafe Graph(int width, int height) {
    texture_ = new UnityEngine.Texture2D(width, height, UnityEngine.TextureFormat.RGBA32, mipChain: false);
    all_black_ = new byte[width * height * sizeof(UnityEngine.Color32)];
    const int offsetof_a = 3;
    for (int i = 0; i < width * height; ++i) {
      all_black_[sizeof(UnityEngine.Color32) * i + offsetof_a] = 255;
    }
  }

  public void LoadPixels(IntPtr rgba32_begin) {
    if (rgba32_begin == IntPtr.Zero) {
      texture_.LoadRawTextureData(all_black_);
    } else {
      texture_.LoadRawTextureData(rgba32_begin, texture_.width * texture_.height * 4);
    }
    texture_.Apply();
    labels_.Clear();
  }

  public void AddLabel(GraphLabel label) {
    labels_.Add(new Label{
        x = label.x, y = label.y, text = label.text,
        colour = Interface.FromRGBA(label.colour),
        anchor = (UnityEngine.TextAnchor)label.anchor
    });
  }

  public void Render() {
    UnityEngine.GUILayout.Box("",
                              UnityEngine.GUILayout.Width(texture_.width),
                              UnityEngine.GUILayout.Height(texture_.height));
    if (UnityEngine.Event.current.type == UnityEngine.EventType.Repaint) {
      var graph_rectangle = UnityEngine.GUILayoutUtility.GetLastRect();
      graph_rectangle = new UnityEngine.Rect(graph_rectangle.xMin,
                                             graph_rectangle.yMax,
                                             graph_rectangle.width,
                                             -graph_rectangle.height);
      UnityEngine.GUI.DrawTexture(graph_rectangle, texture_);
      foreach (var label in labels_) {
        var label_rectangle =
            new UnityEngine.Rect(graph_rectangle.xMin + label.x,
                                 graph_rectangle.yMax - label.y,
                                 Width(2),
                                 Height(1));
        switch (label.anchor) {
          case UnityEngine.TextAnchor.UpperLeft:
          case UnityEngine.TextAnchor.UpperCenter:
          case UnityEngine.TextAnchor.UpperRight:
            break;
          case UnityEngine.TextAnchor.MiddleLeft:
          case UnityEngine.TextAnchor.MiddleCenter:
          case UnityEngine.TextAnchor.MiddleRight:
            label_rectangle.y -= label_rectangle.height / 2;
            break;
          case UnityEngine.TextAnchor.LowerLeft:
          case UnityEngine.TextAnchor.LowerCenter:
          case UnityEngine.TextAnchor.LowerRight:
            label_rectangle.y -= label_rectangle.height;
            break;
        }
        switch (label.anchor) {
          case UnityEngine.TextAnchor.UpperLeft:
          case UnityEngine.TextAnchor.MiddleLeft:
          case UnityEngine.TextAnchor.LowerLeft:
            break;
          case UnityEngine.TextAnchor.UpperCenter:
          case UnityEngine.TextAnchor.MiddleCenter:
          case UnityEngine.TextAnchor.LowerCenter:
            label_rectangle.x -= label_rectangle.width / 2;
            break;
          case UnityEngine.TextAnchor.UpperRight:
          case UnityEngine.TextAnchor.MiddleRight:
          case UnityEngine.TextAnchor.LowerRight:
            label_rectangle.x -= label_rectangle.width;
            break;
        }
        UnityEngine.GUI.Label(label_rectangle,
                              label.text,
                              new UnityEngine.GUIStyle(
                                  UnityEngine.GUI.skin.label) {
                                  focused = {
                                      textColor = label.colour
                                  },
                                  normal = {
                                      textColor = label.colour
                                  },
                                  alignment = label.anchor
                              });
      }
    }
  }

  private struct Label {
    public long x;
    public long y;
    public string text;
    public UnityEngine.Color colour;
    public UnityEngine.TextAnchor anchor;
  };

  private readonly List<Label> labels_ = new List<Label>();
  
  private readonly UnityEngine.Texture2D texture_;
  private readonly byte[] all_black_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
