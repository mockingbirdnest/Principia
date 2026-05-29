using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace principia {
namespace ksp_plugin_adapter {

internal class Graph : ScalingRenderer {
  public Graph(int width, int height) {
    texture_ = new UnityEngine.Texture2D(width, height);
    all_black_ = new UnityEngine.Color[width * height];
    for (int i = 0; i < all_black_.Length; ++i) {
      all_black_[i] = XKCDColors.Black;
    }
  }

  public void PrepareCanvas(Interval x_range, Interval y_range) {
    x_range_ = x_range;
    y_range_ = y_range;
    dirty_ = true;
    texture_.SetPixels(all_black_);
    labels_.Clear();
  }

  public void PlotFunction(Func<double, double> f,
                           Interval x_subrange,
                           UnityEngine.Color colour) {
    for (int i = AbscissaToPixel(x_subrange.min);
         i <= AbscissaToPixel(x_subrange.max);
         ++i) {
      // Honest plotting assuming f is monotone.
      Interval pixel_range = PixelToAbscissa(i).IntersectedWith(x_subrange);
      double f_x_min = f.Invoke(pixel_range.min);
      double f_x_max = f.Invoke(pixel_range.max);
      double f_min;
      double f_max;
      if (f_x_min <= f_x_max) {
        f_min = f_x_min;
        f_max = f_x_max;
      } else {
        f_min = f_x_max;
        f_max = f_x_min;
      }
      for (int j = OrdinateToPixel(f_min); j <= OrdinateToPixel(f_max); ++j) {
        texture_.SetPixel(i, j, colour);
      }
    }
  }

  public void PlotVerticalLine(double x,
                               UnityEngine.Color colour,
                               Interval? y_range = null) {
    if (!x_range_.Contains(x)) {
      return;
    }
    for (int j = y_range == null ? 0 : OrdinateToPixel(y_range.Value.min);
         j <
         (y_range == null
              ? texture_.height
              : OrdinateToPixel(y_range.Value.max));
         ++j) {
      texture_.SetPixel(AbscissaToPixel(x), j, colour);
    }
  }

  public void PlotHorizontalLine(double y, UnityEngine.Color colour) {
    if (!y_range_.Contains(y)) {
      return;
    }
    for (int i = 0; i < texture_.width; ++i) {
      texture_.SetPixel(i, OrdinateToPixel(y), colour);
    }
  }

  public void PlotPoint(double x, double y, UnityEngine.Color colour) {
    dirty_ = true;
    texture_.SetPixel(AbscissaToPixel(x), OrdinateToPixel(y), colour);
  }

  public void AddLabel(double x,
                       double y,
                       string text,
                       UnityEngine.Color colour,
                       UnityEngine.TextAnchor anchor) {
    labels_.Add(new Label{
        x_pixels = AbscissaToPixel(x), y_pixels = OrdinateToPixel(y),
        text = text, colour = colour, anchor = anchor,
    });
  }

  public void Render() {
    if (dirty_) {
      texture_.Apply();
      dirty_ = false;
    }
    UnityEngine.GUILayout.Box("",
                              UnityEngine.GUILayout.Width(texture_.width),
                              UnityEngine.GUILayout.Height(texture_.height));
    if (UnityEngine.Event.current.type == UnityEngine.EventType.Repaint) {
      var graph_rectangle = UnityEngine.GUILayoutUtility.GetLastRect();
      UnityEngine.GUI.DrawTexture(graph_rectangle, texture_);
      foreach (var label in labels_) {
        var label_rectangle =
            new UnityEngine.Rect(graph_rectangle.xMin + label.x_pixels,
                                 graph_rectangle.yMax - label.y_pixels,
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

  private int AbscissaToPixel(double x) {
    return (int)(texture_.width * (x - x_range_.min) / x_range_.measure);
  }

  private Interval PixelToAbscissa(int i) {
    return new Interval{
        min = i * x_range_.measure / texture_.width + x_range_.min,
        max = (i + 1) * x_range_.measure / texture_.width + x_range_.min
    };
  }

  private int OrdinateToPixel(double y) {
    return (int)(texture_.height * (y - y_range_.min) / y_range_.measure);
  }

  private struct Label {
    public int x_pixels;
    public int y_pixels;
    public string text;
    public UnityEngine.Color colour;
    public UnityEngine.TextAnchor anchor;
  };

  private Interval x_range_;
  private Interval y_range_;
  private bool dirty_;
  private List<Label> labels_ = new List<Label>();
  
  private UnityEngine.Texture2D texture_;
  private UnityEngine.Color[] all_black_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
