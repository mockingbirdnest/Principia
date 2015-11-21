using System;
using System.Globalization;

namespace principia {
namespace ksp_plugin_adapter {

internal class DifferentialSlider {
  private float slider_position_ = 0.0f;
  private DateTime last_time_;

  private readonly string label_;
  private readonly string unit_;
  private readonly CultureInfo culture_;

  private readonly double log10_lower_rate_ = -3;
  private readonly double log10_upper_rate_ = 3.5;
  private readonly double min_value_;
  private readonly double max_value_;

  public double value { get; set; }

  public DifferentialSlider(string label,
                            string unit,
                            double log10_lower_rate,
                            double log10_upper_rate)
      : this(label,
             unit,
             log10_lower_rate,
             log10_upper_rate,
             -double.NegativeInfinity,
             -double.PositiveInfinity) {}

  // Rates are in units of |value| per real-time second.
  public DifferentialSlider(string label,
                            string unit,
                            double log10_lower_rate,
                            double log10_upper_rate,
                            double min_value,
                            double max_value) {
    label_ = label;
    unit_ = unit;
    culture_= new CultureInfo("");
    culture_.NumberFormat.NumberGroupSeparator = "'";
    log10_lower_rate_ = log10_lower_rate;
    log10_upper_rate_ = log10_upper_rate;
    min_value_ = min_value;
    max_value_ = max_value;
  }

  // Renders the |DifferentialSlider|.  Returns true if and only if |value|
  // changed.
  public bool Render() {
    UnityEngine.GUI.skin = null;
    UnityEngine.GUILayout.BeginHorizontal();

    UnityEngine.GUILayout.Label(
        text       : label_,
        options    : UnityEngine.GUILayout.Width(100));

    var old_alignment = UnityEngine.GUI.skin.label.alignment;
    UnityEngine.GUI.skin.label.alignment = UnityEngine.TextAnchor.UpperRight;
    UnityEngine.GUILayout.Label(
        text       : value.ToString("#,0.000", culture_),
        options    : UnityEngine.GUILayout.Width(100));
    UnityEngine.GUI.skin.label.alignment = old_alignment;
    UnityEngine.GUILayout.Label(
        text       : unit_,
        options    : UnityEngine.GUILayout.Width(50));
    if (!UnityEngine.Input.GetMouseButton(0)) {
      slider_position_ = 0;
    }
    UnityEngine.GUI.skin.horizontalSlider.fixedHeight = 21;
    UnityEngine.GUI.skin.horizontalSliderThumb.fixedHeight = 21;
    slider_position_ = UnityEngine.GUILayout.HorizontalSlider(
        value      : slider_position_,
        leftValue  : -1,
        rightValue : 1,
        options    : UnityEngine.GUILayout.ExpandWidth(true));
    bool value_changed = false;

    if (UnityEngine.GUILayout.Button("0", UnityEngine.GUILayout.Width(20))) {
      value_changed = true;
      value = 0;
    }
    if (slider_position_ != 0.0) {
      value_changed = true;
      value += Math.Sign(slider_position_) *
               Math.Pow(10, log10_lower_rate_ +
                                (log10_upper_rate_ - log10_lower_rate_) *
                                    Math.Abs(slider_position_)) *
               (System.DateTime.Now - last_time_).TotalSeconds;
    }
    last_time_ = System.DateTime.Now;

    UnityEngine.GUILayout.EndHorizontal();

    UnityEngine.GUI.skin = HighLogic.Skin;
    return value_changed;
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
