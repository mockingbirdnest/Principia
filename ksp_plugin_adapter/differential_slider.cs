using System;
using System.Globalization;

namespace principia {
namespace ksp_plugin_adapter {

internal class DifferentialSlider : ScalingRenderer {
  public delegate string ValueFormatter(double value);
  public delegate bool ValueParser(string s, out double value);

  // Rates are in units of |value| per real-time second.
  public DifferentialSlider(string label,
                            string unit,
                            double log10_lower_rate,
                            double log10_upper_rate,
                            double zero_value = 0,
                            double min_value = double.NegativeInfinity,
                            double max_value = double.PositiveInfinity,
                            ValueFormatter formatter = null,
                            ValueParser parser = null,
                            UnityEngine.Color? text_colour = null) {
    label_ = label;
    unit_ = unit;
    if (formatter == null) {
      formatter_ = v => v.ToString("#,0.000", Culture.culture);
    } else {
      formatter_ = formatter;
    }
    if (parser == null) {
      // As a special exemption we allow a comma as the decimal separator.
      parser_ = (string s, out double value) =>
                    Double.TryParse(s.Replace(',', '.'),
                                    NumberStyles.AllowDecimalPoint |
                                    NumberStyles.AllowLeadingSign |
                                    NumberStyles.AllowLeadingWhite |
                                    NumberStyles.AllowThousands |
                                    NumberStyles.AllowTrailingWhite,
                                    Culture.culture.NumberFormat,
                                    out value);
    } else {
      parser_ = parser;
    }
    log10_lower_rate_ = log10_lower_rate;
    log10_upper_rate_ = log10_upper_rate;
    zero_value_ = zero_value;
    min_value_ = min_value;
    max_value_ = max_value;
    text_colour_ = text_colour;
  }

  public double value {
    get {
      return value_ ?? 0.0;
    }
    set {
      if (!value_.HasValue || value_ != value) {
        value_ = value;
        formatted_value_ = formatter_(value_.Value);
      }
    }
  }

  // Renders the |DifferentialSlider|.  Returns true if and only if |value|
  // changed.
  public bool Render(bool enabled) {
    bool value_changed = false;

    using (new UnityEngine.GUILayout.HorizontalScope()) {
      {
        var style = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.label);
        if (text_colour_.HasValue) {
          style.normal.textColor = text_colour_.Value;
        }
        UnityEngine.GUILayout.Label(text    : label_,
                                    options : GUILayoutWidth(3),
                                    style   : style);
      }

      if (enabled) {
        var style = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.textField);
        style.alignment = UnityEngine.TextAnchor.MiddleRight;

        // If the text is not syntactically correct, inform the user by drawing
        // it in colour.  We don't expect to see the "red" case as we should
        // revert to a parseable value on exit.
        if (!parser_(formatted_value_, out double v1)) {
          style.focused.textColor = XKCDColors.Orange;
          style.normal.textColor = XKCDColors.Red;
        }

        // Draw the text field and give it a name to be able to detect if it has
        // focus.
        String text_field_name = GetHashCode() + ":text_field";
        UnityEngine.GUI.SetNextControlName(text_field_name);
        formatted_value_ = UnityEngine.GUILayout.TextField(
            text    : formatted_value_,
            style   : style,
            options : GUILayoutWidth(5 + (unit_ == null ? 2 : 0)));

        // See if the user typed 'Return' in the field, or moved focus
        // elsewhere, in which case we terminate text entry.
        bool terminate_text_entry = false;
        var current_event = UnityEngine.Event.current;
        if (UnityEngine.Event.current.isKey &&
            UnityEngine.Event.current.keyCode == UnityEngine.KeyCode.Return &&
            UnityEngine.GUI.GetNameOfFocusedControl() == text_field_name) {
          terminate_text_entry = true;
        } else if (UnityEngine.GUI.GetNameOfFocusedControl() !=
                       text_field_name &&
                   formatted_value_ != formatter_(value_.Value)) {
          terminate_text_entry = true;
        }
        if (terminate_text_entry) {
          // Try to parse the input.  If that fails, go back to the previous
          // legal value.
          if (parser_(formatted_value_, out double v2)) {
            value_changed = true;
            value = v2;
            slider_position_ = 0;
            // TODO(phl): If the value computed here is rejected by the C++, we
            // revert to the previous value just as if there was a parsing error
            // and this is not nice.
          } else {
            // Go back to the previous legal value.
            formatted_value_ = formatter_(value_.Value);
          }
        }
      } else {
        var style = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.label);
        style.alignment = UnityEngine.TextAnchor.MiddleRight;
        UnityEngine.GUILayout.Label(
            text    : formatted_value_,
            style   : style,
            options : GUILayoutWidth(5 + (unit_ == null ? 2 : 0)));
      }
      UnityEngine.GUILayout.Label(
          text    : unit_ ?? "",
          options : GUILayoutWidth(unit_ == null ? 0 : 2));

      if (enabled) {
        if (!UnityEngine.Input.GetMouseButton(0)) {
          slider_position_ = 0;
        }
        slider_position_ = UnityEngine.GUILayout.HorizontalSlider(
            value      : slider_position_,
            leftValue  : -1,
            rightValue : 1,
            options    : UnityEngine.GUILayout.ExpandWidth(true));

        if (UnityEngine.GUILayout.Button("0", GUILayoutWidth(1))) {
          value_changed = true;
          // Force a change of value so that any input is discarded.
          value = zero_value_ + 1;
          value = zero_value_;
        }
        if (slider_position_ != 0.0) {
          value_changed = true;
          // Moving the slider doesn't cause a loss of focus so we terminate
          // input if necessary.
          if (parser_(formatted_value_, out double v)) {
            value = v;
          }
          value += Math.Sign(slider_position_) *
                Math.Pow(10, log10_lower_rate_ +
                                (log10_upper_rate_ - log10_lower_rate_) *
                                    Math.Abs(slider_position_)) *
                (DateTime.Now - last_time_).TotalSeconds;
          value = Math.Min(Math.Max(min_value_, value), max_value_);
        }
      } else {
        slider_position_ = 0;
      }
      last_time_ = DateTime.Now;
    }
    return value_changed;
  }

  private readonly string label_;
  private readonly string unit_;

  private readonly double log10_lower_rate_ = -3;
  private readonly double log10_upper_rate_ = 3.5;
  private readonly double zero_value_;
  private readonly double min_value_;
  private readonly double max_value_;

  private readonly ValueFormatter formatter_;
  private readonly ValueParser parser_;
  private readonly UnityEngine.Color? text_colour_;

  private float slider_position_ = 0.0f;
  private DateTime last_time_;
  private double? value_;
  private String formatted_value_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
