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
                            UnityEngine.Color? text_colour = null,
                            int label_width = 3,
                            int field_width = 5) {
    label_ = label;
    label_width_ = label_width;
    field_width_ = field_width;
    unit_ = unit;
    if (formatter == null) {
      // TODO(egg): Is this just "N3"?
      formatter_ = v => v.ToString("#,0.000", Culture.culture);
    } else {
      formatter_ = formatter;
    }
    if (parser == null) {
      // As a special exemption we allow a comma as the decimal separator.
      parser_ = (string s, out double value) => double.TryParse(
          s.Replace(',', '.'),
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

  public double max_value {
    set => max_value_ = value;
  }

  // Set from the UI.
  public double value {
    get => value_ ?? 0.0;
    set {
      value_ = value;
      // Reformat systematically, even if the value has not changed numerically,
      // as it may have been edited all the same (e.g., remove zeroes).
      formatted_value_ = formatter_(value_.Value);
    }
  }

  // Set from programmatic data, e.g., data obtained from C++.  Does nothing if
  // the value has not changed to avoid messing with UI input.
  public double value_if_different {
    set {
      if (!value_.HasValue || value_ != value) {
        value_ = value;
        formatted_value_ = formatter_(value_.Value);
      }
    }
  }

  // TODO(phl): Remove and use value instead.
  public void ResetValue(double new_value) {
    value_ = null;
    value = new_value;
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
                                    options : GUILayoutWidth(label_width_),
                                    style   : style);
      }

      if (enabled) {
        // If the text is not syntactically correct, or it exceeds the upper
        // bound, inform the user by drawing it in the warning style.  Note the
        // fudge factor to account for uncertainty in text/double conversions.
        var style = Style.RightAligned(UnityEngine.GUI.skin.textField);
        if (!parser_(formatted_value_, out double v1) ||
            v1 > max_value_ + 0.1) {
          style = Style.Warning(style);
        }

        var current_event = UnityEngine.Event.current;
        bool text_field_has_focus =
            UnityEngine.GUI.GetNameOfFocusedControl() == text_field_name;

        // Use up the vertical arrow keys before the |TextField| does (it
        // interprets up as home and down as end).
        bool event_was_arrow_key = false;
        if (text_field_has_focus &&
            current_event.type == UnityEngine.EventType.KeyDown &&
            (current_event.keyCode == UnityEngine.KeyCode.UpArrow ||
             current_event.keyCode == UnityEngine.KeyCode.DownArrow)) {
          event_was_arrow_key = true;
          current_event.Use();
        }

        // Draw the text field and give it a name to be able to detect if it has
        // focus.
        UnityEngine.GUI.SetNextControlName(text_field_name);
        formatted_value_ = UnityEngine.GUILayout.TextField(
            text    : formatted_value_,
            style   : style,
            options : GUILayoutWidth(field_width_));
        var text_field = UnityEngine.GUILayoutUtility.GetLastRect();

        DisplayDigitAdjustmentIndicators(text_field, style);

        // Handle digit adjustment using scroll wheel or the arrow keys.
        double increment = GetIncrement(event_was_arrow_key);

        // Check if the user is hovering over a digit.
        if (current_event.type == UnityEngine.EventType.Repaint) {
          scroll_adjustment_ = null;
          if (text_field.Contains(current_event.mousePosition)) {
            int cursor_index = style.GetCursorStringIndex(
                text_field,
                new UnityEngine.GUIContent(formatted_value_),
                current_event.mousePosition);
            if (CanIncrementAt(cursor_index, out double Δ)) {
              scroll_adjustment_ = new DigitAdjustment(cursor_index, Δ);
            }
          }
        }

        // Otherwise, check if the text cursor is before a digit.
        arrows_adjustment_ = null;
        if (scroll_adjustment_ == null && text_field_has_focus) {
          var editor =
              (UnityEngine.TextEditor)UnityEngine.GUIUtility.GetStateObject(
                  typeof(UnityEngine.TextEditor),
                  UnityEngine.GUIUtility.keyboardControl);
          int cursor_index = editor.cursorIndex;
          // If there is a selection, try to edit the last digit of the
          // selection.
          if (cursor_index != editor.selectIndex) {
            cursor_index = Math.Max(cursor_index, editor.selectIndex) - 1;
          }
          if (CanIncrementAt(cursor_index, out double Δ)) {
            arrows_adjustment_ = new DigitAdjustment(cursor_index, Δ);
          }
        }

        // See if the user typed 'Return' in the field, or moved focus
        // elsewhere, or adjusted a digit by scrolling or pressing the arrow
        // keys, in which case we terminate text entry.
        bool terminate_text_entry = false;
        if (current_event.isKey &&
            (current_event.keyCode == UnityEngine.KeyCode.Return ||
             current_event.keyCode == UnityEngine.KeyCode.KeypadEnter) &&
            text_field_has_focus) {
          terminate_text_entry = true;
        } else if (!text_field_has_focus &&
                   formatted_value_ != formatter_(value_.Value)) {
          terminate_text_entry = true;
        } else if (increment != 0) {
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
        if (increment != 0) {
          double incremented_value = value + increment;
          if (incremented_value >= min_value_ &&
              incremented_value <= max_value_) {
            value_changed = true;
            value = incremented_value;
          }
        }
      } else {
        UnityEngine.GUILayout.Label(text    : formatted_value_,
                                    style   : Style.RightAligned(
                                        UnityEngine.GUI.skin.label),
                                    options : GUILayoutWidth(
                                        5 + (unit_ == null ? 2 : 0)));
      }
      UnityEngine.GUILayout.Label(text    : unit_ ?? "",
                                  options : GUILayoutWidth(
                                      unit_ == null ? 0 : 2));

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
          ResetValue(zero_value_);
        }
        if (slider_position_ != 0.0) {
          value_changed = true;
          // Moving the slider doesn't always cause a loss of focus so we
          // terminate input if necessary.
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

  // Adds a marker to hint if a digit can be adjusted by scrolling or pressing
  // the arrow keys.  Note that scrollability takes precedence throughout.
  private void DisplayDigitAdjustmentIndicators(
      UnityEngine.Rect text_field,
      UnityEngine.GUIStyle style) {
    if (scroll_adjustment_.HasValue || arrows_adjustment_.HasValue) {
      UnityEngine.Vector2 indicator_position = style.GetCursorPixelPosition(
          text_field,
          new UnityEngine.GUIContent(formatted_value_),
          (scroll_adjustment_ ?? arrows_adjustment_.Value).index);
      indicator_position.y -= Width(0.08f);
      if (scroll_indicator == null) {
        PrincipiaPluginAdapter.LoadTextureOrDie(
            out scroll_indicator,
            "digit_scroll_indicator.png");
      }
      UnityEngine.GUI.DrawTexture(
          new UnityEngine.Rect(indicator_position,
                                new UnityEngine.Vector2(Width(1.28f),
                                                        Width(1.28f))),
                                scroll_indicator);
    }
  }

  // Returns the increment resulting from any scrolling or arrow keys.
  private double GetIncrement(bool event_was_arrow_key) {
    double increment = 0;
    var current_event = UnityEngine.Event.current;
    if (scroll_adjustment_.HasValue && current_event.isScrollWheel) {
      increment = scroll_adjustment_.Value.increment *
          -UnityEngine.Event.current.delta.normalized.y;
      current_event.Use();
    } else if (arrows_adjustment_.HasValue && event_was_arrow_key) {
      increment = arrows_adjustment_.Value.increment;
      if (current_event.keyCode == UnityEngine.KeyCode.DownArrow) {
        increment = -increment;
      }
    }
    return increment;
  }

  // If |formatted_value_[digit_index]| is a decimal place, returns true and
  // sets |increment| to the value of a unit in that place.  Otherwise, returns
  // false, setting |increment| to 0.
  private bool CanIncrementAt(int digit_index, out double increment) {
    increment = 0;
    // Cannot increment an ill-formed value.
    if (!parser_(formatted_value_, out double base_value)) {
      return false;
    }
    // Can only increment digits.
    if (digit_index < 0 || digit_index >= formatted_value_.Length ||
        !char.IsDigit(formatted_value_[digit_index])) {
      return false;
    }
    // Increment or decrement the digit by 1 (whichever one doesn’t involve
    // carries).  We will take the absolute value of the difference below.
    char[] adjusted_formatted_value = formatted_value_.ToCharArray();
    if (adjusted_formatted_value[digit_index] == '9') {
      --adjusted_formatted_value[digit_index];
    } else {
      ++adjusted_formatted_value[digit_index];
    }
    if (!parser_(new string(adjusted_formatted_value),
                 out double adjusted_value)) {
      return false;
    }
    increment = Math.Abs(adjusted_value - base_value);
    return true;
  }

  private readonly string label_;
  private readonly int label_width_;
  private readonly int field_width_;
  private readonly string unit_;

  private readonly double log10_lower_rate_;
  private readonly double log10_upper_rate_;
  private readonly double zero_value_;
  private readonly double min_value_;

  private readonly ValueFormatter formatter_;
  private readonly ValueParser parser_;
  private readonly UnityEngine.Color? text_colour_;

  private float slider_position_ = 0.0f;
  private DateTime last_time_;
  private double max_value_;
  // It is convenient for the value to be nullable so that we have a way to
  // ensure that an assignment to it will actually have an effect and won't be
  // optimized due to the existing value.  This happens at initialization and
  // during some events handling.
  private double? value_;
  private string formatted_value_;

  // Represents a possible adjustment of the digit at |index| in
  // |formatted_value_|.  The unit in that place is |increment|.
  private struct DigitAdjustment {
    public DigitAdjustment(int index, double increment) {
      this.index = index;
      this.increment = increment;
    }

    public double increment;
    public int index;
  }

  // This field is set if a digit is being hovered over.
  private DigitAdjustment? scroll_adjustment_;
  // This field is set if the text edition cursor is before a digit.
  private DigitAdjustment? arrows_adjustment_;
  private UnityEngine.Texture scroll_indicator;

  private string text_field_name => GetHashCode() + ":text_field";
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
