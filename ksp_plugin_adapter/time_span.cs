using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text.RegularExpressions;

namespace principia {
namespace ksp_plugin_adapter {

// This class starts with "Principia" to avoid confusion with the .Net TimeSpan.
class PrincipiaTimeSpan {
  public PrincipiaTimeSpan(double seconds) {
    seconds_ = seconds;
  }

  public bool Split(out int days,
                    out int hours,
                    out int minutes,
                    out double seconds) {
    days = 0;
    hours = 0;
    minutes = 0;
    seconds = seconds_;
    try {
      seconds = seconds_ % date_time_formatter.Minute;
      minutes = ((int)(seconds_ - seconds) % date_time_formatter.Hour) /
                date_time_formatter.Minute;
      hours =
          ((int)(seconds_ - seconds - minutes * date_time_formatter.Minute) %
           date_time_formatter.Day) /
          date_time_formatter.Hour;
      days = (int)(seconds_ -
                   seconds -
                   minutes * date_time_formatter.Minute -
                   hours * date_time_formatter.Hour) /
             date_time_formatter.Day;
      return true;
    } catch (OverflowException) {
      return false;
    }
  }

  // Formats a duration, optionally omitting leading components if they are 0,
  // and leading 0s on the days; optionally exclude seconds.
  public string Format(bool with_leading_zeroes,
                       bool with_seconds,
                       bool iau_style = false,
                       int fractional_seconds_digits = 1) {
    return seconds_.ToString("+;−") +
           FormatPositive(with_leading_zeroes,
                          with_seconds,
                          iau_style,
                          fractional_seconds_digits);
  }

  public string FormatPositive(
      bool with_leading_zeroes, 
      bool with_seconds,
      bool iau_style = false,
      int fractional_second_digits = 1) {
    if (!Split(out int days,
               out int hours,
               out int minutes,
               out double seconds)) {
      // In case of error, saturate to the largest representable value.
      days = int.MaxValue;
      hours = 0;
      minutes = 0;
      seconds = 0;
    }
    var components = new List<string>();
    if (with_leading_zeroes) {
      components.Add(day_is_short
                         ? days.ToString("0000;0000")
                         : days.ToString("000;000"));
    } else if (days != 0) {
      components.Add(days.ToString("0;0"));
    }
    if (components.Count > 0) {
      components.Add(iau_style ?$"{short_day_symbol}{nbsp}"
                               :$"{nbsp}{day_symbol}{nbsp}");
    }
    if (components.Count > 0 || with_leading_zeroes || hours != 0) {
      components.Add(day_is_short
                         ? hours.ToString("0;0")
                         : hours.ToString("00;00"));
      components.Add(iau_style ? $"ʰ{nbsp}"
                               : $"{nbsp}h{nbsp}");
    }
    if (components.Count > 0 ||
        with_leading_zeroes ||
        minutes != 0 ||
        !with_seconds) {
      components.Add(minutes.ToString("00;00"));
      components.Add(iau_style ? "ᵐ"
                               : $"{nbsp}min");
    }
    if (with_seconds) {
      if (fractional_second_digits > 0) {
        string fractional_format = new string('0', fractional_second_digits);
        string seconds_field =
            Regex.Replace(
                seconds.ToString(
                    $"00.{fractional_format};00.{fractional_format}"),
                @"\d{3}(?=\d)",
                match => match.Value + "'");;
        if (iau_style) {
          components.Add(nbsp + seconds_field.Replace(".", "ˢ."));
        } else {
          components.Add($"{nbsp}{seconds_field}{nbsp}s");
        }
      } else {
        components.Add(iau_style ? $"{nbsp}{seconds:00}ˢ"
                                 : $"{nbsp}{seconds:00}{nbsp}s");
      }
    }
    return string.Join("", components.ToArray());
  }

  public double total_seconds => seconds_;

  public static bool TryParse(string text,
                              out PrincipiaTimeSpan time_span) {
    time_span = new PrincipiaTimeSpan(double.NaN);
    // Using a technology that is customarily used to parse HTML.
    // Wrapping the literal in a Regex constructor and then substituting the day
    // symbols in order to get VS to syntax highlight the regex.
    var regex = new Regex(@"
        ^[+]?\s*
        (?:(?<days>\d+)\s*(?:{day_symbol}|{short_day_symbol})\s*)?
        (?:(?<hours>\d+)\s*[hʰ]\s*)?
        (?:(?<minutes>\d+)\s*(?:[mᵐ]|min)\s*)?
        (?:(?<seconds>[0-9.,']+)\s*[sˢ]\s*|
           (?<integer_seconds>[0-9']+)[sˢ][.,]
                (?<fractional_seconds>[0-9']+))?$",
    RegexOptions.IgnorePatternWhitespace);
    regex = new Regex(
        regex.ToString().Replace("{day_symbol}", day_symbol)
                        .Replace("{short_day_symbol}", short_day_symbol),
        RegexOptions.IgnorePatternWhitespace);
    var match = regex.Match(text);
    if (!match.Success) {
      return false;
    }

    var days_group = match.Groups["days"];
    var hours_group = match.Groups["hours"];
    var minutes_group = match.Groups["minutes"];
    var seconds_group = match.Groups["seconds"];
    var integer_seconds_group = match.Groups["integer_seconds"];
    var fractional_seconds_group = match.Groups["fractional_seconds"];
    string days = days_group.Success ? days_group.Value : "0";
    string hours = hours_group.Success ? hours_group.Value : "0";
    string minutes = minutes_group.Success ? minutes_group.Value : "0";
    string seconds = seconds_group.Success
        ? seconds_group.Value
        : integer_seconds_group.Success
            ? string.Join(".",
                          integer_seconds_group.Value,
                          fractional_seconds_group.Value)
            : "0";
    if (!int.TryParse(days, out int d) ||
        !int.TryParse(hours, out int h) ||
        !int.TryParse(minutes, out int min) ||
        !double.TryParse(seconds.Replace(',', '.').Replace("'", ""),
                         NumberStyles.AllowDecimalPoint |
                         NumberStyles.AllowThousands,
                         Culture.culture.NumberFormat,
                         out double s)) {
      return false;
    }
    time_span = new PrincipiaTimeSpan(
        Convert.ToDouble(d) * date_time_formatter.Day +
        Convert.ToDouble(h) * date_time_formatter.Hour +
        Convert.ToDouble(min) * date_time_formatter.Minute +
        s);
    return true;
  }

  public static int day_duration => date_time_formatter.Day;
  
  public static string day_symbol =>
      hour_divides_day && day_is_short
          ? "d" + (int)(date_time_formatter.Day / date_time_formatter.Hour)
          : "d";

  public static string short_day_symbol =>
      hour_divides_day && day_is_short
          ? "ᵈ" + "⁰¹²³⁴⁵⁶⁷⁸⁹"[date_time_formatter.Day /
                               date_time_formatter.Hour]
          : "ᵈ";

  private static bool day_is_short =>
      date_time_formatter.Day / date_time_formatter.Hour < 10;

  private static bool hour_divides_day =>
      (int)(date_time_formatter.Day / date_time_formatter.Hour) *
      date_time_formatter.Hour ==
      date_time_formatter.Day;

  private static IDateTimeFormatter date_time_formatter =>
      KSPUtil.dateTimeFormatter;

  private readonly double seconds_;
  private const string nbsp = "\xA0";
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
