using System;
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
           date_time_formatter.Day) / date_time_formatter.Hour;
      days = (int)(seconds_ - seconds - minutes * date_time_formatter.Minute -
                   hours * date_time_formatter.Hour) / date_time_formatter.Day;
      return true;
    } catch (OverflowException) {
      return false;
    }
  }

  // Formats a duration, optionally omitting leading components if they are 0,
  // and leading 0s on the days; optionally exclude seconds.
  public string Format(bool with_leading_zeroes, bool with_seconds) {
    return seconds_.ToString("+;-") + FormatPositive(with_leading_zeroes,
                                                     with_seconds);
  }

  public string FormatPositive(bool with_leading_zeroes, bool with_seconds) {
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
    string formatted_days = with_leading_zeroes
                                ? (is_stock_day
                                       ? days.ToString("0000;0000")
                                       : days.ToString("000;000"))
                                : (days == 0 ? "" : days.ToString("0;0"));
    string formatted_hours = with_leading_zeroes || hours != 0
                                 ? (date_time_formatter.Day /
                                    date_time_formatter.Hour < 10
                                        ? hours.ToString("0;0")
                                        : hours.ToString("00;00"))
                                 : "";
    string formatted_minutes = with_leading_zeroes || hours != 0
                                   ? minutes.ToString("00;00")
                                   : "";
    string formatted = formatted_days + $"{nbsp}{day_symbol}{nbsp}" +
                       formatted_hours + $"{nbsp}h{nbsp}" + formatted_minutes +
                       $"{nbsp}min";
    if (with_seconds) {
      formatted += $"{nbsp}" + seconds.ToString("00.0;00.0") + $"{nbsp}s";
    }
    return formatted;
  }

  public double total_seconds => seconds_;

  public static bool TryParse(string text,
                              bool with_seconds,
                              out PrincipiaTimeSpan time_span) {
    time_span = new PrincipiaTimeSpan(double.NaN);
    // Using a technology that is customarily used to parse HTML.
    string pattern = @"^[+]?\s*(\d+)\s*" + day_symbol +
                     @"\s*(\d+)\s*h\s*(\d+)\s*min";
    if (with_seconds) {
      pattern += @"\s*([0-9.,']+)\s*s$";
    } else {
      pattern += @"$";
    }
    var regex = new Regex(pattern);
    var match = regex.Match(text);
    if (!match.Success) {
      return false;
    }
    string days = match.Groups[1].Value;
    string hours = match.Groups[2].Value;
    string minutes = match.Groups[3].Value;
    string seconds = "0";
    if (with_seconds) {
      seconds = match.Groups[4].Value;
    }
    if (!int.TryParse(days, out int d) || !int.TryParse(hours, out int h) ||
        !int.TryParse(minutes, out int min) || !double.TryParse(
            seconds.Replace(',', '.'),
            NumberStyles.AllowDecimalPoint | NumberStyles.AllowThousands,
            Culture.culture.NumberFormat,
            out double s)) {
      return false;
    }
    time_span = new PrincipiaTimeSpan(d * date_time_formatter.Day +
                                      h * date_time_formatter.Hour +
                                      min * date_time_formatter.Minute + s);
    return true;
  }

  public static int day_duration => date_time_formatter.Day;
  public static string day_symbol => is_stock_day ? "d6" : "d";

  private static IDateTimeFormatter date_time_formatter =>
      KSPUtil.dateTimeFormatter;

  private static bool is_stock_day => GameSettings.KERBIN_TIME &&
                                      date_time_formatter.Day == 6 * 60 * 60;

  private readonly double seconds_;
  private const string nbsp = "\xA0";
}

}  // namespace ksp_plugin_adapter
}  // namespace principia