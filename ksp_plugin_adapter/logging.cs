using System;
using System.Runtime.CompilerServices;
using KSP.Localization;

namespace principia {
namespace ksp_plugin_adapter {

internal static class Log {
  internal static string[] severity_names = {
      Localizer.Format("#Principia_Logging_Info"),
      Localizer.Format("#Principia_Logging_Warning"),
      Localizer.Format("#Principia_Logging_Error"),
      Localizer.Format("#Principia_Logging_Fatal")
  };

  internal static void InitGoogleLogging() {
    Interface.InitGoogleLogging();
  }

  internal static void ActivateRecorder(bool activate) {
    Interface.ActivateRecorder(activate);
  }

  internal static void SetBufferedLogging(int max_severity) {
    Interface.SetBufferedLogging(max_severity);
  }

  internal static int GetBufferedLogging() {
    return Interface.GetBufferedLogging();
  }

  internal static void SetBufferDuration(int seconds) {
    Interface.SetBufferDuration(seconds);
  }

  internal static int GetBufferDuration() {
    return Interface.GetBufferDuration();
  }

  internal static void SetSuppressedLogging(int min_severity) {
    Interface.SetSuppressedLogging(min_severity);
  }

  internal static int GetSuppressedLogging() {
    return Interface.GetSuppressedLogging();
  }

  internal static void SetVerboseLogging(int level) {
    Interface.SetVerboseLogging(level);
  }

  internal static int GetVerboseLogging() {
    return Interface.GetVerboseLogging();
  }

  internal static void SetStderrLogging(int min_severity) {
    Interface.SetStderrLogging(min_severity);
  }

  internal static int GetStderrLogging() {
    return Interface.GetStderrLogging();
  }

  internal static void Info(string message,
                            [CallerFilePath] string file = "",
                            [CallerLineNumber] int line = -1) {
    Interface.LogInfo(file, line, message);
  }

  internal static void Warning(string message,
                               [CallerFilePath] string file = "",
                               [CallerLineNumber] int line = -1) {
    Interface.LogWarning(file, line, message);
  }

  internal static void Error(string message,
                             [CallerFilePath] string file = "",
                             [CallerLineNumber] int line = -1) {
    Interface.LogError(file, line, message);
  }

  // Returns an exception so it can be thrown so that the compiler doesn't
  // complain about non-returning code.
  internal static Exception Fatal(string message,
                                  [CallerFilePath] string file = "",
                                  [CallerLineNumber] int line = -1) {
    Interface.LogFatal(file, line, message);
    return new Exception();
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
