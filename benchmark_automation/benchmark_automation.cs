using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace principia {
namespace benchmark_automation {

internal class BenchmarkAutomation {
  private const string benchmark_executable = @".\Release\x64\benchmarks.exe";

  private static void Main(string[] args) {
    var benchmark_directory = new DirectoryInfo(args[0]);
    var mathematica_directory = new DirectoryInfo(args[1]);
    var jenkins_directory = new DirectoryInfo(args[2]);
    DateTime date = DateTime.UtcNow;
    string mathematica_date = date.ToString("{yyyy, M, d, H, m, s.fffffff},");
    string mathematica_output_file =
        Path.Combine(mathematica_directory.FullName,
                     "principia_benchmark_results_" +
                         date.ToString("yyyy_M_d__H_m_s_fffffff.wl"));
    Console.WriteLine("DateList[TimeZone->0]: " + mathematica_date);
    Console.WriteLine("Results will be written to " + mathematica_output_file);
    StreamWriter mathematica_stream =
        new StreamWriter(
            File.Open(mathematica_output_file, FileMode.CreateNew));
    mathematica_stream.WriteLine("{");
    mathematica_stream.WriteLine(mathematica_date);
    FileInfo[] files = benchmark_directory.GetFiles("*.cpp");
    Process list_benchmarks_process = new Process {
      StartInfo = new ProcessStartInfo {
        FileName = benchmark_executable,
        Arguments = "--benchmark_list_tests",
        UseShellExecute = false,
        RedirectStandardOutput = true,
        CreateNoWindow = true
      }
    };
    var seen_benchmarks = new Dictionary<string, bool>();
    list_benchmarks_process.Start();
    while (!list_benchmarks_process.StandardOutput.EndOfStream) {
      seen_benchmarks.Add(list_benchmarks_process.StandardOutput.ReadLine().
          Replace(" / ","/").
          Replace(", ", ",").
          Replace(",\n", ",").
          Replace(",\r\n", ","), false);
    }
    string last_benchmark_name = "";
    foreach (FileInfo file in files) {
      StreamReader stream = file.OpenText();
      while (!stream.EndOfStream) {
        string command_line = stream.ReadLine();
        if (command_line.StartsWith("// " + benchmark_executable)) {
          // Get rid of the // NOLINT comments and of the actual command,
          // leaving only the arguments.
          command_line =
              command_line.Split(
                  separator : new []{benchmark_executable, "//"},
                  options   : StringSplitOptions.None)[2].Trim();
          bool has_repetitions =
              command_line.Contains("--benchmark_repetitions");
          Console.WriteLine(
              "Running benchmarks with arguments from " + file.Name);
          Console.WriteLine(command_line);
          Process process = new Process {
            StartInfo = new ProcessStartInfo {
              FileName = benchmark_executable,
              Arguments = command_line,
              UseShellExecute = false,
              RedirectStandardOutput = true,
              CreateNoWindow = true
            }
          };
          process.Start();
          while (!process.StandardOutput.EndOfStream) {
            // A comma followed by a space or a line break occurs for a
            // templated benchmark that has more than one parameter.  Remove the
            // space or line break as it would confuse word splitting.
            string line = process.StandardOutput.ReadLine().
                              Replace(" / ","/").
                              Replace(", ",",").
                              Replace(",\n", ",").
                              Replace(",\r\n", ",");
            string[] words =
                line.Split(separator : new char[]{' '},
                           options   : StringSplitOptions.RemoveEmptyEntries);
            if (seen_benchmarks.Keys.Contains(words[0])) {
              string benchmark_name = words[0];
              string mathematica = "";
              if (last_benchmark_name != benchmark_name) {
                if (seen_benchmarks[benchmark_name]) {
                  throw new ArgumentException("Benchmark " + benchmark_name +
                                              " matches multiple filters");
                }
                seen_benchmarks[benchmark_name] = true;
                if (last_benchmark_name != "") {
                  mathematica += "}}," + mathematica_stream.NewLine;
                }
                last_benchmark_name = benchmark_name;
                mathematica_stream.Flush();
                mathematica += "{";
                mathematica += "\"" + benchmark_name + "\", {";
              } else {
                mathematica += ", ";
              }
              mathematica += (double.Parse(words[1]) *
                              TimeConversionFactor(words[2]))
                      .ToString().Replace("e", "*^");
              mathematica_stream.Write(mathematica);
              Console.WriteLine(
                  line + " " +
                  mathematica.Replace(mathematica_stream.NewLine,
                                      mathematica_stream.NewLine +
                                          new string(' ', line.Length + 1)));
            }
          }
        }
      }
    }
    Console.WriteLine("}}");
    mathematica_stream.WriteLine("}}");
    mathematica_stream.WriteLine("}");
    mathematica_stream.Close();
    var unseen_benchmarks = (from seen in seen_benchmarks
                             where !seen.Value
                             select seen.Key).ToArray();
    if (unseen_benchmarks.Length > 0 ) {
      throw new ArgumentException(
          "The following benchmarks were not run:\n" +
          string.Join("\n", unseen_benchmarks));
    }
  }

  private static double TimeConversionFactor(string unit) {
    if (unit == "ns") {
      return 1e0;
    } else if (unit == "us") {
      return 1e3;
    } else if (unit == "ms") {
      return 1e6;
    } else if (unit == "s") {
      return 1e9;
    } else {
      throw new ArgumentException("Invalid unit: " + unit);
    }
  }
}

}  // namespace benchmark_automation
}  // namespace principia
