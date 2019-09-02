using System;
using System.Diagnostics;
using System.IO;

namespace principia {
namespace benchmark_automation {

class BenchmarkAutomation {
  private const string benchmark_executable = @".\Release\x64\benchmarks.exe";

  private static void Main(string[] args) {
    DirectoryInfo benchmark_directory = new DirectoryInfo(args[0]);
    DirectoryInfo mathematica_directory = new DirectoryInfo(args[1]);
    DirectoryInfo jenkins_directory = new DirectoryInfo(args[2]);
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
    string csv_benchmark_names = "";
    string csv_means = "";
    mathematica_stream.WriteLine("{");
    mathematica_stream.WriteLine(mathematica_date);
    FileInfo[] files = benchmark_directory.GetFiles("*.cpp");
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
                              Replace(", ",",").
                              Replace(",\n", ",").
                              Replace(",\r\n", ",");
            string[] words =
                line.Split(separator : new char[]{' '},
                           options   : StringSplitOptions.RemoveEmptyEntries);
            const string mean_suffix = "_mean";
            const string median_suffix = "_median";
            const string stddev_suffix = "_stddev";
            if (words[0].StartsWith("BM_")) {
              if (has_repetitions && words[0].EndsWith(mean_suffix)) {
                string benchmark_name =
                    words[0].Substring(
                        startIndex : 0,
                        length     : words[0].Length - mean_suffix.Length);
                Int64 μ = Int64.Parse(words[1]);
                Console.WriteLine(benchmark_name + ": μ = " + μ + " ns");
                CommaSeparatedAppend(
                    ref csv_benchmark_names,
                    "\"" + benchmark_name.Replace("\"", "\"\"") + "\"");
                CommaSeparatedAppend(ref csv_means, μ.ToString());
              } else if (!has_repetitions) {
                string benchmark_name = words[0];
                Int64 μ = Int64.Parse(words[1]);
                CommaSeparatedAppend(
                    ref csv_benchmark_names,
                    "\"" + benchmark_name.Replace("\"", "\"\"") + "\"");
                CommaSeparatedAppend(ref csv_means, μ.ToString());
                Console.WriteLine(benchmark_name + ": μ = " + μ + " ns");
              }
              if (!words[0].EndsWith(mean_suffix) &&
                  !words[0].EndsWith(median_suffix) &&
                  !words[0].EndsWith(stddev_suffix)) {
                string benchmark_name = words[0];
                if (last_benchmark_name != benchmark_name) {
                  if (last_benchmark_name != "") {
                    mathematica_stream.WriteLine("}},");
                  }
                  last_benchmark_name = benchmark_name;
                  mathematica_stream.Flush();
                  mathematica_stream.Write("{");
                  mathematica_stream.Write("\"" + benchmark_name + "\", {");
                } else {
                  mathematica_stream.Write(", ");
                }
                mathematica_stream.Write(Int64.Parse(words[1]));
              }
            }
          }
        }
      }
    }
    mathematica_stream.WriteLine("}}");
    mathematica_stream.WriteLine("}");
    mathematica_stream.Close();
    File.WriteAllText(
        Path.Combine(jenkins_directory.FullName,
                     "jenkins_benchmark_results.csv"),
        csv_benchmark_names + "\n" + csv_means);
  }

  private static void CommaSeparatedAppend(ref string csv, string value) {
    if (csv == "") {
      csv = value;
    } else {
      csv = csv + ", " + value;
    }
  }
}

}  // namespace benchmark_automation
}  // namespace principia
