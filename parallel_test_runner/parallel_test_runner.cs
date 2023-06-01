using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;

namespace principia {
namespace parallel_test_runner {

class ParallelTestRunner {
  enum Granularity {
    Package,
    TestCase,
    Test,
  }

  private static T ParseEnum<T>(string value) {
    return (T)Enum.Parse(typeof(T), value, true);
  }

  private const string vsinstr =
      @"C:\Program Files\Microsoft Visual Studio\2022\Enterprise\Team Tools" +
      @"\Performance Tools\x64\vsinstr.exe";

  static void StubbornStart(Process process) {
    for (int remaining = 10; remaining > 0; --remaining) {
      try {
        process.Start();
        return;
      } catch (Exception e) {
        Console.WriteLine("Exception " +
                          e +
                          process.StartInfo.FileName +
                          " " +
                          process.StartInfo.Arguments);
        if (remaining == 1) {
          throw;
        }
        Thread.Sleep(TimeSpan.FromSeconds(1));
      }
    }
  }

  static Task RunProcessAsync(string file_name, string args) {
    var process = new Process{StartInfo = {FileName = file_name,
                                           Arguments = args,
                                           UseShellExecute = false,
                                           RedirectStandardError = false,
                                           RedirectStandardOutput = true},
                              EnableRaisingEvents = true};
    return new Task(async () => {
      StubbornStart(process);
      while (!process.StandardOutput.EndOfStream) {
        Console.WriteLine(await process.StandardOutput.ReadLineAsync());
      }
      process.WaitForExit();
      process.Close();
    });
  }

  static void Main(string[] args) {
    Granularity? granularity_option = null;
    bool? instrument_option = null;
    bool? also_run_disabled_tests_option = null;
    string filter_option = null;

    var death_test_processes = new List<Process>();
    var processes = new List<Process>();
    int test_process_counter = 0;

    var instrument_tests = new List<Task>();

    foreach (string arg in args) {
      if (arg.StartsWith("--") && arg.Contains(":")) {
        string[] split = arg.Split(new []{"--", ":"}, StringSplitOptions.None);
        string option = split[1];
        string value = split[2];
        if (option == "granularity") {
          granularity_option = ParseEnum<Granularity>(value);
        } else if (option == "instrument") {
          instrument_option = bool.Parse(value);
        } else if (option == "also_run_disabled_tests") {
          also_run_disabled_tests_option = bool.Parse(value);
        } else if (option == "filter") {
          filter_option = value;
        } else {
          Console.WriteLine("Unknown option " + option);
          Environment.Exit(1);
        }
        continue;
      }
      Granularity granularity = granularity_option ?? Granularity.Test;
      bool instrument = instrument_option ?? false;
      bool also_run_disabled_tests = also_run_disabled_tests_option ?? false;
      string filter = filter_option;
      granularity_option = null;
      instrument_option = null;
      filter_option = null;

      if (filter != null && granularity == Granularity.Package) {
        Console.WriteLine(
            "--filter is not supported with --granularity:Package");
        Environment.Exit(1);
      }

      string[] test_binaries = Directory.GetFiles(arg, "*_tests.exe");
      foreach (string test_binary in test_binaries) {
        if (instrument) {
          instrument_tests.Add(
              RunProcessAsync(vsinstr, "/coverage \"" + test_binary + "\""));
        }
        if (granularity == Granularity.Package) {
          var process = new Process{
              StartInfo = {
                  UseShellExecute = false,
                  RedirectStandardOutput = true,
                  RedirectStandardError = true,
                  FileName = test_binary,
                  Arguments = "--gtest_filter=-*DeathTest.*"
              }
          };
          process.StartInfo.Arguments +=
              " --gtest_output=xml:TestResults\\gtest_results_" +
              test_process_counter++ + ".xml";
          processes.Add(process);
          process = new Process{
              StartInfo = {
                  UseShellExecute = false,
                  RedirectStandardOutput = true,
                  RedirectStandardError = true,
                  FileName = test_binary,
                  Arguments = "--gtest_filter=*DeathTest.*"
              }
          };
          process.StartInfo.Arguments +=
              " --gtest_output=xml:TestResults\\gtest_results_" +
              test_process_counter++ + ".xml";
          if (also_run_disabled_tests) {
            process.StartInfo.Arguments += " --gtest_also_run_disabled_tests";
          }
          death_test_processes.Add(process);
          continue;
        }
        var info = new ProcessStartInfo(
            test_binary,
            filter == null
                ? "--gtest_list_tests"
                : $"--gtest_list_tests --gtest_filter={filter}"){
            UseShellExecute = false,
            RedirectStandardOutput = true
        };
        var list_tests = Process.Start(info);
        var output = list_tests.StandardOutput;
        string test_case = null;
        bool? is_death_test = null;
        while(!output.EndOfStream) {
          string line = output.ReadLine();
          if (line[0] != ' ') {
            test_case = line;
            is_death_test = Regex.Match(line, ".*DeathTest").Success;
            if (granularity == Granularity.TestCase) {
              var process = new Process{
                  StartInfo = {
                      UseShellExecute = false,
                      RedirectStandardOutput = true,
                      RedirectStandardError = true,
                      FileName = test_binary,
                      Arguments = "--gtest_filter=" + test_case + "*"
                  }
              };
              process.StartInfo.Arguments +=
                  " --gtest_output=xml:TestResults\\gtest_results_" +
                  test_process_counter++ + ".xml";
              if (also_run_disabled_tests) {
                process.StartInfo.Arguments +=
                    " --gtest_also_run_disabled_tests";
              }
              if (is_death_test.Value) {
                death_test_processes.Add(process);
              } else {
                processes.Add(process);
              }
            }
          } else if (granularity == Granularity.Test) {
            var process = new Process{
                StartInfo = {
                    UseShellExecute = false,
                    RedirectStandardOutput = true,
                    RedirectStandardError = true,
                    FileName = test_binary,
                    Arguments = Encoding.Default.GetString(
                        Encoding.UTF8.GetBytes(
                            "--gtest_filter=" + test_case + line.Split(' ')[2]))
                }
            };
            process.StartInfo.Arguments +=
                " --gtest_output=xml:TestResults\\gtest_results_" +
                test_process_counter++ + ".xml";
            if (also_run_disabled_tests) {
              process.StartInfo.Arguments += " --gtest_also_run_disabled_tests";
            }
            if (is_death_test.Value) {
              death_test_processes.Add(process);
            } else {
              processes.Add(process);
            }
          }
        }
      }
    }
    var stopwatch = new Stopwatch();

    Console.WriteLine("Instrumenting " + instrument_tests.Count +
                      " processes...");
    stopwatch.Start();
    if (instrument_tests.Count > 0) {
      instrument_tests.ForEach(task => task.Start());
      Task.WaitAll(instrument_tests.ToArray());
      Console.WriteLine("Done (" + stopwatch.ElapsedMilliseconds + " ms)");
    }

    Console.WriteLine("Running " + death_test_processes.Count +
                      " death test processes...");
    stopwatch.Restart();
    var errors = new ConcurrentBag<string>();
    foreach (Process process in death_test_processes) {
      process.StartInfo.RedirectStandardOutput = true;
      process.StartInfo.RedirectStandardError = true;
      StubbornStart(process);
      Task standard_output_writer = Task.Run(async () => {
        await HandleOutput("O", process.StandardOutput);
      });
      Task standard_error_writer = Task.Run(async () => {
        await HandleOutput("E", process.StandardError);
      });
      process.WaitForExit();
      Task.WaitAll(new Task[]
                        { standard_output_writer, standard_error_writer });
      if (process.ExitCode != 0) {
        errors.Add("Exit code " + process.ExitCode + " from a death test: " +
                   process.StartInfo.FileName +
                   " " +
                   process.StartInfo.Arguments);
      }
      process.Close();
    }

    Console.WriteLine("Running " + processes.Count + " processes...");
    var process_semaphore = new Semaphore(initialCount: max_parallelism,
                                          maximumCount: max_parallelism);

    Task[] tasks = new Task[processes.Count];
    for (int i = 0; i < processes.Count; ++i) {
      var process = processes[i];
      // We cannot use i in the lambdas, it would be captured by reference.
      int index = i;
      process_semaphore.WaitOne();
      StubbornStart(process);
      tasks[i] = Task.Run(() => {
        Task standard_output_writer = Task.Run(async () => {
          await HandleOutput("O" + index.ToString().PadLeft(4),
                             process.StandardOutput);
        });
        Task standard_error_writer = Task.Run(async () => {
          await HandleOutput("E" + index.ToString().PadLeft(4),
                             process.StandardError);
        });
        process.WaitForExit();
        if (process.ExitCode != 0) {
          errors.Add("Exit code " +
                     process.ExitCode +
                     " from (" +
                     index.ToString() +
                     ") " +
                     process.StartInfo.FileName +
                     " " +
                     process.StartInfo.Arguments);
        }
        Task.WaitAll(new Task[]
                         { standard_output_writer, standard_error_writer });
        process.Close();
        process_semaphore.Release();
      });
    }

    Task.WaitAll(tasks);
    foreach (string error in errors) {
      Console.WriteLine(error);
    }
    Console.WriteLine("Done (" + stopwatch.ElapsedMilliseconds + " ms)");
    Environment.Exit(errors.Count);
  }

  struct Location {
    public string file;
    public string line;
  }

  class Error {
    public List<Location> locations = new List<Location>();
    public string title;
    public string message;

    public void WriteToGitHub() {
      foreach (var location in locations) {
        Console.WriteLine(
            $@"::error file={location.file},line={location.line},title={
                title}::{message.Replace("\n", "%0A")}");
      }
    }
  };

  static async Task HandleOutput(string prefix, StreamReader output) {
    var error_line = new Regex(
        @"^(?:(?<file>{DIR}[^:]*)\((?<line>\d+)\)|unknown file): error: (?<message>.*)");
    error_line = new Regex(error_line.ToString().Replace(
        "{DIR}",
        Regex.Escape(Directory.GetCurrentDirectory())));
    var fatal_line = new Regex(
        @"^F\d{4} \d\d:\d\d:\d\d\.\d{6}\s+\d+ (?<file>[^:]*):(?<line>\d+)\] (?<message>.*)");
    var other_line = new Regex(@"
        ^ (
            (?!(?-x)\[  DEATH   \])  # Miscellaneous gtest output, excluding the
                    \[..........\]   # “actual” side of death tests.
          | (?:(?-x)[IWE]\d{4} \d\d:\d\d:\d\d\.\d{6})  # Non-fatal glog output.
          )",
        RegexOptions.IgnorePatternWhitespace);
    Error error = null;
    while (!output.EndOfStream) {
      string line = await output.ReadLineAsync();
      var error_match = error_line.Match(line);
      var fatal_match = fatal_line.Match(line);
      if (other_line.IsMatch(line) ||
          error_match.Success ||
          fatal_match.Success) {
        error?.WriteToGitHub();
        error = null;
        if (error_match.Success) {
          error = new Error{
              title = "Test failure",
              message = error_match.Groups["message"].Value};
          var location = new Location{
              file = error_match.Groups["file"].Value,
              line = error_match.Groups["line"].Value};
          if (location.file != null && location.line != null) {
            error.locations.Add(location);
          }
        } else if (fatal_match.Success) {
          error = new Error{
              title = "Check failure",
              message = fatal_match.Groups["message"].Value};
        }
      } else if (error is Error e) {
        var stack_path = new Regex(
            @"^(?!\[  DEATH   \]).*\((?<file>{DIR}[^:]*):(?<line>\d+)\)$");
        stack_path = new Regex(stack_path.ToString()
            .Replace("{DIR}", Regex.Escape(Directory.GetCurrentDirectory())));
        var stack_match = stack_path.Match(line);
        if (stack_match.Success) {
          e.locations.Add(new Location{
              file = stack_match.Groups["file"].Value,
              line = stack_match.Groups["line"].Value});
        }
        e.message += "\n";
        e.message += line;
      }
      Console.WriteLine($"{prefix} {line}");
    }
    error?.WriteToGitHub();
  }

  // Maximum number of processes to execute in parallel.
  private static readonly int max_parallelism = 100;
}

}  // namespace parallel_test_runner
}  // namespace principia
