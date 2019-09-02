using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.ExceptionServices;
using System.Text.RegularExpressions;
using Microsoft.VisualStudio.Coverage.Analysis;

namespace principia {
namespace tools {

class CoverageAnalyser {
  private struct CodeLine {
    public Int64 line_number;
    public string file;
  }

  [HandleProcessCorruptedStateExceptions]
  private static void Main(string[] args) {
    Int64 lines_covered = 0;
    Int64 lines_partially_covered = 0;
    Int64 lines_not_covered = 0;
    string breakdown_percent_not_covered = "";
    string breakdown_lines_not_covered = "";
    string unit_names = "";
    FileInfo coverage_file = new FileInfo(args[0]);
    using(CoverageInfo info =
              CoverageInfo.CreateFromFile(coverage_file.FullName)) {
        var lines = new List<BlockLineRange>();
        foreach (ICoverageModule module in info.Modules) {
          Console.WriteLine("Analysing " + module.Name);

          Regex module_regex = new Regex(@"^(.+?)(_tests?)+.exe");
          Match module_match = module_regex.Match(module.Name);
          string tested_unit = module_match.Groups[1].ToString();
          Regex regex;
          if (tested_unit == "ksp_plugin") {
            Console.WriteLine("Covering principia::" + tested_unit +
                              " as well as extern \"C\" interface functions" +
                              " (of the form ::principia__Identifier)");
            regex = new Regex("^principia(::" + tested_unit + "|__)");
          } else {
            Console.WriteLine("Covering principia::" + tested_unit);
            regex = new Regex("^principia::" + tested_unit);
          }
          Regex ignored_files_regex = new Regex(
              @"((_test\.cpp|" + @"\.generated\.cc|" + @"\.generated\.h|" +
              @"\.pb\.h|" + @"\.pb\.cc)$|" + @"\\mock_)");
          var covered = new Dictionary<CodeLine, Dictionary<uint, bool>>();

          byte[] coverage_buffer = module.GetCoverageBuffer(null);
          using(ISymbolReader reader = module.Symbols.CreateReader()) {
            for (;;) {
              string method_name;
              try {
                if (!reader.GetNextMethod(out uint method_id,
                                          out method_name,
                                          out string undecorated_method_name,
                                          out string class_name,
                                          out string namespace_name,
                                          lines)) {
                  break;
                }
              } catch (AccessViolationException e) {
                Console.WriteLine(e.ToString());
                continue;
              }
              if (regex.Match(method_name).Success) {
                foreach (var line in lines) {
                  if (!ignored_files_regex.Match(line.SourceFile).Success) {
                    CoverageStatistics stats = CoverageInfo.GetMethodStatistics(
                        coverage_buffer, new List<BlockLineRange>{line});
                    if (line.StartLine != line.EndLine ||
                        stats.LinesCovered + stats.LinesNotCovered != 1) {
                      Console.WriteLine("in " + method_name);
                      Console.WriteLine(line.SourceFile + ":" + line.StartLine +
                                        "-" + line.EndLine);
                      Console.WriteLine(stats.LinesCovered + "," +
                                        stats.LinesNotCovered);
                      Environment.Exit(1);
                    }
                    bool block_is_covered = stats.LinesCovered == 1;
                    var code_line = new CodeLine{file = line.SourceFile,
                                                 line_number = line.StartLine};
                    if (!covered.ContainsKey(code_line)) {
                      covered.Add(code_line, new Dictionary<uint, bool>());
                    }
                    if (covered[code_line].ContainsKey(line.BlockIndex)) {
                      covered[code_line][line.BlockIndex] |= block_is_covered;
                    } else {
                      covered[code_line].Add(line.BlockIndex, block_is_covered);
                    }
                  }
                }
              }
              lines.Clear();
            }
          }
          Int64 subtotal_lines = covered.Count;
          Int64 subtotal_not_covered = 0;
          foreach (var pair in covered) {
            bool line_is_partially_covered =
                (from block_coverage in pair.Value select block_coverage.Value)
                    .Any(_ => _);
            bool line_is_fully_covered =
                (from block_coverage in pair.Value select block_coverage.Value)
                    .All(_ => _);
            if (line_is_fully_covered) {
              ++lines_covered;
            } else if (line_is_partially_covered) {
              ++lines_partially_covered;
            } else {
              ++subtotal_not_covered;
              ++lines_not_covered;
              Console.WriteLine(pair.Key.file + ":" + pair.Key.line_number);
            }
          }
          CommaSeparatedAppend(ref breakdown_lines_not_covered,
                               subtotal_not_covered.ToString());
          CommaSeparatedAppend(
              ref breakdown_percent_not_covered,
              (((double)subtotal_not_covered / (double)subtotal_lines) * 100.0)
                  .ToString());
          CommaSeparatedAppend(ref unit_names, tested_unit);
      }
    }
    Int64 total = lines_partially_covered + lines_covered + lines_not_covered;
    Console.WriteLine("Covered      : " +
                          ValueAndPercentage(lines_partially_covered +
                                                 lines_covered, total));
    Console.WriteLine("  Partially  : " +
                          ValueAndPercentage(lines_partially_covered, total));
    Console.WriteLine("  Completely : " +
                          ValueAndPercentage(lines_covered, total));
    Console.WriteLine("Not Covered  : " +
                          ValueAndPercentage(lines_not_covered, total));
    File.WriteAllText(
        Path.Combine(coverage_file.Directory.FullName,
                     "jenkins_percent_coverage.csv"),
        "not covered,partially covered,fully covered\n" +
            ((double)lines_not_covered / (double)total) * 100.0 + "," +
            ((double)lines_partially_covered / (double)total) * 100.0 + "," +
            ((double)lines_covered / (double)total) * 100.0);
    File.WriteAllText(Path.Combine(coverage_file.Directory.FullName,
                                   "jenkins_lines_coverage.csv"),
                      "not covered,partially covered,fully covered\n" +
                          lines_not_covered + "," + lines_partially_covered +
                          "," + lines_covered);
    File.WriteAllText(Path.Combine(coverage_file.Directory.FullName,
                                   "jenkins_lines_coverage_breakdown.csv"),
                      unit_names + "\n" + breakdown_lines_not_covered);
    File.WriteAllText(Path.Combine(coverage_file.Directory.FullName,
                                   "jenkins_percent_coverage_breakdown.csv"),
                      unit_names + "\n" + breakdown_percent_not_covered);
  }

  private static string ValueAndPercentage(Int64 value, Int64 total) {
    return $"{value,10} ({(double)value / (double)total,8:P2})";
  }

  private static void CommaSeparatedAppend(ref string csv, string value) {
    if (csv == "") {
      csv = value;
    } else {
      csv = csv + "," + value;
    }
  }
}

}  // namespace tools
}  // namespace principia
