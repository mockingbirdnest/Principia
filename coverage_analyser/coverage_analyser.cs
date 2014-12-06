using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

using Microsoft.VisualStudio.Coverage.Analysis;

namespace coverage_analyser {

class coverage_analyser {

  private struct CodeLine {
    public Int64 line_number;
    public string file;
  }

  private static void Main(string[] args) {
    Int64 lines_covered = 0;
    Int64 lines_partially_covered = 0;
    Int64 lines_not_covered = 0;
    string breakdown_percent_not_covered = "";
    string breakdown_lines_not_covered = "";
    string unit_names = "";
    DirectoryInfo directory = new DirectoryInfo(args[0]);
    Console.WriteLine(directory.FullName);
    FileInfo[] files = directory.GetFiles("*.coverage");
    foreach (FileInfo file in files) {
      Console.WriteLine("Analysing " + file.Name);
      Regex file_regex = new Regex(@"^(.+?)(_tests?)+.exe.coverage");
      Match file_match = file_regex.Match(file.Name);
      string tested_unit = file_match.Groups[1].ToString();
      Console.WriteLine("Covering principia::" + tested_unit);
      Regex regex = new Regex("^principia(::|__)" + tested_unit);
      Regex test_file_regex = new Regex("_test.cpp$");
      var covered = new Dictionary<CodeLine, UInt32>();
      using (CoverageInfo info = CoverageInfo.CreateFromFile(file.FullName)) {
        CoverageDS dataset = info.BuildDataSet();
        foreach (CoverageDSPriv.LinesRow lines in dataset.Lines) {
          if (regex.Match(lines.MethodRow.MethodName).Success &&
              !test_file_regex.Match(
                   dataset.GetSourceFileName(lines)).Success) {
            var code_line = new CodeLine {
                file = dataset.GetSourceFileName(lines),
                line_number = lines.LnStart};
            if (lines.LnStart != lines.LnEnd) {
              Console.WriteLine("lines.LnStart != lines.LnEnd");
              return;
            }
            if (covered.ContainsKey(code_line)) {
              covered[code_line] = Math.Min(covered[code_line], lines.Coverage);
            } else {
              covered.Add(code_line, lines.Coverage);
            }
          }
        }
      }
      Int64 subtotal_lines = covered.Count;
      Int64 subtotal_not_covered = 0;
      foreach (var pair in covered) {
        if (pair.Value > 1) {
          ++subtotal_not_covered;
          ++lines_not_covered;
          Console.WriteLine(pair.Key.file + ":" + pair.Key.line_number);
        } else if (pair.Value == 1) {
          ++lines_partially_covered;
        } else {
          ++lines_covered;
        }
      }
      CommaSeparatedAppend(ref breakdown_lines_not_covered,
                           subtotal_not_covered.ToString());
      CommaSeparatedAppend(
          ref breakdown_percent_not_covered,
          (((double)subtotal_not_covered / (double)subtotal_lines) *
               100.0).ToString());
      CommaSeparatedAppend(ref unit_names, tested_unit);
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
        Path.Combine(directory.FullName, "jenkins_percent_coverage.csv"),
        "not covered, partially covered, fully covered\n" +
            ((double)lines_not_covered / (double)total) * 100.0 + ", " +
            ((double)lines_partially_covered / (double)total) * 100.0 + ", " +
            ((double)lines_covered / (double)total) * 100.0);
    File.WriteAllText(
        Path.Combine(directory.FullName, "jenkins_lines_coverage.csv"),
        "not covered, partially covered, fully covered\n" +
            lines_not_covered + ", " +
            lines_partially_covered + ", " +
            lines_covered);
    File.WriteAllText(
        Path.Combine(directory.FullName,
                     "jenkins_lines_coverage_breakdown.csv"),
        unit_names + "\n" + breakdown_lines_not_covered);
    File.WriteAllText(
        Path.Combine(directory.FullName,
                     "jenkins_percent_coverage_breakdown.csv"),
        unit_names + "\n" + breakdown_percent_not_covered);
  }

  private static string ValueAndPercentage(Int64 value, Int64 total) {
    return String.Format("{0,10} ({1,8:P2})",
                         value,
                         (double)value / (double)total);
  }

  private static void CommaSeparatedAppend(ref string csv, string value) {
    if (csv == "") {
      csv = value;
    } else {
      csv = csv + ", " + value;
    }
  }

}

}  // coverage_analyser
