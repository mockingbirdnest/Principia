using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace principia {
namespace sourcerer {

public class Filenames {
  public static string CorrespondingHeader(FileInfo file_info) {
    return Regex.Replace(
        Path.GetFullPath(file_info.Name, file_info.DirectoryName!),
        @"_body\.hpp$|\.cpp$",
        @".hpp");
  }

  public static string SolutionRelativePath(string solution_full_path,
                                            FileInfo file) {
    string file_full_path = Path.GetFullPath(file.Name, file.DirectoryName!);
    if (file_full_path.StartsWith(solution_full_path)) {
      return file_full_path.Substring(solution_full_path.Length + 1);
    } else {
      return file_full_path;
    }
  }

  public static string ToSlash(string with_backslash) {
    return with_backslash.Replace('\\', '/');
  }

  public static string GetExtension(FileInfo file_info) {
    // If the file is foo.mathematica.h, this returns .mathematica.h.
    return Regex.Replace(file_info.Name, @"^[^.]*\.", ".");
  }

  public static string RemoveExtension(FileInfo file_info) {
    // If the file is foo.mathematica.h, this returns foo.
    return Regex.Replace(file_info.Name, @"\..+$", "");
  }

  public static bool IsBenchmark(FileInfo file_info) {
    return Regex.IsMatch(file_info.Name, @"_benchmark\.cpp$");
  }

  public static bool IsBody(FileInfo file_info, HashSet<string> extra_headers) {
    if (extra_headers.Contains(file_info.Name)) {
      return false;
    } else {
      return Regex.IsMatch(file_info.Name, @"_body\.hpp$|\.cpp$");
    }
  }

  public static bool IsBodyHpp(FileInfo file_info) {
    return Regex.IsMatch(file_info.Name, @"^.*_body\.hpp$");
  }

  public static bool IsGeneratedH(FileInfo file_info) {
    return Regex.IsMatch(file_info.Name,
                         @"^.*\.generated\.h$|^.*\.mathematica\.h$");
  }

  public static bool IsTest(FileInfo file_info) {
    return Regex.IsMatch(file_info.Name, @"_test\.cpp$");
  }
}

} // namespace sourcerer
} // namespace principia
