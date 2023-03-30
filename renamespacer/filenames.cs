using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace principia {
namespace renamespacer {

public class Filenames {

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

  public static bool IsTest(FileInfo file_info) {
    return Regex.IsMatch(file_info.Name, @"_test\.cpp$");
  }
}

} // namespace renamespacer
} // namespace principia
