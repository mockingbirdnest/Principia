using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static principia.sourcerer.Filenames;

namespace principia {
namespace sourcerer {

class IncludeWhatYouUsing {

  public static void Run(string[] args) {
    // Parse the arguments.
    var projects = new List<DirectoryInfo>();
    var excluded = new HashSet<string>();
    bool dry_run = true;
    foreach (string arg in args) {
      if (arg.StartsWith("--") && arg.Contains(":")) {
        string[] split =
            arg.Split(new []{ "--", ":" }, StringSplitOptions.None);
        string option = split[1];
        string value = split[2];
        if (option == "project") {
          projects.Add(new DirectoryInfo(value));
        } else if (option == "dry_run") {
          dry_run = bool.Parse(value);
        } else if (option == "exclude") {
          excluded.Add(value);
        } else {
          throw new ArgumentException("Unknown option " + option);
        }
      }
    }
    if (projects.Count == 0) {
      throw new NullReferenceException();
    }

    // Parse all the files in our projects.
    var file_info_to_file = new Dictionary<FileInfo, Parser.File>();
    foreach (DirectoryInfo project in projects) {
      FileInfo[] hpp_files = project.GetFiles("*.hpp");
      FileInfo[] cpp_files = project.GetFiles("*.cpp");
      FileInfo[] body_hpp_files = Array.FindAll(hpp_files, IsBodyHpp);
      FileInfo[] all_body_files = body_hpp_files.Union(cpp_files).ToArray();
      FileInfo[] all_files = hpp_files.Union(cpp_files).ToArray();
      foreach (FileInfo input_file in all_files) {
        if (excluded.Contains(input_file.Name)) {
          continue;
        }
        Parser.File parser_file =
            Parser.ParseFile(input_file, IsBody(input_file, extra_headers));
        file_info_to_file.Add(input_file, parser_file);
      }
    }
  }
}

}  // namespace sourcerer
}  // namespace principia
