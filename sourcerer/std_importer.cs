using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using static principia.sourcerer.Analyser;
using static principia.sourcerer.Filenames;
using static principia.sourcerer.Parser;
using static principia.sourcerer.Rewriter;


namespace principia {
namespace sourcerer {

class StdImporter {

  public static void Run(string[] args) {
    // Parse the arguments.
    DirectoryInfo? project = null;
    var clients = new List<DirectoryInfo>();
    var excluded = new HashSet<string>();
    bool dry_run = true;
    foreach (string arg in args) {
      if (arg.StartsWith("--") && arg.Contains(":")) {
        string[] split =
            arg.Split(new []{ "--", ":" }, StringSplitOptions.None);
        string option = split[1];
        string value = split[2];
        if (option == "project") {
          project = new DirectoryInfo(value);
        } else if (option == "dry_run") {
          dry_run = bool.Parse(value);
        } else if (option == "exclude") {
          excluded.Add(value);
        } else {
          throw new ArgumentException("Unknown option " + option);
        }
      }
    }
    if (project == null) {
      throw new NullReferenceException();
    }

    // Parse all the files in our project.
    FileInfo[] hpp_files = project.GetFiles("*.hpp");
    FileInfo[] cpp_files = project.GetFiles("*.cpp");
    FileInfo[] body_hpp_files = Array.FindAll(hpp_files, IsBodyHpp);
    FileInfo[] all_body_files = body_hpp_files.Union(cpp_files).ToArray();
    FileInfo[] all_files = hpp_files.Union(cpp_files).ToArray();
    var file_info_to_file = new Dictionary<FileInfo, Parser.File>();
    foreach (FileInfo input_file in all_files) {
      if (excluded.Contains(input_file.Name)) {
        continue;
      }
      Parser.File parser_file =
          Parser.ParseFile(input_file,
                           IsBody(input_file,
                                  extra_headers: new HashSet<string>() ));
      file_info_to_file.Add(input_file, parser_file);
    }

    // Rewrite the files.
    foreach (var kv in file_info_to_file) {
      FileInfo input_file = kv.Key;
      Parser.File parser_file = kv.Value;
      var includes = FindIncludes(parser_file);
      bool has_system_headers = FixSystemHeaders(includes);
      var consecutive_user_includes = FindConsecutiveUserIncludes(includes);
      FixImportStd(parser_file, consecutive_user_includes, has_system_headers);
      RewriteFile(input_file, parser_file, dry_run);
    }
  }

  // Locate the first block of consecutive user includes.  We won't touch
  // subsequent blocks, they probably contain something odd.
  private static List<Include> FindConsecutiveUserIncludes(
      List<Include> includes) {
    var consecutive_user_includes = new List<Include>();
    int? position_in_parent = null;
    foreach (Include inc in includes) {
      if (position_in_parent.HasValue) {
        position_in_parent = position_in_parent.Value + 1;
        if (position_in_parent.Value == inc.position_in_parent) {
          consecutive_user_includes.Add(inc);
        } else {
          break;
        }
      } else {
        position_in_parent = inc.position_in_parent;
        consecutive_user_includes.Add(inc);
      }
    }
    return consecutive_user_includes;
  }

  private static bool FixSystemHeaders(List<Include> includes) {
    var system_includes = includes.FindAll(inc => inc.is_system);
    bool has_system_headers = system_includes.Count > 0;
    int removed = includes.RemoveChildren(system_includes);
    return has_system_headers;
  }

  private static void FixImportStd(Parser.File file,
                                   List<Include> consecutive_user_includes,
                                   bool has_system_headers) {
    // Find the `import std;` statements, if any.
    var import_std = new List<Import>();
    foreach (Node child in file.children) {
      if (child is Import{ name: "std" } import) {
        import_std.Add(import);
      }
    }

    // Remove them all.
    foreach (Import import in import_std) {
      file.RemoveChild(import);
    }

    // Insert an `import std;` if needed, at the end of the first block of
    // consecutive user includes, or at the beginning of the file if there are
    if (has_system_headers) {
      int last_include_position =
          consecutive_user_includes[^1].position_in_parent;
      string import_line = "import std;";
      var import = new Import(import_line, parent: file, name: "std");
      file.InsertChild(last_include_position + 1, import);
    }
  }
}

} // namespace sourcerer
} // namespace principia
