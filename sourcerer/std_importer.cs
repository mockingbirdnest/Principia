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
      var first_ns = FindFirstNamespace(parser_file);
      if (first_ns == null) {
        throw new ArgumentException("File " +
                                    input_file.FullName +
                                    " has no namespace");
      }
      bool has_system_headers = FixSystemHeaders(parser_file, includes);
      FixImportStd(parser_file, first_ns, has_system_headers);
      RewriteFile(input_file, parser_file, dry_run);
    }
  }

  private static Namespace? FindFirstNamespace(
      Parser.File file) {
    foreach (Node child in file.children) {
      if (child is Namespace ns) {
        return ns;
      }
    }
    return null;
  }

  private static bool FixSystemHeaders(Parser.File file,
                                       List<Include> includes) {
    var system_includes = includes.FindAll(inc => inc.is_system);
    bool has_system_headers = system_includes.Count > 0;
    foreach (Include inc in system_includes) {
      int position_in_parent = inc.position_in_parent;
      file.RemoveChild(inc);
      // If removing the system include leaves two consecutive empty lines,
      // remove one of them.
      if (file.children[position_in_parent - 1] is Text t1 &&
          t1.text == "" &&
          file.children[position_in_parent] is Text t2 &&
          t2.text == "") {
        file.RemoveChild(t1);
      }
    }
    return has_system_headers;
  }

  private static void FixImportStd(Parser.File file,
                                   Namespace first_ns,
                                   bool has_system_headers) {
    // If there are no system headers, don't touch anything.  That includes the
    // case where the file has already been converted.
    if (has_system_headers) {
      // Find the `import std;` statements, if any.
      List<Import> std_imports = new List<Import>();
      foreach (Node child in file.children) {
        if (child is Import{ name: "std" } std_import) {
          std_imports.Add(std_import);
        }
      }

      // Remove said import statements.
      foreach (Import std_import in std_imports) {
        file.RemoveChild(std_import);
      }

      // Insert an `import std;` if needed, just before the first namespace.
      int first_ns_position = first_ns.position_in_parent;
      string import_line = "import std;";
      var import = new Import(import_line, parent: null, name: "std");
      file.InsertChild(first_ns_position, import);
      var blank = new Text("", parent: null);
      file.InsertChild(first_ns_position + 1, blank);
    }
  }
}

} // namespace sourcerer
} // namespace principia
