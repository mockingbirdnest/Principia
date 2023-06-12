using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using static principia.sourcerer.Analyser;
using static principia.sourcerer.Filenames;
using static principia.sourcerer.Parser;
using static principia.sourcerer.Rewriter;

namespace principia {
namespace sourcerer {

class IncludeWhatYouUsing {

  public static void Run(string[] args) {
    // Parse the arguments.
    var projects = new List<DirectoryInfo>();
    var excluded = new HashSet<string>();
    var extra_headers = new HashSet<string>();
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
        } else if (option == "extra_header") {
          extra_headers.Add(value);
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
    var corresponding_header = new Dictionary<FileInfo, FileInfo>();
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
        if (Filenames.IsBody(input_file, extra_headers)) {
          var corresponding_header_file =
              new FileInfo(Filenames.CorrespondingHeader(input_file));
          corresponding_header.Add(input_file, corresponding_header_file);
        }
        Parser.File parser_file =
            Parser.ParseFile(input_file, IsBody(input_file, extra_headers));
        file_info_to_file.Add(input_file, parser_file);
      }
      foreach (FileInfo input_file in all_body_files) {
        if (corresponding_header.ContainsKey(input_file)) {
          var parser_file = file_info_to_file[input_file];
          FixRedundantUsingDirectives(parser_file,
                                      file_info_to_file[
                                          corresponding_header[input_file]]);
          RewriteFile(input_file, parser_file, dry_run);
        }
      }
    }
  }

  private static void FixRedundantUsingDirectives(
      Parser.File body,
      Parser.File header) {
    List<UsingDirective> body_using_directives =
        FindUsingDirectives(body, internal_only: true);
    List<UsingDirective> header_using_directives =
        FindUsingDirectives(header, internal_only: true);
    var common_ns =
        (from ud in body_using_directives select ud.ns).Intersect(
            from ud in header_using_directives select ud.ns);
    var unneeded_body_using_directives = from ud in body_using_directives
                                         where common_ns.Contains(ud.ns)
                                         select ud;

    foreach (UsingDirective ud in unneeded_body_using_directives) {
      var parent = ud.parent;
      Debug.Assert(parent is Namespace,
                   "Internal using directive not within a namespace");
      int internal_position_in_parent = ud.position_in_parent;
      var preceding_nodes_in_parent =
          parent.children.Take(internal_position_in_parent).ToList();
      var following_nodes_in_parent = parent.children.
          Skip(internal_position_in_parent + 1).ToList();
      parent.children = preceding_nodes_in_parent;
      parent.AddChildren(following_nodes_in_parent);
    }
  }
}

}  // namespace sourcerer
}  // namespace principia
