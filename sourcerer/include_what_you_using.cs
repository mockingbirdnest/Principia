using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
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

    foreach (DirectoryInfo project in projects) {
      // Parse all the files in this project.
      var file_info_to_file = new Dictionary<FileInfo, Parser.File>();
      var file_name_to_file_info = new Dictionary<string, FileInfo>();
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
        file_name_to_file_info.Add(Path.GetFullPath(input_file.Name,
                                                    input_file.DirectoryName),
                                   input_file);
      }

      // Map the bodies to their headers.
      var corresponding_header = new Dictionary<FileInfo, FileInfo>();
      foreach (FileInfo input_file in all_body_files) {
        if (Filenames.IsTest(input_file)) {
          continue;
        }
        string corresponding_header_file_name =
            Filenames.CorrespondingHeader(input_file);
        var corresponding_header_file =
            file_name_to_file_info[corresponding_header_file_name];
        if (corresponding_header_file != input_file) {
          corresponding_header.Add(input_file, corresponding_header_file);
        }
      }

      // Remove redundant using directives from the bodies.
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
    // Find the using directives that are present in the body and in the header.
    // They are redundant in the body.
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

    // Remove the redundant using directives.
    bool all_body_using_directives_are_unneeded =
        unneeded_body_using_directives.ToList().Count ==
        body_using_directives.Count;
    Node following_node_in_parent = null;
    foreach (UsingDirective ud in unneeded_body_using_directives) {
      var parent = ud.parent;
      Debug.Assert(parent is Namespace,
                   "Internal using directive not within a namespace");
      int ud_position_in_parent = ud.position_in_parent;
      var preceding_nodes_in_parent =
          parent.children.Take(ud_position_in_parent).ToList();
      var following_nodes_in_parent = parent.children.
          Skip(ud_position_in_parent + 1).ToList();
      following_node_in_parent = following_nodes_in_parent[0];
      parent.children = preceding_nodes_in_parent;
      parent.AddChildren(following_nodes_in_parent);
    }

    // If all the using directives are gone, remove any extra blank line that
    // followed them.
    if (all_body_using_directives_are_unneeded &&
        following_node_in_parent is Text{ text: "" } blank_line) {
      var parent = blank_line.parent;
      var blank_line_position_in_parent = blank_line.position_in_parent;
      var preceding_nodes_in_parent =
          parent.children.Take(blank_line_position_in_parent).ToList();
      var following_nodes_in_parent = parent.children.
          Skip(blank_line_position_in_parent + 1).ToList();
      parent.children = preceding_nodes_in_parent;
      parent.AddChildren(following_nodes_in_parent);
    }
  }
}

}  // namespace sourcerer
}  // namespace principia
