using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using static principia.sourcerer.Analyser;
using static principia.sourcerer.Filenames;
using static principia.sourcerer.Parser;
using static principia.sourcerer.Rewriter;

namespace principia {
namespace sourcerer {

public class StringArrayComparer : IComparer<string[]> {
  public int Compare(string[] left, string[] right) {
    int i = 0;
    while (i < left.Length && i < right.Length) {
      int i_compare = string.Compare(left[i], right[i]);
      if (i_compare != 0) {
        return i_compare;
      }
      ++i;
    }
    return left.Length - right.Length;
  }
}

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
        string corresponding_header_file_name =
            Filenames.CorrespondingHeader(input_file);
        // Benchmarks and tests don't have a corresponding header.
        if (file_name_to_file_info.ContainsKey(
                corresponding_header_file_name)) {
          var corresponding_header_file =
              file_name_to_file_info[corresponding_header_file_name];
          if (corresponding_header_file != input_file) {
            corresponding_header.Add(input_file, corresponding_header_file);
          }
        }
      }

      // Remove redundant using directives from the bodies.
      foreach (FileInfo input_file in all_body_files) {
        if (corresponding_header.ContainsKey(input_file)) {
          FixRedundantUsingDirectives(file_info_to_file[input_file],
                                      file_info_to_file[
                                          corresponding_header[input_file]]);
        }
      }

      foreach (FileInfo input_file in all_files) {
          var parser_file = file_info_to_file[input_file];
          FixIncludes(parser_file);
          RewriteFile(input_file, parser_file, dry_run);
      }
    }
  }

  private static void FixIncludes(Parser.File file) {
    var existing_includes =
        (from inc in FindIncludes(file)
        where !inc.is_own_body && !inc.is_own_header
        select inc).ToArray();
    List<UsingDirective> using_directives =
        FindUsingDirectives(file, internal_only: true);

    // Build the sorted set of includes that are needed based on the namespaces
    // mentioned in using directives.
    var using_namespaces = from ud in using_directives select ud.ns;
    var using_paths = new SortedSet<string[]>(new StringArrayComparer());
    foreach (string ns in using_namespaces) {
      var segments = ns.Split("::");
      var using_path = new string[]{};
      foreach (string segment in segments) {
        if (segment == "principia") {
          continue;
        } else if (segment[0] == '_') {
          string header_filename = Regex.Replace(segment, @"^_", "");
          using_path = using_path.Append(header_filename).ToArray();
        } else {
          using_path = using_path.Append(segment).ToArray();
        }
      }
      using_paths.Add(using_path);
    }

    var new_includes = new List<Node>();
    foreach (string[] using_path in using_paths) {
      bool found = false;
      foreach (Include inc in existing_includes) {
        if (inc.path == using_path) {
          // There is already an include for this path, reuse it.
          new_includes.Add(inc);
          found = true;
        }
      }
      // If there is no include for this namespace, add one (in order).
      if (!found) {
        new_includes.Add(new Include(parent: null, using_path, file.file_info));
      }
    }

    int first_include_position = existing_includes[0].position_in_parent;
    int last_include_position = existing_includes[existing_includes.Length - 1].
        position_in_parent;
    var preceding_nodes_in_file =
        file.children.Take(first_include_position).ToList();
    var following_nodes_in_file = file.children.
        Skip(last_include_position + 1).ToList();
    file.children = preceding_nodes_in_file;
    file.AddChildren(new_includes);
    file.AddChildren(following_nodes_in_file);
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

      // Don't process using directives that are not in the namespace for this
      // file.  They probably indicate that we are reopening someone else's
      // namespace, and we wouldn't want to touch that.
      var ns1 = parent as Namespace;
      if (ns1.parent is Namespace ns2 &&
          ns2.name != body.file_namespace_simple_name) {
        continue;
      }

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
