using System.IO;
using static principia.sourcerer.Parser;

namespace principia {
namespace sourcerer {

public class Rewriter {
  public static void RewriteFile(FileInfo input_file,
                                 Parser.File file,
                                 bool dry_run) {
    string input_filename = input_file.FullName;
    string output_filename =
        input_file.DirectoryName + "\\" + input_file.Name + ".new";

    using (StreamWriter writer = System.IO.File.CreateText(output_filename)) {
      RewriteNode(writer, file);
    }
    if (!dry_run) {
      System.IO.File.Move(output_filename, input_filename, overwrite: true);
    }
  }

  private static void RewriteNode(StreamWriter writer, Parser.Node node) {
    foreach (Parser.Node child in node.children) {
      if (child.must_rewrite) {
        writer.WriteLine(child.Cxx());
      } else {
        writer.WriteLine(child.text);
      }
      if (child is Namespace ns) {
        RewriteNode(writer, child);
        if (ns.must_rewrite) {
          writer.WriteLine(ns.ClosingCxx());
        } else {
          writer.WriteLine(ns.closing_text);
        }
      }
    }
  }
}

} // namespace sourcerer
} // namespace principia
