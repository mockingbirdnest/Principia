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
      writer.WriteLine(child.must_rewrite ? child.Cxx() : child.text);
      if (child is Namespace ns) {
        RewriteNode(writer, child);
        writer.WriteLine(ns.must_rewrite ? ns.ClosingCxx() : ns.closing_text);
      }
    }
  }
}

} // namespace sourcerer
} // namespace principia
