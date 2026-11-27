using System.IO;
using static principia.sourcerer.Parser;

namespace principia {
namespace sourcerer {

public class Rewriter {
  public static void RewriteFile(FileInfo input_file,
                                 Parser.File file,
                                 bool dry_run,
                                 bool move_to_ixx = false) {
    string input_filename = input_file.FullName;
    string output_filename =
        input_file.DirectoryName + "\\" + input_file.Name + ".new";

    using (StreamWriter writer = System.IO.File.CreateText(output_filename)) {
      RewriteNode(writer, file);
    }
    if (!dry_run) {
      if (move_to_ixx) {
        string ixx_filename =
            input_file.DirectoryName +
            "\\" +
            Path.GetFileNameWithoutExtension(input_file.Name) +
            ".ixx";
        System.IO.File.Move(output_filename, ixx_filename, overwrite: true);
      } else {
        System.IO.File.Move(output_filename, input_filename, overwrite: true);
      }
    }
  }

  private static void RewriteNode(StreamWriter writer, Parser.Node node) {
    foreach (Parser.Node child in node.children) {
      writer.WriteLine(child.must_rewrite ? child.Cxx() : child.text);
      if (child is Scope s) {
        RewriteNode(writer, child);
        writer.WriteLine(s.must_rewrite ? s.ClosingCxx() : s.closing_text);
      }
    }
  }
}

} // namespace sourcerer
} // namespace principia
