using System;

namespace principia {
namespace sourcerer {

// Usage:
//   sourcerer renamespacer \
//             --project:quantities \
//             --client:base --client:physics \
//             --exclude:macros.hpp --dry_run:false
// This will renamespace quantities and fix the references in the client
// projects.  The files will be overwritten.
internal class Sourcerer {
  static void Main(string[] args) {
    if (args.Length == 0) {
      throw new ArgumentException("Missing command");

    }
    string[] args_after_0 = new string[args.Length - 1];
    Array.Copy(args, 1, args_after_0, 0, args_after_0.Length);
    if (args[0] == "include_what_you_using") {
      IncludeWhatYouUsing.Run(args_after_0);
    } else if (args[0] == "make_editorconfig") {
      MakeEditorconfig.Run(args_after_0);
    } else if (args[0] == "renamespacer") {
      Renamespacer.Run(args_after_0);
    } else {
      throw new ArgumentException("Unknown command " + args[0]);
    }
  }
}

} // namespace sourcerer
} // namespace principia
