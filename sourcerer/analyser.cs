using System.Collections.Generic;
using static principia.sourcerer.Parser;

namespace principia {
namespace sourcerer {

// A class for extracting information from a syntactic tree.  Does not change
// the tree.
public class Analyser {
  // Returns the list of exported (public) declarations in |node| and its
  // descendants.  Note that internal namespace do *not* export their
  // declarations.
  public static List<Declaration> CollectExportedDeclarations(Node node) {
    var exported_declarations = new List<Declaration>();
    foreach (Node child in node.children) {
      if (child is Namespace ns) {
        if (!ns.is_internal) {
          exported_declarations.AddRange(CollectExportedDeclarations(child));
        }
      } else if (child is Declaration decl) {
        exported_declarations.Add(decl);
      }
    }
    return exported_declarations;
  }


  // Given a |using_declaration| and map from (exported) declarations to the
  // file where they are declared, returns the file where the name referenced by
  // the using declaration is declared.
  public static File? FindFileReferencedByUsingDeclaration(
      UsingDeclaration using_declaration,
      Dictionary<Declaration, File> declaration_to_file) {
    var referenced_name = using_declaration.name;
    var referenced_innermost_namespace =
        using_declaration.full_name.Split("::")[^2];
    foreach (var pair in declaration_to_file) {
      var declaration = pair.Key;
      // Poor man's name resolution.
      if (declaration.name == referenced_name &&
          declaration.parent is Namespace ns &&
          ns.name == referenced_innermost_namespace) {
        return pair.Value;
      }
    }
    return null;
  }


  // Returns the list of innermost namespaces found in |node| and its
  // descendants.
  public static List<Namespace> FindInnermostNamespaces(Node node) {
    var innermost_namespaces = new List<Namespace>();
    foreach (Node child in node.children) {
      if (child is Namespace ns) {
        var nested_namespaces = FindInnermostNamespaces(child);
        if (nested_namespaces.Count == 0) {
          innermost_namespaces.Add(ns);
        } else {
          innermost_namespaces.AddRange(nested_namespaces);
        }
      }
    }
    return innermost_namespaces;
  }

  // Returns the last outermost namespace found in |file|.
  // TODO(phl): This should probably be changed to return all the outermost
  // namespaces.
  public static Namespace? FindLastOutermostNamespace(File file) {
    Namespace? last_outermost_namespace = null;
    foreach (Node child in file.children) {
      if (child is Namespace ns) {
        last_outermost_namespace = ns;
      }
    }
    return last_outermost_namespace;
  }

  // Returns the legacy internal namespaces found in |node| and it descendants.
  // A legacy internal namespace is one whose name starts with "internal_".

  public static List<Namespace> FindLegacyInternalNamespaces(Node node) {
    var legacy_internal_namespaces = new List<Namespace>();
    foreach (Node child in node.children) {
      if (child is Namespace internal_namespace &&
          internal_namespace.name.StartsWith("internal_")) {
        legacy_internal_namespaces.Add(internal_namespace);
      } else if (child is Namespace ns) {
        legacy_internal_namespaces.AddRange(
            FindLegacyInternalNamespaces(child));
      }
    }
    return legacy_internal_namespaces;
  }

  // Returns the using declarations found in |node| and its descendants.  If
  // |internal_only| is true, internal namespaces are ignored (that is, only the
  // exported using declarations are returned).
  public static List<UsingDeclaration> FindUsingDeclarations(
      Node node,
      bool internal_only) {
    var internal_using_declarations = new List<UsingDeclaration>();
    foreach (Node child in node.children) {
      if (child is Namespace ns) {
        if (ns.is_internal || !internal_only) {
          foreach (Node grandchild in child.children) {
            if (grandchild is UsingDeclaration ud) {
              internal_using_declarations.Add(ud);
            }
          }
        }
        internal_using_declarations.AddRange(
            FindUsingDeclarations(child, internal_only));
      }
    }
    return internal_using_declarations;
  }

  // Returns the using directives found in |node| and its descendants.  If
  // |internal_only| is true, internal namespaces are ignored (that is, only the
  // exported using directives are returned).
  public static List<UsingDirective> FindUsingDirectives(
      Node node,
      bool internal_only) {
    var internal_using_directives = new List<UsingDirective>();
    foreach (Node child in node.children) {
      if (child is Namespace ns) {
        if (ns.is_internal || !internal_only) {
          foreach (Node grandchild in child.children) {
            if (grandchild is UsingDirective ud) {
              internal_using_directives.Add(ud);
            }
          }
        }
        internal_using_directives.AddRange(
            FindUsingDirectives(child, internal_only));
      }
    }
    return internal_using_directives;
  }
}

} // namespace sourcerer
} // namespace principia
