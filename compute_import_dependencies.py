from collections import defaultdict
import glob;
import os;
import re;

# https://eel.is/c++draft/module.import#nt:module-import-declaration.
# We don’t use module partitions for now.
MODULE_IMPORT_DECLARATION = re.compile(
    r'''
      import \s+ (?P<module_name>[^\s;"<>]+) \s* ;
    | import \s* [<"](?P<header_name>[^"]*)[">] \s* ;
    ''',
    re.VERBOSE)
# https://eel.is/c++draft/module.unit#nt:module-declaration
# No module-partition nor attribute-specifier-seq for now.
MODULE_DECLARATION = re.compile(
    r'''
    (?P<export_keyword> export \s+ )? module \s+ (?P<module_name>[^\s:;\[\]]+) \s* ;
    ''',
    re.VERBOSE)

class Dependency:
  def compilation_dependencies(self):
    return []

  def linkage_dependencies(self):
    return []

class TranslationUnit:
  all               : list["TranslationUnit"] = []
  precompiled_units : list["TranslationUnit"] = []
  compiled_units    : list["TranslationUnit"] = []
  linked_units      : list["TranslationUnit"] = []
  def __init__(self, source_file: str):
    self.source_file = source_file
    self.compiled = f"obj/{self.source_file}.o"
    self.binary = (f"bin/{self.source_file.removesuffix('.cpp')}"
                   if source_file.endswith("_test.cpp") or
                      source_file.rsplit("/", 1)[0] == "main.cpp"
                   else None)
    self.precompiled = None
    self.register()

  def register(self):
    TranslationUnit.all.append(self)
    if self.precompiled:
      TranslationUnit.precompiled_units.append(self)
    if self.compiled:
      TranslationUnit.compiled_units.append(self)
    if self.binary:
      TranslationUnit.linked_units.append(self)

  def __eq__(self, value):
    return isinstance(value, TranslationUnit) and self.source_file == value.source_file

  def __hash__(self):
    return self.source_file.__hash__()

  def __str__(self):
    return f"{self.__class__.__name__} {self.source_file}"


class ModuleInterfaceUnit(TranslationUnit):
  def __init__(self, primary_module_interface_unit: str|None, module_name: str):
    self.module_name = module_name
    self.precompiled = f"pmu/{self.module_name}.pcm"
    self.binary = None
    if primary_module_interface_unit is None:
      # Handled in the manually-written Makefile.
      self.source_file = None
      self.compiled = None
    else:
      self.source_file = primary_module_interface_unit
      self.compiled = f"obj/{self.source_file}.o"
      self.register()

class HeaderUnit(TranslationUnit, Dependency):
  def __init__(self, header_name: str):
    self.header_name = header_name
    self.source_file = header_name
    self.precompiled = f"phu/{self.source_file}.pcm"
    # Header units are precompiled, but not compiled.
    self.compiled = None
    self.binary = None
    self.register()

  def compilation_dependencies(self):
    return [self.precompiled]

class NamedModule(Dependency):
  named_modules : dict[str, "NamedModule"] = {}
  def __init__(
      self,
      module_name: str,
      primary_module_interface_unit: ModuleInterfaceUnit):
    self.module_name = module_name
    self.primary_module_interface_unit = primary_module_interface_unit
    self.module_units : list[TranslationUnit] = [primary_module_interface_unit]
    NamedModule.named_modules[module_name] = self

  def compilation_dependencies(self):
    return [self.primary_module_interface_unit.precompiled]

  def linkage_dependencies(self):
    return [u.compiled for u in self.module_units if u.compiled]

NamedModule("std",
            ModuleInterfaceUnit(None, "std"))

for path in glob.glob("*/*.ixx"):
  path = path.replace("\\", "/")
  with open(path, encoding="utf8") as f:
    for l, line in enumerate(f.readlines(), 1):
      module_declaration = MODULE_DECLARATION.match(line)
      if module_declaration:
        if module_declaration.group("export_keyword"):
          module_name = module_declaration.group("module_name")
          module = NamedModule(
              module_name,
              ModuleInterfaceUnit(path, module_name))
        else:
          raise TypeError(
            f"{path}:{l}: Expected .ixx to be a module interface unit, but it"
            f" is a module implementation unit:\n{line}")
        break
    else:
      raise TypeError(
        f"{path}: Expected .ixx to be a module interface unit, but no module"
        f" declaration was found")


for path in glob.glob("*/*.cpp"):
  path = path.replace("\\", "/")
  translation_unit = TranslationUnit(path)
  with open(path, encoding="utf8") as f:
    for l, line in enumerate(f.readlines(), 1):
      module_declaration = MODULE_DECLARATION.match(line)
      if module_declaration:
        if module_declaration.group("export_keyword"):
          raise TypeError(
            f"{path}:{l}: Expected .cpp to be a module implementation unit, but it"
            f" is a module interface unit:\n{line}")
        else:
          NamedModule.named_modules[
              module_declaration.group("module_name")].module_units.append(
              translation_unit)
        break

transitive_dependencies : dict[
    TranslationUnit, dict[Dependency, None]] = defaultdict(dict)

def add_transitive_imports(translation_unit: TranslationUnit):
  if translation_unit in transitive_dependencies:
    return
  if not os.path.exists(translation_unit.source_file):
    if not isinstance(translation_unit, HeaderUnit):
      raise FileNotFoundError(f"{translation_unit} not found")
    # This must be a third-party header; these don’t have imports at this time.
    return
  with open(translation_unit.source_file, encoding="utf8") as f:
    for line in f.readlines():
      import_declaration = re.match(MODULE_IMPORT_DECLARATION, line)
      if import_declaration:
        module_name = import_declaration.group("module_name")
        if module_name:
          transitive_dependencies[translation_unit][
              NamedModule.named_modules[module_name]] = None
        else:
          imported_unit = HeaderUnit(import_declaration.group("header_name"))
          transitive_dependencies[translation_unit][imported_unit] = None
          if imported_unit.source_file:
            add_transitive_imports(imported_unit)
            transitive_dependencies[translation_unit].update(transitive_dependencies[imported_unit])

for unit in TranslationUnit.all:
  add_transitive_imports(unit)

for unit in sorted(TranslationUnit.precompiled_units, key=lambda u: u.source_file):
  if not os.path.exists(unit.source_file):
    # This must be a third-party dependency, let the makefile provide its path.
    continue
  rule = f"{unit.precompiled}: {unit.source_file}"
  for d in transitive_dependencies[unit]:
    for precompiled in d.compilation_dependencies():
      rule += f" \\\n\t{precompiled}"
  print(rule)

for unit in sorted(TranslationUnit.compiled_units, key=lambda u: u.source_file):
  rule = f"{unit.compiled}: {unit.precompiled or unit.source_file}"
  for d in transitive_dependencies[unit]:
    for precompiled in d.compilation_dependencies():
      rule += f" \\\n\t{precompiled}"
  print(rule)

for unit in sorted(TranslationUnit.linked_units, key=lambda u: u.source_file):
  rule = f"{unit.binary}: {unit.compiled}"
  for d in transitive_dependencies[unit]:
    for compiled in d.linkage_dependencies():
      rule += f" \\\n\t{compiled}"
  print(rule)
