from collections import defaultdict
import glob;
import os;
import re;
import sys;

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

class Importable:
  def precompiled_units(self):
    return []

  def compiled_units(self):
    return []

class TranslationUnit:
  all                   : list["TranslationUnit"] = []
  all_precompiled_units : list["TranslationUnit"] = []
  all_compiled_units    : list["TranslationUnit"] = []
  all_linked_units      : list["TranslationUnit"] = []
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
      TranslationUnit.all_precompiled_units.append(self)
    if self.compiled:
      TranslationUnit.all_compiled_units.append(self)
    if self.binary:
      TranslationUnit.all_linked_units.append(self)

  def __eq__(self, value):
    return isinstance(value, TranslationUnit) and self.source_file == value.source_file

  def __hash__(self):
    return self.source_file.__hash__()

  def __repr__(self):
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

class HeaderUnit(TranslationUnit, Importable):
  def __init__(self, header_name: str):
    self.header_name = header_name
    self.source_file = header_name
    self.precompiled = f"phu/{self.source_file}.pcm"
    # Header units are precompiled, but not compiled.
    self.compiled = None
    self.binary = None
    self.register()

  def precompiled_units(self):
    return [self]

class NamedModule(Importable):
  named_modules : dict[str, "NamedModule"] = {}
  def __init__(
      self,
      module_name: str,
      primary_module_interface_unit: ModuleInterfaceUnit):
    self.module_name = module_name
    self.primary_module_interface_unit = primary_module_interface_unit
    self.module_units : list[TranslationUnit] = [primary_module_interface_unit]
    NamedModule.named_modules[module_name] = self

  def precompiled_units(self):
    return [self.primary_module_interface_unit]

  def compiled_units(self):
    return [u for u in self.module_units if u.compiled]

  def __repr__(self):
    return f"NamedModule {self.module_name}: {self.module_units}"

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

# We do not do preprocessing, so we cannot tell that these headers are not
# actually imported on non-Windows.
SKIPPED_HEADERS = set((
    "intrin.h",
))

# We use dict[T, None] throughout as a set that retains insertion order, so that
# the output Makefile is deterministic.
imports : dict[
    TranslationUnit, dict[Importable, None]] = defaultdict(dict)

def transitively_find_imports(translation_unit: TranslationUnit):
  if translation_unit in imports:
    return
  if not os.path.exists(translation_unit.source_file):
    if not isinstance(translation_unit, HeaderUnit):
      raise FileNotFoundError(f"{translation_unit} not found")
    # This must be a third-party header; these don’t have imports at this time.
    return
  with open(translation_unit.source_file, encoding="utf8") as f:
    for line in f.readlines():
      module_name = None
      module_declaration = MODULE_DECLARATION.match(line)
      if module_declaration and not module_declaration.group("export_keyword"):
        # This is a module implementation unit, it depends on the module
        # interface unit.
        module_name = module_declaration.group("module_name")
      import_declaration = re.match(MODULE_IMPORT_DECLARATION, line)
      if import_declaration:
        module_name = import_declaration.group("module_name")
      if module_name:
        module = NamedModule.named_modules[module_name]
        imports[translation_unit][module] = None
        primary_module_interface_unit = module.primary_module_interface_unit
      elif import_declaration:
        header_name = import_declaration.group("header_name")
        if header_name in SKIPPED_HEADERS:
          continue
        imported_unit = HeaderUnit(header_name)
        imports[translation_unit][imported_unit] = None

for unit in TranslationUnit.all:
  transitively_find_imports(unit)

# The graph of imports yields two distinct dependency graphs, one for
# (pre-)compilation and one for linkage: dependencies on a module need to be
# (pre-)compiled against the precompiled transitive dependencies of the module
# interface units, but need to be linked against the compiled transitive
# dependencies of both the module interface units and the module implementation
# units.  We now compute both transitive closures.

precompiled_dependencies : dict[
    TranslationUnit, dict[Importable, None]] = defaultdict(dict)

compiled_dependencies : dict[
    TranslationUnit, dict[Importable, None]] = defaultdict(dict)

def compute_precompiled_dependencies(translation_unit: TranslationUnit):
  if translation_unit in precompiled_dependencies:
    return
  for direct_dependency in imports[translation_unit]:
    for unit in direct_dependency.precompiled_units():
      precompiled_dependencies[translation_unit][direct_dependency] = None
      compute_precompiled_dependencies(unit)
      precompiled_dependencies[translation_unit].update(
          precompiled_dependencies[unit])

def compute_compiled_dependencies(translation_unit: TranslationUnit):
  if translation_unit in compiled_dependencies:
    return
  for direct_dependency in imports[translation_unit]:
    for unit in direct_dependency.compiled_units():
      compiled_dependencies[translation_unit][direct_dependency] = None
      compute_compiled_dependencies(unit)
      compiled_dependencies[translation_unit].update(
          compiled_dependencies[unit])

for unit in TranslationUnit.all:
  compute_precompiled_dependencies(unit) 
  compute_compiled_dependencies(unit) 

for unit in sorted(TranslationUnit.all_precompiled_units,
		   key=lambda u: u.source_file):
  if not os.path.exists(unit.source_file):
    # This must be a third-party dependency, let the makefile provide its path.
    continue
  rule = f"{unit.precompiled}: {unit.source_file}"
  for d in precompiled_dependencies[unit]:
    for u in d.precompiled_units():
      rule += f" \\\n\t{u.precompiled}"
  print(rule)

for unit in sorted(TranslationUnit.all_compiled_units,
	           key=lambda u: u.source_file):
  rule = f"{unit.compiled}: {unit.precompiled or unit.source_file}"
  for d in precompiled_dependencies[unit]:
    for u in d.precompiled_units():
      rule += f" \\\n\t{u.precompiled}"
  print(rule)

for unit in sorted(TranslationUnit.all_linked_units,
                   key=lambda u: u.source_file):
  rule = f"{unit.binary}: {unit.compiled}"
  for d in compiled_dependencies[unit]:
    for u in d.compiled_units():
      rule += f" \\\n\t{u.compiled}"
  print(rule)
