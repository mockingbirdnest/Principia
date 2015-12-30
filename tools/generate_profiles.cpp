#include "tools/generate_profiles.hpp"

#include <algorithm>
#include <cctype>
#include <experimental/filesystem>  // NOLINT
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "glog/logging.h"
#include "google/protobuf/descriptor.h"
#include "serialization/journal.pb.h"

namespace principia {

using ::google::protobuf::Descriptor;
using ::google::protobuf::FieldDescriptor;
using ::google::protobuf::FieldOptions;
using ::google::protobuf::FileDescriptor;

namespace tools {

namespace {

char const kIn[] = "In";
char const kReturn[] = "Return";
char const kOut[] = "Out";

std::string ToLower(std::string const& s) {
  std::string lower(s.size(), ' ');
  for (int i = 0; i < s.size(); ++i) {
    lower[i] = std::tolower(s[i]);
  }
  return lower;
}

class Generator {
 public:
  void ProcessMethods();

  std::vector<std::string> GetCppMethodTypes() const;

 private:
  void ProcessRepeatedMessageField(FieldDescriptor const* descriptor);

  void ProcessOptionalInt32Field(FieldDescriptor const* descriptor);

  void ProcessRequiredFixed64Field(FieldDescriptor const* descriptor);
  void ProcessRequiredMessageField(FieldDescriptor const* descriptor);
  void ProcessRequiredBoolField(FieldDescriptor const* descriptor);
  void ProcessRequiredDoubleField(FieldDescriptor const* descriptor);
  void ProcessRequiredInt32Field(FieldDescriptor const* descriptor);

  void ProcessSingleStringField(FieldDescriptor const* descriptor);

  void ProcessOptionalField(FieldDescriptor const* descriptor);
  void ProcessRepeatedField(FieldDescriptor const* descriptor);
  void ProcessRequiredField(FieldDescriptor const* descriptor);
  void ProcessSingleField(FieldDescriptor const* descriptor);

  void ProcessField(FieldDescriptor const* descriptor);

  void ProcessInOut(Descriptor const* descriptor,
                    std::vector<FieldDescriptor const*>* field_descriptors);
  void ProcessReturn(Descriptor const* descriptor);

  void ProcessMethodExtension(Descriptor const* descriptor);

 private:
  std::map<FieldDescriptor const*, std::string> serializer_name_;
  std::map<FieldDescriptor const*, std::string> size_field_name_;
  std::set<FieldDescriptor const*> in_out_field_;
  std::map<FieldDescriptor const*, std::string> cpp_field_copy_;
  std::map<FieldDescriptor const*, std::string> cpp_field_setter_;
  std::map<FieldDescriptor const*, std::string> cpp_field_type_;
  std::map<Descriptor const*, std::string> cpp_fill_body_;
  std::map<Descriptor const*, std::string> cpp_method_type_;
  std::map<Descriptor const*, std::string> cpp_nested_type_;
};

void Generator::ProcessMethods() {
  // Get the file containing |Method|.
  Descriptor const* method_descriptor = serialization::Method::descriptor();
  FileDescriptor const* file_descriptor = method_descriptor->file();

  // Process all the messages in that file.
  for (int i = 0; i < file_descriptor->message_type_count(); ++i) {
    Descriptor const* message_descriptor = file_descriptor->message_type(i);
    switch (message_descriptor->extension_count()) {
      case 0: {
        // An auxiliary message that is not an extension.  Nothing to do.
        break;
      }
      case 1: {
        // An extension.  Check that it extends |Method|.
        FieldDescriptor const* extension = message_descriptor->extension(0);
        CHECK(extension->is_extension());
        Descriptor const* containing_type = extension->containing_type();
        CHECK_EQ(method_descriptor, containing_type)
            << message_descriptor->name() << " extends a message other than "
            << method_descriptor->name() << ": " << containing_type->name();
        ProcessMethodExtension(message_descriptor);
        break;
      }
      default: {
        LOG(FATAL) << message_descriptor->name() << " has "
                   << message_descriptor->extension_count() << " extensions";
      }
    }
  }
}

std::vector<std::string> Generator::GetCppMethodTypes() const {
  std::vector<std::string> result;
  for (auto const& pair : cpp_method_type_) {
    result.push_back(pair.second);
  }
  return result;
}

void Generator::ProcessRepeatedMessageField(FieldDescriptor const* descriptor) {
  cpp_field_type_[descriptor] = descriptor->message_type()->name() +
                                " const*";
}

void Generator::ProcessOptionalInt32Field(FieldDescriptor const* descriptor) {
  cpp_field_type_[descriptor] = "int const*";
}

void Generator::ProcessRequiredFixed64Field(FieldDescriptor const* descriptor) {
  FieldOptions const& options = descriptor->options();
  CHECK(options.HasExtension(serialization::pointer_to))
      << descriptor->full_name() << " is missing a pointer_to option";
  std::string const& pointer_to =
      options.GetExtension(serialization::pointer_to);
  cpp_field_type_[descriptor] = pointer_to + "*";
  cpp_field_setter_[descriptor] = "set_" + descriptor->name();
  serializer_name_[descriptor] = "SerializePointer";
}

void Generator::ProcessRequiredMessageField(FieldDescriptor const* descriptor) {
  std::string const& message_type_name = descriptor->message_type()->name();
  cpp_field_type_[descriptor] = message_type_name;
  cpp_field_setter_[descriptor] =
      "mutable_" + descriptor->name() + "()->CopyFrom";
  serializer_name_[descriptor] = "Serialize" + message_type_name;
}

void Generator::ProcessRequiredBoolField(FieldDescriptor const* descriptor) {
  cpp_field_type_[descriptor] = descriptor->cpp_type_name();
  cpp_field_setter_[descriptor] = "set_" + descriptor->name();
}

void Generator::ProcessRequiredDoubleField(FieldDescriptor const* descriptor) {
  cpp_field_type_[descriptor] = descriptor->cpp_type_name();
  cpp_field_setter_[descriptor] = "set_" + descriptor->name();
}

void Generator::ProcessRequiredInt32Field(FieldDescriptor const* descriptor) {
  cpp_field_type_[descriptor] = "int";
  cpp_field_setter_[descriptor] = "set_" + descriptor->name();
}

void Generator::ProcessSingleStringField(FieldDescriptor const* descriptor) {
  cpp_field_type_[descriptor] = "char const*";
}

void Generator::ProcessOptionalField(FieldDescriptor const* descriptor) {
  switch (descriptor->type()) {
    case FieldDescriptor::TYPE_INT32:
      ProcessOptionalInt32Field(descriptor);
      break;
    case FieldDescriptor::TYPE_STRING:
      ProcessSingleStringField(descriptor);
      break;
    default:
      LOG(FATAL) << descriptor->full_name() << " has unexpected type "
                 << descriptor->type_name();
  }
}

void Generator::ProcessRepeatedField(FieldDescriptor const* descriptor) {
  switch (descriptor->type()) {
    case FieldDescriptor::TYPE_MESSAGE:
      ProcessRepeatedMessageField(descriptor);
      break;
    default:
      LOG(FATAL) << descriptor->full_name() << " has unexpected type "
                 << descriptor->type_name();
  }
}

void Generator::ProcessRequiredField(FieldDescriptor const* descriptor) {
  switch (descriptor->type()) {
    case FieldDescriptor::TYPE_BOOL:
      ProcessRequiredBoolField(descriptor);
      break;
    case FieldDescriptor::TYPE_DOUBLE:
      ProcessRequiredDoubleField(descriptor);
      break;
    case FieldDescriptor::TYPE_FIXED64:
      ProcessRequiredFixed64Field(descriptor);
      break;
    case FieldDescriptor::TYPE_INT32:
      ProcessRequiredInt32Field(descriptor);
      break;
    case FieldDescriptor::TYPE_MESSAGE:
      ProcessRequiredMessageField(descriptor);
      break;
    case FieldDescriptor::TYPE_STRING:
      ProcessSingleStringField(descriptor);
      break;
    default:
      LOG(FATAL) << descriptor->full_name() << " has unexpected type "
                 << descriptor->type_name();
  }

  // For in-out fields the data is actually passed with an extra level of
  // indirection.
  if (in_out_field_.find(descriptor) != in_out_field_.end()) {
    cpp_field_type_[descriptor] += "*";
  }
}

void Generator::ProcessSingleField(FieldDescriptor const* descriptor) {
  bool const is_in_out = in_out_field_.find(descriptor) != in_out_field_.end();
  bool const needs_serialization =
      serializer_name_.find(descriptor) != serializer_name_.end();
  cpp_field_copy_[descriptor] = "*to->" + cpp_field_setter_[descriptor] + "(";
  if (needs_serialization) {
    cpp_field_copy_[descriptor] += serializer_name_[descriptor] + "(";
  }
  if (is_in_out) {
    cpp_field_copy_[descriptor] += "*";
  }
  cpp_field_copy_[descriptor] += "from." + descriptor->name();
  if (needs_serialization) {
    cpp_field_copy_[descriptor] += ")";
  }
  cpp_field_copy_[descriptor] += ");";
}

void Generator::ProcessField(FieldDescriptor const* descriptor) {
  switch (descriptor->label()) {
    case FieldDescriptor::LABEL_OPTIONAL:
      ProcessOptionalField(descriptor);
      ProcessSingleField(descriptor);
      break;
    case FieldDescriptor::LABEL_REPEATED:
      ProcessRepeatedField(descriptor);
      break;
    case FieldDescriptor::LABEL_REQUIRED:
      ProcessRequiredField(descriptor);
      ProcessSingleField(descriptor);
      break;
  }
  FieldOptions const& options = descriptor->options();
  if (options.HasExtension(serialization::size)) {
    size_field_name_[descriptor] = options.GetExtension(serialization::size);
  }
}

void Generator::ProcessInOut(
    Descriptor const* descriptor,
    std::vector<FieldDescriptor const*>* field_descriptors) {
  std::string const& name = descriptor->name();
  cpp_fill_body_[descriptor] = "  auto* to = message->mutable_" +
                               ToLower(name) + "();\n";
  cpp_nested_type_[descriptor] = "  struct " + name + " {\n";
  for (int i = 0; i < descriptor->field_count(); ++i) {
    FieldDescriptor const* field_descriptor = descriptor->field(i);
    if (field_descriptors != nullptr) {
      field_descriptors->push_back(field_descriptor);
    }
    ProcessField(field_descriptor);
    cpp_fill_body_[descriptor] += "  " +
                                  cpp_field_copy_[field_descriptor] + "\n";
    cpp_nested_type_[descriptor] += "    " +
                                    cpp_field_type_[field_descriptor] +
                                    " const " +
                                    field_descriptor->name() + ";\n";

    // This this field has a size, generate it now.
    if (size_field_name_.find(field_descriptor) != size_field_name_.end()) {
      cpp_nested_type_[descriptor] +=
          "    int const " + size_field_name_[field_descriptor] + ";\n";
    }
  }
  cpp_nested_type_[descriptor] += "  };\n";
  LOG(ERROR)<<cpp_fill_body_[descriptor];
}

void Generator::ProcessReturn(Descriptor const* descriptor) {
  CHECK_EQ(1, descriptor->field_count())
      << descriptor->full_name() << " must have exactly one field";
  FieldDescriptor const* field_descriptor = descriptor->field(0);
  ProcessField(field_descriptor);
  cpp_nested_type_[descriptor] =
      "  using Return = " + cpp_field_type_[field_descriptor] + ";\n";
}

void Generator::ProcessMethodExtension(Descriptor const* descriptor) {
  bool has_in = false;
  bool has_out = false;
  bool has_return = false;

  // Do a first pass to determine which fields are in-out.  The data produced
  // here will be overwritten by the next pass.
  std::vector<FieldDescriptor const*> field_descriptors;
  for (int i = 0; i < descriptor->nested_type_count(); ++i) {
    Descriptor const* nested_descriptor = descriptor->nested_type(i);
    const std::string& nested_name = nested_descriptor->name();
    if (nested_name == kIn) {
      has_in = true;
      ProcessInOut(nested_descriptor, &field_descriptors);
    } else if (nested_name == kOut) {
      has_out = true;
      ProcessInOut(nested_descriptor, &field_descriptors);
    } else if (nested_name == kReturn) {
      has_return = true;
    } else {
      LOG(FATAL) << "Unexpected nested message "
                 << nested_descriptor->full_name();
    }
  }

  // Now mark the fields that have the same name in In and Out as in-out.
  if (has_in && has_out) {
    std::sort(field_descriptors.begin(),
              field_descriptors.end(),
              [](FieldDescriptor const* left,
                 FieldDescriptor const* right) {
                return left->name() < right->name();
              });
    for (int i = 0; i < field_descriptors.size() - 1; ++i) {
      if (field_descriptors[i]->name() == field_descriptors[i + 1]->name()) {
        in_out_field_.insert(field_descriptors[i]);
        in_out_field_.insert(field_descriptors[i + 1]);
      }
    }
  }

  // The second pass that produces the actual output.
  cpp_method_type_[descriptor] = "struct " + descriptor->name() + " {\n";
  for (int i = 0; i < descriptor->nested_type_count(); ++i) {
    Descriptor const* nested_descriptor = descriptor->nested_type(i);
    const std::string& nested_name = nested_descriptor->name();
    if (nested_name == kIn) {
      ProcessInOut(nested_descriptor, /*field_descriptors=*/nullptr);
    } else if (nested_name == kOut) {
      ProcessInOut(nested_descriptor, /*field_descriptors=*/nullptr);
    } else if (nested_name == kReturn) {
      ProcessReturn(nested_descriptor);
    }
    cpp_method_type_[descriptor] += cpp_nested_type_[nested_descriptor];
  }
  if (has_in || has_out || has_return) {
    cpp_method_type_[descriptor] += "\n";
  }
  cpp_method_type_[descriptor] +=
      "  using Message = serialization::" +
      descriptor->name() + ";\n";
  if (has_in) {
    cpp_method_type_[descriptor] += "  static void Fill(In const& in, "
                                    "not_null<Message*> const message);\n";
  }
  if (has_out) {
    cpp_method_type_[descriptor] += "  static void Fill(Out const& out, "
                                    "not_null<Message*> const message);\n";
  }
  if (has_return) {
    cpp_method_type_[descriptor] += "  static void Fill("
                                    "Return const& result, "
                                    "not_null<Message*> const message);\n";
  }
  cpp_method_type_[descriptor] += "  static void Run("
                                  "Message const& message,\n"
                                  "                  not_null<"
                                  "Player::PointerMap*> const pointer_map);"
                                  "\n";
  cpp_method_type_[descriptor] += "};\n\n";
}

}  // namespace

void GenerateProfiles() {
  Generator generator;
  generator.ProcessMethods();

  // Now write the output.
  std::experimental::filesystem::path const directory =
      SOLUTION_DIR / "journal";
  std::ofstream profiles_hpp(directory / "profiles.hpp");
  CHECK(profiles_hpp.good());

  profiles_hpp << "#pragma once\n";
  profiles_hpp << "\n";
  profiles_hpp << "#include \"base/not_null.hpp\"\n";
  profiles_hpp << "#include \"journal/player.hpp\"\n";
  profiles_hpp << "#include \"ksp_plugin/interface.hpp\"\n";
  profiles_hpp << "#include \"serialization/journal.pb.h\"\n";
  profiles_hpp << "\n";
  profiles_hpp << "namespace principia {\n";
  profiles_hpp << "\n";
  profiles_hpp << "using base::not_null;\n";
  profiles_hpp << "using ksp_plugin::KSPPart;\n";
  profiles_hpp << "using ksp_plugin::LineAndIterator;\n";
  profiles_hpp << "using ksp_plugin::NavigationFrame;\n";
  profiles_hpp << "using ksp_plugin::Plugin;\n";
  profiles_hpp << "using ksp_plugin::QP;\n";
  profiles_hpp << "using ksp_plugin::WXYZ;\n";
  profiles_hpp << "using ksp_plugin::XYZ;\n";
  profiles_hpp << "using ksp_plugin::XYZSegment;\n";
  profiles_hpp << "\n";
  profiles_hpp << "namespace journal {\n";
  profiles_hpp << "\n";

  for (auto const& cpp_method_type : generator.GetCppMethodTypes()) {
    profiles_hpp << cpp_method_type;
  }

  profiles_hpp << "}  // namespace journal\n";
  profiles_hpp << "}  // namespace principia\n";
}

}  // namespace tools
}  // namespace principia
