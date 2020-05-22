
#include "tools/journal_proto_processor.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "base/map_util.hpp"
#include "glog/logging.h"
#include "serialization/journal.pb.h"

namespace principia {
namespace tools {
namespace internal_journal_proto_processor {

using base::Contains;
using ::google::protobuf::MessageOptions;

namespace {

constexpr char method_message_name[] = "Method";
constexpr char in_message_name[] = "In";
constexpr char return_message_name[] = "Return";
constexpr char out_message_name[] = "Out";

std::string Join(std::vector<std::string> const& v, std::string const& joiner) {
  std::string joined;
  bool is_first = true;
  for (auto const& vi : v) {
    if (vi.empty()) {
      continue;
    }
    if (is_first) {
      is_first = false;
      joined = vi;
    } else {
      joined += joiner + vi;
    }
  }
  return joined;
}

std::string ToLower(std::string const& s) {
  std::string lower;
  for (int i = 0; i < s.size(); ++i) {
    if (i > 0 && i < s.size() - 1 &&
        static_cast<bool>(std::isupper(s[i])) &&
        static_cast<bool>(std::islower(s[i + 1]))) {
      lower += "_" + std::string(1, std::tolower(s[i]));
    } else {
      lower += static_cast<char>(std::tolower(s[i]));
    }
  }
  return lower;
}

}  // namespace

void JournalProtoProcessor::ProcessMessages() {
  // Get the file containing |Method|.
  Descriptor const* method_descriptor =
      journal::serialization::Method::descriptor();
  FileDescriptor const* file_descriptor = method_descriptor->file();

  // Process all the messages in that file.
  for (int i = 0; i < file_descriptor->message_type_count(); ++i) {
    Descriptor const* message_descriptor = file_descriptor->message_type(i);
    std::string message_descriptor_name = message_descriptor->name();
    if (message_descriptor->extension_range_count() > 0) {
      // Only the |Method| message should have a range.  Don't generate any code
      // for it.
      CHECK_EQ(method_message_name, message_descriptor_name)
          << message_descriptor_name << " should not have extension ranges";
      continue;
    }
    switch (message_descriptor->extension_count()) {
      case 0: {
        // A message corresponding to a struct interchanged through the
        // interface.
        ProcessInterchangeMessage(message_descriptor);
        break;
      }
      case 1: {
        // An extension.  Check that it extends |Method|.
        FieldDescriptor const* extension = message_descriptor->extension(0);
        CHECK(extension->is_extension());
        Descriptor const* containing_type = extension->containing_type();
        CHECK_EQ(method_descriptor, containing_type)
            << message_descriptor_name << " extends a message other than "
            << method_descriptor->name() << ": " << containing_type->name();
        ProcessMethodExtension(message_descriptor);
        break;
      }
      default: {
        LOG(FATAL) << message_descriptor_name << " has "
                   << message_descriptor->extension_count() << " extensions";
      }
    }
  }
}

std::vector<std::string>
JournalProtoProcessor::GetCsInterfaceMethodDeclarations() const {
  std::vector<std::string> result;
  for (auto const& pair : cs_interface_method_declaration_) {
    result.push_back(pair.second);
  }
  return result;
}

std::vector<std::string>
JournalProtoProcessor::GetCsInterchangeTypeDeclarations() const {
  std::vector<std::string> result;
  for (auto const& pair : cs_interchange_type_declaration_) {
    result.push_back(pair.second);
  }
  return result;
}

std::vector<std::string> JournalProtoProcessor::GetCsCustomMarshalerClasses()
    const {
  std::vector<std::string> result;
  for (auto const& pair : cs_custom_marshaler_class_) {
    result.push_back(pair.second);
  }
  return result;
}

std::vector<std::string>
JournalProtoProcessor::GetCxxInterfaceMethodDeclarations() const {
  std::vector<std::string> result;
  for (auto const& pair : cxx_interface_method_declaration_) {
    result.push_back(pair.second);
  }
  return result;
}

std::vector<std::string>
JournalProtoProcessor::GetCxxInterchangeTypeDeclarations() const {
  std::vector<std::string> result;
  for (auto const& pair : cxx_interchange_type_declaration_) {
    result.push_back(pair.second);
  }
  return result;
}

std::vector<std::string>
JournalProtoProcessor::GetCxxInterchangeImplementations() const {
  std::vector<std::string> result;
  result.push_back("namespace {\n\n");
  for (auto const& pair : cxx_deserialize_definition_) {
    result.push_back(pair.second);
  }
  for (auto const& pair : cxx_serialize_definition_) {
    result.push_back(pair.second);
  }
  result.push_back("}  // namespace\n\n");
  return result;
}

std::vector<std::string>
JournalProtoProcessor::GetCxxMethodImplementations() const {
  std::vector<std::string> result;
  for (auto const& pair : cxx_functions_implementation_) {
    result.push_back(pair.second);
  }
  return result;
}

std::vector<std::string> JournalProtoProcessor::GetCxxMethodTypes() const {
  std::vector<std::string> result;
  for (auto const& pair : cxx_toplevel_type_declaration_) {
    result.push_back(pair.second);
  }
  return result;
}

std::vector<std::string> JournalProtoProcessor::GetCxxPlayStatements() const {
  std::vector<std::string> result;
  result.push_back("{\n  bool ran = false;\n");
  for (auto const& pair : cxx_play_statement_) {
    result.push_back(pair.second);
  }
  result.push_back("  CHECK(ran) << method_in->DebugString() << \"\\n\"\n"
                   "             << method_out_return->DebugString();\n}\n");
  return result;
}

void JournalProtoProcessor::ProcessRepeatedMessageField(
    FieldDescriptor const* descriptor) {
  Descriptor const* message_type = descriptor->message_type();
  std::string const& message_type_name = message_type->name();

  FieldOptions const& options = descriptor->options();
  field_cs_type_[descriptor] = message_type_name + "[]";
  // TODO(phl): Use Contains instead of empty for all the set/map predicates.
  if (cs_custom_marshaler_name_[message_type].empty()) {
    // This wouldn't be hard, we'd need another RepeatedMarshaller that copies
    // structs, but we don't need it yet.
    LOG(FATAL) << "Repeated messages with an element that does not have a "
                  "custom marshaller are not yet implemented.";
  } else {
    field_cs_custom_marshaler_[descriptor] =
        "RepeatedMarshaler<" + message_type_name + ", " +
        cs_custom_marshaler_name_[message_type] + ">";
  }
  field_cxx_type_[descriptor] = message_type_name + " const* const*";

  field_cxx_arguments_fn_[descriptor] =
      [](std::string const& identifier) -> std::vector<std::string> {
        return {"&" + identifier + "[0]"};
      };
  field_cxx_assignment_fn_[descriptor] =
      [this, descriptor, message_type_name](
          std::string const& prefix, std::string const& expr) {
        std::string const& descriptor_name = descriptor->name();
        return "  for (" + field_cxx_type_[descriptor] + " " + descriptor_name +
               " = " + expr + "; " +
               descriptor_name + " != nullptr && "
               "*" + descriptor_name + " != nullptr; "
               "++" + descriptor_name + ") {\n"
               "    *" + prefix + "add_" + descriptor_name + "() = " +
               field_cxx_serializer_fn_[descriptor]("**"+ descriptor_name) +
               ";\n"
               "  }\n";
      };
  std::string const storage_name = descriptor->name() + "_storage";
  field_cxx_deserialization_storage_name_[descriptor] = storage_name;
  field_cxx_deserializer_fn_[descriptor] =
      [message_type_name, storage_name](
          std::string const& expr) {
        // Yes, this lambda generates a lambda.
        return "[&pointer_map, &" + storage_name +
               "](::google::protobuf::RepeatedPtrField<serialization::" +
               message_type_name + "> const& messages) {\n"
               "            for (auto const& message : messages) {\n" +
               "              " + storage_name +
               ".first.push_back(Deserialize" + message_type_name +
               "(message, pointer_map));\n" +
               "              " + storage_name +
               ".second.push_back(&" + storage_name + ".first.back());\n" +
               "            }\n"
               "            return &" + storage_name + ".second[0];\n" +
               "          }(" + expr + ")";
      };
  field_cxx_deserialization_storage_type_[descriptor] =
      "std::pair<"
      "std::vector<" + message_type_name + ">, "
      "std::vector<" + message_type_name + " const*>>";
  field_cxx_serializer_fn_[descriptor] =
      [message_type_name](std::string const& expr) {
        return "Serialize" + message_type_name + "(" + expr + ")";
      };
}

void JournalProtoProcessor::ProcessOptionalNonStringField(
    FieldDescriptor const* descriptor,
    std::string const& cs_boxed_type,
    std::string const& cs_unboxed_type,
    std::string const& cxx_type) {
  FieldOptions const& options = descriptor->options();

  // Build a lambda to construct a marshaler name.
  std::function<std::string(std::string const&)> custom_marshaler_generic_name;
  if (options.HasExtension(journal::serialization::is_produced)) {
    CHECK(options.GetExtension(journal::serialization::is_produced))
        << descriptor->full_name() << " has incorrect (is_produced) option";
    custom_marshaler_generic_name =
        [](std::string const& type) {
          return "OwnershipTransferMarshaler<" + type +
                 ", OptionalMarshaler<" + type + ">>";
        };
  } else {
    custom_marshaler_generic_name =
        [](std::string const& type) {
          return "OptionalMarshaler<" + type + ">";
        };
  }

  // It is not possible to use a custom marshaler on an |T?|, as this raises
  // |System.Runtime.InteropServices.MarshalDirectiveException| with the message
  // "Custom marshalers are only allowed on classes, strings, arrays, and boxed
  // value types.".
  if (Contains(interchange_, descriptor)) {
    // This may be null as we may be called on a scalar field.
    Descriptor const* message_type = descriptor->message_type();
    if (cs_custom_marshaler_name_[message_type].empty()) {
      field_cs_custom_marshaler_[descriptor] =
          custom_marshaler_generic_name(cs_unboxed_type);
    } else {
      // This wouldn't be hard, we'd need another OptionalMarshaller that calls
      // the element's marshaler, but we don't need it yet.
      LOG(FATAL) << "Optional messages with an element that does have a custom "
                    "marshaller are not yet implemented.";
    }
    // For fields of interchange messages we can use a |T?| as the field is not
    // the part being marshaled, it is the entire interchange message.
    field_cs_type_[descriptor] = cs_unboxed_type + "?";
  } else {
    // We could use a boxed |T|, whose type would be |object|, but we would lose
    // static typing.  We use a custom strongly-typed boxed type instead.
    field_cs_custom_marshaler_[descriptor] =
        custom_marshaler_generic_name(cs_unboxed_type);
    field_cs_type_[descriptor] = cs_boxed_type;
  }
  // Unfortunately this pointer type cannot be defaulted to nullptr as it would
  // make the enclosing type non-POD.
  field_cxx_type_[descriptor] = cxx_type + " const*";

  field_cxx_arguments_fn_[descriptor] =
      [](std::string const& identifier) -> std::vector<std::string> {
        return {identifier};
      };
  field_cxx_indirect_member_get_fn_[descriptor] =
      [](std::string const& expr) {
        return "*" + expr;
      };
}

void JournalProtoProcessor::ProcessOptionalScalarField(
    FieldDescriptor const* descriptor,
    std::string const& cxx_type) {
  // What we would want is a reference to a scalar field of a proto, but this
  // is not exposed by the proto API.  Instead, we copy the field in a variable
  // and use the address of that variable.
  std::string const storage_name = descriptor->name() + "_storage";
  field_cxx_deserialization_storage_name_[descriptor] = storage_name;
  field_cxx_deserializer_fn_[descriptor] =
      [cxx_type, storage_name](
          std::string const& expr) {
        // Yes, this lambda generates a lambda.
        return "[&" + storage_name +
               "](" + cxx_type + " value) {\n"
               "            " + storage_name + " = value;\n" +
               "            return &" + storage_name + ";\n"
               "          }(" + expr + ")";
      };
  field_cxx_deserialization_storage_type_[descriptor] = cxx_type;
  field_cxx_optional_pointer_fn_[descriptor] =
      [storage_name](std::string const& condition, std::string const& expr) {
        return condition + " ? " + expr + " : nullptr";
      };
}

void JournalProtoProcessor::ProcessOptionalDoubleField(
    FieldDescriptor const* descriptor) {
  ProcessOptionalNonStringField(
      descriptor,
      /*cs_boxed_type=*/"BoxedDouble",
      /*cs_unboxed_type=*/"double",
      /*cxx_type=*/"double");
  ProcessOptionalScalarField(descriptor, "double");
}

void JournalProtoProcessor::ProcessOptionalInt32Field(
    FieldDescriptor const* descriptor) {
  ProcessOptionalNonStringField(
      descriptor,
      /*cs_boxed_type=*/"BoxedInt32",
      /*cs_unboxed_type=*/"int",
      /*cxx_type=*/"int");
  ProcessOptionalScalarField(descriptor, "int");
}

void JournalProtoProcessor::ProcessOptionalMessageField(
    FieldDescriptor const* descriptor) {
  Descriptor const* message_type = descriptor->message_type();
  std::string const& message_type_name = message_type->name();
  ProcessOptionalNonStringField(
      descriptor,
      /*cs_boxed_type=*/"Boxed" + message_type_name,
      /*cs_unboxed_type=*/message_type_name,
      /*cxx_type=*/message_type_name);

  std::string const storage_name = descriptor->name() + "_storage";
  field_cxx_deserialization_storage_name_[descriptor] = storage_name;
  field_cxx_deserializer_fn_[descriptor] =
      [message_type_name, storage_name](
          std::string const& expr) {
        // Yes, this lambda generates a lambda.
        return "[&pointer_map, &" + storage_name +
               "](serialization::" + message_type_name + " const& message) {\n"
               "            " + storage_name + " = Deserialize" +
               message_type_name + "(message, pointer_map);\n" +
               "            return &" + storage_name + ";\n"
               "          }(" + expr + ")";
      };
  field_cxx_deserialization_storage_type_[descriptor] = message_type_name;
  field_cxx_optional_pointer_fn_[descriptor] =
      [storage_name](std::string const& condition, std::string const& expr) {
        return condition + " ? " + expr + " : nullptr";
      };

  ProcessSingleMessageField(descriptor);
}

void JournalProtoProcessor::ProcessRequiredFixed32Field(
    FieldDescriptor const* descriptor) {
  field_cs_type_[descriptor] = "uint";
  field_cxx_type_[descriptor] = "uint32_t";
}

void JournalProtoProcessor::ProcessRequiredFixed64Field(
    FieldDescriptor const* descriptor) {
  FieldOptions const& options = descriptor->options();
  CHECK(!options.HasExtension(journal::serialization::pointer_to) ||
        !options.HasExtension(journal::serialization::encoding))
      << descriptor->full_name()
      << " cannot have both a (pointer_to) and an (encoding) option";
  std::string pointer_to;
  if (options.HasExtension(journal::serialization::pointer_to)) {
    pointer_to = options.GetExtension(journal::serialization::pointer_to);
  }
  if (options.HasExtension(journal::serialization::encoding)) {
    switch (options.GetExtension(journal::serialization::encoding)) {
      case journal::serialization::UTF_8:
        pointer_to = "char const";
        break;
      case journal::serialization::UTF_16:
        pointer_to = "char16_t const";
        break;
    }
  }

  if (options.HasExtension(journal::serialization::disposable)) {
    CHECK(!options.HasExtension(journal::serialization::is_consumed) ||
          !options.HasExtension(journal::serialization::is_consumed_if))
      << descriptor->full_name() << " must not be consumed to be disposable";
    field_cs_type_[descriptor] =
        options.GetExtension(journal::serialization::disposable);
    field_cs_custom_marshaler_[descriptor] =
        options.GetExtension(journal::serialization::disposable) +
        "Marshaler";
  } else {
    field_cs_type_[descriptor] = "IntPtr";
  }
  if (options.HasExtension(journal::serialization::is_subject)) {
    CHECK(options.GetExtension(journal::serialization::is_subject))
        << descriptor->full_name() << " has incorrect (is_subject) option";
    field_cs_type_[descriptor] = "this " + field_cs_type_[descriptor];
  }
  field_cxx_type_[descriptor] = pointer_to + "*";

  if (Contains(out_, descriptor) && !Contains(in_out_, descriptor)) {
    CHECK(!options.HasExtension(journal::serialization::is_consumed) &&
          !options.HasExtension(journal::serialization::is_consumed_if))
        << "out parameter " << descriptor->full_name() << " cannot be consumed";
  }

  if (options.HasExtension(journal::serialization::is_consumed)) {
    CHECK(options.GetExtension(journal::serialization::is_consumed))
        << descriptor->full_name() << " has incorrect (is_consumed) option";
    field_cxx_deleter_fn_[descriptor] =
        [](std::string const& expr) {
          return "  Delete(" + expr + ", pointer_map);\n";
        };
  }
  if (options.HasExtension(journal::serialization::is_consumed_if)) {
    CHECK(!options.HasExtension(journal::serialization::is_consumed))
        << descriptor->full_name()
        << " has incorrect (is_consumed) and (is_consumed_if) options";
    field_cxx_deleter_fn_[descriptor] =
        [options](std::string const& expr) {
          return "  if (" +
                 options.GetExtension(journal::serialization::is_consumed_if) +
                 ") {\n    Delete(" + expr + ", pointer_map);\n  }\n";
        };
  }
  if (options.HasExtension(journal::serialization::is_produced)) {
    CHECK(options.GetExtension(journal::serialization::is_produced))
        << descriptor->full_name() << " has incorrect (is_produced) option";
    field_cxx_inserter_fn_[descriptor] =
        [](std::string const& expr1, std::string const& expr2) {
          return "  Insert(" + expr1 + ", " + expr2 + ", pointer_map);\n";
        };
  }
  if (options.HasExtension(journal::serialization::is_produced_if)) {
    CHECK(!options.HasExtension(journal::serialization::is_produced))
        << descriptor->full_name()
        << " has incorrect (is_produced) and (is_produced_if) options";
    field_cxx_inserter_fn_[descriptor] =
        [options](std::string const& expr1, std::string const& expr2) {
          return "  if (" +
                 options.GetExtension(journal::serialization::is_produced_if) +
                 ") {\n    Insert(" + expr1 + ", " + expr2 +
                 ", pointer_map);\n  }\n";
        };
  }

  // Special handlings for produced C-style strings: these are seen from the C#
  // as strings, and marshaled with immediate destruction.
  if (options.HasExtension(journal::serialization::encoding) &&
      (options.HasExtension(journal::serialization::is_produced) ||
       options.HasExtension(journal::serialization::is_produced_if))) {
    field_cs_type_[descriptor] = "string";
    switch (options.GetExtension(journal::serialization::encoding)) {
      case journal::serialization::UTF_8:
        field_cs_custom_marshaler_[descriptor] =
            "OwnershipTransferUTF8Marshaler";
        break;
      case journal::serialization::UTF_16:
        field_cs_custom_marshaler_[descriptor] =
            "OwnershipTransferUTF16Marshaler";
        break;
    }
  }

  field_cxx_deserializer_fn_[descriptor] =
      [pointer_to](std::string const& expr) {
        return "DeserializePointer<" + pointer_to + "*>(" + expr +
               ", pointer_map)";
      };
  field_cxx_serializer_fn_[descriptor] =
      [](std::string const& expr) {
        return "SerializePointer(" + expr + ")";
      };
}

void JournalProtoProcessor::ProcessRequiredMessageField(
    FieldDescriptor const* descriptor) {
  FieldOptions const& options = descriptor->options();
  Descriptor const* message_type = descriptor->message_type();
  std::string const& message_type_name = message_type->name();
  field_cs_type_[descriptor] = message_type_name;
  field_cxx_type_[descriptor] = message_type_name;

  if (!cs_custom_marshaler_name_[message_type].empty()) {
    if (Contains(in_, descriptor)) {
      field_cs_custom_marshaler_[descriptor] =
          cs_custom_marshaler_name_[message_type];
      field_cxx_mode_fn_[descriptor] =
          [](std::string const& type) {
            return type + " const&";
          };
      // No need to define field_cxx_indirect_member_get_fn_ here because
      // references don't need a level of indirection.
    }
    if (Contains(return_, descriptor)) {
      if (options.HasExtension(journal::serialization::is_produced)) {
        CHECK(options.GetExtension(journal::serialization::is_produced))
            << descriptor->full_name() << " has incorrect (is_produced) option";
        field_cs_custom_marshaler_[descriptor] =
            "OwnershipTransferMarshaler<" + field_cs_type_[descriptor] + ", " +
            cs_custom_marshaler_name_[message_type] + ">";
      } else {
        field_cs_custom_marshaler_[descriptor] =
            cs_custom_marshaler_name_[message_type];
      }
      field_cxx_mode_fn_[descriptor] =
          [](std::string const& type) {
            return type + "*";
          };
      field_cxx_indirect_member_get_fn_[descriptor] =
          [](std::string const& expr) {
            return "*" + expr;
          };
    }
  }
  std::string const deserialization_storage_arguments =
      cxx_deserialization_storage_arguments_[message_type];
  field_cxx_deserializer_fn_[descriptor] =
      [message_type_name, deserialization_storage_arguments](
          std::string const& expr) {
        return "Deserialize" + message_type_name + "(" + expr +
               ", pointer_map" + deserialization_storage_arguments + ")";
      };

  ProcessSingleMessageField(descriptor);
}

void JournalProtoProcessor::ProcessRequiredBoolField(
    FieldDescriptor const* descriptor) {
  field_cs_type_[descriptor] = "bool";
  field_cs_predefined_marshaler_[descriptor] = "UnmanagedType.I1";
  field_cs_private_type_[descriptor] = "byte";
  field_cs_private_getter_fn_[descriptor] =
      [](std::vector<std::string> const& identifiers) {
        CHECK_EQ(1, identifiers.size());
        return "get { return " + identifiers[0] + " != (byte)0; }";
      };
  field_cs_private_setter_fn_[descriptor] =
      [](std::vector<std::string> const& identifiers) {
        CHECK_EQ(1, identifiers.size());
        return "set { " + identifiers[0] + " = value ? (byte)1 : (byte)0; }";
      };

  field_cxx_type_[descriptor] = descriptor->cpp_type_name();
}

void JournalProtoProcessor::ProcessRequiredBytesField(
  FieldDescriptor const* descriptor) {
  FieldOptions const& options = descriptor->options();
  LOG_IF(FATAL,
    options.HasExtension(journal::serialization::is_produced) ||
    options.HasExtension(journal::serialization::is_produced_if))
    << descriptor->full_name()
    << " is a bytes field and cannot be produced. Use a fixed64 field that "
    << "has the (encoding) option instead.";
  LOG_IF(FATAL,
         !options.HasExtension(journal::serialization::encoding) ||
         options.GetExtension(journal::serialization::encoding) !=
         journal::serialization::UTF_16)
      << descriptor->full_name()
      << " is a bytes field and must have the (encoding) = UTF_16 option.";

  field_cs_predefined_marshaler_[descriptor] = "UnmanagedType.LPWStr";
  field_cs_type_[descriptor] = "string";
  field_cxx_type_[descriptor] = "char16_t const*";
  field_cxx_arguments_fn_[descriptor] =
      [](std::string const& identifier) -> std::vector<std::string> {
        return {identifier + ".c_str()"};
      };
  field_cxx_deserializer_fn_[descriptor] =
      [](std::string const& expr) {
        return "DeserializeUtf16(" + expr + ")";
      };
  field_cxx_serializer_fn_[descriptor] =
      [](std::string const& expr) {
        return "SerializeUtf16(" + expr + ")";
      };
}

void JournalProtoProcessor::ProcessRequiredDoubleField(
    FieldDescriptor const* descriptor) {
  field_cs_type_[descriptor] = "double";
  field_cxx_type_[descriptor] = descriptor->cpp_type_name();
}

void JournalProtoProcessor::ProcessRequiredInt32Field(
    FieldDescriptor const* descriptor) {
  field_cs_type_[descriptor] = "int";
  field_cxx_type_[descriptor] = "int";
}

void JournalProtoProcessor::ProcessRequiredInt64Field(
    FieldDescriptor const* descriptor) {
  field_cs_type_[descriptor] = "long";
  field_cxx_type_[descriptor] = "std::int64_t";
}

void JournalProtoProcessor::ProcessRequiredUint32Field(
    FieldDescriptor const* descriptor) {
  field_cs_type_[descriptor] = "uint";
  field_cxx_type_[descriptor] = "uint32_t";
}

void JournalProtoProcessor::ProcessSingleMessageField(
    FieldDescriptor const* descriptor) {
  Descriptor const* message_type = descriptor->message_type();
  std::string const& message_type_name = message_type->name();

  field_cxx_assignment_fn_[descriptor] =
      [this, descriptor](std::string const& prefix,
                         std::string const& expr) {
        return "  *" + prefix + "mutable_" + descriptor->name() +
               "() = " + field_cxx_serializer_fn_[descriptor](expr) + ";\n";
      };
  field_cxx_serializer_fn_[descriptor] =
      [message_type_name](std::string const& expr) {
        return "Serialize" + message_type_name + "(" + expr + ")";
      };
}

void JournalProtoProcessor::ProcessSingleStringField(
    FieldDescriptor const* descriptor) {
  FieldOptions const& options = descriptor->options();
  LOG_IF(FATAL,
         options.HasExtension(journal::serialization::is_produced) ||
             options.HasExtension(journal::serialization::is_produced_if))
      << descriptor->full_name()
      << " is a string field and cannot be produced. Use a fixed64 field that "
      << "has the (encoding) option instead.";

  field_cs_custom_marshaler_[descriptor] = "NoOwnershipTransferUTF8Marshaler";
  field_cs_type_[descriptor] = "string";
  field_cxx_type_[descriptor] = "char const*";
  field_cxx_deserializer_fn_[descriptor] =
      [](std::string const& expr) {
        return expr + ".c_str()";
      };
}

// Listing all the case values below is not helpful.
#pragma warning(push)
#pragma warning(disable : 4061)
void JournalProtoProcessor::ProcessOptionalField(
    FieldDescriptor const* descriptor) {
  field_cxx_optional_assignment_fn_[descriptor] =
      [](std::string const& expr, std::string const& stmt) {
        return "  if (" + expr + " != nullptr) {\n  " + stmt + "  }\n";
      };
  field_cxx_optional_pointer_fn_[descriptor] =
      [](std::string const& condition, std::string const& expr) {
        return condition + " ? " + expr + " : nullptr";
      };
  switch (descriptor->type()) {
    case FieldDescriptor::TYPE_DOUBLE:
      ProcessOptionalDoubleField(descriptor);
      break;
    case FieldDescriptor::TYPE_INT32:
      ProcessOptionalInt32Field(descriptor);
      break;
    case FieldDescriptor::TYPE_MESSAGE:
      ProcessOptionalMessageField(descriptor);
      break;
    case FieldDescriptor::TYPE_STRING:
      ProcessSingleStringField(descriptor);
      break;
    default:
      LOG(FATAL) << descriptor->full_name() << " has unexpected type "
                 << descriptor->type_name();
  }
}

void JournalProtoProcessor::ProcessRepeatedField(
    FieldDescriptor const* descriptor) {
  switch (descriptor->type()) {
  case FieldDescriptor::TYPE_MESSAGE:
      ProcessRepeatedMessageField(descriptor);
      break;
    default:
      LOG(FATAL) << descriptor->full_name() << " has unexpected type "
                 << descriptor->type_name();
  }
}

void JournalProtoProcessor::ProcessRequiredField(
    FieldDescriptor const* descriptor) {
  switch (descriptor->type()) {
    case FieldDescriptor::TYPE_BOOL:
      ProcessRequiredBoolField(descriptor);
      break;
    case FieldDescriptor::TYPE_BYTES:
      ProcessRequiredBytesField(descriptor);
      break;
    case FieldDescriptor::TYPE_DOUBLE:
      ProcessRequiredDoubleField(descriptor);
      break;
    case FieldDescriptor::TYPE_FIXED32:
      ProcessRequiredFixed32Field(descriptor);
      break;
    case FieldDescriptor::TYPE_FIXED64:
      ProcessRequiredFixed64Field(descriptor);
      break;
    case FieldDescriptor::TYPE_INT32:
      ProcessRequiredInt32Field(descriptor);
      break;
    case FieldDescriptor::TYPE_INT64:
      ProcessRequiredInt64Field(descriptor);
      break;
    case FieldDescriptor::TYPE_MESSAGE:
      ProcessRequiredMessageField(descriptor);
      break;
    case FieldDescriptor::TYPE_STRING:
      ProcessSingleStringField(descriptor);
      break;
    case FieldDescriptor::TYPE_UINT32:
      ProcessRequiredUint32Field(descriptor);
      break;
    default:
      LOG(FATAL) << descriptor->full_name() << " has unexpected type "
                 << descriptor->type_name();
  }

  // For in-out fields the data is actually passed with an extra level of
  // indirection.
  if (Contains(in_out_, descriptor) || Contains(out_, descriptor)) {
    field_cxx_arguments_fn_[descriptor] =
        [](std::string const& identifier) -> std::vector<std::string> {
          return {"&" + identifier};
        };
    field_cxx_indirect_member_get_fn_[descriptor] =
        [](std::string const& expr) {
          return "*" + expr;
        };

    if (Contains(in_out_, descriptor)) {
      field_cs_mode_fn_[descriptor] =
          [](std::string const& type) {
            return "ref " + type;
          };
    } else {
      field_cs_mode_fn_[descriptor] =
          [](std::string const& type) {
            return "out " + type;
          };
    }
    field_cxx_mode_fn_[descriptor] =
        [](std::string const& type) {
          return type + "* const";
        };
  }
}
#pragma warning(pop)

void JournalProtoProcessor::ProcessField(FieldDescriptor const* descriptor) {
  // Useful defaults for the lambdas, which ensure that they are set for all
  // fields.  They will be overwritten by actual processing as needed.
  field_cs_mode_fn_[descriptor] =
      [](std::string const& type) {
        return type;
      };
  field_cxx_arguments_fn_[descriptor] =
      [](std::string const& identifier) -> std::vector<std::string> {
        return {identifier};
      };
  field_cxx_assignment_fn_[descriptor] =
      [this, descriptor](std::string const& prefix, std::string const& expr) {
        return "  " + prefix + "set_" + descriptor->name() + "(" +
               field_cxx_serializer_fn_[descriptor](expr) + ");\n";
      };
  field_cxx_indirect_member_get_fn_[descriptor] =
      [](std::string const& expr) {
        return expr;
      };
  field_cxx_deserializer_fn_[descriptor] =
      [](std::string const& expr) {
        return expr;
      };
  if (Contains(return_, descriptor)) {
    // No const on return types.
    field_cxx_mode_fn_[descriptor] =
        [](std::string const& type) {
          return type;
        };
  } else {
    field_cxx_mode_fn_[descriptor] =
        [](std::string const& type) {
          return type + " const";
        };
  }
  field_cxx_optional_assignment_fn_[descriptor] =
      [](std::string const& expr, std::string const& stmt) {
        return stmt;
      };
  field_cxx_optional_pointer_fn_[descriptor] =
      [](std::string const& condition, std::string const& expr) {
        return expr;
      };
  field_cxx_serializer_fn_[descriptor] =
      [](std::string const& expr) {
        return expr;
      };

  switch (descriptor->label()) {
    case FieldDescriptor::LABEL_OPTIONAL:
      ProcessOptionalField(descriptor);
      break;
    case FieldDescriptor::LABEL_REPEATED:
      ProcessRepeatedField(descriptor);
      break;
    case FieldDescriptor::LABEL_REQUIRED:
      ProcessRequiredField(descriptor);
      break;
    }
}

void JournalProtoProcessor::ProcessInOut(
  Descriptor const* descriptor,
  std::vector<FieldDescriptor const*>* field_descriptors) {
  std::string const& name = descriptor->name();

  std::string cxx_message_prefix;
  {
    std::string const cxx_message_name =
        "message->mutable_" + ToLower(name) + "()";
    // Generate slightly more compact code in the frequent case where the
    // message only has one field.
    if (descriptor->field_count() > 1) {
      cxx_fill_body_[descriptor] =
          "  auto* const m = " + cxx_message_name + ";\n";
      cxx_message_prefix = "m->";
    } else {
      cxx_message_prefix = cxx_message_name + "->";
      cxx_fill_body_[descriptor].clear();
    }
  }

  cs_interface_parameters_[descriptor].clear();
  cxx_interface_parameters_[descriptor].clear();
  cxx_run_body_prolog_[descriptor] =
      "  [[maybe_unused]] auto const& " + ToLower(name) + " = message." +
      ToLower(name) + "();\n";
  cxx_run_arguments_[descriptor].clear();
  cxx_run_body_epilog_[descriptor].clear();

  cxx_nested_type_declaration_[descriptor] = "  struct " + name + " final {\n";
  for (int i = 0; i < descriptor->field_count(); ++i) {
    FieldDescriptor const* field_descriptor = descriptor->field(i);
    std::string const& field_descriptor_name = field_descriptor->name();
    if (field_descriptors != nullptr) {
      field_descriptors->push_back(field_descriptor);
    }
    ProcessField(field_descriptor);

    // For in-out parameters, the code is generated only once, on the in
    // occurrence.
    bool const must_generate_code =
        name == in_message_name || !Contains(in_out_, field_descriptor);

    std::string const cxx_fill_member_name =
        ToLower(name) + "." + field_descriptor_name;
    std::string const cxx_run_field_getter =
        ToLower(name) + "." + field_descriptor_name + "()";
    std::string const run_local_variable = field_descriptor_name;

    cxx_fill_body_[descriptor] +=
        field_cxx_optional_assignment_fn_[field_descriptor](
            cxx_fill_member_name,
            field_cxx_assignment_fn_[field_descriptor](
                cxx_message_prefix,
                field_cxx_indirect_member_get_fn_[field_descriptor](
                    cxx_fill_member_name)));
    std::vector<std::string> const field_arguments =
        field_cxx_arguments_fn_[field_descriptor](run_local_variable);
    if (must_generate_code) {
      std::copy(field_arguments.begin(), field_arguments.end(),
                std::back_inserter(cxx_run_arguments_[descriptor]));

      if (Contains(out_, field_descriptor)) {
        cxx_run_body_prolog_[descriptor] +=
            "  " + field_cxx_type_[field_descriptor] + " " +
            run_local_variable + ";\n";
      } else {
        // If the field is optional, it needs extra storage for deserialization
        // (not passed to Deserialize) irrespective of whether its message type
        // needs extra storage.
        if (Contains(field_cxx_deserialization_storage_name_,
                     field_descriptor)) {
          std::string const cxx_deserialization_storage_declaration =
              "  " + field_cxx_deserialization_storage_type_[field_descriptor] +
              " " + field_cxx_deserialization_storage_name_[field_descriptor] +
              ";\n";
          cxx_run_body_prolog_[descriptor] +=
              cxx_deserialization_storage_declaration;
        }
        // If the field message type needs extra storage for deserialization
        // (passed to Deserialize) generate it now.
        Descriptor const* field_message_type = field_descriptor->message_type();
        if (Contains(cxx_deserialization_storage_declarations_,
                     field_message_type)) {
          cxx_run_body_prolog_[descriptor] +=
              cxx_deserialization_storage_declarations_[field_message_type];
        }
        cxx_run_body_prolog_[descriptor] +=
            "  auto " + run_local_variable + " = " +
            field_cxx_optional_pointer_fn_[field_descriptor](
                ToLower(name) + ".has_" + field_descriptor_name + "()",
                field_cxx_deserializer_fn_[field_descriptor](
                    cxx_run_field_getter)) +
            ";\n";
      }
    }
    if (Contains(field_cxx_deleter_fn_, field_descriptor)) {
      cxx_run_body_epilog_[descriptor] +=
          field_cxx_deleter_fn_[field_descriptor](cxx_run_field_getter);
    }
    if (Contains(field_cxx_inserter_fn_, field_descriptor)) {
      cxx_run_body_epilog_[descriptor] +=
          field_cxx_inserter_fn_[field_descriptor](
              ToLower(name) + "." + field_descriptor_name + "()",
              run_local_variable);
    }

    if (must_generate_code) {
      cs_interface_parameters_[descriptor].push_back(
          "  " + Join({HasMarshaler(field_descriptor)
                           ? "[" + MarshalAs(field_descriptor) + "]"
                           : "",
                       field_cs_mode_fn_[field_descriptor](
                           field_cs_type_[field_descriptor])},
                      /*joiner=*/" ") +
          " " + field_descriptor_name);
      cxx_interface_parameters_[descriptor].push_back(
          field_cxx_mode_fn_[field_descriptor](
              field_cxx_type_[field_descriptor]) +
          " " + field_descriptor_name);
    }
    cxx_nested_type_declaration_[descriptor] +=
        "    " + field_cxx_mode_fn_[field_descriptor](
                     field_cxx_type_[field_descriptor]) +
        " " + field_descriptor_name + ";\n";
  }
  cxx_nested_type_declaration_[descriptor] += "  };\n";
}

void JournalProtoProcessor::ProcessReturn(Descriptor const* descriptor) {
  CHECK_EQ(1, descriptor->field_count())
      << descriptor->full_name() << " must have exactly one field";
  FieldDescriptor const* field_descriptor = descriptor->field(0);
  FieldOptions const& field_options = field_descriptor->options();
  CHECK_EQ(FieldDescriptor::LABEL_REQUIRED, field_descriptor->label())
      << descriptor->full_name() << " must be required";
  return_.insert(field_descriptor);
  ProcessField(field_descriptor);
  cxx_fill_body_[descriptor] =
      field_cxx_assignment_fn_[field_descriptor](
          "message->mutable_return_()->",
          field_cxx_indirect_member_get_fn_[field_descriptor]("result"));
  std::string const cxx_field_getter =
      "message.return_()." + field_descriptor->name() + "()";
  if (Contains(field_cxx_inserter_fn_, field_descriptor)) {
    cxx_run_body_epilog_[descriptor] =
        field_cxx_inserter_fn_[field_descriptor](cxx_field_getter, "result");
  } else if (!field_options.HasExtension(journal::serialization::omit_check)) {
    Descriptor const* field_message_type = field_descriptor->message_type();
    if (Contains(cxx_deserialization_storage_declarations_,
                 field_message_type)) {
      cxx_run_body_epilog_[descriptor] =
          cxx_deserialization_storage_declarations_[field_message_type];
    }
    cxx_run_body_epilog_[descriptor] +=
        "  PRINCIPIA_CHECK_EQ(" +
        field_cxx_deserializer_fn_[field_descriptor](cxx_field_getter) + ", " +
        field_cxx_indirect_member_get_fn_[field_descriptor]("result") + ");\n";
  } else {
    CHECK(field_options.GetExtension(journal::serialization::omit_check))
      << field_descriptor->full_name() << " has incorrect (omit_check) option";
  }
  cs_interface_return_marshal_[descriptor] =
      HasMarshaler(field_descriptor)
          ? "[return : " + MarshalAs(field_descriptor) + "]"
          : "";
  cs_interface_return_type_[descriptor] = field_cs_type_[field_descriptor];
  cxx_interface_return_type_[descriptor] =
      field_cxx_mode_fn_[field_descriptor](field_cxx_type_[field_descriptor]);
  cxx_nested_type_declaration_[descriptor] =
      "  using Return = " + cxx_interface_return_type_[descriptor] + ";\n";
}

void JournalProtoProcessor::ProcessInterchangeMessage(
    Descriptor const* descriptor) {
  std::string const& name = descriptor->name();
  std::string const& parameter_name = ToLower(name);
  MessageOptions const& options = descriptor->options();

  cxx_serialize_definition_[descriptor] =
      "serialization::" + name + " Serialize" + name + "(" + name + " const& " +
      parameter_name + ") {\n  serialization::" + name + " m;\n";

  // Start by processing the fields.  We need to know if any of them has a
  // custom marshaler to decide whether we generate a struct or a class.
  // TODO(phl): Make this const once we generate the marshaler for
  // BodyParameters.
  bool has_custom_marshaler =
      options.HasExtension(journal::serialization::custom_marshaler);
  bool needs_custom_marshaler = false;
  for (int i = 0; i < descriptor->field_count(); ++i) {
    FieldDescriptor const* field_descriptor = descriptor->field(i);
    interchange_.insert(field_descriptor);
    ProcessField(field_descriptor);
    if (!has_custom_marshaler &&
        !field_cs_custom_marshaler_[field_descriptor].empty()) {
      needs_custom_marshaler = true;
    }
  }
  if (has_custom_marshaler) {
    // TODO(phl): Remove this hack once we generate the marshaler for
    // BodyParameters.
    if (options.GetExtension(journal::serialization::custom_marshaler) ==
        "none") {
      has_custom_marshaler = false;
      needs_custom_marshaler = false;
    } else {
      cs_custom_marshaler_name_[descriptor] =
          options.GetExtension(journal::serialization::custom_marshaler);
    }
  } else if (needs_custom_marshaler) {
    cs_custom_marshaler_name_[descriptor] = name + "Marshaler";
  }

  if (has_custom_marshaler || needs_custom_marshaler) {
    cs_interchange_type_declaration_[descriptor] =
        "internal partial class " + name + " {\n";
  } else {
    MessageOptions const& message_options = descriptor->options();
    std::string visibility = "internal";
    if (message_options.HasExtension(journal::serialization::is_public)) {
      CHECK(message_options.GetExtension(journal::serialization::is_public));
      visibility = "public";
    }
    cs_interchange_type_declaration_[descriptor] =
        "[StructLayout(LayoutKind.Sequential)]\n" + visibility +
        " partial struct " + name + " {\n";
  }
  // Produce a class-specific overload of new to be able to safely delete the
  // storage in principia__DeleteVoid.
  cxx_interchange_type_declaration_[descriptor] =
      "extern \"C\"\n"
      "struct " + name + " {\n"
      "  static void* operator new(std::size_t size) {\n"
      "    return ::operator new(size);\n"
      "  };\n";

  // Second pass on the fields to actually generate the code.
  std::vector<std::string> deserialized_expressions;
  for (int i = 0; i < descriptor->field_count(); ++i) {
    FieldDescriptor const* field_descriptor = descriptor->field(i);
    std::string const& field_descriptor_name = field_descriptor->name();

    // If the field needs extra storage for deserialization, generate it now.
    if (Contains(field_cxx_deserialization_storage_name_, field_descriptor)) {
      cxx_deserialization_storage_arguments_[descriptor] +=
          ", " + field_cxx_deserialization_storage_name_[field_descriptor];
      cxx_deserialization_storage_declarations_[descriptor] +=
          "  " + field_cxx_deserialization_storage_type_[field_descriptor] +
          " " + field_cxx_deserialization_storage_name_[field_descriptor] +
          ";\n";
      cxx_deserialization_storage_parameters_[descriptor] +=
          ", " +
          field_cxx_deserialization_storage_type_[field_descriptor] + "& " +
          field_cxx_deserialization_storage_name_[field_descriptor];
    }

    std::string const deserialize_field_checker =
        parameter_name + ".has_" + field_descriptor_name + "()";
    std::string const deserialize_field_getter =
        parameter_name + "." + field_descriptor_name + "()";
    std::string const serialize_member_name =
        parameter_name + "." + field_descriptor_name;

    deserialized_expressions.push_back(
        field_cxx_optional_pointer_fn_[field_descriptor](
            deserialize_field_checker,
            field_cxx_deserializer_fn_[field_descriptor](
                deserialize_field_getter)));
    cxx_serialize_definition_[descriptor] +=
        field_cxx_optional_assignment_fn_[field_descriptor](
            serialize_member_name,
            field_cxx_assignment_fn_[field_descriptor](
                "m.",
                field_cxx_indirect_member_get_fn_[field_descriptor](
                    serialize_member_name)));

    if (field_cs_private_type_[field_descriptor].empty()) {
      cs_interchange_type_declaration_[descriptor] +=
          "  public " + field_cs_type_[field_descriptor] + " " +
          field_descriptor_name + ";\n";
    } else {
      std::string const field_private_member_name = field_descriptor_name + "_";
      std::vector<std::string> fn_arguments = {field_private_member_name};
      cs_interchange_type_declaration_[descriptor] +=
          "  private " + field_cs_private_type_[field_descriptor] + " " +
          field_private_member_name + ";\n";
      cs_interchange_type_declaration_[descriptor] +=
          "  public " + field_cs_type_[field_descriptor] + " " +
          field_descriptor_name + " {\n" +
          "    " + field_cs_private_getter_fn_[field_descriptor](fn_arguments) +
          "\n" +
          "    " + field_cs_private_setter_fn_[field_descriptor](fn_arguments) +
          "\n" +
          "  }\n";
    }
    cxx_interchange_type_declaration_[descriptor] +=
        "  " + field_cxx_type_[field_descriptor] + " " + field_descriptor_name +
        ";\n";

    // If we need to generate a marshaler, do it now.
    if (needs_custom_marshaler) {
      cs_representation_type_declaration_[descriptor] += "    public ";
      if (field_cs_custom_marshaler_[field_descriptor].empty()) {
        cs_representation_type_declaration_[descriptor] +=
            field_cs_type_[field_descriptor] + " " + field_descriptor_name +
            ";\n";
        cs_managed_to_native_definition_[descriptor] +=
            "        " + field_descriptor_name + " = value." +
            field_descriptor_name + ",\n";
        cs_native_to_managed_definition_[descriptor] +=
            "        " + field_descriptor_name + " = representation." +
            field_descriptor_name + ",\n";
      } else {
        cs_representation_type_declaration_[descriptor] +=
            "IntPtr " + field_descriptor_name + ";\n";
        cs_clean_up_native_definition_[descriptor] +=
            "    " + field_cs_custom_marshaler_[field_descriptor] +
            ".GetInstance(null).CleanUpNativeData(representation." +
            field_descriptor_name + ");\n";
        cs_managed_to_native_definition_[descriptor] +=
            "        " + field_descriptor_name + " = " +
            field_cs_custom_marshaler_[field_descriptor] +
            ".GetInstance(null).MarshalManagedToNative(value." +
            field_descriptor_name + "),\n";
        cs_native_to_managed_definition_[descriptor] +=
            "        " + field_descriptor_name + " = " +
            field_cs_custom_marshaler_[field_descriptor] +
            ".GetInstance(null).MarshalNativeToManaged(representation." +
            field_descriptor_name + ") as " + field_cs_type_[field_descriptor] +
            ",\n";
      }
    }
  }
  cxx_serialize_definition_[descriptor] += "  return m;\n}\n\n";

  cxx_deserialize_definition_[descriptor] =
      name + " Deserialize" + name + "(serialization::" + name + " const& " +
      parameter_name + ", Player::PointerMap& pointer_map" +
      cxx_deserialization_storage_parameters_[descriptor] + ") {\n  return {";
  cxx_deserialize_definition_[descriptor] +=
      Join(deserialized_expressions, /*joiner=*/",\n          ") +  // NOLINT
      "};\n}\n\n";

  cs_interchange_type_declaration_[descriptor] += "}\n\n";
  cxx_interchange_type_declaration_[descriptor] +=
      "};\n\nstatic_assert(std::is_pod<" + name +
      ">::value,\n              \"" + name + " is used for interfacing\");\n\n";

  if (needs_custom_marshaler) {
    cs_custom_marshaler_class_[descriptor] =
        "internal class " + cs_custom_marshaler_name_[descriptor] +
        " : MonoMarshaler {\n"
        "  [StructLayout(LayoutKind.Sequential)]\n"
        "  internal struct Representation {\n" +
        cs_representation_type_declaration_[descriptor] +
        "  }\n\n"
        "  public static ICustomMarshaler GetInstance(string s) {\n"
        "    return instance_;\n"
        "  }\n\n"
        "  public override void CleanUpNativeDataImplementation("
        "IntPtr native_data) {\n"
        "    var representation = (Representation)Marshal.PtrToStructure("
        "native_data, typeof(Representation));\n" +
        cs_clean_up_native_definition_[descriptor] +
        "    Marshal.FreeHGlobal(native_data);\n"
        "  }\n\n"
        "  public override IntPtr MarshalManagedToNativeImplementation("
        "object managed_object) {\n"
        "    if (!(managed_object is " + name + " value)) {\n"
        "      throw new NotSupportedException();\n"
        "    }\n"
        "    var representation = new Representation{\n" +
        cs_managed_to_native_definition_[descriptor] +
        "    };\n"
        "    IntPtr buffer = Marshal.AllocHGlobal("
        "Marshal.SizeOf(representation));\n"
        "    Marshal.StructureToPtr("
        "representation, buffer, fDeleteOld: false);\n" +
        "    return buffer;\n"
        "  }\n\n"
        "  public override object MarshalNativeToManaged("
        "IntPtr native_data) {\n"
        "    var representation = (Representation)Marshal.PtrToStructure("
        "native_data, typeof(Representation));\n"
        "    return new " + name + "{\n" +
        cs_native_to_managed_definition_[descriptor] +
        "    };\n"
        "  }\n\n"
        "  private static readonly " + cs_custom_marshaler_name_[descriptor] +
        " instance_ =\n"
        "      new " + cs_custom_marshaler_name_[descriptor] + "();\n"
        "}\n\n";
  }
}

void JournalProtoProcessor::ProcessMethodExtension(
    Descriptor const* descriptor) {
  std::string const& name = descriptor->name();
  bool has_in = false;
  bool has_out = false;
  bool has_return = false;

  // Do a first pass to determine which fields are in-out.  The data produced
  // here will be overwritten by the next pass.
  std::vector<FieldDescriptor const*> field_descriptors;
  for (int i = 0; i < descriptor->nested_type_count(); ++i) {
    Descriptor const* nested_descriptor = descriptor->nested_type(i);
    const std::string& nested_name = nested_descriptor->name();
    if (nested_name == in_message_name) {
      has_in = true;
      std::vector<FieldDescriptor const*> in_field_descriptors;
      ProcessInOut(nested_descriptor, &in_field_descriptors);
      in_.insert(in_field_descriptors.begin(), in_field_descriptors.end());
      std::copy(in_field_descriptors.begin(),
                in_field_descriptors.end(),
                std::back_inserter(field_descriptors));
    } else if (nested_name == out_message_name) {
      has_out = true;
      std::vector<FieldDescriptor const*> out_field_descriptors;
      ProcessInOut(nested_descriptor, &out_field_descriptors);
      out_.insert(out_field_descriptors.begin(), out_field_descriptors.end());
      std::copy(out_field_descriptors.begin(),
                out_field_descriptors.end(),
                std::back_inserter(field_descriptors));
    } else if (nested_name == return_message_name) {
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
              [](FieldDescriptor const* left, FieldDescriptor const* right) {
      return left->name() < right->name();
    });
    for (int i = 0; i < field_descriptors.size() - 1; ++i) {
      if (field_descriptors[i]->name() == field_descriptors[i + 1]->name()) {
        in_out_.insert(field_descriptors[i]);
        in_out_.insert(field_descriptors[i + 1]);
      }
    }
  }

  // The second pass that produces the actual output.
  std::vector<std::string> cs_interface_parameters;
  std::vector<std::string> cxx_interface_parameters;
  std::vector<std::string> cxx_run_arguments;
  std::string cs_interface_return_marshal;
  std::string cs_interface_return_type = "void";
  std::string cxx_interface_return_type = "void";
  std::string cxx_run_prolog;
  std::string cxx_run_epilog;
  cxx_toplevel_type_declaration_[descriptor] =
      "struct " + name + " : not_constructible {\n";
  for (int i = 0; i < descriptor->nested_type_count(); ++i) {
    Descriptor const* nested_descriptor = descriptor->nested_type(i);
    const std::string& nested_name = nested_descriptor->name();
    if (nested_name == in_message_name) {
      ProcessInOut(nested_descriptor, /*field_descriptors=*/nullptr);
      cxx_functions_implementation_[descriptor] +=
          "void " + name + "::Fill(In const& in, "
          "not_null<Message*> const message) {\n" +
          cxx_fill_body_[nested_descriptor] +
          "}\n\n";
      cxx_run_prolog += cxx_run_body_prolog_[nested_descriptor];
      std::copy(cs_interface_parameters_[nested_descriptor].begin(),
                cs_interface_parameters_[nested_descriptor].end(),
                std::back_inserter(cs_interface_parameters));
      std::copy(cxx_interface_parameters_[nested_descriptor].begin(),
                cxx_interface_parameters_[nested_descriptor].end(),
                std::back_inserter(cxx_interface_parameters));
      std::copy(cxx_run_arguments_[nested_descriptor].begin(),
                cxx_run_arguments_[nested_descriptor].end(),
                std::back_inserter(cxx_run_arguments));
    } else if (nested_name == out_message_name) {
      ProcessInOut(nested_descriptor, /*field_descriptors=*/nullptr);
      cxx_functions_implementation_[descriptor] +=
          "void " + name + "::Fill(Out const& out, "
          "not_null<Message*> const message) {\n" +
          cxx_fill_body_[nested_descriptor] +
          "}\n\n";
      cxx_run_prolog += cxx_run_body_prolog_[nested_descriptor];
      std::copy(cs_interface_parameters_[nested_descriptor].begin(),
                cs_interface_parameters_[nested_descriptor].end(),
                std::back_inserter(cs_interface_parameters));
      std::copy(cxx_interface_parameters_[nested_descriptor].begin(),
                cxx_interface_parameters_[nested_descriptor].end(),
                std::back_inserter(cxx_interface_parameters));
      std::copy(cxx_run_arguments_[nested_descriptor].begin(),
                cxx_run_arguments_[nested_descriptor].end(),
                std::back_inserter(cxx_run_arguments));
    } else if (nested_name == return_message_name) {
      ProcessReturn(nested_descriptor);
      cxx_functions_implementation_[descriptor] +=
          "void " + name + "::Fill("
          "Return const& result, "
          "not_null<Message*> const message) {\n" +
          cxx_fill_body_[nested_descriptor] +
          "}\n\n";
      cs_interface_return_marshal =
          cs_interface_return_marshal_[nested_descriptor];
      cs_interface_return_type = cs_interface_return_type_[nested_descriptor];
      cxx_interface_return_type = cxx_interface_return_type_[nested_descriptor];
    }
    cxx_run_epilog += cxx_run_body_epilog_[nested_descriptor];
    cxx_toplevel_type_declaration_[descriptor] +=
        cxx_nested_type_declaration_[nested_descriptor];
  }
  if (has_in || has_out || has_return) {
    cxx_toplevel_type_declaration_[descriptor] += "\n";
  }
  cxx_toplevel_type_declaration_[descriptor] +=
      "  using Message = serialization::" + name + ";\n";
  if (has_in) {
    cxx_toplevel_type_declaration_[descriptor] +=
        "  static void Fill(In const& in, "
        "not_null<Message*> const message);\n";
  }
  if (has_out) {
    cxx_toplevel_type_declaration_[descriptor] +=
        "  static void Fill(Out const& out, "
        "not_null<Message*> const message);\n";
  }
  if (has_return) {
    cxx_toplevel_type_declaration_[descriptor] +=
        "  static void Fill("
        "Return const& result, "
        "not_null<Message*> const message);\n";
  }
  cxx_toplevel_type_declaration_[descriptor] +=
      "  static void Run(Message const& message,\n"
      "                  Player::PointerMap& pointer_map);\n";
  cxx_toplevel_type_declaration_[descriptor] += "};\n\n";

  // The Run method must come after the Fill methods for comparison with manual
  // code.
  cxx_functions_implementation_[descriptor] +=
      "void " + name + "::Run(Message const& message, "
      "Player::PointerMap& pointer_map) {\n" +
      cxx_run_prolog;
  MessageOptions const& options = descriptor->options();
  std::string const run_conditional_compilation_symbol =
      options.HasExtension(
          journal::serialization::run_conditional_compilation_symbol)
          ? options.GetExtension(
                journal::serialization::run_conditional_compilation_symbol)
          : "";
  if (!run_conditional_compilation_symbol.empty()) {
    cxx_functions_implementation_[descriptor] +=
        "#if " + run_conditional_compilation_symbol + "\n";
  }
  if (has_return) {
    cxx_functions_implementation_[descriptor] += "  auto const result = ";
  } else {
    cxx_functions_implementation_[descriptor] += "  ";
  }
  cxx_functions_implementation_[descriptor] +=
      "interface::principia__" + name + "(" +
      Join(cxx_run_arguments, /*joiner=*/", ") + ");\n";
  if (!run_conditional_compilation_symbol.empty()) {
    cxx_functions_implementation_[descriptor] += "#endif\n";
  }
  cxx_functions_implementation_[descriptor] += cxx_run_epilog + "}\n\n";

  cs_interface_method_declaration_[descriptor] = Join(
      {"  [DllImport(dllName           : dll_path,\n"
       "             EntryPoint        = \"principia__" + name + "\",\n"
       "             CallingConvention = CallingConvention.Cdecl)]",
       cs_interface_return_marshal,
       "internal static extern " + cs_interface_return_type + " " + name + "("},
      "\n  ");
  if (!cs_interface_parameters.empty()) {
    cs_interface_method_declaration_[descriptor] +=
        "\n    " + Join(cs_interface_parameters, /*joiner=*/",\n    ");  // NOLINT
  }
  cs_interface_method_declaration_[descriptor] += ");\n\n";

  cxx_interface_method_declaration_[descriptor] =
      "extern \"C\" PRINCIPIA_DLL\n" +
  cxx_interface_return_type + " __cdecl principia__" + name + "(";
  if (!cxx_interface_parameters.empty()) {
    cxx_interface_method_declaration_[descriptor] +=
        "\n    " + Join(cxx_interface_parameters, /*joiner=*/",\n    ");  // NOLINT
  }
  cxx_interface_method_declaration_[descriptor] += ");\n\n";

  cxx_play_statement_[descriptor] =
      "  ran |= RunIfAppropriate<" + name + ">(\n"
      "             *method_in, *method_out_return);\n";
}

bool JournalProtoProcessor::HasMarshaler(
    FieldDescriptor const* descriptor) const {
  auto const it_custom = field_cs_custom_marshaler_.find(descriptor);
  auto const it_predefined = field_cs_predefined_marshaler_.find(descriptor);
  bool const has_custom = it_custom != field_cs_custom_marshaler_.end() &&
                          !it_custom->second.empty();
  bool const has_predefined =
      it_predefined != field_cs_predefined_marshaler_.end() &&
      !it_predefined->second.empty();
  CHECK(!(has_custom && has_predefined)) << descriptor->name();
  return has_custom || has_predefined;
}

std::string JournalProtoProcessor::MarshalAs(
    FieldDescriptor const* descriptor) const {
  auto const it_custom = field_cs_custom_marshaler_.find(descriptor);
  auto const it_predefined = field_cs_predefined_marshaler_.find(descriptor);
  if (it_custom != field_cs_custom_marshaler_.end() &&
      !it_custom->second.empty()) {
    return "MarshalAs(UnmanagedType.CustomMarshaler, "
           "MarshalTypeRef = typeof(" +
           it_custom->second + "))";
  }
  if (it_predefined != field_cs_predefined_marshaler_.end() &&
      !it_predefined->second.empty()) {
    return "MarshalAs(" + it_predefined->second + ")";
  }
  LOG(FATAL) << "Bad marshaler for " << descriptor->name();
}

}  // namespace internal_journal_proto_processor
}  // namespace tools
}  // namespace principia
