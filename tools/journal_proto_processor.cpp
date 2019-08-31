
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

char const method_message_name[] = "Method";
char const in_message_name[] = "In";
char const return_message_name[] = "Return";
char const out_message_name[] = "Out";

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

std::string const& GetSingleElement(std::vector<std::string> const& v) {
  CHECK_EQ(1, v.size());
  return v.front();
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
JournalProtoProcessor::GetCsInterfaceTypeDeclarations() const {
  std::vector<std::string> result;
  for (auto const& pair : cs_interface_type_declaration_) {
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
JournalProtoProcessor::GetCxxInterfaceTypeDeclarations() const {
  std::vector<std::string> result;
  for (auto const& pair : cxx_interface_type_declaration_) {
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
  std::string const& message_type_name = descriptor->message_type()->name();

  FieldOptions const& options = descriptor->options();
  CHECK(options.HasExtension(journal::serialization::size))
      << descriptor->full_name() << " is missing a (size) option";
  size_member_name_[descriptor] =
      options.GetExtension(journal::serialization::size);
  field_cs_type_[descriptor] = message_type_name + "[]";
  field_cxx_type_[descriptor] = message_type_name + " const*";

  field_cxx_arguments_fn_[descriptor] =
      [](std::string const& identifier) -> std::vector<std::string> {
        return {"&" + identifier + "[0]", identifier + ".size()"};
      };
  field_cxx_assignment_fn_[descriptor] =
      [this, descriptor, message_type_name](
          std::string const& prefix, std::string const& expr) {
        std::string const& descriptor_name = descriptor->name();
        // The use of |substr| below is a bit of a cheat because we known the
        // structure of |expr|.
        return "  for (" + message_type_name + " const* " + descriptor_name +
               " = " + expr + "; " + descriptor_name + " < " + expr + " + " +
               expr.substr(0, expr.find('.')) + "." +
               size_member_name_[descriptor] + "; ++" + descriptor_name +
               ") {\n    *" + prefix + "add_" + descriptor_name +
               "() = " +
               field_cxx_serializer_fn_[descriptor]("*"+ descriptor_name) +
               ";\n  }\n";
      };
  std::string const storage_name = descriptor->name() + "_storage";
  field_cxx_deserialization_storage_name_[descriptor] = storage_name;
  field_cxx_deserializer_fn_[descriptor] =
      [message_type_name, storage_name](
          std::string const& expr) -> std::vector<std::string> {
        // Yes, this lambda generates a lambda.
        return {"[&" + storage_name +
                "](::google::protobuf::RepeatedPtrField<serialization::" +
                message_type_name + "> const& messages) {\n"
                "            for (auto const& message : messages) {\n" +
                "              " + storage_name +
                ".push_back(Deserialize" + message_type_name + "(message));\n" +
                "            }\n"
                "            return &" + storage_name + "[0];\n" +
                "          }(" + expr + ")",
                expr + ".size()"};
      };
  field_cxx_deserialization_storage_type_[descriptor] =
      "std::vector<" + message_type_name + ">";
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
  // It is not possible to use a custom marshaler on an |T?|, as this raises
  // |System.Runtime.InteropServices.MarshalDirectiveException| with the message
  // "Custom marshalers are only allowed on classes, strings, arrays, and boxed
  // value types.".  We could use a boxed |T|, whose type would be |object|, but
  // we would lose static typing.  We use a custom strongly-typed boxed type
  // instead.
  field_cs_type_[descriptor] = cs_boxed_type;
  field_cs_marshal_[descriptor] =
      "MarshalAs(UnmanagedType.CustomMarshaler, "
      "MarshalTypeRef = typeof(OptionalMarshaler<" + cs_unboxed_type + ">))";
  field_cxx_type_[descriptor] = cxx_type + " const*";

  field_cxx_arguments_fn_[descriptor] =
      [](std::string const& identifier) -> std::vector<std::string> {
        return {identifier + ".get()"};
      };
  field_cxx_indirect_member_get_fn_[descriptor] =
      [](std::string const& expr) {
        return "*" + expr;
      };
  field_cxx_optional_pointer_fn_[descriptor] = [cxx_type](
      std::string const& condition,
      std::string const& expr) {
    // Tricky.  We need a heap allocation to obtain a pointer to the value.
    return condition + " ? std::make_unique<" + cxx_type + " const>(" + expr +
           ") : nullptr";
  };
}

void JournalProtoProcessor::ProcessOptionalDoubleField(
    FieldDescriptor const* descriptor) {
  ProcessOptionalNonStringField(
      descriptor,
      /*cs_boxed_type=*/"BoxedDouble",
      /*cs_unboxed_type=*/"double",
      /*cxx_type=*/"double");
}

void JournalProtoProcessor::ProcessOptionalInt32Field(
    FieldDescriptor const* descriptor) {
  ProcessOptionalNonStringField(
      descriptor,
      /*cs_boxed_type=*/"BoxedInt32",
      /*cs_unboxed_type=*/"int",
      /*cxx_type=*/"int");
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
  ProcessSingleMessageField(descriptor);
}

void JournalProtoProcessor::ProcessOptionalStringField(
    FieldDescriptor const* descriptor) {
  ProcessSingleStringField(descriptor);
  field_cs_type_[descriptor] = "string?";
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
    field_cs_marshal_[descriptor] =
        "MarshalAs(UnmanagedType.CustomMarshaler, "
        "MarshalTypeRef = typeof(" +
        options.GetExtension(journal::serialization::disposable) +
        "Marshaller))";
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
          return "  Delete(pointer_map, " + expr + ");\n";
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
                 ") {\n    Delete(pointer_map, " + expr + ");\n  }\n";
        };
  }
  if (options.HasExtension(journal::serialization::is_produced)) {
    CHECK(options.GetExtension(journal::serialization::is_produced))
        << descriptor->full_name() << " has incorrect (is_produced) option";
    field_cxx_inserter_fn_[descriptor] =
        [](std::string const& expr1, std::string const& expr2) {
          return "  Insert(pointer_map, " + expr1 + ", " + expr2 + ");\n";
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
                 ") {\n    Insert(pointer_map, " + expr1 + ", " + expr2 +
                 ");\n  }\n";
        };
  }

  // Special handlings for produced C-style strings: these are seen from the C#
  // as strings, and marshalled with immediate destruction.
  if (options.HasExtension(journal::serialization::encoding) &&
      (options.HasExtension(journal::serialization::is_produced) ||
       options.HasExtension(journal::serialization::is_produced_if))) {
    field_cs_type_[descriptor] = "String";
    switch (options.GetExtension(journal::serialization::encoding)) {
      case journal::serialization::UTF_8:
        field_cs_marshal_[descriptor] =
            "MarshalAs(UnmanagedType.CustomMarshaler, "
            "MarshalTypeRef = typeof(OutOwnedUTF8Marshaler))";
        break;
      case journal::serialization::UTF_16:
        field_cs_marshal_[descriptor] =
            "MarshalAs(UnmanagedType.CustomMarshaler, "
            "MarshalTypeRef = typeof(OutOwnedUTF16Marshaler))";
        break;
    }
  }

  field_cxx_deserializer_fn_[descriptor] =
      [pointer_to](std::string const& expr) -> std::vector<std::string> {
        return {"DeserializePointer<" + pointer_to + "*>(pointer_map, " + expr +
                ")"};
      };
  field_cxx_serializer_fn_[descriptor] =
      [](std::string const& expr) {
        return "SerializePointer(" + expr + ")";
      };
}

void JournalProtoProcessor::ProcessRequiredMessageField(
    FieldDescriptor const* descriptor) {
  Descriptor const* message_type = descriptor->message_type();
  std::string const& message_type_name = message_type->name();
  field_cs_type_[descriptor] = message_type_name;
  field_cxx_type_[descriptor] = message_type_name;

  MessageOptions const& message_options = message_type->options();
  if (Contains(in_, descriptor) &&
      message_options.HasExtension(
          journal::serialization::in_custom_marshaler)) {
    field_cs_marshal_[descriptor] =
        "MarshalAs(UnmanagedType.CustomMarshaler, MarshalTypeRef = typeof(" +
        message_options.GetExtension(
            journal::serialization::in_custom_marshaler) + "))";
    field_cxx_mode_fn_[descriptor] =
        [](std::string const& type) {
          return type + " const&";
        };
  }

  ProcessSingleMessageField(descriptor);
}

void JournalProtoProcessor::ProcessRequiredBoolField(
    FieldDescriptor const* descriptor) {
  field_cs_type_[descriptor] = "bool";
  field_cs_marshal_[descriptor] = "MarshalAs(UnmanagedType.I1)";
  field_cs_private_type_[descriptor] = "Byte";
  field_cs_private_getter_fn_[descriptor] =
      [](std::vector<std::string> const& identifiers) {
        CHECK_EQ(1, identifiers.size());
        return "get { return " + identifiers[0] + " != (Byte)0; }";
      };
  field_cs_private_setter_fn_[descriptor] =
      [](std::vector<std::string> const& identifiers) {
        CHECK_EQ(1, identifiers.size());
        return "set { " + identifiers[0] + " = value ? (Byte)1 : (Byte)0; }";
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

  field_cs_marshal_[descriptor] = "MarshalAs(UnmanagedType.LPWStr)";
  field_cs_type_[descriptor] = "String";
  field_cxx_type_[descriptor] = "char16_t const*";
  field_cxx_arguments_fn_[descriptor] =
      [](std::string const& identifier) -> std::vector<std::string> {
        return {identifier + ".c_str()"};
      };
  field_cxx_deserializer_fn_[descriptor] =
      [](std::string const& expr) -> std::vector<std::string> {
        return {"DeserializeUtf16(" + expr + ")"};
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
  field_cs_type_[descriptor] = "Int64";
  field_cxx_type_[descriptor] = "std::int64_t";
}

void JournalProtoProcessor::ProcessRequiredUint32Field(
    FieldDescriptor const* descriptor) {
  field_cs_type_[descriptor] = "uint";
  field_cxx_type_[descriptor] = "uint32_t";
}

void JournalProtoProcessor::ProcessRequiredStringField(
    FieldDescriptor const* descriptor) {
  ProcessSingleStringField(descriptor);
  field_cs_type_[descriptor] = "string";
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
  const std::string deserialization_storage_arguments =
      cxx_deserialization_storage_arguments_[message_type];
  field_cxx_deserializer_fn_[descriptor] =
      [message_type_name, deserialization_storage_arguments](
          std::string const& expr) -> std::vector<std::string> {
        return {"Deserialize" + message_type_name + "(" + expr +
                deserialization_storage_arguments + ")"};
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

  // Note that it is important to use an out marshmallow for return fields,
  // hence the use of the |in_| set here.
  field_cs_marshal_[descriptor] =
      Contains(in_, descriptor) ? "MarshalAs(UnmanagedType.CustomMarshaler, "
                                  "MarshalTypeRef = typeof(InUTF8Marshaler))"
                                : "MarshalAs(UnmanagedType.CustomMarshaler, "
                                  "MarshalTypeRef = typeof(OutUTF8Marshaler))";
  field_cxx_type_[descriptor] = "char const*";
  if (options.HasExtension(journal::serialization::size)) {
    size_member_name_[descriptor] =
        options.GetExtension(journal::serialization::size);

    field_cxx_arguments_fn_[descriptor] =
        [](std::string const& identifier) -> std::vector<std::string> {
          return {identifier + "->c_str()", identifier + "->size()"};
        };
    field_cxx_deserializer_fn_[descriptor] =
        [](std::string const& expr) -> std::vector<std::string> {
          return {"&" + expr};
        };
    field_cxx_indirect_member_get_fn_[descriptor] =
        [this, descriptor](std::string const& expr) {
          return "std::string(" + expr + ", " +
                 expr.substr(0, expr.find('.')) + "." +
                 size_member_name_[descriptor] + ")";
        };
  } else {
    field_cxx_deserializer_fn_[descriptor] =
        [](std::string const& expr) -> std::vector<std::string> {
          return {expr + ".c_str()"};
        };
  }
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
      [](std::string const& expr) -> std::vector<std::string> {
        return {expr};
      };
  field_cxx_mode_fn_[descriptor] =
      [](std::string const& type) {
        return type + " const";
      };
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
      "  auto const& " + ToLower(name) + " = message." +
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
                GetSingleElement(field_cxx_deserializer_fn_[field_descriptor](
                    cxx_run_field_getter))) +
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
          "  " + Join({field_cs_marshal_[field_descriptor].empty()
                           ? ""
                           : "[" + field_cs_marshal_[field_descriptor] + "]",
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

    // If this field has a size, generate it now.
    if (Contains(size_member_name_, field_descriptor)) {
      if (must_generate_code) {
        cs_interface_parameters_[descriptor].push_back(
            "  int " + size_member_name_[field_descriptor]);
        cxx_interface_parameters_[descriptor].push_back(
            "int const " + size_member_name_[field_descriptor]);
      }
      cxx_nested_type_declaration_[descriptor] +=
          "    int const " + size_member_name_[field_descriptor] + ";\n";
    }
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
  ProcessField(field_descriptor);
  cxx_fill_body_[descriptor] =
      field_cxx_assignment_fn_[field_descriptor]("message->mutable_return_()->",
                                             "result");
  std::string const cxx_field_getter =
      "message.return_()." + field_descriptor->name() + "()";
  if (Contains(field_cxx_inserter_fn_, field_descriptor)) {
    cxx_run_body_epilog_[descriptor] =
        field_cxx_inserter_fn_[field_descriptor](cxx_field_getter, "result");
  } else if (!field_options.HasExtension(journal::serialization::omit_check)) {
    cxx_run_body_epilog_[descriptor] =
        "  PRINCIPIA_CHECK_EQ(" +
        GetSingleElement(
            field_cxx_deserializer_fn_[field_descriptor](cxx_field_getter)) +
        ", result);\n";
  } else {
    CHECK(field_options.GetExtension(journal::serialization::omit_check))
      << field_descriptor->full_name() << " has incorrect (omit_check) option";
  }
  cs_interface_return_marshal_[descriptor] =
      field_cs_marshal_[field_descriptor].empty()
          ? ""
          : "[return : " + field_cs_marshal_[field_descriptor] + "]";
  cs_interface_return_type_[descriptor] = field_cs_type_[field_descriptor];
  cxx_interface_return_type_[descriptor] = field_cxx_type_[field_descriptor];
  cxx_nested_type_declaration_[descriptor] =
      "  using Return = " + field_cxx_type_[field_descriptor] + ";\n";
}

void JournalProtoProcessor::ProcessInterchangeMessage(
    Descriptor const* descriptor) {
  std::string const& name = descriptor->name();
  std::string const& parameter_name = ToLower(name);

  cxx_serialize_definition_[descriptor] =
      "serialization::" + name + " Serialize" + name + "(" + name + " const& " +
      parameter_name + ") {\n  serialization::" + name + " m;\n";

  MessageOptions const& options = descriptor->options();
  if (options.HasExtension(journal::serialization::in_custom_marshaler)) {
    cs_interface_type_declaration_[descriptor] =
        "internal partial class " + name + " {\n";
  } else {
    MessageOptions const& message_options = descriptor->options();
    std::string visibility = "internal";
    if (message_options.HasExtension(journal::serialization::is_public)) {
      CHECK(message_options.GetExtension(journal::serialization::is_public));
      visibility = "public";
    }
    cs_interface_type_declaration_[descriptor] =
        "[StructLayout(LayoutKind.Sequential)]\n" + visibility +
        " partial struct " + name + " {\n";
  }
  cxx_interface_type_declaration_[descriptor] =
      "extern \"C\"\nstruct " + name + " {\n";

  std::vector<std::string> deserialized_expressions;
  for (int i = 0; i < descriptor->field_count(); ++i) {
    FieldDescriptor const* field_descriptor = descriptor->field(i);
    std::string const& field_descriptor_name = field_descriptor->name();
    ProcessField(field_descriptor);

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

    // There may be several expressions to deserialize a field because of the
    // size.
    std::vector<std::string> deserialize_fields =
        field_cxx_deserializer_fn_[field_descriptor](deserialize_field_getter);
    for (std::string const& deserialize_field : deserialize_fields) {
      deserialized_expressions.push_back(
          field_cxx_optional_pointer_fn_[field_descriptor](
              deserialize_field_checker,
              deserialize_field));
    }
    cxx_serialize_definition_[descriptor] +=
        field_cxx_optional_assignment_fn_[field_descriptor](
            serialize_member_name,
            field_cxx_assignment_fn_[field_descriptor](
                "m.", serialize_member_name));

    // TODO(phl): field_cs_private_type_ should be set iff field_cs_marshal_ is
    // set.  This is not the case at the moment because of strings being passed
    // in the interchange messages.  This will need fixing if we ever want to
    // pass non-ASCII strings or to return these structs.
    if (field_cs_private_type_[field_descriptor].empty()) {
      cs_interface_type_declaration_[descriptor] +=
          "  public " + field_cs_type_[field_descriptor] + " " +
          field_descriptor_name + ";\n";
    } else {
      std::string const field_private_member_name = field_descriptor_name + "_";
      std::vector<std::string> fn_arguments = {field_private_member_name};
      cs_interface_type_declaration_[descriptor] +=
          "  private " + field_cs_private_type_[field_descriptor] + " " +
          field_private_member_name + ";\n";
      cs_interface_type_declaration_[descriptor] +=
          "  public " + field_cs_type_[field_descriptor] + " " +
          field_descriptor_name + " {\n" +
          "    " + field_cs_private_getter_fn_[field_descriptor](fn_arguments) +
          "\n" +
          "    " + field_cs_private_setter_fn_[field_descriptor](fn_arguments) +
          "\n" +
          "  }\n";
    }
    cxx_interface_type_declaration_[descriptor] +=
        "  " + field_cxx_type_[field_descriptor] + " " + field_descriptor_name +
        ";\n";

    // If this field has a size, generate it now.
    if (Contains(size_member_name_, field_descriptor)) {
      cxx_interface_type_declaration_[descriptor] +=
          "  int " + size_member_name_[field_descriptor] + ";\n";
    }
  }
  cxx_serialize_definition_[descriptor] += "  return m;\n}\n\n";

  cxx_deserialize_definition_[descriptor] =
      name + " Deserialize" + name + "(serialization::" + name + " const& " +
      parameter_name + cxx_deserialization_storage_parameters_[descriptor] +
      ") {\n  return {";
  cxx_deserialize_definition_[descriptor] +=
      Join(deserialized_expressions, /*joiner=*/",\n          ") +  // NOLINT
      "};\n}\n\n";

  cs_interface_type_declaration_[descriptor] += "}\n\n";
  cxx_interface_type_declaration_[descriptor] +=
      "};\n\nstatic_assert(std::is_pod<" + name +
      ">::value,\n              \"" + name + " is used for interfacing\");\n\n";
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
  cxx_interface_return_type + " CDECL principia__" + name + "(";
  if (!cxx_interface_parameters.empty()) {
    cxx_interface_method_declaration_[descriptor] +=
        "\n    " + Join(cxx_interface_parameters, /*joiner=*/",\n    ");  // NOLINT
  }
  cxx_interface_method_declaration_[descriptor] += ");\n\n";

  cxx_play_statement_[descriptor] =
      "  ran |= RunIfAppropriate<" + name + ">(\n"
      "             *method_in, *method_out_return);\n";
}

}  // namespace internal_journal_proto_processor
}  // namespace tools
}  // namespace principia
