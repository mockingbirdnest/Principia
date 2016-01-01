#include "tools/journal_proto_processor.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "glog/logging.h"
#include "serialization/journal.pb.h"

namespace principia {
namespace tools {

namespace {

char const kIn[] = "In";
char const kReturn[] = "Return";
char const kOut[] = "Out";

template<typename Container>
bool Contains(Container const& container,
              typename Container::key_type const key) {
  return container.find(key) != container.end();
}

std::string Join(std::vector<std::string> const& v) {
  std::string joined;
  for (int i = 0; i < v.size(); ++i) {
    if (i == 0) {
      joined = v[i];
    } else {
      joined += ", " + v[i];
    }
  }
  return joined;
}

std::string ToLower(std::string const& s) {
  std::string lower(s.size(), ' ');
  for (int i = 0; i < s.size(); ++i) {
    lower[i] = std::tolower(s[i]);
  }
  return lower;
}

}  // namespace

void JournalProtoProcessor::ProcessMethods() {
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

std::vector<std::string>
JournalProtoProcessor::GetCppMethodImplementations() const {
    std::vector<std::string> result;
  for (auto const& pair : functions_implementation_) {
    result.push_back(pair.second);
  }
  return result;
}

std::vector<std::string> JournalProtoProcessor::GetCppMethodTypes() const {
    std::vector<std::string> result;
  for (auto const& pair : toplevel_type_declaration_) {
    result.push_back(pair.second);
  }
  return result;
}

void JournalProtoProcessor::ProcessRepeatedMessageField(
    FieldDescriptor const* descriptor) {
  std::string const& message_type_name = descriptor->message_type()->name();

  FieldOptions const& options = descriptor->options();
  CHECK(options.HasExtension(serialization::size))
      << descriptor->full_name() << " is missing a (size) option";
  size_member_name_[descriptor] = options.GetExtension(serialization::size);
  field_type_[descriptor] = message_type_name + " const*";

  field_arguments_fn_[descriptor] =
      [](std::string const& identifier) -> std::vector<std::string> {
        return {"&" + identifier + "[0]", identifier + ".size()"};
      };
  field_assignment_fn_[descriptor] =
      [this, descriptor, message_type_name](
          std::string const& identifier, std::string const& expr) {
        std::string const& descriptor_name = descriptor->name();
        // The use of |substr| below is a bit of a cheat because we known the
        // structure of |expr|.
        return "  for (" + message_type_name + " const* " + descriptor_name +
               " = " + expr + "; " + descriptor_name + " < " + expr + " + " +
               expr.substr(0, expr.find('.')) + "." +
               size_member_name_[descriptor] + "; ++" + descriptor_name +
               ") {\n    *" + identifier + "->add_" + descriptor_name +
               "() = " +
               field_serializer_fn_[descriptor]("*"+ descriptor_name) +
               ";\n  }\n";
      };
  field_deserializer_fn_[descriptor] =
      [descriptor, message_type_name](std::string const& expr) {
        std::string const& descriptor_name = descriptor->name();
        // Yes, this lambda generates a lambda.
        return "[](::google::protobuf::RepeatedPtrField<serialization::" +
               message_type_name + "> const& messages) -> std::vector<" +
               message_type_name + "> {\n"
               "      std::vector<" + message_type_name + "> deserialized_" +
               descriptor_name + ";\n" +
               "      for (auto const& message : messages) {\n" +
               "        deserialized_" + descriptor_name +
               ".push_back(Deserialize" + message_type_name + "(message));\n" +
               "      }\n"
               "      return deserialized_" + descriptor_name +
               ";\n    }(" + expr + ")";
      };
  field_serializer_fn_[descriptor] =
      [message_type_name](std::string const& expr) {
        return "Serialize" + message_type_name + "(" + expr + ")";
      };
}

void JournalProtoProcessor::ProcessOptionalInt32Field(
    FieldDescriptor const* descriptor) {
  field_type_[descriptor] = "int const*";

  field_arguments_fn_[descriptor] =
      [](std::string const& identifier) -> std::vector<std::string> {
        return {identifier + ".get()"};
      };
  field_indirect_member_get_fn_[descriptor] =
      [](std::string const& expr) {
        return "*" + expr;
      };
  field_optional_pointer_fn_[descriptor] =
      [this, descriptor](std::string const& condition,
                         std::string const& expr) {
        // Tricky.  We need a heap allocation to obtain a pointer to the value.
        return condition + " ? std::make_unique<int const>(" + expr +
               ") : nullptr";
      };
}

void JournalProtoProcessor::ProcessRequiredFixed64Field(
    FieldDescriptor const* descriptor) {
  FieldOptions const& options = descriptor->options();
  CHECK(options.HasExtension(serialization::pointer_to))
      << descriptor->full_name() << " is missing a (pointer_to) option";
  std::string const& pointer_to =
      options.GetExtension(serialization::pointer_to);
  field_type_[descriptor] = pointer_to + "*";

  if (options.HasExtension(serialization::is_consumed)) {
    CHECK(options.GetExtension(serialization::is_consumed))
        << descriptor->full_name() << " has incorrect (is_consumed) option";
    field_deleter_fn_[descriptor] =
        [](std::string const& expr) {
          return "  Delete(pointer_map, " + expr + ");\n";
        };
  }
  if (options.HasExtension(serialization::is_consumed_if)) {
    CHECK(!options.HasExtension(serialization::is_consumed))
        << descriptor->full_name()
        << " has incorrect )is_consumed) and (is_consumed_if) options";
    field_deleter_fn_[descriptor] =
        [options](std::string const& expr) {
          return "  if (" +
                 options.GetExtension(serialization::is_consumed_if) +
                 ") {\n    Delete(pointer_map, " + expr + ");\n  }\n";
        };
  }
  if (options.HasExtension(serialization::is_produced)) {
    CHECK(options.GetExtension(serialization::is_produced))
        << descriptor->full_name() << " has incorrect (is_produced) option";
    field_inserter_fn_[descriptor] =
        [](std::string const& expr1, std::string const& expr2) {
          return "  Insert(pointer_map, " + expr1 + ", " + expr2 + ");\n";
        };
  }
  if (options.HasExtension(serialization::is_produced_if)) {
    CHECK(!options.HasExtension(serialization::is_produced))
        << descriptor->full_name()
        << " has incorrect (is_produced) and (is_produced_if) options";
    field_inserter_fn_[descriptor] =
        [options](std::string const& expr1, std::string const& expr2) {
          return "  if (" +
                 options.GetExtension(serialization::is_produced_if) +
                 ") {\n    Insert(pointer_map, " + expr1 + ", " + expr2 +
                 ");\n  }\n";
        };
  }

  field_deserializer_fn_[descriptor] =
      [pointer_to](std::string const& expr) {
        return "DeserializePointer<" + pointer_to + "*>(*pointer_map, " + expr +
               ")";
      };
  field_serializer_fn_[descriptor] =
      [](std::string const& expr) {
        return "SerializePointer(" + expr + ")";
      };
}

void JournalProtoProcessor::ProcessRequiredMessageField(
    FieldDescriptor const* descriptor) {
  std::string const& message_type_name = descriptor->message_type()->name();
  field_type_[descriptor] = message_type_name;

  field_assignment_fn_[descriptor] =
      [this, descriptor](std::string const& identifier,
                         std::string const& expr) {
        return "  *" + identifier + "->mutable_" + descriptor->name() +
               "() = " + field_serializer_fn_[descriptor](expr) + ";\n";
      };
  field_deserializer_fn_[descriptor] =
      [message_type_name](std::string const& expr) {
        return "Deserialize" + message_type_name + "(" + expr + ")";
      };
  field_serializer_fn_[descriptor] =
      [message_type_name](std::string const& expr) {
        return "Serialize" + message_type_name + "(" + expr + ")";
      };
}

void JournalProtoProcessor::ProcessRequiredBoolField(
    FieldDescriptor const* descriptor) {
  field_type_[descriptor] = descriptor->cpp_type_name();
}

void JournalProtoProcessor::ProcessRequiredDoubleField(
    FieldDescriptor const* descriptor) {
  field_type_[descriptor] = descriptor->cpp_type_name();
}

void JournalProtoProcessor::ProcessRequiredInt32Field(
    FieldDescriptor const* descriptor) {
  field_type_[descriptor] = "int";
}

void JournalProtoProcessor::ProcessSingleStringField(
    FieldDescriptor const* descriptor) {
  field_type_[descriptor] = "char const*";
  FieldOptions const& options = descriptor->options();
  if (options.HasExtension(serialization::size)) {
    size_member_name_[descriptor] = options.GetExtension(serialization::size);

    field_arguments_fn_[descriptor] =
        [](std::string const& identifier) -> std::vector<std::string> {
          return {identifier + "->c_str()", identifier + "->size()"};
        };
    field_deserializer_fn_[descriptor] =
        [](std::string const& expr) {
          return "&" + expr;
        };
    field_indirect_member_get_fn_[descriptor] =
        [this, descriptor](std::string const& expr) {
          return "std::string(" + expr + ", " +
                 expr.substr(0, expr.find('.')) + "." +
                 size_member_name_[descriptor] + ")";
        };
  } else {
    field_deserializer_fn_[descriptor] =
        [](std::string const& expr) {
          return expr + ".c_str()";
        };
  }
}

void JournalProtoProcessor::ProcessOptionalField(
    FieldDescriptor const* descriptor) {
  field_optional_assignment_fn_[descriptor] =
      [](std::string const& expr, std::string const& stmt) {
        return "  if (" + expr + " != nullptr) {\n  " + stmt + "  }\n";
      };
  field_optional_pointer_fn_[descriptor] =
      [](std::string const& condition, std::string const& expr) {
        return condition + " ? " + expr + " : nullptr";
      };
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
  if (Contains(in_out_, descriptor)) {
    field_type_[descriptor] += "*";

    field_arguments_fn_[descriptor] =
        [](std::string const& identifier) -> std::vector<std::string> {
          return {"&" + identifier};
        };
    field_indirect_member_get_fn_[descriptor] =
        [](std::string const& expr) {
          return "*" + expr;
        };
  }
}

void JournalProtoProcessor::ProcessField(FieldDescriptor const* descriptor) {
  // Useful defaults for the lambdas, which ensure that they are set for all
  // fields.  They will be overwritten by actual processing as needed.
  field_arguments_fn_[descriptor] =
      [](std::string const& identifier) -> std::vector<std::string> {
        return {identifier};
      };
  field_assignment_fn_[descriptor] =
      [this, descriptor](std::string const& identifier,
                         std::string const& expr) {
        return "  " + identifier + "->set_" + descriptor->name() + "(" +
               field_serializer_fn_[descriptor](expr) + ");\n";
      };
  field_indirect_member_get_fn_[descriptor] =
      [](std::string const& expr) {
        return expr;
      };
  field_deserializer_fn_[descriptor] =
      [](std::string const& expr) {
        return expr;
      };
  field_optional_assignment_fn_[descriptor] =
      [](std::string const& expr, std::string const& stmt) {
        return stmt;
      };
  field_optional_pointer_fn_[descriptor] =
      [](std::string const& condition, std::string const& expr) {
        return expr;
      };
  field_serializer_fn_[descriptor] =
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

  // Generate slightly more compact code in the frequent case where the message
  // only has one field.
  std::string cpp_message_name = "message->mutable_" + ToLower(name) + "()";
  if (descriptor->field_count() > 1) {
    fill_body_[descriptor] = "  auto* const m = " + cpp_message_name + ";\n";
    cpp_message_name = "m";
  } else {
    fill_body_[descriptor].clear();
  }
  run_body_prolog_[descriptor] =
      "  auto const& " + ToLower(name) + " = message." +
      ToLower(name) + "();\n";
  run_arguments_[descriptor].clear();
  run_body_epilog_[descriptor].clear();

  nested_type_declaration_[descriptor] = "  struct " + name + " {\n";
  for (int i = 0; i < descriptor->field_count(); ++i) {
    FieldDescriptor const* field_descriptor = descriptor->field(i);
    std::string const& field_descriptor_name = field_descriptor->name();
    if (field_descriptors != nullptr) {
      field_descriptors->push_back(field_descriptor);
    }
    ProcessField(field_descriptor);

    std::string const fill_member_name =
        ToLower(name) + "." + field_descriptor_name;
    std::string const run_field_getter =
        ToLower(name) + "." + field_descriptor_name + "()";
    std::string const run_local_variable = field_descriptor_name;

    fill_body_[descriptor] +=
        field_optional_assignment_fn_[field_descriptor](
            fill_member_name,
            field_assignment_fn_[field_descriptor](
                cpp_message_name,
                field_indirect_member_get_fn_[field_descriptor](
                    fill_member_name)));
    std::vector<std::string> const field_arguments =
        field_arguments_fn_[field_descriptor](run_local_variable);
    std::copy(field_arguments.begin(), field_arguments.end(),
              std::back_inserter(run_arguments_[descriptor]));
    run_body_prolog_[descriptor] +=
        "  auto " + run_local_variable + " = " +
        field_optional_pointer_fn_[field_descriptor](
            ToLower(name) + ".has_" + field_descriptor_name + "()",
            field_deserializer_fn_[field_descriptor](run_field_getter)) +
        ";\n";
    if (Contains(field_deleter_fn_, field_descriptor)) {
      run_body_epilog_[descriptor] +=
          field_deleter_fn_[field_descriptor](run_field_getter);
    }
    if (Contains(field_inserter_fn_, field_descriptor)) {
      // The reference to |message| below avoids having to generate names for
      // the |out| fields (we wouldn't use them anywhere else).  This works
      // because we know that we insert only out parameters.
      run_body_epilog_[descriptor] +=
          field_inserter_fn_[field_descriptor](
              "message." + ToLower(name) + "()." + field_descriptor_name + "()",
              run_local_variable);
    }
    nested_type_declaration_[descriptor] += "    " +
                                    field_type_[field_descriptor] +
                                    " const " +
                                    field_descriptor_name + ";\n";

    // If this field has a size, generate it now.
    if (Contains(size_member_name_, field_descriptor)) {
      nested_type_declaration_[descriptor] +=
          "    int const " + size_member_name_[field_descriptor] + ";\n";
    }
  }
  nested_type_declaration_[descriptor] += "  };\n";
}

void JournalProtoProcessor::ProcessReturn(Descriptor const* descriptor) {
  CHECK_EQ(1, descriptor->field_count())
      << descriptor->full_name() << " must have exactly one field";
  FieldDescriptor const* field_descriptor = descriptor->field(0);
  CHECK_EQ(FieldDescriptor::LABEL_REQUIRED, field_descriptor->label())
      << descriptor->full_name() << " must be required";
  ProcessField(field_descriptor);
  fill_body_[descriptor] =
      field_assignment_fn_[field_descriptor]("message->mutable_return_()",
                                             "result");
  std::string const field_name =
      "message.return_()." + field_descriptor->name() + "()";
  if (Contains(field_inserter_fn_, field_descriptor)) {
    run_body_epilog_[descriptor] =
        field_inserter_fn_[field_descriptor](field_name, "result");
  } else {
    run_body_epilog_[descriptor] =
        "  CHECK(" + field_deserializer_fn_[field_descriptor](field_name) +
        " == result);\n";
  }
  nested_type_declaration_[descriptor] =
      "  using Return = " + field_type_[field_descriptor] + ";\n";
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
        in_out_.insert(field_descriptors[i]);
        in_out_.insert(field_descriptors[i + 1]);
      }
    }
  }

  // The second pass that produces the actual output.
  std::string cpp_run_arguments;
  std::string cpp_run_prolog;
  std::string cpp_run_epilog;
  toplevel_type_declaration_[descriptor] = "struct " + name + " {\n";
  for (int i = 0; i < descriptor->nested_type_count(); ++i) {
    Descriptor const* nested_descriptor = descriptor->nested_type(i);
    const std::string& nested_name = nested_descriptor->name();
    if (nested_name == kIn) {
      ProcessInOut(nested_descriptor, /*field_descriptors=*/nullptr);
      functions_implementation_[descriptor] +=
          "void " + name + "::Fill(In const& in, "
          "not_null<Message*> const message) {\n" +
          fill_body_[nested_descriptor] +
          "}\n\n";
      cpp_run_arguments += Join(run_arguments_[nested_descriptor]);
      cpp_run_prolog += run_body_prolog_[nested_descriptor];
    } else if (nested_name == kOut) {
      ProcessInOut(nested_descriptor, /*field_descriptors=*/nullptr);
      functions_implementation_[descriptor] +=
          "void " + name + "::Fill(Out const& out, "
          "not_null<Message*> const message) {\n" +
          fill_body_[nested_descriptor] +
          "}\n\n";
    } else if (nested_name == kReturn) {
      ProcessReturn(nested_descriptor);
      functions_implementation_[descriptor] +=
          "void " + name + "::Fill("
          "Return const& result, "
          "not_null<Message*> const message) {\n" +
          fill_body_[nested_descriptor] +
          "}\n\n";
    }
    cpp_run_epilog += run_body_epilog_[nested_descriptor];
    toplevel_type_declaration_[descriptor] +=
        nested_type_declaration_[nested_descriptor];
  }
  if (has_in || has_out || has_return) {
    toplevel_type_declaration_[descriptor] += "\n";
  }
  toplevel_type_declaration_[descriptor] +=
      "  using Message = serialization::" + name + ";\n";
  if (has_in) {
    toplevel_type_declaration_[descriptor] +=
        "  static void Fill(In const& in, "
        "not_null<Message*> const message);\n";
  }
  if (has_out) {
    toplevel_type_declaration_[descriptor] +=
        "  static void Fill(Out const& out, "
        "not_null<Message*> const message);\n";
  }
  if (has_return) {
    toplevel_type_declaration_[descriptor] +=
        "  static void Fill("
        "Return const& result, "
        "not_null<Message*> const message);\n";
  }
  toplevel_type_declaration_[descriptor] +=
      "  static void Run("
      "Message const& message,\n"
      "                  not_null<"
      "Player::PointerMap*> const pointer_map);"
      "\n";
  toplevel_type_declaration_[descriptor] += "};\n\n";

  // Must come after the Fill methods for comparison with manual code.
  functions_implementation_[descriptor] +=
      "void " + name + "::Run(Message const& message, "
      "not_null<Player::PointerMap*> const pointer_map) {\n" +
      cpp_run_prolog;
  if (has_return) {
    functions_implementation_[descriptor] += "  auto const result = ";
  } else {
    functions_implementation_[descriptor] += "  ";
  }
  functions_implementation_[descriptor] +=
      "ksp_plugin::principia__" + name + "(" + cpp_run_arguments + ");\n";
  functions_implementation_[descriptor] += cpp_run_epilog + "}\n\n";
}

}  // namespace tools
}  // namespace principia
