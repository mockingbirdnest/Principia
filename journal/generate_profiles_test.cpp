//TODO(phl): This should totally no be a test.

#include <vector>

#include "glog/logging.h"
#include "google/protobuf/descriptor.h"
#include "gtest/gtest.h"
#include "serialization/journal.pb.h"

namespace principia {

using ::google::protobuf::Descriptor;
using ::google::protobuf::FieldDescriptor;
using ::google::protobuf::FieldOptions;
using ::google::protobuf::FileDescriptor;

namespace journal {

namespace {
char const kIn[] = "In";
char const kReturn[] = "Return";
char const kOut[] = "Out";
}  // namespace

class GeneratorTest : public testing::Test {
 protected:
  ~GeneratorTest() {
    for (auto const& pair : cpp_method_type_) {
      std::cout<<pair.second;
    }
  }

   void ProcessOptionalInt32Field(FieldDescriptor const* descriptor) {
     cpp_field_type_[descriptor] = "int const*";
   }

   void ProcessRequiredFixed64Field(FieldDescriptor const* descriptor) {
    FieldOptions const& options = descriptor->options();
    CHECK(options.HasExtension(serialization::pointer_to))
        << descriptor->full_name() << " is missing a pointer_to option";
    std::string const& pointer_to =
        options.GetExtension(serialization::pointer_to);
    cpp_field_type_[descriptor] = pointer_to + "*";
  }

  void ProcessRequiredMessageField(FieldDescriptor const* descriptor) {
    cpp_field_type_[descriptor] = descriptor->message_type()->name();
  }

  void ProcessRequiredBoolField(FieldDescriptor const* descriptor) {
    cpp_field_type_[descriptor] = descriptor->cpp_type_name();
  }

  void ProcessRequiredDoubleField(FieldDescriptor const* descriptor) {
    cpp_field_type_[descriptor] = descriptor->cpp_type_name();
  }

  void ProcessRequiredInt32Field(FieldDescriptor const* descriptor) {
    cpp_field_type_[descriptor] = "int";
  }

  void ProcessSingleStringField(FieldDescriptor const* descriptor) {
    cpp_field_type_[descriptor] = "char const*";
  }

  void ProcessOptionalField(FieldDescriptor const* descriptor) {
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

  void ProcessRepeatedField(FieldDescriptor const* descriptor) {}

  void ProcessRequiredField(FieldDescriptor const* descriptor) {
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

  void ProcessField(FieldDescriptor const* descriptor) {
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

  void ProcessInOut(Descriptor const* descriptor,
                    std::vector<FieldDescriptor const*>* field_descriptors) {
    cpp_nested_type_[descriptor] = "  struct " + descriptor->name() + " {\n";
    for (int i = 0; i < descriptor->field_count(); ++i) {
      FieldDescriptor const* field_descriptor = descriptor->field(i);
      if (field_descriptors != nullptr) {
        field_descriptors->push_back(field_descriptor);
      }
      ProcessField(field_descriptor);
      cpp_nested_type_[descriptor] += "    " +
                                      cpp_field_type_[field_descriptor] +
                                      " const " +
                                      field_descriptor->name() + ";\n";
    }
    cpp_nested_type_[descriptor] += "  };\n";
  }

  void ProcessReturn(Descriptor const* descriptor) {
    CHECK_EQ(1, descriptor->field_count())
        << descriptor->full_name() << " must have exactly one field";
    FieldDescriptor const* field_descriptor = descriptor->field(0);
    ProcessField(field_descriptor);
    cpp_nested_type_[descriptor] =
        "  using Return = " + cpp_field_type_[field_descriptor] + ";\n";
  }

  void ProcessMethodExtension(Descriptor const* descriptor) {
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
    if (has_in || has_out) {
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
        std::string("  using Message = serialization::") +
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

 private:
  std::set<FieldDescriptor const*> in_out_field_;
  std::map<Descriptor const*, std::string> cpp_method_type_;
  std::map<Descriptor const*, std::string> cpp_nested_type_;
  std::map<FieldDescriptor const*, std::string> cpp_field_type_;
};

TEST_F(GeneratorTest, JustDoIt) {
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

}  // namespace journal
}  // namespace principia
