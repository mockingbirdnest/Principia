//TODO(phl): This should totally no be a test.

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
  void ProcessOptionalField(const FieldDescriptor* descriptor) {
    switch (descriptor->type()) {
    case FieldDescriptor::TYPE_INT32:
      break;
    case FieldDescriptor::TYPE_STRING:
      break;
    default:
      LOG(FATAL) << descriptor->full_name() << " has unexpected type "
        << descriptor->type_name();
    }
  }

  void ProcessRepeatedField(const FieldDescriptor* descriptor) {}

  void ProcessRequiredFixed64Field(const FieldDescriptor* descriptor) {
    FieldOptions const& options = descriptor->options();
    CHECK(options.HasExtension(serialization::pointer_to))
        << descriptor->full_name() << " is missing a pointer_to option";
    std::string const& pointer_to =
        options.GetExtension(serialization::pointer_to);
  }

  void ProcessRequiredField(const FieldDescriptor* descriptor) {
    switch (descriptor->type()) {
      case FieldDescriptor::TYPE_BOOL:
        break;
      case FieldDescriptor::TYPE_DOUBLE:
        break;
      case FieldDescriptor::TYPE_FIXED64:
        ProcessRequiredFixed64Field(descriptor);
        break;
      case FieldDescriptor::TYPE_INT32:
        break;
      case FieldDescriptor::TYPE_MESSAGE:
        break;
      case FieldDescriptor::TYPE_STRING:
        break;
      default:
        LOG(FATAL) << descriptor->full_name() << " has unexpected type "
                   << descriptor->type_name();
    }
  }

  void ProcessField(const FieldDescriptor* descriptor) {
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

  void ProcessIn(const Descriptor* descriptor) {
    for (int i = 0; i < descriptor->field_count(); ++i) {
      const FieldDescriptor* field_descriptor = descriptor->field(i);
      ProcessField(field_descriptor);
    }
  }

  void ProcessOut(const Descriptor* descriptor) {
    for (int i = 0; i < descriptor->field_count(); ++i) {
      const FieldDescriptor* field_descriptor = descriptor->field(i);
      ProcessField(field_descriptor);
    }
  }

  void ProcessReturn(const Descriptor* descriptor) {
    CHECK_EQ(1, descriptor->field_count())
        << descriptor->full_name() << " must have exactly one field";
    ProcessField(descriptor->field(0));
  }

  void ProcessMethodExtension(const Descriptor* descriptor) {
    for (int i = 0; i < descriptor->nested_type_count(); ++i) {
      const Descriptor* nested_descriptor = descriptor->nested_type(i);
      const std::string& nested_name = nested_descriptor->name();
      if (nested_name == kIn) {
        ProcessIn(nested_descriptor);
      } else if (nested_name == kOut) {
        ProcessOut(nested_descriptor);
      } else if (nested_name == kReturn) {
        ProcessReturn(nested_descriptor);
      } else {
        LOG(FATAL) << "Unexpected nested message "
                   << nested_descriptor->full_name();
      }
    }
  }
};

TEST_F(GeneratorTest, JustDoIt) {
  // Get the file containing |Method|.
  const Descriptor* method_descriptor = serialization::Method::descriptor();
  const FileDescriptor* file_descriptor = method_descriptor->file();

  // Process all the messages in that file.
  for (int i = 0; i < file_descriptor->message_type_count(); ++i) {
    const Descriptor* message_descriptor = file_descriptor->message_type(i);
    switch (message_descriptor->extension_count()) {
      case 0: {
        // An auxiliary message that is not an extension.  Nothing to do.
        break;
      }
      case 1: {
        // An extension.  Check that it extends |Method|.
        const FieldDescriptor* extension = message_descriptor->extension(0);
        CHECK(extension->is_extension());
        const Descriptor* containing_type = extension->containing_type();
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
