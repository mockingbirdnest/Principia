#pragma once

#include <functional>
#include <map>

#include "google/protobuf/descriptor.h"

namespace principia {

using ::google::protobuf::Descriptor;
using ::google::protobuf::FieldDescriptor;
using ::google::protobuf::FieldOptions;
using ::google::protobuf::FileDescriptor;

namespace tools {

class JournalProtoProcessor {
 public:
  void ProcessMethods();

  std::vector<std::string> GetCppMethodImplementations() const;
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

  void ProcessField(FieldDescriptor const* descriptor);

  void ProcessInOut(Descriptor const* descriptor,
                    std::vector<FieldDescriptor const*>* field_descriptors);
  void ProcessReturn(Descriptor const* descriptor);

  void ProcessMethodExtension(Descriptor const* descriptor);

  // As the recursive methods above traverse the protocol buffer type
  // declarations, they enter in the following maps (and set) various pieces of
  // information to help in generating C++ code.  For the simplest use cases
  // (mostly, the generation of the .hpp file), the values are merely one or
  // several strings for C++ code snippets.  For more complex use cases (the
  // generation of the .cpp file) the values are lambdas which transform one or
  // two C++ code snippets by wrapping them in a more complex structure.

  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& name,
                                     std::string const& expr)>>
           field_copy_wrapper_;
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
           field_deleter_wrapper_;
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
           field_deserializer_wrapper_;
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr1,
                                     std::string const& expr2)>>
           field_inserter_wrapper_;
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
           field_serializer_wrapper_;
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
           indirect_field_wrapper_;
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& condition,
                                     std::string const& expr)>>
           optional_field_get_wrapper_;
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr,
                                     std::string const& stmt)>>
           optional_field_wrapper_;
  std::map<FieldDescriptor const*,
           std::function<std::vector<std::string>(std::string const& name)>>
           field_arguments_wrapper_;

  std::map<FieldDescriptor const*, std::string> size_field_name_;
  std::set<FieldDescriptor const*> in_out_field_;
  std::map<FieldDescriptor const*, std::string> cpp_field_type_;

  // The entire sequence of statements for the body of a Fill function.  The key
  // is a descriptor for an In or Out message.
  std::map<Descriptor const*, std::string> cpp_fill_body_;

  // A code snippet that goes before the call to the interface in the
  // body of the Run function.  The key is a descriptor for an In or Out
  // message.  Produced but not used for Out messages.
  std::map<Descriptor const*, std::string> cpp_run_prolog_;

  // A list of code snippets for arguments to be passed to the interface in the
  // body of the Run function.
  std::map<Descriptor const*, std::vector<std::string>> cpp_run_arguments_;

  // A code snippet that goes after the call to the interface in the
  // body of the Run function.  The key is a descriptor for an In or Out
  // message.
  std::map<Descriptor const*, std::string> cpp_run_epilog_;
  std::map<Descriptor const*, std::string> cpp_method_impl_;
  std::map<Descriptor const*, std::string> cpp_method_type_;
  std::map<Descriptor const*, std::string> cpp_nested_type_;
};

}  // namespace tools
}  // namespace principia
