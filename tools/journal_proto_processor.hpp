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

  // The fields that are in-out, i.e. for which fields of the same name exist in
  // both the In and the Out message.  Note that both fields are present in this
  // set.  Those fields are transmitted through the interface with an extra
  // level of indirection.
  std::set<FieldDescriptor const*> in_out_;

  // For fields that have a (size) option, the name of the size member variable
  // in the In or Out struct.  Special processing is required when filling those
  // fields from the struct members.  No data for other fields.
  std::map<FieldDescriptor const*, std::string> size_member_name_;

  // For all fields, a lambda that takes the name of a local variable containing
  // data extracted (and deserialized) from the field and returns a list of
  // expressions to be passed to the interface.  Deals with passing address and
  // size for fields that have a size member, and with passing by reference for
  // fields that are in-out or optional.
  std::map<FieldDescriptor const*,
           std::function<std::vector<std::string>(
                             std::string const& identifier)>>
      field_arguments_fn_;

  // For all fields, a lambda that takes a serialized expression |expr| and a
  // protocol buffer object |identifier| and returns a statement to assign
  // |expr| to the proper field of |identifier|.  |identifier| must be a
  // pointer.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& name,
                                     std::string const& expr)>>
      field_assignment_fn_;

  // For fields that have a (is_deleted) or (is_deleted_if) option, a lambda
  // producing a statement to call Delete() to remove the appropriate entry from
  // the pointer_map.  |expr| is an uint64 expression for the entry to be
  // removed (typically something like |message.in().bar()|).  No data for other
  // fields.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
      field_deleter_fn_;

  // For all fields, a lambda that takes an expression for reading a protobuf
  // field (typically something like |message.in().bar()|) and returns an
  // expression for the deserialized form of |expr| suitable for storing in a
  // local variable (typically a call to some Deserialize function, but other
  // transformations are possible).
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
      field_deserializer_fn_;

  // For fields that have a (is_inserted) or (is_inserted_if) option, a lambda
  // producing a statement to call Insert() to enter the appropriate entry into
  // the pointer_map.  |expr1| is an uint64 expression for the serialized value
  // of the pointer (typically something like |message.in().bar()|), |expr2| is
  // a pointer expression for the current value of the pointer (typically the
  // name of a local variable).  No data for other fields.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr1,
                                     std::string const& expr2)>>
      field_inserter_fn_;

  // For all fields, a lambda that takes an expression for reading a local
  // variable (possibly with dereferencing) and returns a protocol buffer
  // expression suitable for assigning to some field either using set_bar() or
  // mutable_bar() (typically the result is a call to some Serialize function).
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
           field_serializer_fn_;

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
