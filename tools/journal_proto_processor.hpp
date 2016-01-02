#pragma once

#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "google/protobuf/descriptor.h"

namespace principia {

using ::google::protobuf::Descriptor;
using ::google::protobuf::FieldDescriptor;
using ::google::protobuf::FieldOptions;
using ::google::protobuf::FileDescriptor;

namespace tools {

class JournalProtoProcessor {
 public:
  void ProcessMessages();

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
  void ProcessRequiredUint32Field(FieldDescriptor const* descriptor);

  void ProcessSingleStringField(FieldDescriptor const* descriptor);

  void ProcessOptionalField(FieldDescriptor const* descriptor);
  void ProcessRepeatedField(FieldDescriptor const* descriptor);
  void ProcessRequiredField(FieldDescriptor const* descriptor);

  void ProcessField(FieldDescriptor const* descriptor);

  void ProcessInOut(Descriptor const* descriptor,
                    std::vector<FieldDescriptor const*>* field_descriptors);
  void ProcessReturn(Descriptor const* descriptor);

  void ProcessInterchangeMessage(Descriptor const* descriptor);
  void ProcessMethodExtension(Descriptor const* descriptor);

  // As the recursive methods above traverse the protocol buffer type
  // declarations, they enter in the following maps (and set) various pieces of
  // information to help in generating C++ code.  For the simplest use cases
  // (mostly, the generation of the .hpp file), the values are merely one or
  // several strings for C++ code snippets.  For more complex use cases (the
  // generation of the .cpp file) the values are lambdas which transform one or
  // two C++ code snippets by wrapping them in a more complex structure.

  // The fields that are in-out, i.e. for which fields of the same name exist in
  // both the In and the Out messages.  Note that both fields are present in
  // this set.  Those fields are transmitted through the interface with an extra
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
  // pointer.  The lambda calls |field_serializer_fn_| to serialize expressions
  // as necessary; thus, |expr| must *not* be serialized.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& identifier,
                                     std::string const& expr)>>
      field_assignment_fn_;

  // For fields that have an (is_consumed) or (is_consumed_if) option, a lambda
  // producing a statement to call Delete() to remove the appropriate entry from
  // the pointer_map.  |expr| is a uint64 expression for the entry to be
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

  // For all fields, a lambda that takes an expression for a struct member and
  // returns an expression that dereferences it if the field uses a level of
  // indirection (e.g., is optional or in-out).
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
      field_indirect_member_get_fn_;

  // For fields that have an (is_produced) or (is_produced_if) option, a lambda
  // producing a statement to call Insert() to enter the appropriate entry into
  // the pointer_map.  |expr1| is an uint64 expression for the serialized value
  // of the pointer (typically something like |message.in().bar()|), |expr2| is
  // a pointer expression for the current value of the pointer (typically the
  // name of a local variable).  No data for other fields.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr1,
                                     std::string const& expr2)>>
      field_inserter_fn_;

  // For all fields, a lambda that takes a pointer expression and a statement
  // generated by |field_assignment_fn_|.  If the field is optional, returns an
  // if statement that only executes |stmt| if |expr| in nonnull.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr,
                                     std::string const& stmt)>>
      field_optional_assignment_fn_;

  // For all fields, a lambda that takes a condition to take for the presence
  // of an optional field (typically something like |message.in().has_bar()|)
  // and a deserialized expression for reading the field (typically the result
  // of |field_deserializer_fn_|) and returns a conditional expression for
  // either a pointer to the deserialized value or nullptr.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& condition,
                                     std::string const& expr)>>
      field_optional_pointer_fn_;

  // For all fields, a lambda that takes an expression for reading a local
  // variable (possibly with dereferencing) and returns a protocol buffer
  // expression suitable for assigning to some field either using set_bar() or
  // mutable_bar() (typically the result is a call to some Serialize function).
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
      field_serializer_fn_;

  // The C++ type for a field, suitable for use in a member or parameter
  // declaration, in a typedef, etc.
  std::map<FieldDescriptor const*, std::string> field_type_;

  //TODO(phl):comment
  std::map<Descriptor const*, std::string> deserialize_body_;
  std::map<Descriptor const*, std::string> serialize_body_;

  // The entire sequence of statements for the body of a Fill function.  The key
  // is a descriptor for an In or Out message.
  std::map<Descriptor const*, std::string> fill_body_;

  // A code snippet that goes before the call to the interface in the
  // body of the Run function.  The key is a descriptor for an In or Out
  // message.  Produced but not used for Out messages.
  std::map<Descriptor const*, std::string> run_body_prolog_;

  // A list of code snippets for arguments to be passed to the interface in the
  // body of the Run function.
  std::map<Descriptor const*, std::vector<std::string>> run_arguments_;

  // A code snippet that goes after the call to the interface in the
  // body of the Run function.  The key is a descriptor for an In or Out
  // message.
  std::map<Descriptor const*, std::string> run_body_epilog_;

  // A code snippet for the implementation of the Fill and Run functions.  The
  // key is a descriptor for a method message.
  std::map<Descriptor const*, std::string> functions_implementation_;

  // A code snippet for the declaration of the top-level struct for a method.
  // The key is a descriptor for a method message.
  std::map<Descriptor const*, std::string> toplevel_type_declaration_;

  // A code snippet for the declaration of a nested In or Out struct or a Return
  // typedef.  The key is a descriptor for an In, Out or Return message.
  std::map<Descriptor const*, std::string> nested_type_declaration_;
};

}  // namespace tools
}  // namespace principia
