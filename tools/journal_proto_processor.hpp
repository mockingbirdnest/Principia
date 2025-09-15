#pragma once

#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "google/protobuf/descriptor.h"

namespace principia {
namespace tools {
namespace _journal_proto_processor {
namespace internal {

using ::google::protobuf::Descriptor;
using ::google::protobuf::FieldDescriptor;
using ::google::protobuf::FieldOptions;
using ::google::protobuf::FileDescriptor;

class JournalProtoProcessor final {
 public:
  void ProcessMessages();

  // ksp_plugin_adapter/interface.generated.cs
  std::vector<std::string> GetCsInterfaceMethodDeclarations() const;
  std::vector<std::string> GetCsInterfaceSymbolDeclarations() const;
  std::vector<std::string> GetCsInterchangeTypeDeclarations() const;

  // ksp_plugin_adapter/marshalers.generated.cs
  std::vector<std::string> GetCsCustomMarshalerClasses() const;

  // ksp_plugin/interface.generated.h
  std::vector<std::string> GetCxxInterfaceMethodDeclarations() const;
  std::vector<std::string> GetCxxInterchangeTypeDeclarations() const;

  // journal/profiles.generated.{h,cc}
  std::vector<std::string> GetCxxInterchangeImplementations() const;
  std::vector<std::string> GetCxxMethodImplementations() const;
  std::vector<std::string> GetCxxMethodTypes() const;

  // journal/player.generated.cc
  std::vector<std::string> GetCxxPlayStatements() const;

 private:
  void ProcessRepeatedNonStringField(FieldDescriptor const* descriptor,
                                     std::string const& cs_boxed_type,
                                     std::string const& cs_unboxed_type,
                                     std::string const& cxx_type);
  void ProcessRepeatedScalarField(FieldDescriptor const* descriptor,
                                  std::string const& cxx_type);
  void ProcessRepeatedDoubleField(FieldDescriptor const* descriptor);
  void ProcessRepeatedInt32Field(FieldDescriptor const* descriptor);
  void ProcessRepeatedInt64Field(FieldDescriptor const* descriptor);
  void ProcessRepeatedUint32Field(FieldDescriptor const* descriptor);
  void ProcessRepeatedMessageField(FieldDescriptor const* descriptor);
  void ProcessRepeatedStringField(FieldDescriptor const* descriptor);

  void ProcessOptionalNonStringField(FieldDescriptor const* descriptor,
                                     std::string const& cs_boxed_type,
                                     std::string const& cs_unboxed_type,
                                     std::string const& cxx_type);
  void ProcessOptionalScalarField(FieldDescriptor const* descriptor,
                                  std::string const& cxx_type);
  void ProcessOptionalDoubleField(FieldDescriptor const* descriptor);
  void ProcessOptionalInt32Field(FieldDescriptor const* descriptor);
  void ProcessOptionalInt64Field(FieldDescriptor const* descriptor);
  void ProcessOptionalUint32Field(FieldDescriptor const* descriptor);
  void ProcessOptionalMessageField(FieldDescriptor const* descriptor);

  void ProcessRequiredFixed32Field(FieldDescriptor const* descriptor);
  void ProcessRequiredFixed64Field(FieldDescriptor const* descriptor);
  void ProcessRequiredMessageField(FieldDescriptor const* descriptor);
  void ProcessRequiredBoolField(FieldDescriptor const* descriptor);
  void ProcessRequiredBytesField(FieldDescriptor const* descriptor);
  void ProcessRequiredDoubleField(FieldDescriptor const* descriptor);
  void ProcessRequiredInt32Field(FieldDescriptor const* descriptor);
  void ProcessRequiredInt64Field(FieldDescriptor const* descriptor);
  void ProcessRequiredUint32Field(FieldDescriptor const* descriptor);

  void ProcessSingleMessageField(FieldDescriptor const* descriptor);
  void ProcessSingleStringField(FieldDescriptor const* descriptor);

  void ProcessOptionalField(FieldDescriptor const* descriptor);
  void ProcessRepeatedField(FieldDescriptor const* descriptor);
  void ProcessRequiredField(FieldDescriptor const* descriptor);

  void ProcessField(FieldDescriptor const* descriptor);

  void ProcessAddressOfSizeOf(Descriptor const* descriptor);

  void ProcessInOut(Descriptor const* descriptor,
                    std::vector<FieldDescriptor const*>* field_descriptors);
  void ProcessReturn(Descriptor const* descriptor);

  void ProcessInterchangeMessage(Descriptor const* descriptor);
  void ProcessMethodExtension(Descriptor const* descriptor);

  bool HasMarshaler(FieldDescriptor const* descriptor) const;
  std::string MarshalAs(FieldDescriptor const* descriptor) const;

  // As the recursive methods above traverse the protocol buffer type
  // declarations, they enter in the following maps (and set) various pieces of
  // information to help in generating C++ and C# code.  For the simplest use
  // cases (mostly, the generation of declarations), the values are merely one
  // or several strings for C++ or C# code snippets.  For more complex use cases
  // (the generation of implementation code) the values are lambdas which
  // transform one or two code snippets by wrapping them in a more complex
  // structure.
  // We use cxx to designate C++ code and cs to designate C# code.

  // The fields that are in.  Note that the out fields present in `in_out_` are
  // not in `in_`.
  std::set<FieldDescriptor const*> in_;

  // The fields that are in-out, i.e. for which fields of the same name exist in
  // both the In and the Out messages.  Note that both fields are present in
  // this set.  Those fields are transmitted through the interface with an extra
  // level of indirection.
  std::set<FieldDescriptor const*> in_out_;

  // The fields that are out.  Those fields are transmitted through the
  // interface with an extra level of indirection.  Note that the in fields
  // present in `in_out_` are not in `out_`.
  std::set<FieldDescriptor const*> out_;

  // The fields that are returned.
  std::set<FieldDescriptor const*> return_;

  // The fields that are part of interchange messages.
  std::set<FieldDescriptor const*> interchange_;

  // For a field that has an (address_of) option, field_cxx_address_of_ has that
  // field as a key and the field designated by the option as its value.
  // field_cxx_address_ is the inverse map.
  std::map<FieldDescriptor const*, FieldDescriptor const*>
      field_cxx_address_of_;
  std::map<FieldDescriptor const*, FieldDescriptor const*> field_cxx_address_;

  // For a field that has a (size_of) option, field_cxx_size_of_ has that field
  // as a key and the field designated by the option as its value.
  // field_cxx_size_ is the inverse map.
  std::map<FieldDescriptor const*, FieldDescriptor const*> field_cxx_size_of_;
  std::map<FieldDescriptor const*, FieldDescriptor const*> field_cxx_size_;

  // For all fields, a lambda that takes the name of a local variable containing
  // data extracted (and deserialized) from the field and returns a list of
  // expressions to be passed to the interface.  Deals with passing by reference
  // for fields that are in-out or optional.
  std::map<FieldDescriptor const*,
           std::function<std::vector<std::string>(
                             std::string const& identifier)>>
      field_cxx_arguments_fn_;

  // For all fields, a lambda that takes a serialized expression `expr` and a
  // protocol buffer denoted by `prefix` and returns a statement to assign
  // `expr` to the proper field of `prefix`.  `prefix` must be suitable as a
  // prefix of a call, i.e., it must be a pointer followed by "->" or a
  // reference followed by ".".  The lambda calls `field_cxx_serializer_fn_` to
  // serialize expressions as necessary; thus, `expr` must *not* be serialized.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& prefix,
                                     std::string const& expr)>>
      field_cxx_assignment_fn_;

  // For fields that have an (is_consumed) or (is_consumed_if) option, a lambda
  // producing a statement to call Delete() to remove the appropriate entry from
  // the pointer_map.  `expr` is a uint64 expression for the entry to be
  // removed (typically something like `message.in().bar()`).  No data for other
  // fields.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
      field_cxx_deleter_fn_;

  // For all fields, a lambda that takes an expression for reading a protobuf
  // field (typically something like `message.in().bar()`) and returns an
  // expression for the deserialized form of `expr` suitable for storing in a
  // local variable (typically a call to some Deserialize function, but other
  // transformations are possible).  Deals with arrays of pointers for repeated
  // fields.
  // Note that for a field that is designated by a (size_of) field, the
  // (size_of) field is passed to this lambda instead of the field itself.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
      field_cxx_deserializer_fn_;

  // For all fields, a lambda that takes an expression for a struct member and
  // returns an expression that dereferences it if the field uses a level of
  // indirection (e.g., is optional or in-out).
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
      field_cxx_indirect_member_get_fn_;

  // For fields that have an (is_produced) or (is_produced_if) option, a lambda
  // producing a statement to call Insert() to enter the appropriate entry into
  // the pointer_map.  `expr1` is an uint64 expression for the serialized value
  // of the pointer (typically something like `message.in().bar()`), `expr2` is
  // a pointer expression for the current value of the pointer (typically the
  // name of a local variable).  No data for other fields.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr1,
                                     std::string const& expr2)>>
      field_cxx_inserter_fn_;

  // For all fields, a lambda that takes a C# parameter type as stored in
  // `field_cs_type_`, and adds a mode to it.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& type)>>
      field_cs_mode_fn_;
  // For all fields, a lambda that takes a C++ parameter or member type as
  // stored in `field_cxx_type_`, and adds a mode to it.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& type)>>
      field_cxx_mode_fn_;

  // For all fields, a lambda that takes a C# parameter type as stored in
  // `field_cs_type_`, and adds `this` to it if `is_subject`; the result can be
  // used to declare an parameter in an extension method.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& type)>>
      field_cs_extension_method_parameter_fn_;

  // For all fields, a lambda that takes a pointer expression and a statement
  // generated by `field_cxx_assignment_fn_`.  If the field is optional, returns
  // an if statement that only executes `stmt` if `expr` in nonnull.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr,
                                     std::string const& stmt)>>
      field_cxx_optional_assignment_fn_;

  // For all fields, a lambda that takes a condition to check for the presence
  // of an optional field (typically something like `message.in().has_bar()`)
  // and a deserialized expression for reading the field (typically the result
  // of `field_cxx_deserializer_fn_`) and returns a conditional expression for
  // either a pointer to the deserialized value or nullptr.
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& condition,
                                     std::string const& expr)>>
      field_cxx_optional_pointer_fn_;

  // For all fields, a lambda that takes an expression for reading a local
  // variable (possibly with dereferencing) and returns a protocol buffer
  // expression suitable for assigning to some field either using set_bar() or
  // mutable_bar() (typically the result is a call to some Serialize function).
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
      field_cxx_serializer_fn_;

  // For fields that are implemented using private members, a lambda that takes
  // the names of the members and returns the C# implementation of the getter
  // and setter.
  std::map<FieldDescriptor const*,
           std::function<
              std::string(std::vector<std::string> const& identifiers)>>
      field_cs_private_getter_fn_;
  std::map<FieldDescriptor const*,
           std::function<
              std::string(std::vector<std::string> const& identifiers)>>
      field_cs_private_setter_fn_;

  // The C# class for marshalling a field.  This may be either a custom
  // marshaler (derived from ICustomMarshaler) or one predefined by .Net.  At
  // most one is set for a given field.
  std::map<FieldDescriptor const*, std::string> field_cs_custom_marshaler_;
  std::map<FieldDescriptor const*, std::string> field_cs_predefined_marshaler_;

  // The fields that must be marshaled by simply copying their fields.  This is
  // useful for classes-within-classes when we don't need a level of
  // indirection.
  std::set<FieldDescriptor const*> field_cs_marshal_by_copy_;

  // The C# type for a field, suitable for use in a private member when the
  // actual data cannot be exposed directly (think bool).
  std::map<FieldDescriptor const*, std::string> field_cs_private_type_;

  // For fields that require deserialization storage, the C++ type and name of
  // the storage object.
  std::map<FieldDescriptor const*, std::string>
      field_cxx_deserialization_storage_name_;
  std::map<FieldDescriptor const*, std::string>
      field_cxx_deserialization_storage_type_;

  // The C#/C++ type for a field, suitable for use in a member or parameter
  // declaration, in a typedef, etc.
  std::map<FieldDescriptor const*, std::string> field_cs_type_;
  std::map<FieldDescriptor const*, std::string> field_cxx_type_;

  // The C#/C++ declaration of an interface method corresponding to a method
  // message.  The key is a descriptor for a method message.
  std::map<Descriptor const*, std::string> cs_interface_method_declaration_;
  std::map<Descriptor const*, std::string> cxx_interface_method_declaration_;

  // A partial declaration of the C# class Interface.Symbols containing the
  // declaration of the delegate corresponding to a method message.  The key is
  // a descriptor for a method message.
  std::map<Descriptor const*, std::string> cs_interface_symbol_declaration_;

  // A list of C#/C++ parameters for an interface method.  The key is a
  // descriptor for an In or Out message.  Produced but not used for Out
  // messages.
  std::map<Descriptor const*, std::vector<std::string>>
      cs_interface_parameters_;
  std::map<Descriptor const*, std::vector<std::string>>
      cxx_interface_parameters_;
  // Same as `cs_interface_parameters_`, but with `MarshalAs` attributes.
  std::map<Descriptor const*, std::vector<std::string>>
      cs_interface_marshaled_parameters_;
  // A list of arguments for a call to the interface method from C#, including
  // modes `out` or `ref`, with the same identifiers as the parameter names.
  std::map<Descriptor const*, std::vector<std::string>> cs_interface_arguments_;

  // The C#/C++ return type of an interface method.  The key is a descriptor for
  // a Return message.
  std::map<Descriptor const*, std::string> cs_interface_return_type_;
  std::map<Descriptor const*, std::string> cxx_interface_return_type_;

  // The C# attribute for marshalling the return value of an interface method.
  // The key is a descriptor for a Return message.
  std::map<Descriptor const*, std::string> cs_interface_return_marshal_;

  // The interchange messages that are represented by a class (as opposed to a
  // struct) in the C# code.
  std::set<Descriptor const*> cs_interchange_classes_;

  // The C#/C++ definition of a type corresponding to an interchange message.
  // The key is a descriptor for an interchange message.
  std::map<Descriptor const*, std::string> cs_interchange_type_declaration_;
  std::map<Descriptor const*, std::string> cxx_interchange_type_declaration_;

  // The full name of the C# class that implements a custom marshaler for an
  // interchange message.
  std::map<Descriptor const*, std::string> cs_custom_marshaler_name_;

  // The definition of the C# class that implements a custom marshaler for an
  // interchange message.
  std::map<Descriptor const*, std::string> cs_custom_marshaler_class_;

  // The C# declarations of fields in the Representation struct of a custom
  // marshaler.
  std::map<Descriptor const*, std::string> cs_representation_type_declaration_;

  // The C# statements in the CleanUpNative, MarshalManagedToNative and
  // MarshalNativeToManaged functions of a custom marshaller.
  std::map<Descriptor const*, std::string> cs_clean_up_native_definition_;
  std::map<Descriptor const*, std::string> cs_managed_to_native_definition_;
  std::map<Descriptor const*, std::string> cs_native_to_managed_definition_;

  // The definitions of the Serialize, Deserialize and (if needed) Insert
  // functions for interchange messages.  The key is a descriptor for an
  // interchange message.
  std::map<Descriptor const*, std::string> cxx_deserialize_definition_;
  std::map<Descriptor const*, std::string> cxx_serialize_definition_;
  std::map<Descriptor const*, std::string> cxx_insert_definition_;

  // For interchange messages that require deserialization storage, the
  // arguments for the storage in the call to the Deserialize function.  Starts
  // with a comma, arguments are comma-separated.
  std::map<Descriptor const*, std::string>
      cxx_deserialization_storage_arguments_;

  // For interchange messages that require deserialization storage, the
  // declaration of the storage in the caller of the Deserialize function.
  std::map<Descriptor const*, std::string>
      cxx_deserialization_storage_declarations_;

  // For interchange messages that require deserialization storage, the
  // parameters for the storage in the Deserialize function.
  std::map<Descriptor const*, std::string>
      cxx_deserialization_storage_parameters_;

  // The statements to be included in the body of the Play function.
  std::map<Descriptor const*, std::string> cxx_play_statement_;

  // The entire sequence of statements for the body of a Fill function.  The key
  // is a descriptor for an In or Out message.
  std::map<Descriptor const*, std::string> cxx_fill_body_;

  // A code snippet that goes before the call to the interface in the
  // body of the Run function.  The key is a descriptor for an In or Out
  // message.  Produced but not used for Out messages.
  std::map<Descriptor const*, std::string> cxx_run_body_prolog_;

  // A list of code snippets for arguments to be passed to the interface in the
  // body of the Run function.
  std::map<Descriptor const*, std::vector<std::string>> cxx_run_arguments_;

  // A code snippet that goes after the call to the interface in the
  // body of the Run function.  The key is a descriptor for an In or Out
  // message.
  std::map<Descriptor const*, std::string> cxx_run_body_epilog_;

  // A code snippet for the implementation of the Fill and Run functions.  The
  // key is a descriptor for a method message.
  std::map<Descriptor const*, std::string> cxx_functions_implementation_;

  // A code snippet for the declaration of the top-level struct for a method.
  // The key is a descriptor for a method message.
  std::map<Descriptor const*, std::string> cxx_toplevel_type_declaration_;

  // A code snippet for the declaration of a nested In or Out struct or a Return
  // typedef.  The key is a descriptor for an In, Out or Return message.
  std::map<Descriptor const*, std::string> cxx_nested_type_declaration_;
};

}  // namespace internal

using internal::JournalProtoProcessor;

}  // namespace _journal_proto_processor
}  // namespace tools
}  // namespace principia
