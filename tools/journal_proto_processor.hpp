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

private:
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& name,
                                     std::string const& expr)>>
           field_copy_wrapper_;
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
           field_serializer_wrapper_;
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr)>>
           indirect_field_wrapper_;
  std::map<FieldDescriptor const*,
           std::function<std::string(std::string const& expr,
                                     std::string const& stmt)>>
           optional_field_wrapper_;

  std::map<FieldDescriptor const*, std::string> size_field_name_;
  std::set<FieldDescriptor const*> in_out_field_;
  std::map<FieldDescriptor const*, std::string> cpp_field_type_;
  std::map<Descriptor const*, std::string> cpp_fill_body_;
  std::map<Descriptor const*, std::string> cpp_method_impl_;
  std::map<Descriptor const*, std::string> cpp_method_type_;
  std::map<Descriptor const*, std::string> cpp_nested_type_;
};

}  // namespace tools
}  // namespace principia
