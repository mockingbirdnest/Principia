
#pragma once

#include "geometry/frame.hpp"

#include <string>

#include "base/fingerprint2011.hpp"
#include "geometry/named_quantities.hpp"
#include "google/protobuf/descriptor.h"

namespace principia {
namespace geometry {
namespace internal_frame {

using base::Fingerprint2011;

// Utility for fingerprinting.
inline uint32_t Fingerprint(std::string const& s) {
  return Fingerprint2011(s.c_str(), s.size()) & 0xFFFFFFFF;
}

template<typename FrameTag, FrameTag tag_,
         FrameMotion motion_, Handedness handedness_>
void Frame<FrameTag, tag_, motion_, handedness_>::WriteToMessage(
  not_null<serialization::Frame*> const message) {
  if constexpr (google::protobuf::is_proto_enum<FrameTag>::value) {
    std::string const& tag_type_full_name =
        google::protobuf::GetEnumDescriptor<Tag>()->full_name();

    message->set_tag_type_fingerprint(Fingerprint(tag_type_full_name));
    message->set_tag(static_cast<int>(tag));
    message->set_is_inertial(is_inertial);
  } else {
    LOG(FATAL) << "Non serializable frame with tag " << static_cast<int>(tag_)
               << ", motion " << motion_ << ", handedness "
               << static_cast<int>(handedness_);
  }
}

template<typename FrameTag, FrameTag tag_,
         FrameMotion motion_, Handedness handedness_>
template<typename T, typename>
void Frame<FrameTag, tag_, motion_, handedness_>::ReadFromMessage(
  serialization::Frame const& message) {
  std::string const& tag_type_full_name =
      google::protobuf::GetEnumDescriptor<Tag>()->full_name();

  CHECK_EQ(Fingerprint(tag_type_full_name), message.tag_type_fingerprint())
      << tag_type_full_name;
  CHECK_EQ(static_cast<int>(tag), message.tag());
  CHECK_EQ(is_inertial, message.is_inertial());
}

// Default-initialized to {0, 0, 0}.
template<typename FrameTag, FrameTag tag_,
         FrameMotion motion_, Handedness handedness_>
Position<Frame<FrameTag, tag_, motion_, handedness_>> const
Frame<FrameTag, tag_, motion_, handedness_>::origin;
template<typename FrameTag, FrameTag tag_,
         FrameMotion motion_, Handedness handedness_>
Velocity<Frame<FrameTag, tag_, motion_, handedness_>> const
Frame<FrameTag, tag_, motion_, handedness_>::unmoving;

inline void ReadFrameFromMessage(
    serialization::Frame const& message,
    google::protobuf::EnumValueDescriptor const*& enum_value_descriptor,
    bool& is_inertial) {
  // Look at the enumeration types nested in serialization::Frame for one that
  // matches our fingerprint.
  const google::protobuf::Descriptor* frame_descriptor =
      serialization::Frame::descriptor();
  enum_value_descriptor = nullptr;
  for (int i = 0; i < frame_descriptor->enum_type_count(); ++i) {
    const google::protobuf::EnumDescriptor* enum_type_descriptor =
        frame_descriptor->enum_type(i);
    std::string const& enum_type_full_name = enum_type_descriptor->full_name();
    if (Fingerprint(enum_type_full_name) == message.tag_type_fingerprint()) {
      enum_value_descriptor =
          enum_type_descriptor->FindValueByNumber(message.tag());
      break;
    }
  }
  CHECK_NOTNULL(enum_value_descriptor);
  is_inertial = message.is_inertial();
}

}  // namespace internal_frame
}  // namespace geometry
}  // namespace principia
