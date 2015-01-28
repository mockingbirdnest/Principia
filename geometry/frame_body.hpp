#pragma once

#include "geometry/frame.hpp"

#include <string>

#include "base/fingerprint2011.hpp"
#include "geometry/named_quantities.hpp"
#include "google/protobuf/descriptor.h"

namespace principia {
namespace geometry {

template<typename Tag, Tag tag, bool frame_is_inertial>
void Frame<Tag, tag, frame_is_inertial>::WriteToMessage(
    not_null<serialization::Frame*> const message) {
  string const& tag_type_full_name =
      google::protobuf::GetEnumDescriptor<Tag>()->full_name();

  message->set_tag_type_fingerprint(
      Fingerprint2011(tag_type_full_name.c_str(), tag_type_full_name.size()));
  message->set_tag(tag);
  message->set_is_inertial(frame_is_inertial);
}

template<typename Tag, Tag tag, bool frame_is_inertial>
void Frame<Tag, tag, frame_is_inertial>::ReadFromMessage(
    serialization::Frame const& message) {
  string const& tag_type_full_name =
      google::protobuf::GetEnumDescriptor<Tag>()->full_name();

  CHECK_EQ(Fingerprint2011(tag_type_full_name.c_str(),
                           tag_type_full_name.size()),
           message.tag_type_fingerprint()) << tag_type_full_name;
  CHECK_EQ(tag, message.tag());
  CHECK_EQ(frame_is_inertial, message.is_inertial());
}

// Default-initialized to {0, 0, 0}.
template<typename Tag, Tag tag, bool frame_is_inertial>
Position<Frame<Tag, tag, frame_is_inertial>> const
Frame<Tag, tag, frame_is_inertial>::origin;

}  // namespace geometry
}  // namespace principia
