#pragma once

#include "physics/solar_system_factory.hpp"

#include <fstream>

#include "glog/logging.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/text_format.h"
#include "serialization/astronomy.pb.h"

namespace principia {
namespace physics {

void SolarSystemFactory::Initialize(std::string const& initial_state_filename,
                                    std::string const& gravity_model_filename) {
  serialization::SolarSystemFile initial_state;
  std::ifstream initial_state_ifstream(initial_state_filename);
  CHECK(initial_state_ifstream.good());
  google::protobuf::io::IstreamInputStream initial_state_zcs(
                                               &initial_state_ifstream);
  CHECK(google::protobuf::TextFormat::Parse(&initial_state_zcs,
                                            &initial_state));

  serialization::SolarSystemFile gravity_model;
  std::ifstream gravity_model_ifstream(gravity_model_filename);
  CHECK(gravity_model_ifstream.good());
  google::protobuf::io::IstreamInputStream gravity_model_zcs(
                                               &gravity_model_ifstream);
  CHECK(google::protobuf::TextFormat::Parse(&gravity_model_zcs,
                                            &gravity_model));
}

}  // namespace physics
}  // namespace principia
