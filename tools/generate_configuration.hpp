
#pragma once

#include <string>

#include "geometry/named_quantities.hpp"

namespace principia {
namespace tools {
namespace internal_generate_configuration {

using geometry::Instant;

void GenerateConfiguration(Instant const& game_epoch,
                           std::string const& gravity_model_stem,
                           std::string const& initial_state_stem);

}  // namespace internal_generate_configuration

using internal_generate_configuration::GenerateConfiguration;

}  // namespace tools
}  // namespace principia
