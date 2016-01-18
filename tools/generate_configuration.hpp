
#pragma once

#include <string>

#include "geometry/named_quantities.hpp"

namespace principia {

using geometry::Instant;

namespace tools {

void GenerateConfiguration(Instant const& game_epoch,
                           std::string const& gravity_model_stem,
                           std::string const& initial_state_stem);

}  // namespace tools
}  // namespace principia
