
#pragma once

#include <string>

namespace principia {
namespace tools {

void GenerateConfiguration(std::string const& game_epoch,
                           std::string const& gravity_model_stem,
                           std::string const& initial_state_stem,
                           std::string const& numerics_blueprint_stem,
                           std::string const& needs);

}  // namespace tools
}  // namespace principia
