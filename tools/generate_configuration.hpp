#pragma once

#include <string>

namespace principia {
namespace tools {
namespace _generate_configuration {
namespace internal {

void GenerateConfiguration(std::string const& game_epoch,
                           std::string const& gravity_model_stem,
                           std::string const& initial_state_stem,
                           std::string const& numerics_blueprint_stem,
                           std::string const& needs);

}  // namespace internal

using internal::GenerateConfiguration;

}  // namespace _generate_configuration
}  // namespace tools
}  // namespace principia
