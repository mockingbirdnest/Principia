
#include <iostream>
#include <string>

#include "geometry/epoch.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "quantities/parser.hpp"
#include "tools/generate_configuration.hpp"
#include "tools/generate_profiles.hpp"

int main(int argc, char const* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::LogToStderr();
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " command [arguments...]\n";
    return 1;
  }
  std::string command = argv[1];
  if (command == "generate_configuration") {
    if (argc != 5) {
      // tools.exe generate_configuration \
      //     2433647.5 gravity_model initial_state_jd_2433282_500000000
      std::cerr << "Usage: " << argv[0] << " " << argv[1] << " "
                << "game_epoch_jd gravity_model_stem initial_state_stem\n";
      return 2;
    }
    principia::geometry::Instant game_epoch =
        principia::geometry::JulianDate(
            principia::quantities::ParseQuantity<double>(argv[1]));
    std::string const gravity_model_stem = argv[2];
    std::string const initial_state_stem = argv[3];
    principia::tools::GenerateConfiguration(game_epoch,
                                            gravity_model_stem,
                                            initial_state_stem);
    return 0;
  } else if (command == "generate_profiles") {
    if (argc != 2) {
      // tools.exe generate_profiles
      std::cerr << "Usage: " << argv[0] << " " << argv[1] << "\n";
      return 3;
    }
    principia::tools::GenerateProfiles();
    return 0;
  } else {
    std::cerr << "Usage: " << argv[0]
              << " generate_configuration|generate_profiles\n";
    return 4;
  }
}
