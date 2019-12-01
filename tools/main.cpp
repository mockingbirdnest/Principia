
#include <iostream>
#include <string>

#include "astronomy/epoch.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "quantities/parser.hpp"
#include "tools/generate_configuration.hpp"
#include "tools/generate_kopernicus.hpp"
#include "tools/generate_profiles.hpp"

int __cdecl main(int argc, char const* argv[]) {
  google::SetLogFilenameExtension(".log");
  google::InitGoogleLogging(argv[0]);
  google::LogToStderr();
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " command [arguments...]\n";
    return 1;
  }
  std::string command = argv[1];
  if (command == "generate_configuration") {
    if (argc != 7) {
      // tools.exe generate_configuration \
      //     JD2433647.5 \
      //     sol_gravity_model \
      //     sol_initial_state_jd_2433282_500000000 \
      //     sol_numerics_blueprint \
      //     RealSolarSystem
      // tools.exe generate_configuration \
      //     JD2457000.0 \
      //     trappist_gravity_model \
      //     trappist_initial_state_jd_2457000_000000000 \
      //     trappist_numerics_blueprint \
      //     aSLIPPIST-1
      std::cerr << "Usage: " << argv[0] << " " << argv[1] << " "
                << "game_epoch_jd "
                << "gravity_model_stem "
                << "initial_state_stem "
                << "numerics_blueprint_stem "
                << "needs\n";
      return 2;
    }
    std::string const game_epoch = argv[2];
    std::string const gravity_model_stem = argv[3];
    std::string const initial_state_stem = argv[4];
    std::string const numerics_blueprint_stem = argv[5];
    std::string const needs = argv[6];
    principia::tools::GenerateConfiguration(game_epoch,
                                            gravity_model_stem,
                                            initial_state_stem,
                                            numerics_blueprint_stem,
                                            needs);
    return 0;
  } else if (command == "generate_kopernicus") {
    if (argc != 4) {
      // tools.exe generate_kopernicus \
      //     trappist_gravity_model \
      //     trappist_initial_state_jd_2457000_000000000
      std::cerr << "Usage: " << argv[0] << " " << argv[1] << " "
                << "gravity_model_stem "
                << "initial_state_stem\n";
      return 5;
    }
    std::string const gravity_model_stem = argv[2];
    std::string const initial_state_stem = argv[3];
    principia::tools::GenerateKopernicusForSlippist1(gravity_model_stem,
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
