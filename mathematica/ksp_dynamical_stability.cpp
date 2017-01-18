
#include "mathematica/ksp_dynamical_stability.hpp"

#include <fstream>
#include <memory>

#include "base/array.hpp"
#include "base/get_line.hpp"
#include "base/hexadecimal.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/hierarchical_system.hpp"
#include "quantities/astronomy.hpp"

namespace principia {

using base::GetLine;
using base::HexadecimalDecode;
using base::UniqueBytes;
using geometry::Instant;
using geometry::Position;
using ksp_plugin::Barycentric;
using physics::DegreesOfFreedom;
using physics::Ephemeris;
using physics::HierarchicalSystem;
using physics::MassiveBody;
using quantities::astronomy::JulianYear;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;

namespace mathematica {

namespace {

constexpr Instant ksp_epoch;

template<typename Message>
std::unique_ptr<Message> Read(std::ifstream& file) {
  std::string const line = GetLine(file);
  if (line.empty()) {
    return nullptr;
  }
  std::uint8_t const* const hexadecimal =
      reinterpret_cast<std::uint8_t const*>(line.data());
  int const hexadecimal_size = line.size();
  UniqueBytes bytes(hexadecimal_size >> 1);
  HexadecimalDecode({hexadecimal, hexadecimal_size},
                    {bytes.data.get(), bytes.size});
  auto message = std::make_unique<Message>();
  CHECK(
      message->ParseFromArray(bytes.data.get(), static_cast<int>(bytes.size)));
  return std::move(message);
}

}

void SimulateStockSystem() {
  std::ifstream file("ksp_stock_system.proto.hex", std::ios::in);
  CHECK(!file.fail());
  HierarchicalSystem<Barycentric>::BarycentricSystem system;

  for (auto body = Read<serialization::MassiveBody>(file);
       body != nullptr;
       body = Read<serialization::MassiveBody>(file)) {
    auto const degrees_of_freedom = Read<serialization::Pair>(file);
    CHECK(degrees_of_freedom != nullptr);
    system.bodies.push_back(MassiveBody::ReadFromMessage(*body));
    system.degrees_of_freedom.push_back(
        DegreesOfFreedom<Barycentric>::ReadFromMessage(*degrees_of_freedom));
    LOG(ERROR) << system.bodies.back()->gravitational_parameter();
    LOG(ERROR) << system.degrees_of_freedom.back();
  }
  file.close();
  Ephemeris<Barycentric> ephemeris(
      std::move(system.bodies),
      system.degrees_of_freedom,
      ksp_epoch,
      /*fitting_tolerance=*/1 * Milli(Metre),
      Ephemeris<Barycentric>::FixedStepParameters(
          integrators::QuinlanTremaine1990Order12<Position<Barycentric>>(),
          /*step=*/8 * Minute));
  ephemeris.Prolong(ksp_epoch + 100 * JulianYear);
}

void SimulateFixedSystem() {
}

}  // namespace mathematica
}  // namespace principia
