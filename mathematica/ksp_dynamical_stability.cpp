
#include "mathematica/ksp_dynamical_stability.hpp"

#include <fstream>
#include <memory>

#include "base/array.hpp"
#include "base/get_line.hpp"
#include "base/hexadecimal.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/hierarchical_system.hpp"

namespace principia {

using base::GetLine;
using base::HexadecimalDecode;
using base::UniqueBytes;
using ksp_plugin::Barycentric;
using physics::DegreesOfFreedom;
using physics::HierarchicalSystem;
using physics::MassiveBody;

namespace mathematica {

namespace {

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
}

void SimulateFixedSystem() {
}

}  // namespace mathematica
}  // namespace principia
