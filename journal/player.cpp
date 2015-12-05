#include "journal/player.hpp"

#include <string>

#include "base/array.hpp"
#include "base/get_line.hpp"
#include "base/hexadecimal.hpp"
#include "journal/profiles.hpp"
#include "glog/logging.h"

namespace principia {

using base::GetLine;
using base::HexadecimalDecode;
using base::UniqueBytes;

namespace journal {

Player::Player(std::experimental::filesystem::path const& path)
    : stream_(path, std::ios::in) {
  CHECK(!stream_.fail());
}

bool Player::Play() {
  std::unique_ptr<serialization::Method> method = Read();
  if (method == nullptr) {
    LOG(ERROR)<<"null";
    return true;//false;
  }
  LOG(ERROR)<<method->DebugString();

  bool ran = false;
  ran |= RunIfAppropriate<ActivateRecorder>(*method);
  ran |= RunIfAppropriate<AddVesselToNextPhysicsBubble>(*method);
  ran |= RunIfAppropriate<AdvanceTime>(*method);
  ran |= RunIfAppropriate<AtEnd>(*method);
  ran |= RunIfAppropriate<BubbleDisplacementCorrection>(*method);
  ran |= RunIfAppropriate<BubbleVelocityCorrection>(*method);
  ran |= RunIfAppropriate<CelestialFromParent>(*method);
  ran |= RunIfAppropriate<CurrentTime>(*method);
  ran |= RunIfAppropriate<DeleteLineAndIterator>(*method);
  ran |= RunIfAppropriate<DeleteNavigationFrame>(*method);
  ran |= RunIfAppropriate<DeletePlugin>(*method);
  ran |= RunIfAppropriate<DeletePluginSerialization>(*method);
  ran |= RunIfAppropriate<DeserializePlugin>(*method);
  ran |= RunIfAppropriate<DirectlyInsertCelestial>(*method);
  ran |= RunIfAppropriate<EndInitialization>(*method);
  ran |= RunIfAppropriate<FetchAndIncrement>(*method);
  ran |= RunIfAppropriate<FlightPlanSize>(*method);
  ran |= RunIfAppropriate<ForgetAllHistoriesBefore>(*method);
  ran |= RunIfAppropriate<GetBufferDuration>(*method);
  ran |= RunIfAppropriate<GetBufferedLogging>(*method);
  ran |= RunIfAppropriate<GetStderrLogging>(*method);
  ran |= RunIfAppropriate<GetSuppressedLogging>(*method);
  ran |= RunIfAppropriate<GetVerboseLogging>(*method);
  ran |= RunIfAppropriate<HasPrediction>(*method);
  ran |= RunIfAppropriate<HasVessel>(*method);
  ran |= RunIfAppropriate<InitGoogleLogging>(*method);
  ran |= RunIfAppropriate<InsertCelestial>(*method);
  ran |= RunIfAppropriate<InsertOrKeepVessel>(*method);
  ran |= RunIfAppropriate<InsertSun>(*method);
  ran |= RunIfAppropriate<LogError>(*method);
  ran |= RunIfAppropriate<LogFatal>(*method);
  ran |= RunIfAppropriate<LogInfo>(*method);
  ran |= RunIfAppropriate<LogWarning>(*method);
  ran |= RunIfAppropriate<NavballOrientation>(*method);
  ran |= RunIfAppropriate<NewBarycentricRotatingNavigationFrame>(*method);
  ran |= RunIfAppropriate<NewBodyCentredNonRotatingNavigationFrame>(*method);
  ran |= RunIfAppropriate<NewPlugin>(*method);
  ran |= RunIfAppropriate<NumberOfSegments>(*method);
  ran |= RunIfAppropriate<PhysicsBubbleIsEmpty>(*method);
  ran |= RunIfAppropriate<RenderedFlightPlan>(*method);
  ran |= RunIfAppropriate<RenderedPrediction>(*method);
  ran |= RunIfAppropriate<RenderedVesselTrajectory>(*method);
  ran |= RunIfAppropriate<SayHello>(*method);
  ran |= RunIfAppropriate<SerializePlugin>(*method);
  ran |= RunIfAppropriate<SetBufferDuration>(*method);
  ran |= RunIfAppropriate<SetBufferedLogging>(*method);
  ran |= RunIfAppropriate<SetPredictionLength>(*method);
  ran |= RunIfAppropriate<SetPredictionLengthTolerance>(*method);
  ran |= RunIfAppropriate<SetPredictionSpeedTolerance>(*method);
  ran |= RunIfAppropriate<SetStderrLogging>(*method);
  ran |= RunIfAppropriate<SetSuppressedLogging>(*method);
  ran |= RunIfAppropriate<SetVerboseLogging>(*method);
  ran |= RunIfAppropriate<SetVesselStateOffset>(*method);
  ran |= RunIfAppropriate<UpdateCelestialHierarchy>(*method);
  ran |= RunIfAppropriate<UpdatePrediction>(*method);
  ran |= RunIfAppropriate<VesselFromParent>(*method);
  ran |= RunIfAppropriate<VesselBinormal>(*method);
  ran |= RunIfAppropriate<VesselNormal>(*method);
  ran |= RunIfAppropriate<VesselTangent>(*method);
  CHECK(ran) << method->DebugString();

  return true;
}

std::unique_ptr<serialization::Method> Player::Read() {
  std::string const line = GetLine(&stream_);
  if (line.empty()) {
    return nullptr;
  }

  uint8_t const* const hexadecimal =
      reinterpret_cast<uint8_t const*>(line.c_str());
  int const hexadecimal_size = strlen(line.c_str());
  UniqueBytes bytes(hexadecimal_size >> 1);
  HexadecimalDecode({hexadecimal, hexadecimal_size},
                    {bytes.data.get(), bytes.size});
  auto method = std::make_unique<serialization::Method>();
  CHECK(method->ParseFromArray(bytes.data.get(),
                               static_cast<int>(bytes.size)));

  return method;
}

}  // namespace journal
}  // namespace principia
