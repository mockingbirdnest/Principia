#pragma once

#include "integrators/ordinary_differential_equations.hpp"

namespace principia {
namespace integrators {

template<typename Position>
void
SpecialSecondOrderDifferentialEquation<Position>::SystemState::WriteToMessage(
        not_null<serialization::SystemState*> const message) const {
  for (auto const& position : positions) {
    position.WriteToMessage(message->add_position());
  }
  for (auto const& velocity : velocities) {
    velocity.WriteToMessage(message->add_velocity());
  }
  time.WriteToMessage(message->mutable_time());
}

template<typename Position>
typename SpecialSecondOrderDifferentialEquation<Position>::SystemState
SpecialSecondOrderDifferentialEquation<Position>::SystemState::ReadFromMessage(
        serialization::SystemState const& message) {
  for (auto const p : message.position()) {
    positions.push_back(DoublePrecision<Position>::ReadFromMessage(p));
  }
  for (auto const v : message.velocity()) {
    velocities.push_back(DoublePrecision<Velocity>::ReadFromMessage(v));
  }
  time = Instant::ReadFromMessage(message.time());
}

template<typename DifferentialEquation>
FixedStepSizeIntegrator<DifferentialEquation>::FixedStepSizeIntegrator(
    serialization::FixedStepSizeIntegrator::Kind const kind) : kind_(kind) {}

template<typename DifferentialEquation>
void FixedStepSizeIntegrator<DifferentialEquation>::WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> const message) const {
  message->set_kind(kind_);
}

template<typename DifferentialEquation>
FixedStepSizeIntegrator<DifferentialEquation> const&
FixedStepSizeIntegrator<DifferentialEquation>::ReadFromMessage(
      serialization::FixedStepSizeIntegrator const& message) {
  using FSSI = serialization::FixedStepSizeIntegrator;
  switch (message.kind()) {
    case FSSI::BLANES_MOAN_2002_SRKN_6B:
      return BlanesMoan2002SRKN6B();
    case FSSI::BLANES_MOAN_2002_SRKN_11B:
      return BlanesMoan2002SRKN11B();
    case FSSI::BLANES_MOAN_2002_SRKN_14A:
      return BlanesMoan2002SRKN14A();
    case FSSI::MCLACHLAN_1995_SB3A_4:
      return McLachlan1995SB3A4();
    case FSSI::MCLACHLAN_1995_SB3A_5:
      return McLachlan1995SB3A5();
    case FSSI::MCLACHLAN_ATELA_1992_ORDER_4_OPTIMAL:
      return McLachlanAtela1992Order4Optimal();
    case FSSI::MCLACHLAN_ATELA_1992_ORDER_5_OPTIMAL:
      return McLachlanAtela1992Order5Optimal();
    case FSSI::OKUNBOR_SKEEL_1994_ORDER_6_METHOD_13:
      return OkunborSkeel1994Order6Method13();
  }
}

template<typename DifferentialEquation>
AdaptiveStepSizeIntegrator<DifferentialEquation>::AdaptiveStepSizeIntegrator(
    serialization::AdaptiveStepSizeIntegrator::Kind const kind) : kind_(kind) {}

template<typename DifferentialEquation>
void AdaptiveStepSizeIntegrator<DifferentialEquation>::WriteToMessage(
    not_null<serialization::AdaptiveStepSizeIntegrator*> const message) const {
  message->set_kind(kind_);
}

template<typename DifferentialEquation>
AdaptiveStepSizeIntegrator<DifferentialEquation> const&
AdaptiveStepSizeIntegrator<DifferentialEquation>::ReadFromMessage(
    serialization::AdaptiveStepSizeIntegrator const& message) {
  using ASSI = serialization::AdaptiveStepSizeIntegrator;
  switch (message.kind()) {
    case ASSI::DORMAND_ELMIKKAWY_PRINCE_1986_RKN_434FM:
      return DormandElMikkawyPrince1986RKN434FM();
  }
}

}  // namespace integrators
}  // namespace principia
