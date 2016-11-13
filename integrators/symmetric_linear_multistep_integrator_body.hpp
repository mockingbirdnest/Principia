#pragma once

#include "integrators/symmetric_linear_multistep_integrator.hpp"

namespace principia {
namespace integrators {

namespace {

template<int size>
constexpr FixedVector<double, size> Divide(
    FixedVector<double, size> const& numerators,
    double const denominator) {
  FixedVector<double, size> result = numerators;
  for (int i = 0; i < size; ++i) {
    numerators /= denominator;
  }
  return result;
}

}  // namespace

template<typename Position, int order_>
SymmetricLinearMultistepIntegrator<Position, order_>::
SymmetricLinearMultistepIntegrator(
    serialization::FixedStepSizeIntegrator::Kind const kind,
    FixedVector<double, half_order_> const & ɑ,
    FixedVector<double, half_order_> const & β_numerators,
    double const β_denominator)
    : ɑ_(ɑ), β_(Divide(β_numerators, β_denominator)) {}

template<typename Position, int order_>
void SymmetricLinearMultistepIntegrator<Position, order_>::Solve(
    IntegrationProblem<ODE> const& problem,
    Time const& step) const {}

template <typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const& Quinlan1999Order8A() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_1999_ORDER_8A,
      {1, -2, 2, -2, 2},
      {0, 22081, -29418, 75183, -75212},
      15120);
  return integrator;
}

template <typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const& Quinlan1999Order8B() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_1999_ORDER_8B,
      {1, 0, 0, -1 / 2, -1},
      {0, 192481, 6582, 816783, -156812},
      120960);
  return integrator;
}

template <typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const&
QuinlanTremaine1990Order8() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_8,
      {1, -2, 2, -1, 0},
      {0, 17671, -23622, 61449, 50516},
      12096);
  return integrator;
}

template <typename Position>
SymmetricLinearMultistepIntegrator<Position, 10> const&
QuinlanTremaine1990Order10() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_10,
      {1, -1, 1, -1, 1, 2},
      {0, 399187, -485156, 2391436, -2816732, 4651330},
      241920);
  return integrator;
}

template <typename Position>
SymmetricLinearMultistepIntegrator<Position, 12> const&
QuinlanTremaine1990Order12() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_12,
      {1, -2, 2, -1, 0, 0, 0},
      {0,
       90987349,
       -229596838,
       812627169,
       -1628539944,
       2714971338,
       -3041896548},
      53222400);
  return integrator;
}

template <typename Position>
SymmetricLinearMultistepIntegrator<Position, 14> const&
QuinlanTremaine1990Order14() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_14,
      {1, -2, 2, -1, 0, 0, 0, 0},
      {433489274083,
       -1364031998256,
       5583113380398,
       -14154444148720,
       28630585332045,
       -42056933842656,
       48471792742212},
      237758976000);
  return integrator;
}

}  // namespace integrators
}  // namespace principia
