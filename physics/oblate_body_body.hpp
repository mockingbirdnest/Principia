
#pragma once

#include "physics/oblate_body.hpp"

#include <algorithm>
#include <set>
#include <vector>

#include "astronomy/epoch.hpp"
#include "numerics/legendre_normalization_factor.mathematica.h"
#include "quantities/constants.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_oblate_body {

using astronomy::J2000;
using geometry::AngularVelocity;
using geometry::Instant;
using numerics::LegendreNormalizationFactor;
using quantities::Angle;
using quantities::SIUnit;
using quantities::si::Radian;
using quantities::si::Second;

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(double const j2,
                                          Length const& reference_radius)
    : reference_radius_(reference_radius),
      j2_(j2),
      j2_over_μ_(j2 * reference_radius * reference_radius),
      cos_(typename OblateBody<Frame>::GeopotentialCoefficients()),
      sin_(typename OblateBody<Frame>::GeopotentialCoefficients()),
      degree_(2),
      is_zonal_(true) {
  CHECK_LT(0.0, j2) << "Oblate body must have positive j2";
  cos_[2][0] = -j2 / LegendreNormalizationFactor[2][0];
}

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(Length const& reference_radius)
    : reference_radius_(reference_radius),
      cos_(typename OblateBody<Frame>::GeopotentialCoefficients()),
      sin_(typename OblateBody<Frame>::GeopotentialCoefficients()),
      degree_(0),
      is_zonal_(false) {}

template<typename Frame>
typename OblateBody<Frame>::Parameters
OblateBody<Frame>::Parameters::ReadFromMessage(
    serialization::OblateBody::Geopotential const& message,
    Length const& reference_radius) {
  Parameters parameters(reference_radius);
  std::set<int> degrees_seen;
  for (auto const& row : message.row()) {
    const int n = row.degree();
    CHECK_LE(n, OblateBody<Frame>::max_geopotential_degree);
    CHECK(degrees_seen.insert(n).second)
        << "Degree " << n << " specified multiple times";
    CHECK_LE(row.column_size(), n + 1)
        << "Degree " << n << " has " << row.column_size() << " coefficients";
    std::set<int> orders_seen;
    for (auto const& column : row.column()) {
      const int m = column.order();
      CHECK_LE(m, n);
      CHECK(orders_seen.insert(m).second)
          << "Degree " << n << " order " << m << " specified multiple times";

      // If j was specified, check that it is legit and compute cos.
      double cos = column.cos();
      if (m == 0) {
        if (column.has_j()) {
          CHECK(!column.has_cos()) << "Cos and J specified for degree " << n;
          cos = -column.j() / LegendreNormalizationFactor[n][0];
        } else {
          CHECK(column.has_cos())
              << "Cos missing for degree " << n << " order " << m;
        }
      } else {
        CHECK(!column.has_j())
            << "J specified for degree " << n << " and nonzero order " << m;
        CHECK(column.has_cos())
            << "Cos missing for degree " << n << " order " << m;
      }
      parameters.cos_[n][m] = cos;
      parameters.sin_[n][m] = column.sin();
    }
  }
  parameters.degree_ = *degrees_seen.crbegin();

  // Unnormalization.
  parameters.j2_ = -parameters.cos_[2][0] * LegendreNormalizationFactor[2][0];
  parameters.j2_over_μ_ = -parameters.cos_[2][0] *
                          LegendreNormalizationFactor[2][0] * reference_radius *
                          reference_radius;

  // Zonalness.
  parameters.is_zonal_ = true;
  for (int n = 0; n <= parameters.degree_; ++n) {
    for (int m = 1; m <= n; ++m) {
      if (parameters.cos_[n][m] != 0 || parameters.sin_[n][m] != 0) {
        parameters.is_zonal_ = false;
        break;
      }
    }
  }

  return parameters;
}

template<typename Frame>
void OblateBody<Frame>::Parameters::WriteToMessage(
    not_null<serialization::OblateBody::Geopotential*> const message) const {
  for (int n = 0; n <= degree_; ++n) {
    auto const row = message->add_row();
    row->set_degree(n);
    for (int m = 0; m <= n; ++m) {
      auto const column = row->add_column();
      column->set_order(m);
      column->set_cos(cos_[n][m]);
      column->set_sin(sin_[n][m]);
    }
  }
}

template<typename Frame>
OblateBody<Frame>::OblateBody(
    MassiveBody::Parameters const& massive_body_parameters,
    typename RotatingBody<Frame>::Parameters const& rotating_body_parameters,
    Parameters const& parameters)
    : RotatingBody<Frame>(massive_body_parameters, rotating_body_parameters),
      parameters_(parameters) {}

template<typename Frame>
double OblateBody<Frame>::j2() const {
  return parameters_.j2_;
}

template<typename Frame>
Quotient<Degree2SphericalHarmonicCoefficient,
         GravitationalParameter> const& OblateBody<Frame>::j2_over_μ() const {
  return parameters_.j2_over_μ_;
}

template<typename Frame>
typename OblateBody<Frame>::GeopotentialCoefficients const&
OblateBody<Frame>::cos() const {
  return parameters_.cos_;
}

template<typename Frame>
typename OblateBody<Frame>::GeopotentialCoefficients const&
OblateBody<Frame>::sin() const {
  return parameters_.sin_;
}

template<typename Frame>
int OblateBody<Frame>::geopotential_degree() const {
  return parameters_.degree_;
}

template<typename Frame>
bool OblateBody<Frame>::is_zonal() const {
  return parameters_.is_zonal_;
}

template<typename Frame>
Length const& OblateBody<Frame>::reference_radius() const {
  return parameters_.reference_radius_;
}

template<typename Frame>
bool OblateBody<Frame>::is_massless() const {
  return false;
}

template<typename Frame>
bool OblateBody<Frame>::is_oblate() const {
  return true;
}

template<typename Frame>
void OblateBody<Frame>::WriteToMessage(
    not_null<serialization::Body*> const message) const {
  WriteToMessage(message->mutable_massive_body());
}

template<typename Frame>
void OblateBody<Frame>::WriteToMessage(
    not_null<serialization::MassiveBody*> const message) const {
  RotatingBody<Frame>::WriteToMessage(message);
  not_null<serialization::OblateBody*> const oblate_body =
      message->MutableExtension(serialization::RotatingBody::extension)->
               MutableExtension(serialization::OblateBody::extension);
  parameters_.reference_radius_.WriteToMessage(
      oblate_body->mutable_reference_radius());
  parameters_.WriteToMessage(oblate_body->mutable_geopotential());
}

template<typename Frame>
not_null<std::unique_ptr<OblateBody<Frame>>> OblateBody<Frame>::ReadFromMessage(
      serialization::OblateBody const& message,
      MassiveBody::Parameters const& massive_body_parameters,
      typename RotatingBody<Frame>::Parameters const&
          rotating_body_parameters) {
  std::unique_ptr<Parameters> parameters;
  switch (message.oblateness_case()) {
    case serialization::OblateBody::OblatenessCase::kPreDiophantosJ2: {
      // In the legacy case we didn't record the reference radius, so we use a
      // dummy value to achieve the right effect.
      CHECK(!message.has_reference_radius()) << message.DebugString();
      Length const reference_radius = SIUnit<Length>();
      parameters = std::make_unique<Parameters>(
          Degree2SphericalHarmonicCoefficient::ReadFromMessage(
              message.pre_diophantos_j2()) /
              (massive_body_parameters.gravitational_parameter() *
               reference_radius * reference_radius),
          reference_radius);
      break;
    }
    case serialization::OblateBody::OblatenessCase::kGeopotential:
      CHECK(message.has_reference_radius()) << message.DebugString();
      parameters = std::make_unique<Parameters>(Parameters::ReadFromMessage(
          message.geopotential(),
          Length::ReadFromMessage(message.reference_radius())));
      break;
    case serialization::OblateBody::OblatenessCase::OBLATENESS_NOT_SET:
      LOG(FATAL) << message.DebugString();
      base::noreturn();
  }
  return std::make_unique<OblateBody<Frame>>(massive_body_parameters,
                                             rotating_body_parameters,
                                             *parameters);
}

}  // namespace internal_oblate_body
}  // namespace physics
}  // namespace principia
