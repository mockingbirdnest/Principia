#pragma once

#include "geometry/grassmann.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

using principia::geometry::Bivector;
using principia::geometry::Vector;
using principia::quantities::Acceleration;
using principia::quantities::AngularFrequency;
using principia::quantities::Force;
using principia::quantities::Time;

namespace principia {
namespace physics {

// The tag is just a way to unambiguously distinguish frames.
template<int tag>
class Frame {
  // The construction parameters describe the |acceleration| and |rotation| of
  // this frame with respect to the |ReferenceFrame|.
  //TODO(phl): ReferenceFrame is unclear.
  template<typename ReferenceFrame>
  Frame(Vector<Acceleration> (*acceleration)(Time const& time),
        Bivector<AngularFrequency> (*rotation)(Time const& time));

  template<typename ReferenceFrame>
  Vector<Force> FictitiousForce(Time const& time) const;
};

// The frame with no fictitious.
Frame<0> Inertial();

}  // namespace physics
}  // namespace principia
