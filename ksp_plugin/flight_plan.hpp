#pragma once

#include <vector>

#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {
namespace ksp_plugin {

class FlightPlan {
 public:
  FlightPlan(not_null<DiscreteTrajectory<Barycentric>*> root);
  ~FlightPlan() = default;

  int size() const;
  Manœuvre<Barycentric, Navigation> const& Get(int const index) const;
  void Delete(int const index);

  // The following three functions have no effect and return false if the
  // manœuvre is invalid.
  bool InsertBefore(int const index,
                    Manœuvre<Barycentric, Navigation> const& manœuvre);
  bool InsertAfter(int const index,
                   Manœuvre<Barycentric, Navigation> const& manœuvre);
  bool Replace(int const index,
               Manœuvre<Barycentric, Navigation> const& manœuvre);


 private:
  std
  std::vector<not_null<DiscreteTrajectory<Barycentric>*>> trajectories_;
  std::vector<Manœuvre<Barycentric, Navigation>> manœuvres_;
};

}  // namespace ksp_plugin
}  // namespace principia
