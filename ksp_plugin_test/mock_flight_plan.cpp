#include "ksp_plugin_test/mock_flight_plan.hpp"

namespace principia {
namespace ksp_plugin {

bool MockFlightPlan::Append(Burn burn) {
  return AppendConstRef(burn);
}

bool MockFlightPlan::ReplaceLast(Burn burn) {
  return ReplaceLastConstRef(burn);
}

}  // namespace ksp_plugin
}  // namespace principia
