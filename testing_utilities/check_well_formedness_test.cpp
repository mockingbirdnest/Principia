#include "testing_utilities/check_well_formedness.hpp"

#include <string>
#include <utility>

PRINCIPIA_CHECK_ILL_FORMED_WITH_TYPES(T::npos, (typename T = int));
PRINCIPIA_CHECK_WELL_FORMED_WITH_TYPES(T::npos, (typename T = std::string));
PRINCIPIA_CHECK_ILL_FORMED(s + 3, with_variable<std::string> s);
PRINCIPIA_CHECK_ILL_FORMED(s + "3", with_variable<std::pair<double, double>> s);
PRINCIPIA_CHECK_WELL_FORMED(s + t,
                            with_variable<std::string> s,
                            with_variable<std::string> t);
PRINCIPIA_CHECK_ILL_FORMED(s - t,
                           with_variable<std::string> s,
                           with_variable<std::string> t);
