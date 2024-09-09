#include "base/traits.hpp"

namespace principia {
namespace base {

using namespace principia::base::_traits;

template<typename T>
struct T1;

template<typename T>
struct T2;

template<typename T>
using AliasT1 = T1<T>;

static_assert(is_instance_of_v<T1, T1<int>>);
static_assert(!is_instance_of_v<T2, T1<int>>);

// Note that `is_instance_of` does not work on alias templates.
static_assert(!is_instance_of_v<AliasT1, T1<int>>);


}  // namespace base
}  // namespace principia
