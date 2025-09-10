#include "base/traits.hpp"

namespace principia {
namespace base {

using namespace principia::base::_traits;

template<typename T>
struct T1;

template<typename T>
struct T2;

template<typename T1, typename T2>
struct T3;

template<typename T>
using AliasT1 = T1<T>;

static_assert(internal::can_be_instantiated<AliasT1, double>::value);
static_assert(!internal::can_be_instantiated<AliasT1, double, double>::value);

static_assert(
    internal::matches_instantiation<AliasT1<double>, true, AliasT1, double>::
        value);
static_assert(
    internal::matches_instantiation<T1<double>, true, AliasT1, double>::value);
static_assert(
    internal::matches_instantiation<AliasT1<double>, true, T1, double>::value);
static_assert(
    internal::matches_instantiation<T1<double>, true, T1, double>::value);

static_assert(is_instance_of_v<T1, T1<int>>);
static_assert(is_instance_of_v<AliasT1, T1<int>>);
static_assert(!is_instance_of_v<AliasT1, double>);
static_assert(!is_instance_of_v<T2, T1<int>>);
static_assert(!is_instance_of_v<AliasT1, T3<double, double>>);

}  // namespace base
}  // namespace principia
