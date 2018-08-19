
#pragma once

namespace principia {
namespace quantities {
namespace internal_generators {

// These structs have a |Type| member that is a |Quantity| suitable for
// the result of the operation applied to argument(s) of the |Quantity| types
// given as template parameter(s).

template<typename Q, int n>
struct ExponentiationGenerator;

// Only legal if |n| divides the dimensions of |Q|.
template<typename Q, int n, typename = void>
struct NthRootGenerator;

template<typename Left, typename Right>
struct ProductGenerator;

template<typename Left, typename Right>
struct QuotientGenerator;

// This struct has a |Type| member which is a tuple obtained by applying
// |Transform| to each element type in |Type| (which must be a tuple or an array
// or a pair).
template<typename Qs,
         template<typename> class Transform,
         typename = std::make_integer_sequence<int, std::tuple_size_v<Qs>>>
struct TupleGenerator;

}  // namespace internal_generators
}  // namespace quantities
}  // namespace principia

// Because of circular dependencies, this file doesn't include
// generators_body.hpp.  This will be done by quantities.hpp.
