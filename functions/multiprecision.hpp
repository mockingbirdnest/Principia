#pragma once

#include "boost/multiprecision/gmp.hpp"

namespace principia {
namespace functions {
namespace _multiprecision {
namespace internal {

using namespace boost::multiprecision;

gmp_float<20> Sin(mpq_rational const& angle);

}  // namespace internal

using internal::Sin;

}  // namespace _multiprecision
}  // namespace functions
}  // namespace principia
