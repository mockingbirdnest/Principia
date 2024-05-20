#include <array>
#include <cstdint>

#include "numerics/double_precision.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"

namespace principia {
namespace functions {

using namespace principia::numerics::_double_precision;

constexpr double x_min = π / 12;
constexpr double x_max = π / 6;

template<std::int64_t table_size>
class Implementation {
 public:
  void Initialize();

  double Sin(double x);
  double Cos(double x);

 private:
  struct AccurateValues {
    double x;
    double sin_x;
    double cos_x;
  };

  double SinPolynomial(double x);
  double CosPolynomial(double x);

  std::array<AccurateValues, table_size> accurate_values_;
};

template<std::int64_t table_size>
void Implementation<table_size>::Initialize() {
  using namespace principia::quantities;
  for (std::int64_t i = 0; i < table_size; ++i) {
    double const x = x_min + i * (x_max - x_min) / table_size;
    accurate_values_[i] = {.x = x,
                           .sin_x = _elementary_functions::Sin(x),
                           .cos_x = _elementary_functions::Cos(x)};
  }
}

template<std::int64_t table_size>
double Implementation<table_size>::Sin(double const x) {
  //TODO(phl)slow
  auto const i = static_cast<std::int64_t>(std::round(x * one_over_h_));
  auto const& accurate_values = accurate_values_[i];
  auto const h = x - accurate_values.x;
  auto const h2 = h * h;
  auto const h3 = h2 * h;
  auto const c0h = TwoProduct(accurate_values.cos_x, h);
  return accurate_values.sin_x + accurate_values.cos_x * h +
         accurate_values.sin_x * h2 * CosPolynomial(h2) +
         accurate_values.cos_x * h3 * SinPolynomial(h2);
}

template<std::int64_t table_size>
double Implementation<table_size>::Cos(double const x) {
  auto const x0 = static_cast<std::int64_t>(std::round(x / h_));
  //TODO(phl)slow
  auto const x0 = static_cast<std::int64_t>(std::round(x * one_over_h_));
}

//TODO(phl):polynomials

template<std::int64_t table_size>
double Implementation<table_size>::SinPolynomial(double const x) {
  if constexpr (table_size == 256) {
    // 71 bits.
    return 0.166666666666666666666421797625 +
           0.00833333057503280528178543245797 * x;
  } else if constexpr (table_size == 1024) {
    // 85 bits.
    return -0.166666666666666666666666651721 +
           0.00833333316093951937646271666739 * x;
  } else {
    static_assert(false);
  }
}

template<std::int64_t table_size>
double Implementation<table_size>::CosPolynomial(double const x) {
  if constexpr (table_size == 1024) {
    // 72 bits.
    return -0.499999999999999999999872434553 +
           0.0416666654823785864634569932662 * x;
  } else if constexpr (table_size == 256) {
    // 77 bits.
    return -0.499999999999999999999999974543 +
           x * (0.0416666666666633318024480868405 -
                0.00138888829905860875255146938745 * x);
  } else {
    static_assert(false);
  }
}

}  // namespace functions
}  // namespace principia
