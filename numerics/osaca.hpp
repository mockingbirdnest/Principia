#pragma once

#define PRINCIPIA_USE_OSACA 0

// The macros OSACA_FUNCTION_BEGIN and OSACA_RETURN are used to analyse the
// latency of a double -> double function, as measured, e.g., by the
// nanobenchmarks; this notionally corresponds to the duration of an iteration
// of a loop `x = f(x)`.
// The latency-critical path of the function is reported as the loop-carried
// dependency by OSACA, and as the critical path by IACA in throughput analysis
// mode.
// OSACA and IACA are only suitable for benchmarking a single path.  Any if
// statements, including in inline functions called by the function being
// analysed, should be replaced with OSACA_IF.  The path should be determined by
// providing any necessary constexpr expressions in UNDER_OSACA_HYPOTHESES.
// If OSACA_EVALUATE_CONDITIONS is set to 1, code will be generated to evaluate
// the conditions for the branches taken (outside the loop-carried path, as they
// would be with correct branch prediction).
// OSACA_EVALUATE_CONDITIONS can be set to 0 to get a more streamlined
// dependency graph when the conditions are irrelevant.  Note that this can
// result in artificially low port pressures.
#if 0
// Usage:
  double f(double const x) {
    double const abs_x = std::abs(x);
    if (abs_x < 0.1) {
      return x;
    } else if (x < 0) {
      return -f_positive(abs_x);
    } else {
      return f_positive(abs_x);
    }
  }
  double f_positive(double const abs_x) {
    if (abs_x > 3) {
      return 1;
    } else {
      // …
    }
  }
// becomes:
  double f(double x) {  // The argument cannot be const.
    OSACA_FUNCTION_BEGIN(x);
    double const abs_x = std::abs(x);
    OSACA_IF(abs_x < 0.1) {
      OSACA_RETURN(x);
    } OSACA_ELSE_IF(x < 0) {  // `else OSACA_IF` works, but upsets the linter.
      OSACA_RETURN(-f_positive(abs_x));
    } else {
      OSACA_RETURN(f_positive(abs_x));
    }
  }
  // Other functions can have const arguments and return normally, but need to
  // use OSACA_IF:
  double f_positive(double const abs_x) {
    OSACA_IF(abs_x > 3) {
      return 1;
    } else {
      // …
    }
  }
// To analyse it near x = 5:
#  define OSACA_ANALYSED_FUNCTION f
#  define OSACA_ANALYSED_FUNCTION_NAMESPACE
#  define OSACA_ANALYSED_FUNCTION_TEMPLATE_PARAMETERS
#  define UNDER_OSACA_HYPOTHESES(expression)                                 \
    [&] {                                                                    \
      constexpr double x = 5;                                                \
      /* To avoid inconsistent definitions, constexpr code can be used as */ \
      /* needed.                                                          */ \
      constexpr double abs_x = x > 0 ? x : -x;                               \
      return (expression);                                                   \
    }()

// If multiple functions principia::ns1::f and principia::ns2::f are marked up
// for analysis in a translation unit, use
#  define OSACA_ANALYSED_FUNCTION_NAMESPACE ns1::
// If f is templatized as, e.g, `template<RoundingMode mode>`, use
      OSACA_FUNCTION_BEGIN(x, <mode>)
// and
#  define OSACA_ANALYSED_FUNCTION_TEMPLATE_PARAMETERS \
    <RoundingMode::NearestTiesToEven>
#endif
// In principle, the loop-carried dependency for function latency analysis is
// achieved by wrapping the body of the function in an infinite loop, assigning
// the result to the argument at each iteration, and adding IACA markers outside
// the loop.  There are some subtleties:
// — We need to trick the compiler into believing the loop is finite, so that it
//   doesn’t optimize away the end marker or even the function.  This is
//   achieved by exiting based on the value of `OSACA_loop_terminator`.
// — Return statements may be in if statements, and there may be several of
//   them, so they cannot be the end of a loop started unconditionally.  Instead
//   we loop with goto.
// — We need to prevent the compiler from moving the start and end markers into
//   the middle of register saving and restoring code, which would mess up the
//   dependency analysis.  This is done with additional conditional gotos.
// — Some volatile reads and writes are used to clarify identity of the
//   registers in the generated code (where the names of `OSACA_result` and, if
//   `OSACA_CARRY_LOOP_THROUGH_REGISTER` is set to 0, `OSACA_loop_carry` appear
//   in movsd instructions).
//
// Carrying the loop dependency through a memory load and store can make the
// dependency graph easier to understand, as it forces any usage of the input to
// depend on the initial movsd, with the loop carried by a single backward edge
// to that initial movsd.
// If the loop is carried through a register, multiple usages of the input may
// result in multiple back edges from the final instruction that computed the
// result.  Carrying the loop through the memory could also potentially prevent
// the compiler from reusing intermediate values in the next iteration, e.g., if
// the the computation of f(x) depends on -x and produces -f(x) before f(x), as
// in an even function defined in terms of its positive half, the compiler might
// reuse -f(x₀)=-x₁ instead of computing -x₁ from x₁=f(x₀).  However:
// — it adds a spurious move to the latency;
// — some tools (IACA) cannot see the dependency through memory.
// Set OSACA_CARRY_LOOP_THROUGH_REGISTER to 0 to carry the loop through memory.

#define OSACA_EVALUATE_CONDITIONS 1
#define OSACA_CARRY_LOOP_THROUGH_REGISTER 1

#if PRINCIPIA_USE_OSACA

#include "intel/iacaMarks.h"

#define OSACA_QUALIFIED_ANALYSED_FUNCTION                   \
  OSACA_ANALYSED_FUNCTION_NAMESPACE OSACA_ANALYSED_FUNCTION \
      OSACA_ANALYSED_FUNCTION_TEMPLATE_PARAMETERS

static bool volatile OSACA_loop_terminator = false;

#define OSACA_FUNCTION_BEGIN(arg, ...)                              \
  double OSACA_LOOP_CARRY_QUALIFIER OSACA_loop_carry = arg;         \
  OSACA_outer_loop:                                                 \
  constexpr auto* OSACA_analysed_function_with_current_parameters = \
      &OSACA_ANALYSED_FUNCTION __VA_ARGS__;                         \
  if constexpr (std::string_view(__func__) ==                       \
                    STRINGIFY_EXPANSION(OSACA_ANALYSED_FUNCTION) && \
                OSACA_analysed_function_with_current_parameters ==  \
                    &OSACA_QUALIFIED_ANALYSED_FUNCTION) {           \
    IACA_VC64_START;                                                \
  }                                                                 \
  _Pragma("warning(push)");                                         \
  _Pragma("warning(disable : 4102)");                               \
  OSACA_loop:                                                       \
  _Pragma("warning(pop)");                                          \
  arg = OSACA_loop_carry

#define OSACA_RETURN(result)                                                \
  do {                                                                      \
    if constexpr (std::string_view(__func__) ==                             \
                      STRINGIFY_EXPANSION(OSACA_ANALYSED_FUNCTION) &&       \
                  OSACA_analysed_function_with_current_parameters ==        \
                      &OSACA_QUALIFIED_ANALYSED_FUNCTION) {                 \
      OSACA_loop_carry = (result);                                          \
      if (!OSACA_loop_terminator) {                                         \
        goto OSACA_loop;                                                    \
      }                                                                     \
      double volatile OSACA_result = OSACA_loop_carry;                      \
      IACA_VC64_END;                                                        \
      /* The outer loop prevents the the start and end marker from being */ \
      /* interleaved with register saving and restoring moves.           */ \
      if (!OSACA_loop_terminator) {                                         \
        goto OSACA_outer_loop;                                              \
      }                                                                     \
      return OSACA_result;                                                  \
    } else {                                                                \
      return (result);                                                      \
    }                                                                       \
  } while (false)

#if OSACA_CARRY_LOOP_THROUGH_REGISTER
#define OSACA_LOOP_CARRY_QUALIFIER
#else
#define OSACA_LOOP_CARRY_QUALIFIER volatile
#endif

// The branch not taken, determined by evaluating the condition
// `UNDER_OSACA_HYPOTHESES`, is eliminated by `if constexpr`; the condition is
// also compiled normally.  Whether this results in any generated code depends
// on `OSACA_EVALUATE_CONDITIONS`.  Note that, with `OSACA_EVALUATE_CONDITIONS`,
// in `OSACA_IF(p) { } OSACA_ELSE_IF(q) { }`, if  `p` holds
// `UNDER_OSACA_HYPOTHESES`, code is generated to evaluate `p`, but  not `q`.

#define OSACA_IF(condition)                                             \
  if constexpr (int const OSACA_evaluate =                              \
                    UNDER_OSACA_HYPOTHESES(condition)                   \
                        ? ((condition) ? 0 : OSACA_branch_not_taken())  \
                        : ((condition) ? OSACA_branch_not_taken() : 0); \
                UNDER_OSACA_HYPOTHESES(condition))

#if OSACA_EVALUATE_CONDITIONS
[[noreturn]] inline int OSACA_branch_not_taken() {
  std::abort();
}
#else
inline int OSACA_branch_not_taken() {
  return 0;
}
#endif

#else  // if !PRINCIPIA_USE_OSACA

#define OSACA_FUNCTION_BEGIN(...) ((void) 0)
#define OSACA_RETURN(result) return (result)
#define OSACA_IF(condition) if (condition)

#endif  // PRINCIPIA_USE_OSACA

#define OSACA_ELSE_IF else OSACA_IF  // NOLINT
