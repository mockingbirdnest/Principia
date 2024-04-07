#pragma once
#include <type_traits>

#define PRINCIPIA_EXPAND(...) __VA_ARGS__

#define PRINCIPIA_CONCATENATE(s1, s2) s1##s2
#define PRINCIPIA_CONCATENATE_EXPANSION(s1, s2) PRINCIPIA_CONCATENATE(s1, s2)

#define PRINCIPIA_COMMA_IF_NONEMPTY(...) __VA_OPT__(,)  // NOLINT

#define PRINCIPIA_CHECK_WELL_FORMEDNESS(unique_concept_name,                  \
                                        expected_well_formedness,             \
                                        expected_well_formedness_description, \
                                        expression,                           \
                                        types,                                \
                                        ...)                                  \
  template<template<typename> typename with_variable = std::type_identity_t   \
               PRINCIPIA_COMMA_IF_NONEMPTY types PRINCIPIA_EXPAND types>      \
  concept unique_concept_name = requires(__VA_ARGS__) { (expression); };      \
  static_assert(expected_well_formedness unique_concept_name<>,               \
                "Expected\n  " #expression                                    \
                "\n" expected_well_formedness_description)
#define PRINCIPIA_CHECK_ILL_FORMED(expression, ...)                \
  PRINCIPIA_CHECK_WELL_FORMEDNESS(                                 \
      PRINCIPIA_CONCATENATE_EXPANSION(test_concept_, __COUNTER__), \
      !,                                                           \
      "not to compile",                                            \
      expression,                                                  \
      (),                                                          \
      __VA_ARGS__)
#define PRINCIPIA_CHECK_WELL_FORMED(expression, ...)               \
  PRINCIPIA_CHECK_WELL_FORMEDNESS(                                 \
      PRINCIPIA_CONCATENATE_EXPANSION(test_concept_, __COUNTER__), \
      ,                                                            \
      "to compile",                                                \
      expression,                                                  \
      (),                                                          \
      __VA_ARGS__)
#define PRINCIPIA_CHECK_ILL_FORMED_WITH_TYPES(expression, types, ...) \
  PRINCIPIA_CHECK_WELL_FORMEDNESS(                                    \
      PRINCIPIA_CONCATENATE_EXPANSION(test_concept_, __COUNTER__),    \
      !,                                                              \
      "not to compile",                                               \
      expression,                                                     \
      types,                                                          \
      __VA_ARGS__)
#define PRINCIPIA_CHECK_WELL_FORMED_WITH_TYPES(expression, types, ...) \
  PRINCIPIA_CHECK_WELL_FORMEDNESS(                                     \
      PRINCIPIA_CONCATENATE_EXPANSION(test_concept_, __COUNTER__),     \
      ,                                                                \
      "to compile",                                                    \
      expression,                                                      \
      types,                                                           \
      __VA_ARGS__)
