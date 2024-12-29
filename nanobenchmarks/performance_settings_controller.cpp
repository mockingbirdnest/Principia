#include "nanobenchmarks/performance_settings_controller.hpp"

#include <iostream>
#include <print>
#include <utility>
#include <string_view>

#include <windows.h>
#include <powersetting.h>
#include <powrprof.h>

#include "absl/flags/flag.h"
#include "glog/logging.h"

ABSL_FLAG(bool,
          keep_perf_boost,
          false,
          "Whether to retain the processor performance boost mode during "
          "benchmark execution");

namespace principia {
namespace nanobenchmarks {
namespace _performance_settings_controller {
namespace internal {

class WindowsPerformanceSettingsController
    : public PerformanceSettingsController {
 public:
  WindowsPerformanceSettingsController();
  ~WindowsPerformanceSettingsController() override;

 private:
  void NotifyPowerSetting(ULONG type, PVOID setting);
  std::pair<DWORD, DWORD> ReadAndPrintPerfBoostModeACDC() const;
  static std::string_view PerfBoostModeToString(DWORD mode);

  DWORD perf_boost_mode_ac_;
  DWORD perf_boost_mode_dc_;
  GUID* active_power_scheme_;
  BYTE ac_line_status_;
  HPOWERNOTIFY power_setting_notification_;
  bool settings_changed_ = false;
};

not_null<std::unique_ptr<PerformanceSettingsController>>
PerformanceSettingsController::Make() {
  return make_not_null_unique<WindowsPerformanceSettingsController>();
}

WindowsPerformanceSettingsController::WindowsPerformanceSettingsController() {
  SYSTEM_POWER_STATUS power_status;
  CHECK(GetSystemPowerStatus(&power_status));
  ac_line_status_ = power_status.ACLineStatus;
  std::println(std::cout,
               "ACLineStatus={} ({})",
               ac_line_status_,
               ac_line_status_ == 0   ? "Offline"
               : ac_line_status_ == 1 ? "Online"
                                      : "Unknown");

  DEVICE_NOTIFY_SUBSCRIBE_PARAMETERS parameters{
      [](PVOID context, ULONG type, PVOID setting) -> ULONG {
        static_cast<WindowsPerformanceSettingsController*>(context)
            ->NotifyPowerSetting(type, setting);
        return ERROR_SUCCESS;
      },
      this};
  CHECK_EQ(PowerSettingRegisterNotification(&GUID_ACDC_POWER_SOURCE,
                                            DEVICE_NOTIFY_CALLBACK,
                                            &parameters,
                                            &power_setting_notification_),
           ERROR_SUCCESS);

  CHECK_EQ(PowerGetActiveScheme(nullptr, &active_power_scheme_), ERROR_SUCCESS);
  std::tie(perf_boost_mode_ac_, perf_boost_mode_dc_) =
      ReadAndPrintPerfBoostModeACDC();
  if (!absl::GetFlag(FLAGS_keep_perf_boost)) {
    std::println("Disabling perf boost mode…");
    std::println(R"(If interrupted, restore with
      POWERCFG /SETACVALUEINDEX SCHEME_CURRENT SUB_PROCESSOR PERFBOOSTMODE {}
      POWERCFG /SETDCVALUEINDEX SCHEME_CURRENT SUB_PROCESSOR PERFBOOSTMODE {})",
                 perf_boost_mode_ac_,
                 perf_boost_mode_dc_);
    std::println(R"(Check with
       POWERCFG /QUERY SCHEME_CURRENT SUB_PROCESSOR PERFBOOSTMODE)");
    CHECK_EQ(PowerWriteACValueIndex(nullptr,
                                    active_power_scheme_,
                                    &GUID_PROCESSOR_SETTINGS_SUBGROUP,
                                    &GUID_PROCESSOR_PERF_BOOST_MODE,
                                    PROCESSOR_PERF_BOOST_MODE_DISABLED),
             ERROR_SUCCESS);
    CHECK_EQ(PowerWriteDCValueIndex(nullptr,
                                    active_power_scheme_,
                                    &GUID_PROCESSOR_SETTINGS_SUBGROUP,
                                    &GUID_PROCESSOR_PERF_BOOST_MODE,
                                    PROCESSOR_PERF_BOOST_MODE_DISABLED),
             ERROR_SUCCESS);
    auto const [updated_perf_boost_mode_ac, updated_perf_boost_mode_dc] =
        ReadAndPrintPerfBoostModeACDC();
    CHECK_EQ(updated_perf_boost_mode_ac, PROCESSOR_PERF_BOOST_MODE_DISABLED);
    CHECK_EQ(updated_perf_boost_mode_dc, PROCESSOR_PERF_BOOST_MODE_DISABLED);
  }
 }

 WindowsPerformanceSettingsController::~WindowsPerformanceSettingsController() {
  CHECK(UnregisterPowerSettingNotification(power_setting_notification_));
  if (settings_changed_) {
    std::println("!!! Power settings changed during benchmarking.");
  }
  if (!absl::GetFlag(FLAGS_keep_perf_boost)) {
    std::println("Restoring perf boost mode…");
    CHECK_EQ(PowerWriteACValueIndex(nullptr,
                                    active_power_scheme_,
                                    &GUID_PROCESSOR_SETTINGS_SUBGROUP,
                                    &GUID_PROCESSOR_PERF_BOOST_MODE,
                                    perf_boost_mode_ac_),
             ERROR_SUCCESS);
    CHECK_EQ(PowerWriteDCValueIndex(nullptr,
                                    active_power_scheme_,
                                    &GUID_PROCESSOR_SETTINGS_SUBGROUP,
                                    &GUID_PROCESSOR_PERF_BOOST_MODE,
                                    perf_boost_mode_dc_),
             ERROR_SUCCESS);
    ReadAndPrintPerfBoostModeACDC();
  }
}

void WindowsPerformanceSettingsController::NotifyPowerSetting(DWORD type,
                                                              PVOID setting) {
  SYSTEM_POWER_STATUS power_status;
  CHECK(GetSystemPowerStatus(&power_status));
  if (ac_line_status_ != power_status.ACLineStatus) {
    settings_changed_ = true;
    ac_line_status_ = power_status.ACLineStatus;
    std::println(std::cout,
                 "!!! ACLineStatus changed to {} ({})",
                 ac_line_status_,
                 ac_line_status_ == 0   ? "Offline"
                 : ac_line_status_ == 1 ? "Online"
                                        : "Unknown");
  }
}

#define PRINCIPIA_PROCESSOR_PERF_BOOST_MODE_CASE(value) \
  case PROCESSOR_PERF_BOOST_MODE_##value:               \
    return #value
std::string_view WindowsPerformanceSettingsController::PerfBoostModeToString(
    DWORD mode) {
  switch (mode) {
    PRINCIPIA_PROCESSOR_PERF_BOOST_MODE_CASE(DISABLED);
    PRINCIPIA_PROCESSOR_PERF_BOOST_MODE_CASE(ENABLED);
    PRINCIPIA_PROCESSOR_PERF_BOOST_MODE_CASE(AGGRESSIVE);
    PRINCIPIA_PROCESSOR_PERF_BOOST_MODE_CASE(EFFICIENT_ENABLED);
    PRINCIPIA_PROCESSOR_PERF_BOOST_MODE_CASE(EFFICIENT_AGGRESSIVE);
    PRINCIPIA_PROCESSOR_PERF_BOOST_MODE_CASE(AGGRESSIVE_AT_GUARANTEED);
    PRINCIPIA_PROCESSOR_PERF_BOOST_MODE_CASE(EFFICIENT_AGGRESSIVE_AT_GUARANTEED);
    default:
      return "Unknown";
  }
};
#undef PRINCIPIA_PROCESSOR_PERF_BOOST_MODE_CASE

std::pair<DWORD, DWORD>
WindowsPerformanceSettingsController::ReadAndPrintPerfBoostModeACDC() const {
  DWORD ac;
  DWORD dc;
  DWORD perf_boost_mode_size = sizeof(ac);
  CHECK_EQ(PowerReadACValue(nullptr,
                            active_power_scheme_,
                            &GUID_PROCESSOR_SETTINGS_SUBGROUP,
                            &GUID_PROCESSOR_PERF_BOOST_MODE,
                            nullptr,
                            reinterpret_cast<LPBYTE>(&ac),
                            &perf_boost_mode_size),
           ERROR_SUCCESS)
      << perf_boost_mode_size;
  CHECK_EQ(PowerReadDCValue(nullptr,
                            active_power_scheme_,
                            &GUID_PROCESSOR_SETTINGS_SUBGROUP,
                            &GUID_PROCESSOR_PERF_BOOST_MODE,
                            nullptr,
                            reinterpret_cast<LPBYTE>(&dc),
                            &perf_boost_mode_size),
           ERROR_SUCCESS)
      << perf_boost_mode_size;
  std::println("PERF_BOOST_MODE AC={} ({})", ac, PerfBoostModeToString(ac));
  std::println("PERF_BOOST_MODE DC={} ({})", dc, PerfBoostModeToString(dc));
  return {ac, dc};
}

}  // namespace internal
}  // namespace _performance_settings_controller
}  // namespace nanobenchmarks
}  // namespace principia
