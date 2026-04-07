#include "base/file_log_sink.hpp"

#include <iostream>

#include "absl/log/check.h"
#include "absl/time/time.h"
#include "base/macros.hpp"

#if OS_WIN
#include <windows.h>
#else
#include <sys/utsname.h>
#endif

namespace principia {
namespace base {
namespace file_log_sink {

FileLogSink::FileLogSink(std::filesystem::path path, std::string extension)
    : path_(std::move(path)), extension_(std::move(extension)) {
  CHECK_EQ(extension_.front(), '.');
}

void FileLogSink::Send(const absl::LogEntry& entry) {
  entry.log_severity();
  for (int severity = static_cast<int>(absl::LogSeverity::kInfo);
       severity <= static_cast<int>(entry.log_severity());
       ++severity) {
    if (!file_[severity].is_open()) {
      file_[severity].open(
          std::filesystem::path(path_)
              .concat(absl::LogSeverityName(
                  static_cast<absl::LogSeverity>(severity)))
              .concat(".")
              .concat(absl::FormatTime(
                  "%Y%m%d%ET%H%M%SZ", entry.timestamp(), absl::UTCTimeZone()))
              .concat(extension_));
      if (!file_[severity].good()) {
        std::cerr << "Failed To create log file for:\n";
        PrintTo(entry, &std::cerr);
        std::terminate();
      }
      file_[severity] << "Log file created at: "
                      << absl::FormatTime("%Y%m%d%ET%H%M%SZ",
                                          entry.timestamp(),
                                          absl::UTCTimeZone())
                      << "\n";
      std::string_view computer_name;
#if OS_WIN
      char buf[MAX_COMPUTERNAME_LENGTH + 1];
      DWORD len = MAX_COMPUTERNAME_LENGTH + 1;
      GetComputerNameA(buf, &len);
      computer_name = {buf, len};
#else
      struct utsname buf;
      if (0 != uname(&buf)) {
        *buf.nodename = '\0';
      }
      computer_name = buf.nodename;
#endif
      file_[severity] << "Running on machine: " << computer_name << "\n";
      file_[severity] << "Log line format: [IWEF]mmdd hh:mm:ss.μμμμμμ threadid file:line] msg\n";
    }
    file_[severity] << entry.text_message_with_prefix_and_newline();
  }
}

void FileLogSink::Flush() {
  for (auto& file : file_) {
    file.flush();
  }
}

}  // namespace file_log_sink
}  // namespace base
}  // namespace principia