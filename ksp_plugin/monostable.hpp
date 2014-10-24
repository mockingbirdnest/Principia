#pragma once

namespace principia {
namespace ksp_plugin {

// A boolean which is constructed true and becomes false when |Flop| is called.
class Monostable {
 public:
  Monostable() = default;
  ~Monostable() = default;

  Monostable(Monostable const&) = delete;
  Monostable(Monostable&&) = delete;
  Monostable& operator=(Monostable const&) = delete;
  Monostable& operator=(Monostable&&) = delete;

  void Flop();
  operator bool() const;

 private:
  bool transient_ = true;
};

}  // namespace ksp_plugin
}  // namespace principia
