
#pragma once

template<typename T>
std::ostream& operator<<(std::ostream& out,
                         std::experimental::optional<T> const& optional) {
  if (optional) {
    return out << *optional;
  } else {
    return out << "nullopt";
  }
}
