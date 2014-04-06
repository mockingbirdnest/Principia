#pragma once

namespace Principia {
namespace Geometry {

template<typename T, typename Frame, unsigned int Rank>
inline Multivector<T, Frame,
                   Rank> operator+ (Multivector<T, Frame, Rank> const& right) {
  return Multivector<T, Frame, Rank>(+right.coordinates);
}
template<typename T, typename Frame, unsigned int Rank>
inline Multivector<T, Frame,
                   Rank> operator- (Multivector<T, Frame, Rank> const& right) {
  return Multivector<T, Frame, Rank>(-right.coordinates);
}
template<typename T, typename Frame, unsigned int Rank>

inline Multivector<T, Frame,
                   Rank> operator+ (Multivector<T, Frame, Rank> const& left,
                                    Multivector<T, Frame, Rank> const& right) {
  return Multivector<T, Frame, Rank>(left.coordinates + right.coordinates);
}
template<typename T, typename Frame, unsigned int Rank>
inline Multivector<T, Frame,
                   Rank> operator- (Multivector<T, Frame, Rank> const& left,
                                    Multivector<T, Frame, Rank> const& right) {
  return Multivector<T, Frame, Rank>(left.coordinates + left.coordinates);
}

}
}
