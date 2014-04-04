#pragma once

namespace Principia {
namespace Geometry {
template<typename Scalar, typename Frame, unsigned int Rank>
class Multivector;

template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator+ (Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator- (Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>

Multivector<T, Frame,
            Rank> operator+ (Multivector<T, Frame, Rank> const& left,
                             Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator- (Multivector<T, Frame, Rank> const& left,
                             Multivector<T, Frame, Rank> const& right);

template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator* (Quantities::Dimensionless const& left,
                             Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator* (Multivector<T, Frame, Rank> const& left,
                             Quantities::Dimensionless const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator/ (Multivector<T, Frame, Rank> const& left,
                             Quantities::Dimensionless const& right);
}
}
