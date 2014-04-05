#pragma once

#include "..\Quantities\Dimensionless.hpp"
#include "..\Quantities\Quantities.hpp"

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

template<typename T, typename U, typename Frame, unsigned int Rank>
Multivector<Quantities::Product<U, T> , Frame,
            Rank> operator* (U const& left,
                             Multivector<T, Frame, Rank> const& right);
template<typename T, template U, typename Frame, unsigned int Rank>
Multivector<Quantities::Product<T, U>, Frame,
            Rank> operator* (Multivector<T, Frame, Rank> const& left,
                             Quantities::Dimensionless const& right);
template<typename T, typename U, typename Frame, unsigned int Rank>
Multivector<Quotient<T, U>, Frame,
            Rank> operator/ (Multivector<T, Frame, Rank> const& left,
                             Quantities::Dimensionless const& right);
}
}
