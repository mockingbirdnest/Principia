#pragma once

namespace Principia {
namespace Geometry {
template<typename Scalar, typename Space, unsigned int Rank>
class Multivector;

template<typename T, typename Space, unsigned int Rank>
Multivector<T, Space,
            Rank> operator+ (Multivector<T, Space, Rank> const& right);
template<typename T, typename Space, unsigned int Rank>
Multivector<T, Space,
            Rank> operator- (Multivector<T, Space, Rank> const& right);
template<typename T, typename Space, unsigned int Rank>

Multivector<T, Space,
            Rank> operator+ (Multivector<T, Space, Rank> const& left,
                             Multivector<T, Space, Rank> const& right);
template<typename T, typename Space, unsigned int Rank>
Multivector<T, Space,
            Rank> operator- (Multivector<T, Space, Rank> const& left,
                             Multivector<T, Space, Rank> const& right);

template<typename T, typename Space, unsigned int Rank>
Multivector<T, Space,
            Rank> operator* (Quantities::Dimensionless const& left,
                             Multivector<T, Space, Rank> const& right);
template<typename T, typename Space, unsigned int Rank>
Multivector<T, Space,
            Rank> operator* (Multivector<T, Space, Rank> const& left,
                             Quantities::Dimensionless const& right);
template<typename T, typename Space, unsigned int Rank>
Multivector<T, Space,
            Rank> operator/ (Multivector<T, Space, Rank> const& left,
                             Quantities::Dimensionless const& right);

template<typename T, typename Space, unsigned int Rank, typename U>
R3Element<Quantities::Product<U, T>> operator* (U const& left,
                                                Multivector<T, Space,
                                                            Rank> const& right);
template<typename T, typename Space, unsigned int Rank, typename U>
R3Element<Quantities::Product<T, U>> operator* (Multivector<T, Space,
                                                            Rank> const& left,
                                                U const& right);
template<typename T, typename Space, unsigned int Rank, typename U>
R3Element<Quantities::Quotient<T, U>> operator/ (Multivector<T, Space,
                                                             Rank> const& left,
                                                 U const& right);

template<typename T, typename Space, unsigned int Rank>
Multivector<T, Space,
            Rank> operator+= (Multivector<T, Space, Rank>& left,
                              Multivector<T, Space, Rank> const& right);
template<typename T, typename Space, unsigned int Rank>
Multivector<T, Space,
            Rank> operator-= (Multivector<T, Space, Rank>& left,
                              Multivector<T, Space, Rank> const& right);

template<typename T, typename Space, unsigned int Rank>
Multivector<T, Space,
            Rank> operator*= (Multivector<T, Space, Rank>& left,
                              Quantities::Dimensionless const& right);
template<typename T, typename Space, unsigned int Rank>
Multivector<T, Space,
            Rank> operator/= (Multivector<T, Space, Rank>& left,
                              Quantities::Dimensionless const& right);

}
}
