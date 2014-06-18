namespace principia {
namespace integrators {

template<typename Position, typename Momentum>
template<typename Scalar>
inline SymplecticIntegrator<Position, Momentum>::
DoublePrecision<Scalar>::DoublePrecision(Scalar const& value)
    : value(value) {}

}  // namespace integrators
}  // namespace principia
