namespace principia {
namespace integrators {

template<typename Position, typename Momentum>
template<typename Scalar>
inline SymplecticIntegrator<Position, Momentum>::ValueAndError<Scalar>::ValueAndError(
    Scalar const& value) : value(value) {}

}  // namespace integrators
}  // namespace principia
