namespace principia {
namespace integrators {

template<typename Position, typename Momentum>
SymplecticIntegrator<Position, Momentum>::Parameters::Parameters()
    : p_error(nullptr),
      q_error(nullptr) {}

}  // namespace integrators
}  // namespace principia
