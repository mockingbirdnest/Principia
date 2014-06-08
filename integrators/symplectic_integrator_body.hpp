namespace principia {
namespace integrators {

template<typename Position, typename Momentum>
SymplecticIntegrator<Position, Momentum>::Parameters::Parameters()
    : p_error(nullptr),
      q_error(nullptr),
      t_error(0) {}

}  // namespace integrators
}  // namespace principia
