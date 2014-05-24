namespace principia {
namespace integrators {

Integrator::Integrator() {}

Integrator::~Integrator() {}

Integrator::Parameters::Parameters() 
    : p_error(nullptr),
      q_error(nullptr),
      t_error(0) {}

}  // namespace integrators
}  // namespace principia
