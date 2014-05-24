namespace principia {
namespace integrators {

Integrator::Integrator() {}

Integrator::~Integrator() {}

Integrator::Parameters::Parameters() 
    : p_error(nullptr),
      q_error(nullptr) {}

}  // namespace integrators
}  // namespace principia
