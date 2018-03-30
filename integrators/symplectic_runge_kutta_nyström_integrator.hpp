
// The files containing the tree of of child classes of |Integrator| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
//  parent.
#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#include "integrators/integrators.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_SYMPLECTIC_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_
#define PRINCIPIA_INTEGRATORS_SYMPLECTIC_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_

#include "base/status.hpp"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace integrators {
namespace internal_symplectic_runge_kutta_nyström_integrator {

using base::not_null;
using base::Status;
using geometry::Instant;
using numerics::FixedVector;
using quantities::Time;

// This class solves ordinary differential equations of following forms using a
// symplectic Runge-Kutta-Nyström method:
// (1).  (q, p)′ = X(q, p, t), with X = A(q, p) + B(q, p, t) and known
//       evolutions exp hA and exp hB, where [B, [B, [B, A]]] = 0;
// (2).  The above case, where (exp hB)(q, p, t) = (q, p) + h B(q, p, t),
//       (exp hA)(q, p) = (q, p) + h A(q, p), and A and B are known;
// (3).  q″ = f(q, t), with known f = -M⁻¹ ∇q V(q, t);
// Only (3) is implemented at this time.
// See the documentation for a proof that (3) is a special case of (2), and thus
// of (1), and for an explanation of the relation to Hamiltonian mechanics.

// Each step of size h is computed using the composition of evolutions
//   exp(aᵣ₋₁ h A) exp(bᵣ₋₁ h B) ... exp(a₀ h A) exp(b₀ h B);
// the integrator thus is a composition method.  If the appropriate coefficients
// vanish, the above can be reformulated as either
//   exp(aᵣ₋₁ h A) exp(bᵣ₋₁ h B) ... exp(b₁ h B) exp(a₀ h A)  or
//   exp(bᵣ₋₁ h B) exp(aᵣ₋₂ h A) ... exp(a₀ h A) exp(b₀ h B).
// The former is called type ABA, the latter type BAB, following the conventions
// used in Blanes, Casas and Ros (2001),
// New Families of Symplectic Runge-Kutta-Nyström Integration Methods,
// http://www.gicas.uji.es/Fernando/Proceedings/2000NAA.pdf.
// In the implementation, we call |stages_| the integer r above.  The number of
// |evaluations| is r-1 in the ABA and BAB cases, and r otherwise.
// See the documentation for an explanation of how types ABA and BAB reduce the
// number of evaluations required, especially in cases (2) and (3).

// As above, we follow the most widespread notation, calling the position
// weights a and the momentum weights b.  Note that our indices are 0-based.
// The following notations appear in the literature:
//   (a, b) in most treatments of the subject;
//   (d, c) in Ruth, Yoshida, as well as Forest and Ruth;
//   (B, b) in Sofroniou and Spaletta;
//   (<unnamed>, B) in Okunbor and Skeel;
//   (<unnamed>, b) in Calvo and Sanz-Serna.
// Moreover, we follow the convention of Sofroniou and Spaletta (2002),
// Symplectic Methods for Separable Hamiltonian Systems, in calling c the
// nodes used for the time argument of the evolution of B, with
//   c₀ = 0, cᵢ = cᵢ₋₁ + aᵢ₋₁ for i > 0.
// The notation γ is used for these nodes in Calvo and Sanz-Serna.
// See the documentation for a description of the correspondence between
// these coefficients and those of a general Runge-Kutta-Nyström method.

template<typename Method, typename Position>
class SymplecticRungeKuttaNyströmIntegrator
    : public FixedStepSizeIntegrator<
                 SpecialSecondOrderDifferentialEquation<Position>> {
 public:
  using ODE = SpecialSecondOrderDifferentialEquation<Position>;
  using AppendState = typename Integrator<ODE>::AppendState;

  static constexpr auto order = Method::order;
  static constexpr auto time_reversible = Method::time_reversible;
  static constexpr auto evaluations = Method::evaluations;
  static constexpr auto composition = Method::composition;

  class Instance : public FixedStepSizeIntegrator<ODE>::Instance {
   public:
    Status Solve(Instant const& t_final) override;
    SymplecticRungeKuttaNyströmIntegrator const& integrator() const override;
    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;

   private:
    Instance(IntegrationProblem<ODE> const& problem,
             AppendState const& append_state,
             Time const& step,
             SymplecticRungeKuttaNyströmIntegrator const& integrator);

    SymplecticRungeKuttaNyströmIntegrator const& integrator_;
    friend class SymplecticRungeKuttaNyströmIntegrator;
  };

  SymplecticRungeKuttaNyströmIntegrator();

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      IntegrationProblem<ODE> const& problem,
      AppendState const& append_state,
      Time const& step) const override;

  void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> message) const override;

 private:
  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> ReadFromMessage(
      serialization::FixedStepSizeIntegratorInstance const& message,
      IntegrationProblem<ODE> const& problem,
      AppendState const& append_state,
      Time const& step) const override;

  static constexpr auto BA = serialization::FixedStepSizeIntegrator::BA;
  static constexpr auto ABA = serialization::FixedStepSizeIntegrator::ABA;
  static constexpr auto BAB = serialization::FixedStepSizeIntegrator::BAB;

  static constexpr auto stages_ = Method::stages;
  static constexpr auto a_ = Method::a;
  static constexpr auto b_ = Method::b;

  FixedVector<double, Method::stages> c_;
};

}  // namespace internal_symplectic_runge_kutta_nyström_integrator

template<typename Method, typename Position>
internal_symplectic_runge_kutta_nyström_integrator::
    SymplecticRungeKuttaNyströmIntegrator<Method, Position> const&
SymplecticRungeKuttaNyströmIntegrator();

template<typename Method,
         methods::SymplecticRungeKuttaNyström::CompositionMethod composition,
         typename Position>
internal_symplectic_runge_kutta_nyström_integrator::
    SymplecticRungeKuttaNyströmIntegrator<
        typename methods::AsSymplecticRungeKuttaNyström<Method,
                                                        composition>::Method,
        Position> const&
SymplecticRungeKuttaNyströmIntegrator();

}  // namespace integrators
}  // namespace principia

#include "integrators/symplectic_runge_kutta_nyström_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_SYMPLECTIC_RUNGE_KUTTA_NYSTRÖM_INTEGRATOR_HPP_
#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
