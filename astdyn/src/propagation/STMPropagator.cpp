/**
 * @file STMPropagator.cpp
 * @brief Implementation of State Transition Matrix propagator
 */

#include "astdyn/propagation/STMPropagator.hpp"
#include <cmath>

namespace astdyn::propagation {

STMPropagator::STMPropagator(
    std::unique_ptr<Integrator> integrator,
    ForceFunction force_func,
    JacobianFunction jac_func
)
    : integrator_(std::move(integrator))
    , force_func_(std::move(force_func))
    , jac_func_(std::move(jac_func))
{}

Eigen::Matrix<double, 6, 6> STMPropagator::numerical_jacobian(
    double t,
    const Eigen::Vector<double, 6>& x
) const {
    Eigen::Matrix<double, 6, 6> J;
    
    // Compute f(t, x)
    auto f0 = force_func_(t, x);
    
    // Finite difference for each component
    for (int j = 0; j < 6; ++j) {
        Eigen::Vector<double, 6> x_pert = x;
        double h = jac_epsilon_ * std::max(std::abs(x(j)), 1.0);
        x_pert(j) += h;
        
        auto f_pert = force_func_(t, x_pert);
        J.col(j) = (f_pert - f0) / h;
    }
    
    return J;
}

Eigen::VectorXd STMPropagator::combined_derivative(
    double t,
    const Eigen::VectorXd& y
) const {
    // Unpack state and STM
    Eigen::Vector<double, 6> x;
    Eigen::Matrix<double, 6, 6> Phi;
    unpack(y, x, Phi);
    
    // Compute state derivative: dx/dt = f(t, x)
    auto dxdt = force_func_(t, x);
    
    // Compute Jacobian: J = ∂f/∂x
    Eigen::Matrix<double, 6, 6> J;
    if (jac_func_) {
        J = jac_func_(t, x);
    } else {
        J = numerical_jacobian(t, x);
    }
    
    // Compute STM derivative: dΦ/dt = J Φ
    Eigen::Matrix<double, 6, 6> dPhidt = J * Phi;
    
    // Pack result
    return pack(dxdt, dPhidt);
}

Eigen::VectorXd STMPropagator::pack(
    const Eigen::Vector<double, 6>& state,
    const Eigen::Matrix<double, 6, 6>& stm
) {
    Eigen::VectorXd y(42);
    
    // First 6 components: state
    y.head<6>() = state;
    
    // Next 36 components: STM (column-major)
    for (int j = 0; j < 6; ++j) {
        for (int i = 0; i < 6; ++i) {
            y(6 + j * 6 + i) = stm(i, j);
        }
    }
    
    return y;
}

void STMPropagator::unpack(
    const Eigen::VectorXd& y,
    Eigen::Vector<double, 6>& state,
    Eigen::Matrix<double, 6, 6>& stm
) {
    // Extract state
    state = y.head<6>();
    
    // Extract STM
    for (int j = 0; j < 6; ++j) {
        for (int i = 0; i < 6; ++i) {
            stm(i, j) = y(6 + j * 6 + i);
        }
    }
}

STMPropagator::PropagationResult STMPropagator::propagate(
    const Eigen::Vector<double, 6>& x0,
    double t0,
    double tf,
    const Eigen::Matrix<double, 6, 6>& Phi0
) {
    // Pack initial conditions
    auto y0 = pack(x0, Phi0);
    
    // Create derivative function for integrator
    auto derivative = [this](double t, const Eigen::VectorXd& y) {
        return this->combined_derivative(t, y);
    };
    
    // Integrate
    auto y_final = integrator_->integrate(derivative, y0, t0, tf);
    
    // Unpack result
    PropagationResult result;
    unpack(y_final, result.state, result.stm);
    result.stats = integrator_->statistics();
    
    return result;
}

std::vector<STMPropagator::PropagationResult> STMPropagator::propagate_to_epochs(
    const Eigen::Vector<double, 6>& x0,
    double t0,
    const std::vector<double>& epochs
) {
    std::vector<PropagationResult> results;
    results.reserve(epochs.size());
    
    // Current state and STM
    auto x_current = x0;
    Eigen::Matrix<double, 6, 6> Phi_current;
    Phi_current.setIdentity();  // Initialize as identity matrix
    double t_current = t0;
    
    for (double t_target : epochs) {
        if (t_target < t_current) {
            throw std::invalid_argument("Epochs must be in increasing order");
        }
        
        // Propagate from current to target
        auto result = propagate(x_current, t_current, t_target, Phi_current);
        results.push_back(result);
        
        // Update current state
        x_current = result.state;
        Phi_current.noalias() = result.stm;  // Use noalias() for Eigen assignment
        t_current = t_target;
    }
    
    return results;
}

} // namespace astdyn::propagation
