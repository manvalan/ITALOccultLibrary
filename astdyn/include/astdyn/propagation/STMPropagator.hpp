/**
 * @file STMPropagator.hpp
 * @brief State Transition Matrix propagator for orbit determination
 * @author AstDyn Team
 * @date 2025-12-09
 * 
 * Propagates both the state vector and the State Transition Matrix (STM)
 * simultaneously, enabling sensitivity analysis and orbit determination.
 * 
 * The STM Φ(t,t₀) describes how small variations in initial conditions
 * propagate forward in time:
 * 
 *   δx(t) = Φ(t,t₀) δx(t₀)
 * 
 * where x = [r, v]ᵀ is the 6D state vector.
 */

#ifndef ASTDYN_STM_PROPAGATOR_HPP
#define ASTDYN_STM_PROPAGATOR_HPP

#include "Integrator.hpp"
#include <Eigen/Dense>
#include <functional>

namespace astdyn::propagation {

/**
 * @brief State Transition Matrix propagator
 * 
 * Propagates state (6D) and STM (6×6 = 36D) simultaneously.
 * Total system dimension: 42D
 */
class STMPropagator {
public:
    /**
     * @brief Force function signature: f(t, x) → dx/dt
     */
    using ForceFunction = std::function<Eigen::Vector<double, 6>(double, const Eigen::Vector<double, 6>&)>;
    
    /**
     * @brief Jacobian function signature: J(t, x) → ∂f/∂x
     * If not provided, computed numerically
     */
    using JacobianFunction = std::function<Eigen::Matrix<double, 6, 6>(double, const Eigen::Vector<double, 6>&)>;
    
    /**
     * @brief Construct STM propagator
     * 
     * @param integrator Numerical integrator (RKF78 recommended)
     * @param force_func Force function f(t, x)
     * @param jac_func Jacobian function (optional, computed numerically if nullptr)
     */
    explicit STMPropagator(
        std::unique_ptr<Integrator> integrator,
        ForceFunction force_func,
        JacobianFunction jac_func = nullptr
    );
    
    /**
     * @brief Propagate state and STM
     * 
     * @param x0 Initial state [r, v]ᵀ (6D)
     * @param t0 Initial time
     * @param tf Final time
     * @param Phi0 Initial STM (default: identity)
     * @return Final state (6D) and STM (6×6)
     */
    struct PropagationResult {
        Eigen::Vector<double, 6> state;
        Eigen::Matrix<double, 6, 6> stm;
        IntegrationStatistics stats;
    };
    
    PropagationResult propagate(
        const Eigen::Vector<double, 6>& x0,
        double t0,
        double tf,
        const Eigen::Matrix<double, 6, 6>& Phi0 = Eigen::Matrix<double, 6, 6>::Identity()
    );
    
    /**
     * @brief Propagate to multiple epochs
     * 
     * Useful for computing residuals at observation times
     */
    std::vector<PropagationResult> propagate_to_epochs(
        const Eigen::Vector<double, 6>& x0,
        double t0,
        const std::vector<double>& epochs
    );
    
    /**
     * @brief Set numerical Jacobian epsilon
     */
    void set_jacobian_epsilon(double eps) { jac_epsilon_ = eps; }
    
    /**
     * @brief Get integration statistics
     */
    const IntegrationStatistics& statistics() const { return stats_; }
    
private:
    std::unique_ptr<Integrator> integrator_;
    ForceFunction force_func_;
    JacobianFunction jac_func_;
    double jac_epsilon_ = 1e-8;
    IntegrationStatistics stats_;
    
    /**
     * @brief Compute numerical Jacobian ∂f/∂x
     */
    Eigen::Matrix<double, 6, 6> numerical_jacobian(
        double t,
        const Eigen::Vector<double, 6>& x
    ) const;
    
    /**
     * @brief Combined derivative function for [state, STM]
     * 
     * Computes d/dt[x, vec(Φ)] where vec() stacks columns
     * 
     * Equations:
     *   dx/dt = f(t, x)
     *   dΦ/dt = J(t, x) Φ
     */
    Eigen::VectorXd combined_derivative(
        double t,
        const Eigen::VectorXd& y
    ) const;
    
    /**
     * @brief Pack state and STM into single vector
     */
    static Eigen::VectorXd pack(
        const Eigen::Vector<double, 6>& state,
        const Eigen::Matrix<double, 6, 6>& stm
    );
    
    /**
     * @brief Unpack state and STM from single vector
     */
    static void unpack(
        const Eigen::VectorXd& y,
        Eigen::Vector<double, 6>& state,
        Eigen::Matrix<double, 6, 6>& stm
    );
};

} // namespace astdyn::propagation

#endif // ASTDYN_STM_PROPAGATOR_HPP
