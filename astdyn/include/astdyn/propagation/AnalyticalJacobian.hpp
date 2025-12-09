/**
 * @file AnalyticalJacobian.hpp
 * @brief Analytical Jacobian for orbital dynamics
 * @author AstDyn Team
 * @date 2025-12-09
 * 
 * Provides analytical computation of ∂f/∂x for various force models,
 * which is much faster and more accurate than numerical differentiation.
 */

#ifndef ASTDYN_ANALYTICAL_JACOBIAN_HPP
#define ASTDYN_ANALYTICAL_JACOBIAN_HPP

#include <Eigen/Dense>

namespace astdyn::propagation {

/**
 * @brief Analytical Jacobian calculator for orbital dynamics
 */
class AnalyticalJacobian {
public:
    /**
     * @brief Compute Jacobian for 2-body problem
     * 
     * For dx/dt = f(x) where x = [r, v]ᵀ:
     * 
     *   f = [v, -μr/r³]ᵀ
     * 
     * Jacobian:
     *   J = [  0₃ₓ₃    I₃ₓ₃  ]
     *       [ ∂a/∂r   0₃ₓ₃  ]
     * 
     * where ∂a/∂r = -μ/r³ (I - 3rr^T/r²)
     * 
     * @param x State [r, v]ᵀ
     * @param mu Gravitational parameter
     * @return 6×6 Jacobian matrix
     */
    static Eigen::Matrix<double, 6, 6> two_body(
        const Eigen::Vector<double, 6>& x,
        double mu
    );
    
    /**
     * @brief Compute Jacobian for N-body problem
     * 
     * Includes perturbations from multiple bodies.
     * 
     * @param x State [r, v]ᵀ
     * @param mu_central Central body GM
     * @param perturber_positions Positions of perturbing bodies
     * @param perturber_mus GMs of perturbing bodies
     * @return 6×6 Jacobian matrix
     */
    static Eigen::Matrix<double, 6, 6> n_body(
        const Eigen::Vector<double, 6>& x,
        double mu_central,
        const std::vector<Eigen::Vector3d>& perturber_positions,
        const std::vector<double>& perturber_mus
    );
    
    /**
     * @brief Compute Jacobian for J2 perturbation
     * 
     * Adds oblateness effect (J2 term).
     * 
     * @param x State [r, v]ᵀ
     * @param mu Central body GM
     * @param J2 J2 coefficient
     * @param R_eq Equatorial radius
     * @return 6×6 Jacobian matrix
     */
    static Eigen::Matrix<double, 6, 6> with_j2(
        const Eigen::Vector<double, 6>& x,
        double mu,
        double J2,
        double R_eq
    );
    
private:
    /**
     * @brief Compute ∂a/∂r for point mass
     * 
     * For acceleration a = -μr/r³:
     * 
     *   ∂a/∂r = -μ/r³ (I - 3rr^T/r²)
     */
    static Eigen::Matrix3d acceleration_gradient(
        const Eigen::Vector3d& r,
        double mu
    );
};

} // namespace astdyn::propagation

#endif // ASTDYN_ANALYTICAL_JACOBIAN_HPP
