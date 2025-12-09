/**
 * @file AnalyticalJacobian.cpp
 * @brief Implementation of analytical Jacobian
 */

#include "astdyn/propagation/AnalyticalJacobian.hpp"
#include <cmath>

namespace astdyn::propagation {

Eigen::Matrix3d AnalyticalJacobian::acceleration_gradient(
    const Eigen::Vector3d& r,
    double mu
) {
    double r_norm = r.norm();
    double r2 = r_norm * r_norm;
    double r3 = r_norm * r2;  // FIX: was r_norm * r3 (circular!)
    double r5 = r2 * r3;
    
    // ∂a/∂r = -μ/r³ (I - 3rr^T/r²)
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d rr_T = r * r.transpose();
    
    return -mu / r3 * (I - 3.0 * rr_T / r2);
}

Eigen::Matrix<double, 6, 6> AnalyticalJacobian::two_body(
    const Eigen::Vector<double, 6>& x,
    double mu
) {
    Eigen::Vector3d r = x.head<3>();
    
    Eigen::Matrix<double, 6, 6> J = Eigen::Matrix<double, 6, 6>::Zero();
    
    // Upper-right block: ∂v/∂v = I
    J.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
    
    // Lower-left block: ∂a/∂r
    J.block<3, 3>(3, 0) = acceleration_gradient(r, mu);
    
    // Other blocks are zero
    
    return J;
}

Eigen::Matrix<double, 6, 6> AnalyticalJacobian::n_body(
    const Eigen::Vector<double, 6>& x,
    double mu_central,
    const std::vector<Eigen::Vector3d>& perturber_positions,
    const std::vector<double>& perturber_mus
) {
    Eigen::Vector3d r = x.head<3>();
    
    // Start with 2-body Jacobian
    Eigen::Matrix<double, 6, 6> J = two_body(x, mu_central);
    
    // Add perturbations
    for (size_t i = 0; i < perturber_positions.size(); ++i) {
        Eigen::Vector3d r_pert = perturber_positions[i];
        double mu_pert = perturber_mus[i];
        
        // Relative position: asteroid - perturber
        Eigen::Vector3d d = r - r_pert;
        
        // Direct term: ∂/∂r[-μ_p d/d³]
        J.block<3, 3>(3, 0) += acceleration_gradient(d, mu_pert);
        
        // Indirect term is constant (∂/∂r[μ_p r_p/r_p³] = 0)
    }
    
    return J;
}

Eigen::Matrix<double, 6, 6> AnalyticalJacobian::with_j2(
    const Eigen::Vector<double, 6>& x,
    double mu,
    double J2,
    double R_eq
) {
    Eigen::Vector3d r = x.head<3>();
    
    // Start with 2-body
    Eigen::Matrix<double, 6, 6> J = two_body(x, mu);
    
    // Add J2 perturbation gradient
    double r_norm = r.norm();
    double r2 = r_norm * r_norm;
    double z2 = r(2) * r(2);
    
    double factor = 1.5 * J2 * mu * R_eq * R_eq / std::pow(r_norm, 5);
    
    // Simplified J2 gradient (full derivation is complex)
    // This is an approximation - full version requires tensor calculus
    Eigen::Matrix3d J2_grad = Eigen::Matrix3d::Zero();
    
    // Diagonal terms
    J2_grad(0, 0) = factor * (1.0 - 5.0 * z2 / r2);
    J2_grad(1, 1) = factor * (1.0 - 5.0 * z2 / r2);
    J2_grad(2, 2) = factor * (3.0 - 5.0 * z2 / r2);
    
    // Off-diagonal terms
    J2_grad(0, 2) = J2_grad(2, 0) = -5.0 * factor * r(0) * r(2) / r2;
    J2_grad(1, 2) = J2_grad(2, 1) = -5.0 * factor * r(1) * r(2) / r2;
    
    J.block<3, 3>(3, 0) += J2_grad;
    
    return J;
}

} // namespace astdyn::propagation
