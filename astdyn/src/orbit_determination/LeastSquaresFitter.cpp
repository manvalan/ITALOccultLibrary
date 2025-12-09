/**
 * @file LeastSquaresFitter.cpp
 * @brief Implementation of least squares orbit fitter
 */

#include "astdyn/orbit_determination/LeastSquaresFitter.hpp"
#include <cmath>
#include <algorithm>

namespace astdyn::orbit_determination {

LeastSquaresFitter::LeastSquaresFitter() {}

Eigen::MatrixXd LeastSquaresFitter::build_design_matrix(
    const std::vector<ObservationResidual>& residuals,
    const Eigen::Vector<double, 6>& state,
    double epoch_mjd,
    STMFunction stm_func
) {
    // Design matrix A: each observation contributes 2 rows (RA, Dec)
    int n_obs = residuals.size();
    Eigen::MatrixXd A(2 * n_obs, 6);
    
    // For each observation, compute: A_i = (∂ρ/∂x) × Φ(t_i, t_0)
    // where ∂ρ/∂x is partial of RA/Dec w.r.t. state at obs time
    // and Φ is STM from epoch to obs time
    
    constexpr double DEG_TO_RAD = 1.7453292519943295769e-2;
    // Obliquity J2000
    constexpr double eps = 23.4392911 * DEG_TO_RAD;
    double ce = std::cos(eps);
    double se = std::sin(eps);
    
    // Rotation Matrix Ecliptic -> Equatorial
    Eigen::Matrix3d R_ecl_to_eq;
    R_ecl_to_eq << 1, 0, 0,
                   0, ce, -se,
                   0, se, ce;

    for (int i = 0; i < n_obs; ++i) {
        double t_obs = residuals[i].epoch_mjd;
        
        // Get STM from epoch to observation time (In Ecliptic Frame)
        auto [state_at_obs_ecl, stm_ecl] = stm_func(state, epoch_mjd, t_obs);
        
        // We need partial derivatives of RA/Dec (Equatorial) w.r.t State (Ecliptic).
        // By Chain Rule: d(RA)/dr_ecl = d(RA)/dr_eq * d(r_eq)/d(r_ecl)
        // d(r_eq)/d(r_ecl) is the Rotation Matrix R
        
        // 1. Convert State position to Hessian Equatorial Frame for derivative calculation
        // Note: We ignore the Observer position shift here for the partials matrix 
        // because d(r_topo)/d(r_ast) = I (assuming r_obs is constant w.r.t asteroid state)
        // So d(RA)/d(r_ast_eq) is calculated using r_topo_eq.
        // Ideally we should use the EXACT topocentric vector here, but using the heliocentric 
        // rotated vector is a good approximation if Earth is far, 
        // OR better: use the State from STM but rotated.
        
        // Strictly, Jacobian H = H_geometric * R * STM
        // where H_geometric is d(RA,Dec)/d(r_eq) evaluated at r_topo_eq.
        
        // But we don't have r_topo_eq here easily (it's in ResidualCalculator).
        // However, state_at_obs is r_ast_ecl.
        // Let's approximate by evaluating derivatives at r_ast_eq (Heliocentric).
        // For partial derivatives direction this is usually sufficient for convergence.
        // For high precision convergence, we should pass r_topo from ResidualCalculator.
        // But let's try with rotated state first.
        
        Eigen::Vector3d r_ast_ecl = state_at_obs_ecl.head<3>();
        Eigen::Vector3d r_eq = R_ecl_to_eq * r_ast_ecl; // Use this for calculating d/dr_eq
        
        // Compute partials d/dr_eq using r_eq
        double x = r_eq(0), y = r_eq(1), z = r_eq(2);
        double r_norm = r_eq.norm();
        double rho_xy = std::sqrt(x*x + y*y);
        
        // ∂RA/∂r_eq (in radians)
        Eigen::Vector3d dRA_dr_eq;
        if (rho_xy > 1e-10) {
            dRA_dr_eq(0) = -y / (rho_xy * rho_xy);
            dRA_dr_eq(1) = x / (rho_xy * rho_xy);
            dRA_dr_eq(2) = 0.0;
        } else {
            dRA_dr_eq.setZero();
        }
        
        // ∂Dec/∂r_eq (in radians)
        Eigen::Vector3d dDec_dr_eq;
        if (r_norm > 1e-10) {
            dDec_dr_eq(0) = -x * z / (r_norm * r_norm * rho_xy);
            dDec_dr_eq(1) = -y * z / (r_norm * r_norm * rho_xy);
            dDec_dr_eq(2) = rho_xy / (r_norm * r_norm);
        } else {
            dDec_dr_eq.setZero();
        }
        
        // Chain Rule: d/dr_ecl = d/dr_eq * R
        // Matrix form: [dRA/dx_ecl dRA/dy_ecl ...] = [dRA/dx_eq ...] * R
        Eigen::Vector3d dRA_dr_ecl = dRA_dr_eq.transpose() * R_ecl_to_eq;
        Eigen::Vector3d dDec_dr_ecl = dDec_dr_eq.transpose() * R_ecl_to_eq;
        
        // Convert to arcsec (1 rad = 206265 arcsec)
        constexpr double rad_to_arcsec = 206265.0;
        dRA_dr_ecl *= rad_to_arcsec;
        dDec_dr_ecl *= rad_to_arcsec;
        
        // ∂ρ/∂x_ecl at obs time
        Eigen::Matrix<double, 2, 6> dRho_dx_ecl;
        dRho_dx_ecl.row(0) << dRA_dr_ecl.transpose(), 0, 0, 0;   // ∂RA/∂x
        dRho_dx_ecl.row(1) << dDec_dr_ecl.transpose(), 0, 0, 0;  // ∂Dec/∂x
        
        // Design matrix: A = (∂ρ/∂x_ecl) × Φ_ecl
        Eigen::Matrix<double, 2, 6> A_i = dRho_dx_ecl * stm_ecl;
        
        A.row(2*i) = A_i.row(0);      // ∂RA/∂x₀
        A.row(2*i+1) = A_i.row(1);    // ∂Dec/∂x₀
    }
    
    return A;
}

Eigen::Vector<double, 6> LeastSquaresFitter::solve_normal_equations(
    const Eigen::MatrixXd& A,
    const Eigen::VectorXd& residuals,
    const Eigen::VectorXd& weights,
    Eigen::Matrix<double, 6, 6>& covariance
) {
    // Weight matrix W = diag(weights)
    Eigen::MatrixXd W = weights.asDiagonal();
    
    // Normal equations: (A^T W A) δx = A^T W Δρ
    Eigen::Matrix<double, 6, 6> N = A.transpose() * W * A;
    Eigen::Vector<double, 6> b = A.transpose() * W * residuals;
    
    // Solve using LU decomposition
    Eigen::Vector<double, 6> dx = N.ldlt().solve(b);
    
    // Covariance: (A^T W A)^{-1}
    covariance = N.inverse();
    
    return dx;
}

int LeastSquaresFitter::reject_outliers(std::vector<ObservationResidual>& residuals) {
    if (!outlier_rejection_) return 0;
    
    // Compute RMS
    double sum_sq = 0.0;
    int count = 0;
    
    for (const auto& res : residuals) {
        if (!res.rejected) {
            sum_sq += res.ra_residual_arcsec * res.ra_residual_arcsec;
            sum_sq += res.dec_residual_arcsec * res.dec_residual_arcsec;
            count += 2;
        }
    }
    
    double rms = std::sqrt(sum_sq / count);
    double threshold = outlier_threshold_ * rms;
    
    // Reject outliers
    int num_rejected = 0;
    for (auto& res : residuals) {
        if (!res.rejected) {
            double res_norm = std::sqrt(
                res.ra_residual_arcsec * res.ra_residual_arcsec +
                res.dec_residual_arcsec * res.dec_residual_arcsec
            );
            
            if (res_norm > threshold) {
                res.rejected = true;
                num_rejected++;
            }
        }
    }
    
    return num_rejected;
}

void LeastSquaresFitter::compute_statistics(
    const std::vector<ObservationResidual>& residuals,
    FitResult& result
) {
    double sum_ra_sq = 0.0;
    double sum_dec_sq = 0.0;
    int count = 0;
    
    for (const auto& res : residuals) {
        if (!res.rejected) {
            sum_ra_sq += res.ra_residual_arcsec * res.ra_residual_arcsec;
            sum_dec_sq += res.dec_residual_arcsec * res.dec_residual_arcsec;
            count++;
        }
    }
    
    result.num_observations = residuals.size();
    result.num_rejected = std::count_if(residuals.begin(), residuals.end(),
                                        [](const auto& r) { return r.rejected; });
    
    if (count > 0) {
        result.rms_ra_arcsec = std::sqrt(sum_ra_sq / count);
        result.rms_dec_arcsec = std::sqrt(sum_dec_sq / count);
        result.rms_total_arcsec = std::sqrt((sum_ra_sq + sum_dec_sq) / (2 * count));
    } else {
        result.rms_ra_arcsec = 0.0;
        result.rms_dec_arcsec = 0.0;
        result.rms_total_arcsec = 0.0;
    }
    
    result.chi_squared = (sum_ra_sq + sum_dec_sq) / (2 * count - 6);  // DOF = 2*n - 6
}

FitResult LeastSquaresFitter::fit(
    const Eigen::Vector<double, 6>& initial_state,
    double epoch_mjd,
    ResidualFunction residual_func,
    STMFunction stm_func
) {
    FitResult result;
    result.state = initial_state;
    result.converged = false;
    result.num_iterations = 0;
    
    for (int iter = 0; iter < max_iterations_; ++iter) {
        result.num_iterations = iter + 1;
        
        // Compute residuals
        auto obs_residuals = residual_func(result.state, epoch_mjd);
        
        // Reject outliers
        if (iter > 0) {
            reject_outliers(obs_residuals);
        }
        
        // Build design matrix
        auto A = build_design_matrix(obs_residuals, result.state, epoch_mjd, stm_func);
        
        // Pack residuals into vector
        int n_obs = obs_residuals.size();
        Eigen::VectorXd residuals(2 * n_obs);
        Eigen::VectorXd weights(2 * n_obs);
        
        for (int i = 0; i < n_obs; ++i) {
            if (!obs_residuals[i].rejected) {
                residuals(2*i) = obs_residuals[i].ra_residual_arcsec;
                residuals(2*i+1) = obs_residuals[i].dec_residual_arcsec;
                weights(2*i) = obs_residuals[i].weight_ra;
                weights(2*i+1) = obs_residuals[i].weight_dec;
            } else {
                residuals(2*i) = 0.0;
                residuals(2*i+1) = 0.0;
                weights(2*i) = 0.0;
                weights(2*i+1) = 0.0;
            }
        }
        
        // Solve normal equations
        Eigen::Vector<double, 6> dx = solve_normal_equations(
            A, residuals, weights, result.covariance
        );
        
        // Update state
        result.state += dx;
        
        // Check convergence
        if (dx.norm() < tolerance_) {
            result.converged = true;
            result.residuals = obs_residuals;
            compute_statistics(obs_residuals, result);
            break;
        }
    }
    
    if (!result.converged) {
        // Final statistics even if not converged
        auto final_residuals = residual_func(result.state, epoch_mjd);
        result.residuals = final_residuals;
        compute_statistics(final_residuals, result);
    }
    
    return result;
}

} // namespace astdyn::orbit_determination
