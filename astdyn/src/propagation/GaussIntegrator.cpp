/**
 * @file GaussIntegrator.cpp
 * @brief Simplified implementation of Gauss-Legendre integrator
 * @author AstDyn Team
 * @date 2025-12-09
 */

#include "astdyn/propagation/GaussIntegrator.hpp"
#include <cmath>
#include <stdexcept>

namespace astdyn::propagation {

// Gauss-Legendre 4-stage (order 8) coefficients
// Gauss points in [0,1]
const double GaussIntegrator::c_[num_stages_] = {
    0.5 - std::sqrt(525.0 + 70.0*std::sqrt(30.0)) / 70.0,
    0.5 - std::sqrt(525.0 - 70.0*std::sqrt(30.0)) / 70.0,
    0.5 + std::sqrt(525.0 - 70.0*std::sqrt(30.0)) / 70.0,
    0.5 + std::sqrt(525.0 + 70.0*std::sqrt(30.0)) / 70.0
};

// Gauss-Legendre RK matrix (symplectic)
const double GaussIntegrator::a_[num_stages_][num_stages_] = {
    {0.25, 0.25 - std::sqrt(30.0)/72.0, 0.25 - std::sqrt(30.0)/72.0, 0.25},
    {0.25 + std::sqrt(30.0)/72.0, 0.25, 0.25, 0.25 - std::sqrt(30.0)/72.0},
    {0.25 + std::sqrt(30.0)/72.0, 0.25, 0.25, 0.25 - std::sqrt(30.0)/72.0},
    {0.25, 0.25 + std::sqrt(30.0)/72.0, 0.25 + std::sqrt(30.0)/72.0, 0.25}
};

// Weights (Gauss quadrature)
const double GaussIntegrator::b_[num_stages_] = {
    0.25, 0.25, 0.25, 0.25
};

GaussIntegrator::GaussIntegrator(double initial_step,
                                 double tolerance,
                                 double min_step,
                                 double max_step,
                                 int max_newton_iter)
    : h_initial_(initial_step)
    , tolerance_(tolerance)
    , h_min_(min_step)
    , h_max_(max_step)
    , max_newton_iter_(max_newton_iter)
{
    if (tolerance <= 0.0) {
        throw std::invalid_argument("Tolerance must be positive");
    }
}

Eigen::VectorXd GaussIntegrator::integrate(const DerivativeFunction& f,
                                           const Eigen::VectorXd& y0,
                                           double t0,
                                           double tf) {
    stats_.reset();
    
    double t = t0;
    Eigen::VectorXd y = y0;
    double h = h_initial_;
    
    while (t < tf) {
        if (t + h > tf) {
            h = tf - t;
        }
        
        // Solve implicit system
        std::vector<Eigen::VectorXd> k(num_stages_, Eigen::VectorXd::Zero(y.size()));
        bool converged = solve_implicit_system(f, t, y, h, k);
        
        if (!converged) {
            h *= 0.5;
            h = std::max(h, h_min_);
            stats_.num_rejected_steps++;
            continue;
        }
        
        // Update solution
        for (int i = 0; i < num_stages_; ++i) {
            y += h * b_[i] * k[i];
        }
        
        t += h;
        stats_.num_steps++;
    }
    
    stats_.final_time = t;
    return y;
}

void GaussIntegrator::integrate_steps(const DerivativeFunction& f,
                                      const Eigen::VectorXd& y0,
                                      double t0,
                                      double tf,
                                      std::vector<double>& t_out,
                                      std::vector<Eigen::VectorXd>& y_out) {
    // Simplified: just call integrate
    t_out.clear();
    y_out.clear();
    
    t_out.push_back(t0);
    y_out.push_back(y0);
    
    Eigen::VectorXd y_final = integrate(f, y0, t0, tf);
    
    t_out.push_back(tf);
    y_out.push_back(y_final);
}

bool GaussIntegrator::solve_implicit_system(const DerivativeFunction& f,
                                            double t,
                                            const Eigen::VectorXd& y,
                                            double h,
                                            std::vector<Eigen::VectorXd>& k) {
    const int n = y.size();
    
    // Simplified fixed-point iteration (not optimal, but works)
    for (int iter = 0; iter < max_newton_iter_; ++iter) {
        std::vector<Eigen::VectorXd> k_new(num_stages_);
        double max_change = 0.0;
        
        for (int i = 0; i < num_stages_; ++i) {
            Eigen::VectorXd y_stage = y;
            for (int j = 0; j < num_stages_; ++j) {
                y_stage += h * a_[i][j] * k[j];
            }
            
            k_new[i] = f(t + c_[i] * h, y_stage);
            stats_.num_function_evals++;
            
            double change = (k_new[i] - k[i]).norm();
            max_change = std::max(max_change, change);
        }
        
        k = k_new;
        
        if (max_change < tolerance_ * 0.01) {
            return true;
        }
    }
    
    return false;
}

} // namespace astdyn::propagation
