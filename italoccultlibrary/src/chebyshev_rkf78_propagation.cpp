/**
 * @file chebyshev_rkf78_propagation.cpp
 * @brief Implementazione ChebyshevRKF78Propagator con pattern PIMPL
 * @author ITALOccultLibrary Development Team
 * @date 4 December 2025
 */

#include "chebyshev_rkf78_propagation.h"
#include "astdyn_wrapper.h"
#include <iostream>
#include <stdexcept>

namespace ioccultcalc {

// ============================================================================
// PIMPL Implementation Class
// ============================================================================

class ChebyshevRKF78Propagator::Impl {
public:
    AstDynWrapper wrapper_;
    RKF78PropagationConfig config_;
    
    explicit Impl(const std::string& eq1_file)
        : wrapper_(PropagationSettings::highAccuracy()),
          config_() {
        
        if (!wrapper_.loadFromEQ1File(eq1_file)) {
            throw std::runtime_error("Impossibile caricare file .eq1: " + eq1_file);
        }
    }
};

// ============================================================================
// ChebyshevRKF78Propagator - Costruttore/Distruttore/Move
// ============================================================================

ChebyshevRKF78Propagator::ChebyshevRKF78Propagator(const std::string& eq1_file)
    : pImpl_(std::make_unique<Impl>(eq1_file)) {
}

ChebyshevRKF78Propagator::~ChebyshevRKF78Propagator() = default;

ChebyshevRKF78Propagator::ChebyshevRKF78Propagator(ChebyshevRKF78Propagator&&) noexcept = default;
ChebyshevRKF78Propagator& ChebyshevRKF78Propagator::operator=(ChebyshevRKF78Propagator&&) noexcept = default;

// ============================================================================
// propagateForChebyshev
// ============================================================================

std::vector<Eigen::Vector3d> ChebyshevRKF78Propagator::propagateForChebyshev(
    double start_epoch,
    double end_epoch,
    size_t num_points) {
    
    if (num_points < 3) {
        throw std::runtime_error("num_points deve essere >= 3");
    }
    
    if (start_epoch >= end_epoch) {
        throw std::runtime_error("start_epoch deve essere < end_epoch");
    }
    
    std::vector<Eigen::Vector3d> positions;
    positions.reserve(num_points);
    
    for (size_t i = 0; i < num_points; ++i) {
        double t = start_epoch + i * (end_epoch - start_epoch) / (num_points - 1);
        
        try {
            CartesianStateICRF state = pImpl_->wrapper_.propagateToEpoch(t);
            positions.push_back(state.position);
        } catch (const std::exception& e) {
            std::cerr << "Errore nella propagazione all'epoca " << t 
                      << ": " << e.what() << std::endl;
            throw;
        }
    }
    
    return positions;
}

// ============================================================================
// propagateWithVelocities
// ============================================================================

std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>>
ChebyshevRKF78Propagator::propagateWithVelocities(
    double start_epoch,
    double end_epoch,
    size_t num_points) {
    
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Vector3d> velocities;
    positions.reserve(num_points);
    velocities.reserve(num_points);
    
    for (size_t i = 0; i < num_points; ++i) {
        double t = start_epoch + i * (end_epoch - start_epoch) / (num_points - 1);
        
        CartesianStateICRF state = pImpl_->wrapper_.propagateToEpoch(t);
        positions.push_back(state.position);
        velocities.push_back(state.velocity);
    }
    
    return {positions, velocities};
}

// ============================================================================
// getConfig / setConfig
// ============================================================================

const RKF78PropagationConfig& ChebyshevRKF78Propagator::getConfig() const {
    return pImpl_->config_;
}

void ChebyshevRKF78Propagator::setConfig(const RKF78PropagationConfig& cfg) {
    pImpl_->config_ = cfg;
}

// ============================================================================
// Factory Functions
// ============================================================================

ChebyshevRKF78Propagator createChebyshevPropagatorFullCorrections(
    const std::string& eq1_file) {
    return ChebyshevRKF78Propagator(eq1_file);
}

std::pair<size_t, size_t> getRecommendedChebyshevParameters(
    double interval_days,
    const std::string& desired_accuracy) {
    
    if (desired_accuracy == "quick") {
        if (interval_days <= 1.0) return {10, 4};
        if (interval_days <= 7.0) return {30, 5};
        if (interval_days <= 14.0) return {50, 6};
        return {80, 7};
    }
    else if (desired_accuracy == "standard") {
        if (interval_days <= 1.0) return {20, 5};
        if (interval_days <= 7.0) return {50, 6};
        if (interval_days <= 14.0) return {100, 8};
        return {150, 9};
    }
    else if (desired_accuracy == "high") {
        if (interval_days <= 1.0) return {30, 6};
        if (interval_days <= 7.0) return {80, 8};
        if (interval_days <= 14.0) return {150, 10};
        return {200, 12};
    }
    else if (desired_accuracy == "ultra") {
        if (interval_days <= 1.0) return {50, 8};
        if (interval_days <= 7.0) return {120, 10};
        if (interval_days <= 14.0) return {200, 12};
        return {300, 15};
    }
    
    return {100, 8};  // Default: standard
}

} // namespace ioccultcalc
