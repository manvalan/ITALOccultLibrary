/**
 * @file astdyn_wrapper.cpp
 * @brief Implementazione wrapper AstDyn
 * @author IOccultCalc Integration Team
 * @date 1 Dicembre 2025
 */

#include "astdyn_wrapper.h"
#include <chrono>
#include <stdexcept>
#include <sstream>

namespace ioccultcalc {

// ============================================================================
// Costruttore / Distruttore
// ============================================================================

AstDynWrapper::AstDynWrapper(const AstDynConfig& config)
    : config_(config), propagator_(nullptr), stats_()
{
    if (!config_.isValid()) {
        throw std::invalid_argument("Invalid AstDynConfig");
    }
    initialize();
}

AstDynWrapper::~AstDynWrapper() = default;

AstDynWrapper::AstDynWrapper(AstDynWrapper&&) noexcept = default;
AstDynWrapper& AstDynWrapper::operator=(AstDynWrapper&&) noexcept = default;

// ============================================================================
// Inizializzazione
// ============================================================================

void AstDynWrapper::initialize() {
    try {
        // Crea propagatore con tolleranza
        propagator_ = std::make_unique<astdyn::AstDynPropagator>(
            config_.tolerance,
            config_.max_step_days
        );
        
        // Configura perturbazioni
        configurePerturbations();
        
        // Reset statistiche
        stats_ = Statistics();
        
    } catch (const std::exception& e) {
        std::ostringstream oss;
        oss << "Failed to initialize AstDynPropagator: " << e.what();
        throw std::runtime_error(oss.str());
    }
}

void AstDynWrapper::configurePerturbations() {
    if (!propagator_) {
        throw std::runtime_error("Propagator not initialized");
    }
    
    stats_.num_perturbations = 0;
    
    // 1. Perturbazioni planetarie
    if (config_.enable_planets) {
        for (int i = 1; i <= config_.num_planets; ++i) {
            propagator_->enablePlanet(i);
            stats_.num_perturbations++;
        }
    }
    
    // 2. Relatività generale (Schwarzschild)
    if (config_.enable_relativity) {
        propagator_->enableRelativity();
        stats_.num_perturbations++;
    }
    
    // 3. Perturbazioni asteroidi (AST17: Cerere + Vesta)
    if (config_.enable_asteroids) {
        propagator_->enableAsteroidPerturbations();
        stats_.num_perturbations += 2;  // Cerere + Vesta
    }
    
    // 4. Carica ephemeris se specificato
    if (!config_.ephemeris_path.empty()) {
        propagator_->loadEphemeris(config_.ephemeris_path);
    }
}

// ============================================================================
// Propagazione
// ============================================================================

PropagationResult AstDynWrapper::propagate(
    const Eigen::Vector3d& r0,
    const Eigen::Vector3d& v0,
    double t0,
    double tf)
{
    PropagationResult result;
    result.epoch_jd = tf;
    
    // Validazione input
    std::string error;
    if (!validateInput(r0, v0, t0, tf, error)) {
        result.success = false;
        result.error_message = error;
        return result;
    }
    
    // Check inizializzazione
    if (!propagator_) {
        result.success = false;
        result.error_message = "Propagator not initialized";
        return result;
    }
    
    try {
        // Timing
        auto start = std::chrono::high_resolution_clock::now();
        
        // Crea stato iniziale AstDyn
        astdyn::StateVector state_in;
        state_in.position = r0;
        state_in.velocity = v0;
        state_in.epoch_jd = t0;
        
        // Propaga
        astdyn::StateVector state_out = propagator_->propagate(
            state_in, 
            tf
        );
        
        // Timing fine
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            end - start);
        
        // Popola risultato
        result.position = state_out.position;
        result.velocity = state_out.velocity;
        result.epoch_jd = state_out.epoch_jd;
        result.num_steps = propagator_->getLastStepCount();
        result.computation_time_ms = duration.count() / 1000.0;
        result.success = true;
        
        // Aggiorna statistiche
        stats_.total_steps += result.num_steps;
        stats_.total_time_ms += result.computation_time_ms;
        if (result.num_steps > 0) {
            stats_.avg_step_size_days = 
                std::abs(tf - t0) / static_cast<double>(result.num_steps);
        }
        
    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = std::string("Propagation failed: ") + e.what();
    }
    
    return result;
}

std::vector<PropagationResult> AstDynWrapper::propagateMultiple(
    const Eigen::Vector3d& r0,
    const Eigen::Vector3d& v0,
    double t0,
    const std::vector<double>& output_times)
{
    std::vector<PropagationResult> results;
    results.reserve(output_times.size());
    
    if (output_times.empty()) {
        return results;
    }
    
    // Propaga sequenzialmente (AstDyn è stabile per propagazioni continue)
    Eigen::Vector3d r_current = r0;
    Eigen::Vector3d v_current = v0;
    double t_current = t0;
    
    for (double t_target : output_times) {
        auto result = propagate(r_current, v_current, t_current, t_target);
        
        if (!result.success) {
            // Propagazione fallita - ritorna risultati parziali
            results.push_back(result);
            break;
        }
        
        results.push_back(result);
        
        // Aggiorna stato per prossima propagazione
        r_current = result.position;
        v_current = result.velocity;
        t_current = result.epoch_jd;
    }
    
    return results;
}

// ============================================================================
// Configurazione
// ============================================================================

void AstDynWrapper::reconfigure(const AstDynConfig& config) {
    if (!config.isValid()) {
        throw std::invalid_argument("Invalid AstDynConfig");
    }
    
    config_ = config;
    reset();
}

void AstDynWrapper::reset() {
    propagator_.reset();
    stats_ = Statistics();
    initialize();
}

// ============================================================================
// Validazione
// ============================================================================

bool AstDynWrapper::validateInput(
    const Eigen::Vector3d& r0,
    const Eigen::Vector3d& v0,
    double t0, double tf,
    std::string& error) const
{
    // 1. Check valori finiti
    if (!std::isfinite(r0.norm()) || !std::isfinite(v0.norm())) {
        error = "Non-finite position or velocity";
        return false;
    }
    
    // 2. Check range posizione [0.1, 100 AU]
    double r_mag = r0.norm();
    if (r_mag < 0.1 || r_mag > 100.0) {
        error = "Position out of range [0.1, 100] AU";
        return false;
    }
    
    // 3. Check velocità ragionevole [< 100 AU/day ≈ 173 km/s]
    double v_mag = v0.norm();
    if (v_mag > 100.0) {
        error = "Velocity too large (> 100 AU/day)";
        return false;
    }
    
    // 4. Check epoche valide [1900-2100 circa]
    if (t0 < 2415020.5 || t0 > 2488069.5) {  // JD 1900-2100
        error = "Initial epoch out of range [1900-2100]";
        return false;
    }
    if (tf < 2415020.5 || tf > 2488069.5) {
        error = "Final epoch out of range [1900-2100]";
        return false;
    }
    
    // 5. Check intervallo ragionevole [< 1000 giorni]
    double dt = std::abs(tf - t0);
    if (dt > 1000.0) {
        error = "Time span too large (> 1000 days)";
        return false;
    }
    
    // 6. Check t0 != tf
    if (std::abs(tf - t0) < 1e-9) {
        error = "Initial and final epochs are identical";
        return false;
    }
    
    return true;
}

} // namespace ioccultcalc
