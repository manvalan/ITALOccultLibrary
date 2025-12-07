/**
 * @file astdyn_wrapper.h
 * @brief Wrapper semplificato per AstDyn library
 * @author IOccultCalc Integration Team
 * @date 1 Dicembre 2025
 * 
 * Wrapper che incapsula l'uso di AstDyn Propagator per la propagazione
 * orbitale di asteroidi con perturbazioni complete.
 * 
 * API Reale AstDyn:
 * - astdyn::propagation::Propagator (propagatore principale)
 * - astdyn::propagation::RKF78Integrator (integratore Runge-Kutta-Fehlberg 7/8)
 * - astdyn::ephemeris::PlanetaryEphemeris (effemeridi planetarie)
 * - astdyn::io::parsers::OrbFitEQ1Parser (parser file .eq1)
 */

#ifndef IOCCULTCALC_ASTDYN_WRAPPER_H
#define IOCCULTCALC_ASTDYN_WRAPPER_H

#include <astdyn/AstDyn.hpp>
#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/coordinates/CartesianState.hpp>
#include <astdyn/coordinates/KeplerianElements.hpp>
#include <astdyn/io/parsers/OrbFitEQ1Parser.hpp>
#include <Eigen/Dense>
#include <memory>
#include <string>

namespace ioccultcalc {

/**
 * @brief Configurazione per la propagazione con AstDyn
 */
struct AstDynConfig {
    /// Tolleranza assoluta per RKF78 (default: 1e-12 AU)
    double tolerance = 1e-12;
    
    /// Step massimo in giorni (default: 10 giorni)
    double max_step_days = 10.0;
    
    /// Abilita perturbazioni planetarie (default: true)
    bool enable_planets = true;
    
    /// Abilita relatività generale (Schwarzschild) (default: true)
    bool enable_relativity = true;
    
    /// Abilita perturbazioni asteroidi (AST17) (default: true)
    bool enable_asteroids = true;
    
    /// Numero di pianeti da considerare [1-8] (default: 8 - tutti)
    int num_planets = 8;
    
    /// Ephemeris file path (default: vuoto - usa embedded DE440)
    std::string ephemeris_path = "";
    
    /**
     * @brief Verifica validità configurazione
     */
    bool isValid() const {
        return tolerance > 0.0 && 
               tolerance < 1e-6 &&
               max_step_days > 0.0 &&
               max_step_days <= 100.0 &&
               num_planets >= 1 &&
               num_planets <= 8;
    }
    
    /**
     * @brief Configurazione JPL-compliant (massima accuratezza)
     */
    static AstDynConfig jplCompliant() {
        AstDynConfig cfg;
        cfg.tolerance = 1e-12;
        cfg.max_step_days = 5.0;
        cfg.enable_planets = true;
        cfg.enable_relativity = true;
        cfg.enable_asteroids = true;
        cfg.num_planets = 8;
        return cfg;
    }
    
    /**
     * @brief Configurazione veloce (screening iniziale)
     */
    static AstDynConfig fast() {
        AstDynConfig cfg;
        cfg.tolerance = 1e-9;
        cfg.max_step_days = 20.0;
        cfg.enable_planets = true;
        cfg.enable_relativity = false;
        cfg.enable_asteroids = false;
        cfg.num_planets = 4;  // Solo pianeti interni
        return cfg;
    }
    
    /**
     * @brief Configurazione bilanciata (uso generale)
     */
    static AstDynConfig balanced() {
        AstDynConfig cfg;
        cfg.tolerance = 1e-11;
        cfg.max_step_days = 10.0;
        cfg.enable_planets = true;
        cfg.enable_relativity = true;
        cfg.enable_asteroids = true;
        cfg.num_planets = 8;
        return cfg;
    }
};

/**
 * @brief Risultato della propagazione
 */
struct PropagationResult {
    /// Stato finale (ICRF, AU, AU/day)
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    
    /// Epoca finale (JD TDB)
    double epoch_jd;
    
    /// Numero di step eseguiti
    int num_steps = 0;
    
    /// Tempo computazionale (ms)
    double computation_time_ms = 0.0;
    
    /// Success flag
    bool success = false;
    
    /// Messaggio di errore (se !success)
    std::string error_message;
    
    /**
     * @brief Verifica validità risultato
     */
    bool isValid() const {
        return success &&
               std::isfinite(position.norm()) &&
               std::isfinite(velocity.norm()) &&
               position.norm() > 0.1 &&
               position.norm() < 100.0;  // Range ragionevole [AU]
    }
};

/**
 * @brief Wrapper per AstDynPropagator con configurazione semplificata
 * 
 * Esempio d'uso:
 * @code
 * AstDynWrapper wrapper(AstDynConfig::jplCompliant());
 * 
 * Eigen::Vector3d r0(1.0, 0.0, 0.0);
 * Eigen::Vector3d v0(0.0, 0.017, 0.0);
 * double t0 = 2460000.5;
 * double tf = 2460100.5;
 * 
 * auto result = wrapper.propagate(r0, v0, t0, tf);
 * if (result.success) {
 *     std::cout << "Position: " << result.position.transpose() << std::endl;
 * }
 * @endcode
 */
class AstDynWrapper {
public:
    /**
     * @brief Costruttore con configurazione
     * @param config Configurazione propagatore
     */
    explicit AstDynWrapper(const AstDynConfig& config = AstDynConfig::balanced());
    
    /**
     * @brief Distruttore
     */
    ~AstDynWrapper();
    
    // Non copiabile (contiene unique_ptr)
    AstDynWrapper(const AstDynWrapper&) = delete;
    AstDynWrapper& operator=(const AstDynWrapper&) = delete;
    
    // Movibile
    AstDynWrapper(AstDynWrapper&&) noexcept;
    AstDynWrapper& operator=(AstDynWrapper&&) noexcept;
    
    /**
     * @brief Propaga stato da t0 a tf
     * @param r0 Posizione iniziale ICRF [AU]
     * @param v0 Velocità iniziale ICRF [AU/day]
     * @param t0 Epoca iniziale [JD TDB]
     * @param tf Epoca finale [JD TDB]
     * @return Risultato propagazione
     */
    PropagationResult propagate(
        const Eigen::Vector3d& r0,
        const Eigen::Vector3d& v0,
        double t0,
        double tf
    );
    
    /**
     * @brief Propaga stato con output intermedi
     * @param r0 Posizione iniziale ICRF [AU]
     * @param v0 Velocità iniziale ICRF [AU/day]
     * @param t0 Epoca iniziale [JD TDB]
     * @param output_times Epoche desiderate [JD TDB]
     * @return Vettore di risultati (uno per ogni epoca)
     */
    std::vector<PropagationResult> propagateMultiple(
        const Eigen::Vector3d& r0,
        const Eigen::Vector3d& v0,
        double t0,
        const std::vector<double>& output_times
    );
    
    /**
     * @brief Riconfigura propagatore
     * @param config Nuova configurazione
     */
    void reconfigure(const AstDynConfig& config);
    
    /**
     * @brief Ottieni configurazione corrente
     */
    const AstDynConfig& getConfig() const { return config_; }
    
    /**
     * @brief Verifica se propagatore è inizializzato
     */
    bool isInitialized() const { return propagator_ != nullptr; }
    
    /**
     * @brief Reset propagatore (forza reinizializzazione)
     */
    void reset();
    
    /**
     * @brief Ottieni statistiche ultimo run
     */
    struct Statistics {
        int total_steps = 0;
        double total_time_ms = 0.0;
        double avg_step_size_days = 0.0;
        int num_perturbations = 0;
    };
    Statistics getStatistics() const { return stats_; }

private:
    /// Configurazione corrente
    AstDynConfig config_;
    
    /// Propagatore AstDyn (pimpl)
    std::unique_ptr<astdyn::AstDynPropagator> propagator_;
    
    /// Statistiche
    Statistics stats_;
    
    /**
     * @brief Inizializza propagatore con configurazione corrente
     */
    void initialize();
    
    /**
     * @brief Configura perturbazioni
     */
    void configurePerturbations();
    
    /**
     * @brief Valida input propagazione
     */
    bool validateInput(const Eigen::Vector3d& r0,
                       const Eigen::Vector3d& v0,
                       double t0, double tf,
                       std::string& error) const;
};

/**
 * @brief Factory per creare wrapper preconfigurati
 */
class AstDynWrapperFactory {
public:
    /**
     * @brief Crea wrapper per occultazioni (massima accuratezza)
     */
    static std::unique_ptr<AstDynWrapper> forOccultations() {
        return std::make_unique<AstDynWrapper>(AstDynConfig::jplCompliant());
    }
    
    /**
     * @brief Crea wrapper per screening veloce
     */
    static std::unique_ptr<AstDynWrapper> forScreening() {
        return std::make_unique<AstDynWrapper>(AstDynConfig::fast());
    }
    
    /**
     * @brief Crea wrapper bilanciato
     */
    static std::unique_ptr<AstDynWrapper> forGeneral() {
        return std::make_unique<AstDynWrapper>(AstDynConfig::balanced());
    }
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_ASTDYN_WRAPPER_H
