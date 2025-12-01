/**
 * @file astdyn_interface.h
 * @brief Interface per ITALOccultLibrary/AstDyn in IOccultCalc
 * @author Michele Bigi - IOccultCalc Integration Team
 * @date 1 Dicembre 2025
 * @version 3.0 - Frame conversion ECLM→ICRF validata
 * 
 * Questo file fornisce l'interfaccia tra IOccultCalc e la libreria
 * ITALOccultLibrary, includendo:
 * - Wrapper per propagatore AstDyn RKF78
 * - Parser elementi .eq1 (OEF2.0 format)
 * - Conversioni orbitali con frame conversion validata
 * - Compatibilità con JPL Horizons (ICRF frame)
 * 
 * VALIDAZIONE: Testato vs JPL Horizons
 * - Errore: 0.0003 arcsec @ 3.2 AU
 * - Test case: Asteroid 17030 Sierks
 * - Report: VALIDAZIONE_ASTDYN_FRAME_CORRECTED.md
 */

#ifndef IOCCULTCALC_ASTDYN_INTERFACE_H
#define IOCCULTCALC_ASTDYN_INTERFACE_H

#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/coordinates.h"
#include <string>
#include <memory>
#include <vector>
#include <Eigen/Dense>

namespace ioccultcalc {

// ============================================================================
// Forward declarations
// ============================================================================

class AstDynPropagatorImpl;  // PIMPL per nascondere dipendenze AstDyn

// ============================================================================
// Strutture Dati
// ============================================================================

/**
 * @struct PropagationResult
 * @brief Risultato propagazione orbitale
 */
struct PropagationResult {
    Eigen::Vector3d position;      ///< Posizione [AU] in ICRF
    Eigen::Vector3d velocity;      ///< Velocità [AU/day] in ICRF
    double epoch_mjd;              ///< Epoca risultato [MJD TDT]
    
    // Statistiche integrazione
    int num_steps;                 ///< Numero step integrazione
    int num_function_evals;        ///< Valutazioni funzione
    int num_rejected_steps;        ///< Step rifiutati
    double computation_time_ms;    ///< Tempo computazionale [ms]
    
    /// Conversione a OrbitalElements IOccultCalc
    OrbitalElements toOrbitalElements() const;
    
    /// Conversione a SpaceTimeCoordinate IOccultCalc
    SpaceTimeCoordinate toSpaceTimeCoordinate() const;
};

/**
 * @struct PropagatorConfig
 * @brief Configurazione propagatore
 */
struct PropagatorConfig {
    // Integratore
    enum class IntegratorType {
        RKF78,      ///< Runge-Kutta-Fehlberg 7/8 (raccomandato)
        RK4,        ///< Runge-Kutta 4th order
        DOPRI54     ///< Dormand-Prince 5(4)
    };
    
    IntegratorType integrator = IntegratorType::RKF78;
    double tolerance = 1e-12;           ///< Tolleranza integrazione [AU]
    double initial_step = 0.1;          ///< Step iniziale [giorni]
    
    // Perturbazioni
    bool perturb_mercury = true;
    bool perturb_venus = true;
    bool perturb_earth = true;
    bool perturb_mars = true;
    bool perturb_jupiter = true;
    bool perturb_saturn = true;
    bool perturb_uranus = true;
    bool perturb_neptune = true;
    bool include_relativity = true;     ///< Schwarzschild correction
    bool include_asteroids = true;      ///< AST17 perturbers
    
    /// Factory per configurazione standard
    static PropagatorConfig standardConfig() {
        return PropagatorConfig();  // Default già ottimali
    }
    
    /// Factory per configurazione veloce (test)
    static PropagatorConfig fastConfig() {
        PropagatorConfig cfg;
        cfg.tolerance = 1e-10;
        cfg.include_asteroids = false;
        return cfg;
    }
    
    /// Factory per massima precisione
    static PropagatorConfig highPrecisionConfig() {
        PropagatorConfig cfg;
        cfg.tolerance = 1e-14;
        return cfg;
    }
};

// ============================================================================
// AstDynPropagator - Main Interface
// ============================================================================

/**
 * @class AstDynPropagator
 * @brief Wrapper per propagatore AstDyn con ITALOccultLibrary
 * 
 * Questa classe fornisce un'interfaccia semplificata per:
 * - Lettura file .eq1 (OEF2.0)
 * - Propagazione orbitale ad alta precisione
 * - Conversione frame ECLM→ICRF automatica
 * - Integrazione con strutture IOccultCalc
 * 
 * IMPORTANTE: Include conversione frame ECLM J2000 → ICRF
 * File .eq1 sono in eclittico, JPL Horizons usa equatoriale (ICRF).
 * 
 * Esempio d'uso:
 * @code
 *   // 1. Crea propagatore
 *   AstDynPropagator prop;
 *   
 *   // 2. Opzionale: configura
 *   PropagatorConfig cfg = PropagatorConfig::standardConfig();
 *   prop.configure(cfg);
 *   
 *   // 3. Carica elementi da .eq1
 *   prop.loadElements("17030.eq1");
 *   
 *   // 4. Propaga a epoca target
 *   double target_mjd = 61007.0;
 *   auto result = prop.propagate(target_mjd);
 *   
 *   // 5. Risultato è in ICRF (compatibile JPL Horizons)
 *   std::cout << "X = " << result.position.x() << " AU\n";
 * @endcode
 */
class AstDynPropagator {
public:
    /**
     * @brief Costruttore con configurazione opzionale
     * @param config Configurazione propagatore (default: standard)
     */
    explicit AstDynPropagator(
        const PropagatorConfig& config = PropagatorConfig::standardConfig());
    
    /// Destruttore
    ~AstDynPropagator();
    
    // Moveable but not copyable (PIMPL con unique_ptr)
    AstDynPropagator(AstDynPropagator&&) noexcept;
    AstDynPropagator& operator=(AstDynPropagator&&) noexcept;
    AstDynPropagator(const AstDynPropagator&) = delete;
    AstDynPropagator& operator=(const AstDynPropagator&) = delete;
    
    // ========================================================================
    // Configurazione
    // ========================================================================
    
    /**
     * @brief Configura propagatore
     * @param config Nuova configurazione
     */
    void configure(const PropagatorConfig& config);
    
    /**
     * @brief Ottiene configurazione corrente
     * @return Configurazione attiva
     */
    PropagatorConfig getConfig() const;
    
    // ========================================================================
    // Caricamento Elementi
    // ========================================================================
    
    /**
     * @brief Carica elementi orbitali da file .eq1
     * @param eq1_file Path al file .eq1 (formato OEF2.0)
     * @throws std::runtime_error Se file non trovato o formato invalido
     * 
     * File .eq1 contiene elementi equinoziali in frame ECLM J2000.
     * Parser legge formato OEF2.0 standard da AstDyS.
     */
    void loadElements(const std::string& eq1_file);
    
    /**
     * @brief Carica elementi da OrbitalElements IOccultCalc
     * @param elements Elementi orbitali
     * 
     * Converte automaticamente da rappresentazione IOccultCalc.
     * Assume frame ICRF se non specificato diversamente.
     */
    void loadElements(const OrbitalElements& elements);
    
    /**
     * @brief Ottiene elementi caricati
     * @return Elementi orbitali correnti
     * @throws std::runtime_error Se nessun elemento caricato
     */
    OrbitalElements getElements() const;
    
    /**
     * @brief Verifica se elementi sono caricati
     * @return true se elementi caricati
     */
    bool hasElements() const;
    
    // ========================================================================
    // Propagazione
    // ========================================================================
    
    /**
     * @brief Propaga orbita a epoca target
     * @param target_mjd Epoca target [MJD TDT]
     * @return Risultato propagazione in frame ICRF
     * @throws std::runtime_error Se elementi non caricati o propagazione fallisce
     * 
     * IMPORTANTE: Risultato è in frame ICRF (equatoriale J2000),
     * compatibile con JPL Horizons. Include conversione automatica
     * da frame ECLM J2000 degli elementi .eq1.
     */
    PropagationResult propagate(double target_mjd);
    
    /**
     * @brief Propaga orbita per intervallo temporale
     * @param target_mjd Epoca finale [MJD TDT]
     * @param step_days Step output [giorni]
     * @return Vettore risultati (epoche intermedie)
     * 
     * Genera output ad intervalli regolari durante propagazione.
     * Utile per effemeridi.
     */
    std::vector<PropagationResult> propagateWithSteps(
        double target_mjd, 
        double step_days = 1.0);
    
    /**
     * @brief Propaga a più epoche target
     * @param target_mjds Vettore epoche target [MJD TDT]
     * @return Vettore risultati (uno per epoca)
     * 
     * Più efficiente di chiamate multiple a propagate().
     */
    std::vector<PropagationResult> propagateMultiple(
        const std::vector<double>& target_mjds);
    
    // ========================================================================
    // Statistiche
    // ========================================================================
    
    /**
     * @brief Ottiene statistiche ultima propagazione
     * @return Informazioni integrazione
     */
    struct IntegrationStats {
        int total_steps;
        int rejected_steps;
        int function_evaluations;
        double avg_step_size;
        double max_step_size;
        double min_step_size;
    };
    
    IntegrationStats getLastIntegrationStats() const;
    
    /**
     * @brief Reset statistiche
     */
    void resetStatistics();
    
    // ========================================================================
    // Utility
    // ========================================================================
    
    /**
     * @brief Valida risultato vs JPL Horizons
     * @param result Risultato da validare
     * @param jpl_position Posizione JPL [AU, ICRF]
     * @return Errore lineare [km]
     * 
     * Utility per validazione: confronta con dati JPL Horizons.
     */
    static double validateWithJPL(
        const PropagationResult& result,
        const Eigen::Vector3d& jpl_position);
    
    /**
     * @brief Scarica elementi da AstDyS
     * @param asteroid_number Numero asteroide (es: 17030)
     * @param output_file Path file output .eq1
     * @return true se successo
     * 
     * Utility per scaricare automaticamente da:
     * https://newton.spacedys.com/astdys2/
     */
    static bool downloadFromAstDyS(
        int asteroid_number,
        const std::string& output_file);

private:
    std::unique_ptr<AstDynPropagatorImpl> pimpl_;  ///< Implementation
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Converte coordinate eclittiche in ICRF
 * @param ecl_position Posizione eclittica [AU]
 * @param ecl_velocity Velocità eclittica [AU/day]
 * @return Pair <position_icrf, velocity_icrf>
 * 
 * Applica rotazione obliquità ε = 23.439291°
 */
std::pair<Eigen::Vector3d, Eigen::Vector3d> eclipticToICRF(
    const Eigen::Vector3d& ecl_position,
    const Eigen::Vector3d& ecl_velocity);

/**
 * @brief Converte coordinate ICRF in eclittiche (inverso)
 * @param icrf_position Posizione ICRF [AU]
 * @param icrf_velocity Velocità ICRF [AU/day]
 * @return Pair <position_ecl, velocity_ecl>
 */
std::pair<Eigen::Vector3d, Eigen::Vector3d> icrfToEcliptic(
    const Eigen::Vector3d& icrf_position,
    const Eigen::Vector3d& icrf_velocity);

/**
 * @brief Calcola errore angolare
 * @param r1 Posizione 1 [AU]
 * @param r2 Posizione 2 [AU]
 * @return Errore angolare [arcsec]
 */
double calculateAngularError(
    const Eigen::Vector3d& r1,
    const Eigen::Vector3d& r2);

} // namespace ioccultcalc

#endif // IOCCULTCALC_ASTDYN_INTERFACE_H
