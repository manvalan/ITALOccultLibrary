/**
 * @file astdyn_wrapper.cpp
 * @brief Implementazione del wrapper AstDyn
 * @author IOccultCalc Integration Team
 * @date 1 Dicembre 2025
 */

#include "astdyn_wrapper.h"
#include <sstream>
#include <chrono>
#include <cmath>
#include <fstream>

namespace ioccultcalc {

// Helper per estrarre il nome dell'oggetto dal file .eq1
static std::string extractObjectNameFromEQ1(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        return "";
    }
    
    bool header_ended = false;
    std::string line;
    while (std::getline(file, line)) {
        // Cerca fine dell'header
        if (line.find("END_OF_HEADER") != std::string::npos) {
            header_ended = true;
            continue;
        }
        
        // Dopo END_OF_HEADER, la prima linea non-commento è il nome/numero
        if (header_ended) {
            // Salta righe vuote e commenti
            if (line.empty() || line[0] == '!') {
                continue;
            }
            
            // Rimuovi spazi bianchi iniziali e finali
            size_t start = line.find_first_not_of(" \t\r\n");
            size_t end = line.find_last_not_of(" \t\r\n");
            
            if (start != std::string::npos && end != std::string::npos) {
                return line.substr(start, end - start + 1);
            }
        }
    }
    
    return "";
}

AstDynWrapper::AstDynWrapper(const PropagationSettings& settings)
    : settings_(settings)
    , current_epoch_mjd_(0.0)
    , initialized_(false)
{
    // Crea effemeridi planetarie (usa DE440 embedded di default)
    ephemeris_ = std::make_shared<astdyn::ephemeris::PlanetaryEphemeris>();
    
    // Inizializza propagatore
    initializePropagator();
}

AstDynWrapper::~AstDynWrapper() = default;

void AstDynWrapper::initializePropagator() {
    // Crea integratore RKF78
    auto integrator = std::make_unique<astdyn::propagation::RKF78Integrator>(
        settings_.initial_step,
        settings_.tolerance
    );
    
    // Configura perturbazioni
    astdyn::propagation::PropagatorSettings prop_settings;
    prop_settings.include_planets = settings_.include_planets;
    prop_settings.include_relativity = settings_.include_relativity;
    prop_settings.include_asteroids = settings_.include_asteroids;
    prop_settings.perturb_mercury = settings_.perturb_mercury;
    prop_settings.perturb_venus = settings_.perturb_venus;
    prop_settings.perturb_earth = settings_.perturb_earth;
    prop_settings.perturb_mars = settings_.perturb_mars;
    prop_settings.perturb_jupiter = settings_.perturb_jupiter;
    prop_settings.perturb_saturn = settings_.perturb_saturn;
    prop_settings.perturb_uranus = settings_.perturb_uranus;
    prop_settings.perturb_neptune = settings_.perturb_neptune;
    
    // Crea propagatore
    propagator_ = std::make_unique<astdyn::propagation::Propagator>(
        std::move(integrator),
        ephemeris_,
        prop_settings
    );
}

bool AstDynWrapper::loadFromEQ1File(const std::string& filepath) {
    try {
        // Usa parser AstDyn
        astdyn::io::parsers::OrbFitEQ1Parser parser;
        current_elements_ = parser.parse(filepath);
        
        // Salva epoca
        current_epoch_mjd_ = current_elements_.epoch_mjd_tdb;
        
        // Il parser AstDyn non estrae object_name dal file .eq1
        // Lo leggiamo manualmente dalla prima linea non-commento
        object_name_ = extractObjectNameFromEQ1(filepath);
        
        // Se non trovato, usa un nome vuoto
        if (object_name_.empty()) {
            object_name_ = "Unknown";
        }
        
        initialized_ = true;
        return true;
        
    } catch (const std::exception& e) {
        initialized_ = false;
        last_stats_ = std::string("Errore caricamento: ") + e.what();
        return false;
    }
}

void AstDynWrapper::setKeplerianElements(
    double a, double e, double i,
    double Omega, double omega, double M,
    double epoch_mjd_tdb,
    const std::string& name)
{
    // Salva elementi per uso futuro (usa parser format)
    current_elements_.semi_major_axis = a;
    current_elements_.eccentricity = e;
    current_elements_.inclination = i;
    current_elements_.longitude_asc_node = Omega;
    current_elements_.argument_perihelion = omega;
    current_elements_.mean_anomaly = M;
    current_elements_.epoch_mjd_tdb = epoch_mjd_tdb;
    current_elements_.object_name = name;
    
    current_epoch_mjd_ = epoch_mjd_tdb;
    object_name_ = name;
    initialized_ = true;
}

CartesianStateICRF AstDynWrapper::propagateToEpoch(double target_mjd_tdb) {
    if (!initialized_) {
        throw std::runtime_error("AstDynWrapper: elementi non inizializzati. "
                                 "Chiamare loadFromEQ1File() o setKeplerianElements() prima.");
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Crea KeplerianElements per propagazione (struct in propagation namespace)
    astdyn::propagation::KeplerianElements kep_initial;
    kep_initial.epoch_mjd_tdb = current_elements_.epoch_mjd_tdb;
    kep_initial.semi_major_axis = current_elements_.semi_major_axis;
    kep_initial.eccentricity = current_elements_.eccentricity;
    kep_initial.inclination = current_elements_.inclination;
    kep_initial.longitude_ascending_node = current_elements_.longitude_asc_node;
    kep_initial.argument_perihelion = current_elements_.argument_perihelion;
    kep_initial.mean_anomaly = current_elements_.mean_anomaly;
    kep_initial.gravitational_parameter = astdyn::constants::GMS;  // AU³/day²
    
    // Propagazione con AstDyn
    auto kep_final = propagator_->propagate_keplerian(kep_initial, target_mjd_tdb);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Statistiche
    std::ostringstream oss;
    oss << "Propagazione: " << current_epoch_mjd_ << " → " << target_mjd_tdb 
        << " MJD (" << (target_mjd_tdb - current_epoch_mjd_) << " giorni)\n"
        << "Tempo: " << duration.count() << " ms";
    last_stats_ = oss.str();
    
    // Converti in cartesiano (frame ECLM J2000)
    auto cart_ecl = astdyn::propagation::keplerian_to_cartesian(kep_final);
    
    // Converti da ECLM J2000 a ICRF
    return eclipticToICRF(cart_ecl, target_mjd_tdb);
}

CartesianStateICRF AstDynWrapper::eclipticToICRF(
    const astdyn::propagation::CartesianElements& ecl_state,
    double epoch_mjd_tdb)
{
    // Obliquità eclittica J2000
    constexpr double epsilon = 23.4393 * M_PI / 180.0;
    const double cos_eps = std::cos(epsilon);
    const double sin_eps = std::sin(epsilon);
    
    // Rotazione posizione
    const double x_ecl = ecl_state.position.x();
    const double y_ecl = ecl_state.position.y();
    const double z_ecl = ecl_state.position.z();
    
    Eigen::Vector3d pos_icrf;
    pos_icrf.x() = x_ecl;
    pos_icrf.y() = y_ecl * cos_eps - z_ecl * sin_eps;
    pos_icrf.z() = y_ecl * sin_eps + z_ecl * cos_eps;
    
    // Rotazione velocità
    const double vx_ecl = ecl_state.velocity.x();
    const double vy_ecl = ecl_state.velocity.y();
    const double vz_ecl = ecl_state.velocity.z();
    
    Eigen::Vector3d vel_icrf;
    vel_icrf.x() = vx_ecl;
    vel_icrf.y() = vy_ecl * cos_eps - vz_ecl * sin_eps;
    vel_icrf.z() = vy_ecl * sin_eps + vz_ecl * cos_eps;
    
    // Crea risultato
    CartesianStateICRF result;
    result.position = pos_icrf;
    result.velocity = vel_icrf;
    result.epoch_mjd_tdb = epoch_mjd_tdb;
    
    return result;
}

} // namespace ioccultcalc
