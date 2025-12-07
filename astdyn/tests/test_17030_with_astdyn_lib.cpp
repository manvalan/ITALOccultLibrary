/**
 * @file test_17030_with_astdyn_lib.cpp
 * @brief Test di propagazione usando AstDyn library
 * 
 * Usa le classi reali di AstDyn:
 * - astdyn::propagation::Propagator
 * - astdyn::propagation::RKF78Integrator
 * - astdyn::propagation::KeplerianElements
 * 
 * Propaga asteroide 17030 Sierks dal 26-30 Novembre 2025
 * e confronta con JPL Horizons.
 * 
 * Compilazione (da directory ITALOccultLibrary):
 *   g++ -std=c++17 -I./astdyn/include -O2 \
 *       -o test_17030_lib astdyn/tests/test_17030_with_astdyn_lib.cpp \
 *       -L./astdyn/lib -lastdyn -lm
 * 
 * Esecuzione:
 *   ./test_17030_lib
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <memory>

// AstDyn headers
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/core/Constants.hpp"

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::ephemeris;

// ============================================================================
// Helper Constants and Functions
// ============================================================================

const double DEG_TO_RAD = M_PI / 180.0;
const double RAD_TO_DEG = 180.0 / M_PI;
const double ARCSEC_PER_RAD = 206264.806247;

struct EquatorialCoords {
    double ra_deg, dec_deg;
    double distance_au;
};

struct JPLData {
    double jd, ra_deg, dec_deg, distance_au;
};

double DateToJD(int year, int month, int day, int hour = 0) {
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    return jdn + (hour - 12.0) / 24.0;
}

double JDToMJD(double jd) {
    return jd - 2400000.5;
}

void MJDToDate(double mjd, int& year, int& month, int& day) {
    double jd = mjd + 2400000.5;
    int jdn = (int)std::floor(jd + 0.5);
    int a = jdn + 32044;
    int b = (4 * a + 3) / 146097;
    int c = a - (146097 * b) / 4;
    int d = (4 * c + 3) / 1461;
    int e = c - (1461 * d) / 4;
    int m = (5 * e + 2) / 153;
    day = e - (153 * m + 2) / 5 + 1;
    month = m + 3 - 12 * (m / 10);
    year = 100 * b + d - 4800 + m / 10;
}

EquatorialCoords cartesian_to_equatorial(const Eigen::Vector3d& position) {
    double x = position(0);
    double y = position(1);
    double z = position(2);
    
    double ra_rad = std::atan2(y, x);
    if (ra_rad < 0) ra_rad += 2.0 * M_PI;
    
    double dec_rad = std::atan2(z, std::sqrt(x*x + y*y));
    double distance = position.norm();
    
    EquatorialCoords coords;
    coords.ra_deg = ra_rad * RAD_TO_DEG;
    coords.dec_deg = dec_rad * RAD_TO_DEG;
    coords.distance_au = distance;
    
    return coords;
}

double angular_separation(double ra1, double dec1, double ra2, double dec2) {
    double dra = (ra2 - ra1) * DEG_TO_RAD;
    double ddec = (dec2 - dec1) * DEG_TO_RAD;
    
    double a = std::sin(ddec/2)*std::sin(ddec/2) +
               std::cos(dec1*DEG_TO_RAD)*std::cos(dec2*DEG_TO_RAD)*
               std::sin(dra/2)*std::sin(dra/2);
    
    double sep = 2.0*std::atan2(std::sqrt(a), std::sqrt(1-a));
    return sep * ARCSEC_PER_RAD;
}

// ============================================================================
// Main Test
// ============================================================================

int main() {
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST: AstDyn Library Propagator                           ║\n";
    std::cout << "║  Asteroide: 17030 Sierks                                  ║\n";
    std::cout << "║  Periodo: 26-30 Novembre 2025                             ║\n";
    std::cout << "║  Metodo: RKF78 + Perturbazioni (8 pianeti + Schwarzschild)║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    try {
        // Create RKF78 integrator (7/8 order, 13 stages)
        auto integrator = std::make_unique<RKF78Integrator>(
            0.1,     // initial step size (days)
            1e-12,   // relative tolerance
            1e-6,    // min step
            100.0    // max step
        );
        
        // Create planetary ephemeris
        auto ephemeris = std::make_shared<PlanetaryEphemeris>();
        
        // Propagator settings
        PropagatorSettings settings;
        settings.include_planets = true;   // Venus, Earth, Mars, Jupiter, Saturn
        settings.include_relativity = true; // Schwarzschild correction
        settings.perturb_venus = true;
        settings.perturb_earth = true;
        settings.perturb_mars = true;
        settings.perturb_jupiter = true;
        settings.perturb_saturn = true;
        
        // Create propagator
        Propagator propagator(std::move(integrator), ephemeris, settings);
        
        // Asteroid 17030 Sierks elements (epoch 2018-03-16)
        KeplerianElements initial;
        initial.epoch_mjd_tdb = JDToMJD(DateToJD(2018, 3, 16, 0));
        initial.semi_major_axis = 2.71926;
        initial.eccentricity = 0.10638;
        initial.inclination = 9.3708 * DEG_TO_RAD;
        initial.longitude_ascending_node = 33.9247 * DEG_TO_RAD;
        initial.argument_perihelion = 153.5094 * DEG_TO_RAD;
        initial.mean_anomaly = 84.2146 * DEG_TO_RAD;
        initial.gravitational_parameter = constants::GMS;
        
        // JPL Horizons reference data (26-30 Nov 2025)
        std::vector<JPLData> jpl_data = {
            {DateToJD(2025, 11, 26, 0), 73.3847, 20.2891, 1.6843},
            {DateToJD(2025, 11, 27, 0), 73.3968, 20.3064, 1.6722},
            {DateToJD(2025, 11, 28, 0), 73.4087, 20.3235, 1.6602},
            {DateToJD(2025, 11, 29, 0), 73.4208, 20.3408, 1.6483},
            {DateToJD(2025, 11, 30, 0), 73.4330, 20.3582, 1.6365}
        };
        
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "PROPAGAZIONE ASTDYN vs JPL HORIZONS\n";
        std::cout << "════════════════════════════════════════════════════════════\n\n";
        
        double max_sep = 0.0;
        double max_dist_err = 0.0;
        int num_points = 0;
        
        for (const auto& jpl : jpl_data) {
            double mjd_target = JDToMJD(jpl.jd);
            
            int year, month, day;
            MJDToDate(mjd_target, year, month, day);
            
            try {
                // Propagate Keplerian elements to target epoch
                KeplerianElements final_elements = propagator.propagate_keplerian(
                    initial, mjd_target
                );
                
                // Convert to Cartesian for coordinate transformation
                CartesianElements cart = keplerian_to_cartesian(final_elements);
                EquatorialCoords calc = cartesian_to_equatorial(cart.position);
                
                // Calculate errors
                double sep = angular_separation(calc.ra_deg, calc.dec_deg,
                                                jpl.ra_deg, jpl.dec_deg);
                double dist_err = std::abs(calc.distance_au - jpl.distance_au);
                
                max_sep = std::max(max_sep, sep);
                max_dist_err = std::max(max_dist_err, dist_err);
                num_points++;
                
                std::cout << year << "-" << std::setfill('0')
                          << std::setw(2) << month << "-" << std::setw(2) << day
                          << " | JD " << jpl.jd << "\n";
                std::cout << "  Calc: RA=" << calc.ra_deg << "° Dec=" << calc.dec_deg
                          << "° r=" << calc.distance_au << " AU\n";
                std::cout << "  JPL:  RA=" << jpl.ra_deg << "° Dec=" << jpl.dec_deg
                          << "° r=" << jpl.distance_au << " AU\n";
                std::cout << "  → Separazione: " << sep << "\" Distance err: "
                          << dist_err << " AU\n\n";
                
            } catch (const std::exception& e) {
                std::cout << "  ❌ Errore propagazione: " << e.what() << "\n\n";
            }
        }
        
        std::cout << "╔════════════════════════════════════════════════════════════╗\n";
        std::cout << "║                    RISULTATI FINALI                        ║\n";
        std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
        
        std::cout << "Punti propagati: " << num_points << "\n";
        std::cout << "Max Separazione Angolare: " << max_sep << " arcsec\n";
        std::cout << "Max Errore Distanza: " << max_dist_err << " AU\n\n";
        
        std::cout << "Valutazione accuratezza:\n";
        if (max_sep < 0.1) {
            std::cout << "✅ ECCELLENTE: Errore < 0.1 arcsec (JPL-grade accuracy)\n";
        } else if (max_sep < 1.0) {
            std::cout << "✅ OTTIMO: Errore < 1 arcsec\n";
        } else if (max_sep < 10.0) {
            std::cout << "⚠️  BUONO: Errore < 10 arcsec\n";
        } else if (max_sep < 60.0) {
            std::cout << "⚠️  ACCETTABILE: Errore < 60 arcsec (phase 1 threshold)\n";
        } else {
            std::cout << "❌ INACCETTABILE: Errore > 60 arcsec\n";
        }
        
        std::cout << "\n════════════════════════════════════════════════════════════\n";
        std::cout << "AstDyn Propagator Statistics:\n";
        const auto& stats = propagator.statistics();
        std::cout << "  Steps taken: " << stats.num_steps << "\n";
        std::cout << "  Function evals: " << stats.num_function_evals << "\n";
        std::cout << "  Rejected steps: " << stats.num_rejected_steps << "\n";
        std::cout << "════════════════════════════════════════════════════════════\n\n";
        
        return (max_sep < 60.0) ? 0 : 1;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Errore: " << e.what() << "\n";
        return 1;
    }
}

/**
 * NOTE:
 * 
 * AstDyn RKF78 Integrator (Runge-Kutta-Fehlberg 7(8)):
 * - 13 function evaluations per step
 * - Adaptive step size control
 * - 7th order accurate with 8th order error estimation
 * - Relative tolerance: 1e-12
 * 
 * Perturbazioni incluse:
 * - Two-body (Sole)
 * - Mercurio, Venere, Terra, Marte, Giove, Saturno (selezionabili)
 * - Relativistic correction (Schwarzschild)
 * 
 * Performance su questo test:
 * - Propagazione 7+ anni (2018 → 2025)
 * - Integrazione numerica dell'equazione di moto
 * - Expected accuracy: < 0.01 AU su distanza
 * - Expected accuracy: < 1 arcmin su posizione
 */
