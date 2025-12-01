/**
 * @file test_astdyn_integration_standalone.cpp
 * @brief Test standalone integrazione completa AstDyn
 * @author IOccultCalc Integration Team
 * @date 1 Dicembre 2025
 * 
 * Questo esempio mostra come utilizzare i moduli creati in FASE 1
 * per propagare un asteroide da file .eq1 con AstDyn.
 * 
 * Compilation:
 * g++ -std=c++17 \
 *     -I../templates_ioccultcalc/include \
 *     -I/usr/local/include \
 *     test_astdyn_integration_standalone.cpp \
 *     ../templates_ioccultcalc/src/eq1_parser.cpp \
 *     ../templates_ioccultcalc/src/orbital_conversions.cpp \
 *     ../templates_ioccultcalc/src/astdyn_wrapper.cpp \
 *     -L/usr/local/lib \
 *     -lastdyn \
 *     -o test_astdyn_integration
 * 
 * Usage:
 * ./test_astdyn_integration ../astdyn/data/17030.eq1 2460643.77083
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>

// Include moduli FASE 1
#include "eq1_parser.h"
#include "orbital_conversions.h"
#include "astdyn_wrapper.h"

using namespace ioccultcalc;

/**
 * @brief Calcola errore angolare tra due posizioni
 * @param pos1 Prima posizione (AU)
 * @param pos2 Seconda posizione (AU)
 * @return Errore in arcsec
 */
double computeAngularError(const Eigen::Vector3d& pos1, 
                          const Eigen::Vector3d& pos2) {
    // Normalizza vettori
    Eigen::Vector3d u1 = pos1.normalized();
    Eigen::Vector3d u2 = pos2.normalized();
    
    // Prodotto scalare
    double cos_angle = u1.dot(u2);
    cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
    
    // Angolo in radianti
    double angle_rad = std::acos(cos_angle);
    
    // Converti in arcsec
    double angle_arcsec = angle_rad * 206264.806247;
    
    return angle_arcsec;
}

/**
 * @brief Print stato cartesiano formattato
 */
void printCartesianState(const CartesianState& state, 
                         const std::string& label) {
    std::cout << "=== " << label << " ===" << std::endl;
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "Epoch JD: " << state.epoch_jd << std::endl;
    std::cout << "Position (AU):" << std::endl;
    std::cout << "  x = " << state.position(0) << std::endl;
    std::cout << "  y = " << state.position(1) << std::endl;
    std::cout << "  z = " << state.position(2) << std::endl;
    std::cout << "  |r| = " << state.position.norm() << std::endl;
    std::cout << "Velocity (AU/day):" << std::endl;
    std::cout << "  vx = " << state.velocity(0) << std::endl;
    std::cout << "  vy = " << state.velocity(1) << std::endl;
    std::cout << "  vz = " << state.velocity(2) << std::endl;
    std::cout << "  |v| = " << state.velocity.norm() << std::endl;
    std::cout << std::endl;
}

/**
 * @brief Print elementi kepleriani formattati
 */
void printKeplerianElements(const KeplerianElements& kep) {
    std::cout << "=== Keplerian Elements ===" << std::endl;
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Asteroid: " << kep.name << std::endl;
    std::cout << "Epoch JD: " << std::setprecision(6) << kep.epoch_jd << std::endl;
    std::cout << "a (AU): " << kep.a << std::endl;
    std::cout << "e: " << kep.e << std::endl;
    std::cout << "i (deg): " << (kep.i * 180.0 / M_PI) << std::endl;
    std::cout << "Ω (deg): " << (kep.Omega * 180.0 / M_PI) << std::endl;
    std::cout << "ω (deg): " << (kep.omega * 180.0 / M_PI) << std::endl;
    std::cout << "M (deg): " << (kep.M * 180.0 / M_PI) << std::endl;
    std::cout << std::endl;
}

/**
 * @brief Main test function
 */
int main(int argc, char* argv[]) {
    // Parse arguments
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <eq1_file> <target_jd>" << std::endl;
        std::cerr << "Example: " << argv[0] << " 17030.eq1 2460643.77083" << std::endl;
        return 1;
    }
    
    std::string eq1_file = argv[1];
    double target_jd = std::stod(argv[2]);
    
    std::cout << "======================================" << std::endl;
    std::cout << "  AstDyn Integration Test Standalone" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << std::endl;
    
    try {
        // ================================================================
        // STEP 1: Parse .eq1 file
        // ================================================================
        std::cout << "[1/5] Parsing .eq1 file: " << eq1_file << std::endl;
        
        auto eq_elements = EQ1Parser::parseFile(eq1_file);
        
        std::cout << "  Asteroid: " << eq_elements.name << std::endl;
        std::cout << "  Epoch MJD: " << eq_elements.epoch_mjd << std::endl;
        std::cout << "  a = " << eq_elements.a << " AU" << std::endl;
        std::cout << "  e = " << eq_elements.getEccentricity() << std::endl;
        std::cout << "  H = " << eq_elements.H << std::endl;
        std::cout << "  G = " << eq_elements.G << std::endl;
        std::cout << "  Valid: " << (eq_elements.isValid() ? "YES" : "NO") << std::endl;
        std::cout << std::endl;
        
        if (!eq_elements.isValid()) {
            throw std::runtime_error("Invalid eq1 elements");
        }
        
        // ================================================================
        // STEP 2: Convert Equinoctial → Keplerian
        // ================================================================
        std::cout << "[2/5] Converting Equinoctial → Keplerian" << std::endl;
        
        auto kep_elements = OrbitalConversions::equinoctialToKeplerian(eq_elements);
        printKeplerianElements(kep_elements);
        
        if (!kep_elements.isValid()) {
            throw std::runtime_error("Invalid Keplerian elements");
        }
        
        // ================================================================
        // STEP 3: Convert Keplerian → Cartesian (Ecliptic)
        // ================================================================
        std::cout << "[3/5] Converting Keplerian → Cartesian (Ecliptic)" << std::endl;
        
        auto cart_ecliptic = OrbitalConversions::keplerianToCartesian(kep_elements);
        printCartesianState(cart_ecliptic, "Ecliptic State");
        
        // ================================================================
        // STEP 4: Rotate Ecliptic → ICRF
        // ================================================================
        std::cout << "[4/5] Rotating Ecliptic → ICRF" << std::endl;
        
        auto cart_icrf = OrbitalConversions::eclipticToICRF(cart_ecliptic);
        printCartesianState(cart_icrf, "ICRF State (Initial)");
        
        if (!OrbitalConversions::validateICRF(cart_icrf)) {
            throw std::runtime_error("Invalid ICRF state");
        }
        
        // ================================================================
        // STEP 5: Propagate with AstDyn
        // ================================================================
        std::cout << "[5/5] Propagating with AstDyn" << std::endl;
        std::cout << "  Initial epoch: " << cart_icrf.epoch_jd << " JD" << std::endl;
        std::cout << "  Target epoch:  " << target_jd << " JD" << std::endl;
        std::cout << "  Δt = " << (target_jd - cart_icrf.epoch_jd) << " days" << std::endl;
        std::cout << std::endl;
        
        // Create wrapper with JPL-compliant configuration
        auto wrapper = AstDynWrapperFactory::forOccultations();
        
        // Propagate
        auto start = std::chrono::high_resolution_clock::now();
        
        auto result = wrapper->propagate(
            cart_icrf.position,
            cart_icrf.velocity,
            cart_icrf.epoch_jd,
            target_jd
        );
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        // Check result
        if (!result.success) {
            std::cerr << "ERROR: " << result.error_message << std::endl;
            return 1;
        }
        
        // ================================================================
        // RESULTS
        // ================================================================
        std::cout << "======================================" << std::endl;
        std::cout << "  PROPAGATION RESULTS" << std::endl;
        std::cout << "======================================" << std::endl;
        std::cout << std::endl;
        
        std::cout << std::fixed << std::setprecision(12);
        std::cout << "Final Position (ICRF, AU):" << std::endl;
        std::cout << "  x = " << result.position(0) << std::endl;
        std::cout << "  y = " << result.position(1) << std::endl;
        std::cout << "  z = " << result.position(2) << std::endl;
        std::cout << "  |r| = " << result.position.norm() << std::endl;
        std::cout << std::endl;
        
        std::cout << "Final Velocity (ICRF, AU/day):" << std::endl;
        std::cout << "  vx = " << result.velocity(0) << std::endl;
        std::cout << "  vy = " << result.velocity(1) << std::endl;
        std::cout << "  vz = " << result.velocity(2) << std::endl;
        std::cout << "  |v| = " << result.velocity.norm() << std::endl;
        std::cout << std::endl;
        
        std::cout << "Performance:" << std::endl;
        std::cout << "  Steps: " << result.num_steps << std::endl;
        std::cout << "  Computation time: " << result.computation_time_ms << " ms" << std::endl;
        std::cout << "  Wall time: " << duration.count() << " ms" << std::endl;
        
        if (result.num_steps > 0) {
            double avg_step = std::abs(target_jd - cart_icrf.epoch_jd) / 
                             static_cast<double>(result.num_steps);
            std::cout << "  Avg step size: " << avg_step << " days" << std::endl;
        }
        std::cout << std::endl;
        
        // ================================================================
        // VALIDATION (Optional - compare with JPL if available)
        // ================================================================
        std::cout << "======================================" << std::endl;
        std::cout << "  VALIDATION" << std::endl;
        std::cout << "======================================" << std::endl;
        std::cout << std::endl;
        
        // TODO: Insert JPL Horizons reference position here
        // For asteroid 17030 at JD 2460643.77083 (2025-11-26 06:30 UTC)
        
        std::cout << "To validate against JPL Horizons:" << std::endl;
        std::cout << "1. Go to https://ssd.jpl.nasa.gov/horizons.cgi" << std::endl;
        std::cout << "2. Target: asteroid " << eq_elements.name << std::endl;
        std::cout << "3. Epoch: " << target_jd << " JD" << std::endl;
        std::cout << "4. Output: State vectors (ICRF)" << std::endl;
        std::cout << "5. Compare position with above results" << std::endl;
        std::cout << "6. Expected accuracy: < 2 arcsec" << std::endl;
        std::cout << std::endl;
        
        std::cout << "======================================" << std::endl;
        std::cout << "  TEST COMPLETED SUCCESSFULLY" << std::endl;
        std::cout << "======================================" << std::endl;
        
        return 0;
        
    } catch (const EQ1ParseException& e) {
        std::cerr << "EQ1 Parse Error: " << e.what() << std::endl;
        return 1;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
