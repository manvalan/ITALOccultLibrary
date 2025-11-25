/**
 * Test diagnostico per verificare la propagazione
 */

#include <iostream>
#include <iomanip>
#include "orbfit/propagation/OrbitalElements.hpp"
#include "orbfit/propagation/Propagator.hpp"
#include "orbfit/propagation/Integrator.hpp"
#include "orbfit/ephemeris/PlanetaryEphemeris.hpp"
#include "orbfit/core/Constants.hpp"

using namespace orbfit;
using namespace orbfit::propagation;
using namespace orbfit::constants;

int main() {
    std::cout << std::setprecision(12);
    
    // Orbita circolare semplice (quasi Terra)
    KeplerianElements kep;
    kep.epoch_mjd_tdb = MJD2000;
    kep.semi_major_axis = 1.0;          // 1 AU
    kep.eccentricity = 0.0;             // Circolare
    kep.inclination = 0.0;
    kep.longitude_ascending_node = 0.0;
    kep.argument_perihelion = 0.0;
    kep.mean_anomaly = 0.0;             // Al perihelion
    kep.gravitational_parameter = GMS;   // Sole
    
    std::cout << "=== TEST DIAGNOSTICO PROPAGAZIONE ===\n\n";
    std::cout << "Orbita iniziale:\n";
    std::cout << "  a = " << kep.semi_major_axis << " AU\n";
    std::cout << "  e = " << kep.eccentricity << "\n";
    std::cout << "  GM = " << kep.gravitational_parameter << " AU³/day²\n";
    std::cout << "  Periodo = " << kep.period() << " giorni\n\n";
    
    // Converti a Cartesiano
    CartesianElements cart = keplerian_to_cartesian(kep);
    std::cout << "Stato Cartesiano iniziale:\n";
    std::cout << "  r = [" << cart.position(0) << ", " 
              << cart.position(1) << ", " << cart.position(2) << "] AU\n";
    std::cout << "  v = [" << cart.velocity(0) << ", " 
              << cart.velocity(1) << ", " << cart.velocity(2) << "] AU/day\n";
    std::cout << "  |r| = " << cart.position.norm() << " AU\n";
    std::cout << "  |v| = " << cart.velocity.norm() << " AU/day\n\n";
    
    // Verifica energia
    double r_mag = cart.position.norm();
    double v_mag = cart.velocity.norm();
    double energy = 0.5 * v_mag * v_mag - GMS / r_mag;
    double expected_energy = -GMS / (2.0 * kep.semi_major_axis);
    std::cout << "Energia specifica:\n";
    std::cout << "  Calcolata: " << energy << " AU²/day²\n";
    std::cout << "  Attesa: " << expected_energy << " AU²/day²\n";
    std::cout << "  Errore: " << std::abs(energy - expected_energy) << "\n\n";
    
    // Propaga per 1 giorno
    double dt = 1.0;
    double target_mjd = kep.epoch_mjd_tdb + dt;
    
    std::cout << "Propagazione analitica (2-body):\n";
    TwoBodyPropagator prop_analytical;
    KeplerianElements kep_analytical = prop_analytical.propagate(kep, target_mjd);
    std::cout << "  a_final = " << kep_analytical.semi_major_axis << " AU\n";
    std::cout << "  e_final = " << kep_analytical.eccentricity << "\n";
    std::cout << "  M_final = " << kep_analytical.mean_anomaly * RAD_TO_DEG << " deg\n";
    std::cout << "  Delta_a = " << (kep_analytical.semi_major_axis - kep.semi_major_axis) << " AU\n\n";
    
    // Propaga numericamente con RK4
    std::cout << "Propagazione numerica (RK4, h=0.1 giorni):\n";
    auto integrator_rk4 = std::make_unique<RK4Integrator>(0.1);
    auto ephemeris = std::make_shared<ephemeris::PlanetaryEphemeris>();
    PropagatorSettings settings;
    settings.include_planets = false;  // Solo 2-body
    settings.central_body_gm = GMS;
    
    std::cout << "  Settings: central_body_gm = " << settings.central_body_gm << " AU³/day²\n";
    
    Propagator prop_rk4(std::move(integrator_rk4), ephemeris, settings);
    
    try {
        KeplerianElements kep_rk4 = prop_rk4.propagate_keplerian(kep, target_mjd);
        std::cout << "  a_final = " << kep_rk4.semi_major_axis << " AU\n";
        std::cout << "  e_final = " << kep_rk4.eccentricity << "\n";
        std::cout << "  M_final = " << kep_rk4.mean_anomaly * RAD_TO_DEG << " deg\n";
        std::cout << "  Delta_a = " << (kep_rk4.semi_major_axis - kep.semi_major_axis) << " AU\n";
        std::cout << "  Delta_e = " << (kep_rk4.eccentricity - kep.eccentricity) << "\n";
        
        // Confronta con analitico
        double error_a = std::abs(kep_rk4.semi_major_axis - kep_analytical.semi_major_axis);
        double error_e = std::abs(kep_rk4.eccentricity - kep_analytical.eccentricity);
        std::cout << "\n  Errore vs analitico:\n";
        std::cout << "    |Delta_a| = " << error_a << " AU\n";
        std::cout << "    |Delta_e| = " << error_e << "\n";
        
    } catch (const std::exception& e) {
        std::cout << "  ERRORE: " << e.what() << "\n";
    }
    
    // ========================================================================
    // TEST 2: Orbita asteroide (come nel benchmark)
    // ========================================================================
    
    std::cout << "\n\n=== TEST ORBITA ASTEROIDE ===\n\n";
    
    KeplerianElements ast;
    ast.epoch_mjd_tdb = MJD2000;
    ast.semi_major_axis = 2.5;
    ast.eccentricity = 0.15;
    ast.inclination = 10.0 * DEG_TO_RAD;
    ast.longitude_ascending_node = 80.0 * DEG_TO_RAD;
    ast.argument_perihelion = 73.0 * DEG_TO_RAD;
    ast.mean_anomaly = 45.0 * DEG_TO_RAD;
    ast.gravitational_parameter = GMS;
    
    std::cout << "Orbita asteroide:\n";
    std::cout << "  a = " << ast.semi_major_axis << " AU\n";
    std::cout << "  e = " << ast.eccentricity << "\n";
    std::cout << "  i = " << ast.inclination * RAD_TO_DEG << " deg\n";
    std::cout << "  Ω = " << ast.longitude_ascending_node * RAD_TO_DEG << " deg\n";
    std::cout << "  ω = " << ast.argument_perihelion * RAD_TO_DEG << " deg\n";
    std::cout << "  M = " << ast.mean_anomaly * RAD_TO_DEG << " deg\n";
    std::cout << "  Periodo = " << ast.period() << " giorni\n\n";
    
    // Test round-trip
    CartesianElements ast_cart = keplerian_to_cartesian(ast);
    std::cout << "Conversione a Cartesiano:\n";
    std::cout << "  r = [" << ast_cart.position(0) << ", " 
              << ast_cart.position(1) << ", " << ast_cart.position(2) << "] AU\n";
    std::cout << "  v = [" << ast_cart.velocity(0) << ", " 
              << ast_cart.velocity(1) << ", " << ast_cart.velocity(2) << "] AU/day\n\n";
    
    // Verifica energia prima della conversione inversa
    double r_ast_mag = ast_cart.position.norm();
    double v_ast_mag = ast_cart.velocity.norm();
    double energy_ast = 0.5 * v_ast_mag * v_ast_mag - GMS / r_ast_mag;
    double expected_a_ast = -GMS / (2.0 * energy_ast);
    
    std::cout << "Verifica energia prima di Cart->Kep:\n";
    std::cout << "  |r| = " << r_ast_mag << " AU\n";
    std::cout << "  |v| = " << v_ast_mag << " AU/day\n";
    std::cout << "  Energia = " << energy_ast << " AU²/day²\n";
    std::cout << "  a da energia = " << expected_a_ast << " AU\n";
    std::cout << "  a atteso = " << ast.semi_major_axis << " AU\n";
    std::cout << "  GM usato = " << ast_cart.gravitational_parameter << " AU³/day²\n\n";
    
    KeplerianElements ast_back = cartesian_to_keplerian(ast_cart);
    std::cout << "Round-trip Cartesiano->Kepleriano:\n";
    std::cout << "  a_back = " << ast_back.semi_major_axis << " AU (Delta = " 
              << (ast_back.semi_major_axis - ast.semi_major_axis) << ")\n";
    std::cout << "  e_back = " << ast_back.eccentricity << " (Delta = " 
              << (ast_back.eccentricity - ast.eccentricity) << ")\n";
    std::cout << "  i_back = " << ast_back.inclination * RAD_TO_DEG << " deg (Delta = " 
              << (ast_back.inclination - ast.inclination) * RAD_TO_DEG << ")\n";
    std::cout << "  Ω_back = " << ast_back.longitude_ascending_node * RAD_TO_DEG << " deg (Delta = " 
              << (ast_back.longitude_ascending_node - ast.longitude_ascending_node) * RAD_TO_DEG << ")\n";
    std::cout << "  ω_back = " << ast_back.argument_perihelion * RAD_TO_DEG << " deg (Delta = " 
              << (ast_back.argument_perihelion - ast.argument_perihelion) * RAD_TO_DEG << ")\n";
    std::cout << "  M_back = " << ast_back.mean_anomaly * RAD_TO_DEG << " deg (Delta = " 
              << (ast_back.mean_anomaly - ast.mean_anomaly) * RAD_TO_DEG << ")\n\n";
    
    // Propagazione analitica
    double dt_ast = 10.0;
    double target_mjd_ast = ast.epoch_mjd_tdb + dt_ast;
    
    KeplerianElements ast_analytical = TwoBodyPropagator::propagate(ast, target_mjd_ast);
    std::cout << "Propagazione analitica (" << dt_ast << " giorni):\n";
    std::cout << "  a_final = " << ast_analytical.semi_major_axis << " AU\n";
    std::cout << "  e_final = " << ast_analytical.eccentricity << "\n";
    std::cout << "  M_final = " << ast_analytical.mean_anomaly * RAD_TO_DEG << " deg\n";
    std::cout << "  Delta_a = " << (ast_analytical.semi_major_axis - ast.semi_major_axis) << " AU\n\n";
    
    // Propagazione numerica
    std::cout << "Propagazione numerica (RK4, h=0.1 giorni):\n";
    auto integrator_ast = std::make_unique<RK4Integrator>(0.1);
    Propagator prop_ast(std::move(integrator_ast), ephemeris, settings);
    
    try {
        KeplerianElements ast_num = prop_ast.propagate_keplerian(ast, target_mjd_ast);
        std::cout << "  a_final = " << ast_num.semi_major_axis << " AU\n";
        std::cout << "  e_final = " << ast_num.eccentricity << "\n";
        std::cout << "  M_final = " << ast_num.mean_anomaly * RAD_TO_DEG << " deg\n";
        std::cout << "  Delta_a = " << (ast_num.semi_major_axis - ast.semi_major_axis) << " AU\n";
        std::cout << "  Delta_e = " << (ast_num.eccentricity - ast.eccentricity) << "\n";
        
        double error_a_ast = std::abs(ast_num.semi_major_axis - ast_analytical.semi_major_axis);
        double error_e_ast = std::abs(ast_num.eccentricity - ast_analytical.eccentricity);
        std::cout << "\n  Errore vs analitico:\n";
        std::cout << "    |Delta_a| = " << error_a_ast << " AU\n";
        std::cout << "    |Delta_e| = " << error_e_ast << "\n";
        
    } catch (const std::exception& e) {
        std::cout << "  ERRORE: " << e.what() << "\n";
    }
    
    return 0;
}
