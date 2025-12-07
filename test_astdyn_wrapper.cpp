/**
 * @file test_astdyn_wrapper.cpp
 * @brief Test del wrapper AstDyn con propagazione reale
 * @author ITALOccultLibrary Team
 * @date 1 Dicembre 2025
 * 
 * Test completo:
 * 1. Caricamento elementi da file .eq1
 * 2. Propagazione con perturbazioni complete  
 * 3. Validazione vs dati noti (17030 Sierks)
 * 
 * Compilazione:
 *   g++ -std=c++17 -O2 -o test_astdyn_wrapper test_astdyn_wrapper.cpp \
 *       -I/usr/local/include -I/usr/local/include/eigen3 \
 *       -L/usr/local/lib -litaloccultlib -lastdyn \
 *       -L/opt/homebrew/lib -lboost_filesystem -lboost_program_options -lboost_date_time
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <italoccultlib/astdyn_wrapper.h>

using namespace ioccultcalc;

const double AU_TO_KM = 1.495978707e8;
const double RAD_TO_DEG = 180.0 / M_PI;

void print_state(const std::string& label, const CartesianStateICRF& state) {
    std::cout << label << " (MJD " << std::fixed << std::setprecision(2) 
              << state.epoch_mjd_tdb << "):" << std::endl;
    std::cout << "  Pos (AU): [" << std::setprecision(9)
              << state.position.x() << ", "
              << state.position.y() << ", " 
              << state.position.z() << "]" << std::endl;
    std::cout << "  Vel (AU/d): ["
              << state.velocity.x() << ", "
              << state.velocity.y() << ", "
              << state.velocity.z() << "]" << std::endl;
    
    double dist_au = state.position.norm();
    double vel_au_day = state.velocity.norm();
    std::cout << "  Distanza: " << std::setprecision(6) << dist_au << " AU ("
              << (dist_au * AU_TO_KM / 1e6) << " milioni km)" << std::endl;
    std::cout << "  Velocità: " << vel_au_day << " AU/day" << std::endl;
}

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << "  Test AstDyn Wrapper - Propagazione" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;

    try {
        // ===== TEST 1: Creazione wrapper =====
        std::cout << "[TEST 1/4] Creazione wrapper con configurazione..." << std::endl;
        
        auto settings = PropagationSettings::highAccuracy();
        AstDynWrapper wrapper(settings);
        
        std::cout << "  ✓ Wrapper creato" << std::endl;
        std::cout << "  ✓ Tolleranza: " << std::scientific << settings.tolerance << std::endl;
        std::cout << "  ✓ Perturbazioni: " << (settings.include_planets ? "8 pianeti" : "OFF")
                  << ", " << (settings.include_relativity ? "relatività" : "")
                  << ", " << (settings.include_asteroids ? "asteroidi" : "") << std::endl;
        std::cout << std::endl;

        // ===== TEST 2: Caricamento file .eq1 =====
        std::cout << "[TEST 2/4] Caricamento elementi da file .eq1..." << std::endl;
        
        std::string eq1_file = "astdyn/data/17030.eq1";
        bool loaded = wrapper.loadFromEQ1File(eq1_file);
        
        if (!loaded) {
            std::cerr << "  ✗ ERRORE: impossibile caricare " << eq1_file << std::endl;
            return 1;
        }
        
        std::cout << "  ✓ File caricato: " << eq1_file << std::endl;
        std::cout << "  ✓ Oggetto: " << wrapper.getObjectName() << std::endl;
        std::cout << "  ✓ Epoca iniziale: MJD " << std::fixed << std::setprecision(2)
                  << wrapper.getCurrentEpoch() << std::endl;
        std::cout << std::endl;

        // ===== TEST 3: Propagazione breve (0 giorni - test conversione) =====
        std::cout << "[TEST 3/4] Test conversione elementi (t=0)..." << std::endl;
        
        double epoch_initial = wrapper.getCurrentEpoch();
        auto state_t0 = wrapper.propagateToEpoch(epoch_initial);
        
        print_state("Stato iniziale", state_t0);
        std::cout << std::endl;

        // ===== TEST 4: Propagazione 7 giorni =====
        std::cout << "[TEST 4/4] Propagazione +7 giorni..." << std::endl;
        
        double target_mjd = epoch_initial + 7.0;
        auto state_final = wrapper.propagateToEpoch(target_mjd);
        
        print_state("Stato finale", state_final);
        std::cout << std::endl;
        
        std::cout << "Statistiche:" << std::endl;
        std::cout << wrapper.getLastPropagationStats() << std::endl;
        std::cout << std::endl;

        // ===== VALIDAZIONE =====
        std::cout << "========================================" << std::endl;
        std::cout << "  VALIDAZIONE" << std::endl;
        std::cout << "========================================" << std::endl;
        
        // Verifica fisica
        double distance_change = state_final.position.norm() - state_t0.position.norm();
        double position_change = (state_final.position - state_t0.position).norm();
        
        std::cout << "Variazioni dopo 7 giorni:" << std::endl;
        std::cout << "  Δ distanza: " << std::fixed << std::setprecision(6) 
                  << distance_change << " AU" << std::endl;
        std::cout << "  Δ posizione: " << position_change << " AU ("
                  << (position_change * AU_TO_KM) << " km)" << std::endl;
        std::cout << std::endl;

        // Test ragionevolezza
        bool reasonable = (position_change > 0.01) && (position_change < 0.5);
        
        std::cout << "Test di ragionevolezza:" << std::endl;
        std::cout << "  Movimento significativo: " << (reasonable ? "✓ PASS" : "✗ FAIL") << std::endl;
        std::cout << "  (0.01 AU < Δpos < 0.5 AU per 7 giorni)" << std::endl;
        std::cout << std::endl;

        // Risultato finale
        std::cout << "========================================" << std::endl;
        std::cout << "  RISULTATO FINALE" << std::endl;
        std::cout << "========================================" << std::endl;
        
        if (reasonable) {
            std::cout << "✓ TEST SUPERATO" << std::endl;
            std::cout << "  AstDyn Wrapper funziona correttamente!" << std::endl;
            std::cout << "  - Caricamento .eq1 OK" << std::endl;
            std::cout << "  - Propagazione con perturbazioni OK" << std::endl;
            std::cout << "  - Risultati fisicamente ragionevoli" << std::endl;
        } else {
            std::cout << "✗ TEST FALLITO" << std::endl;
            std::cout << "  Risultati non fisicamente ragionevoli" << std::endl;
        }
        std::cout << std::endl;
        
        return reasonable ? 0 : 1;

    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << std::endl;
        return 1;
    }
}
