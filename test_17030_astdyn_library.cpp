/**
 * @file test_17030_astdyn_library.cpp
 * @brief Test di AstDynPropagator (RKF78 + Perturbazioni)
 * 
 * Utilizza la classe AstDynPropagator da astdyn_propagator.cpp
 * per propagare asteroide 17030 Sierks dal 26-30 Novembre 2025.
 * 
 * Integratore: RKF78 (Runge-Kutta-Fehlberg 7/8 order)
 * Perturbazioni: 8 pianeti + asteroidali (AST17) + Schwarzschild
 * 
 * Compilazione:
 *   g++ -std=c++17 -O2 -o test_17030_lib test_17030_astdyn_library.cpp -lm
 * 
 * Esecuzione:
 *   ./test_17030_lib
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

// ============================================================================
// Include AstDynPropagator (disabilita main tramite pre-compilazione)
// ============================================================================

// Leggi il file astdyn_propagator.cpp fino al main
// Nota: dovremo compilare solo le parti necessarie
// Per ora, usiamo le classi di base

// Includiamo direttamente le parti che ci servono:
// - Costanti astronomiche
// - Strutture dati
// - RKF78Integrator
// - AstDynPropagator

// Leggi il file fino a prima del main e copialo qui
#include "astdyn/tools/astdyn_propagator.cpp"

using namespace astdyn;

// ============================================================================
// Dati e funzioni helper
// ============================================================================

struct JPLPoint {
    double jd, ra_deg, dec_deg, distance_au;
};

double angular_distance(double ra1, double dec1, double ra2, double dec2) {
    const double DEG2RAD = M_PI / 180.0;
    double ra1_rad = ra1 * DEG2RAD;
    double dec1_rad = dec1 * DEG2RAD;
    double ra2_rad = ra2 * DEG2RAD;
    double dec2_rad = dec2 * DEG2RAD;
    
    double cos_dist = std::sin(dec1_rad)*std::sin(dec2_rad) + 
                      std::cos(dec1_rad)*std::cos(dec2_rad)*std::cos(ra1_rad - ra2_rad);
    return std::acos(std::max(-1.0, std::min(1.0, cos_dist))) * 180.0 / M_PI;
}

double DateToJD(int year, int month, int day, int hour = 0) {
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    return jdn + (hour - 12.0) / 24.0;
}

// ============================================================================
// MAIN TEST
// ============================================================================

int main() {
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST: AstDynPropagator (RKF78 + Perturbazioni)            ║\n";
    std::cout << "║  Asteroide: 17030 Sierks                                  ║\n";
    std::cout << "║  Periodo: 26-30 Novembre 2025                             ║\n";
    std::cout << "║  Riferimento: JPL Horizons                                ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    try {
        // Elementi AstDys per 17030 Sierks (epoca 2018-03-16)
        double epoch_jd = DateToJD(2018, 3, 16, 0);
        
        OrbitalElements elem;
        elem.name = "17030";
        elem.epoch = epoch_jd;
        elem.a = 2.71926;           // AU
        elem.e = 0.10638;
        elem.i = 9.3708;            // gradi
        elem.Omega = 33.9247;       // gradi (nodo ascendente)
        elem.omega = 153.5094;      // gradi (argomento perielio)
        elem.M = 84.2146;           // gradi (anomalia media)
        
        std::cout << "Elementi orbitali (epoca 2018-03-16):\n";
        std::cout << std::fixed << std::setprecision(8);
        std::cout << "  a = " << elem.a << " AU\n";
        std::cout << "  e = " << elem.e << "\n";
        std::cout << "  i = " << elem.i << "°\n";
        std::cout << "  Ω = " << elem.Omega << "°\n";
        std::cout << "  ω = " << elem.omega << "°\n";
        std::cout << "  M = " << elem.M << "°\n\n";
        
        // Crea propagatore AstDyn
        std::cout << "Creazione propagatore AstDynPropagator (RKF78, tol=1e-12)...\n";
        AstDynPropagator propagator(1e-12);
        propagator.usePlanets(true);
        propagator.useAST17(false);
        propagator.useRelativity(true);
        std::cout << "  ✓ Perturbazioni planetarie: ON\n";
        std::cout << "  ✓ Correzione relativistica: ON\n";
        std::cout << "  ✓ Tolleranza RKF78: 1e-12\n\n";
        
        // Dati JPL Horizons (26-30 Nov 2025)
        std::vector<JPLPoint> jpl_data = {
            {DateToJD(2025, 11, 26, 0), 73.3847, 20.2891, 1.6843},
            {DateToJD(2025, 11, 27, 0), 73.3968, 20.3064, 1.6722},
            {DateToJD(2025, 11, 28, 0), 73.4087, 20.3235, 1.6602},
            {DateToJD(2025, 11, 29, 0), 73.4208, 20.3408, 1.6483},
            {DateToJD(2025, 11, 30, 0), 73.4330, 20.3582, 1.6365}
        };
        
        std::cout << "PROPAGAZIONE CON AstDynPropagator\n";
        std::cout << "════════════════════════════════════════════════════════════\n\n";
        
        std::cout << "Data         RA (AstDyn)  Dec (AstDyn) RA (JPL)    Dec (JPL)   \n";
        std::cout << "             Distance(AU) ΔRA (arcsec) ΔDec (arcsec) Separazione\n";
        std::cout << "─────────────────────────────────────────────────────────────────\n";
        
        double max_sep = 0.0;
        double max_dist_err = 0.0;
        double max_ra_err = 0.0;
        double max_dec_err = 0.0;
        
        for (const auto& jpl : jpl_data) {
            // Propaga con AstDynPropagator
            EquatorialCoords coords = propagator.propagateElements(elem, jpl.jd);
            
            // Converti da radianti a gradi
            double ra_deg = coords.ra * 180.0 / M_PI;
            double dec_deg = coords.dec * 180.0 / M_PI;
            
            // Errori
            double delta_ra = (ra_deg - jpl.ra_deg) * 3600.0;  // arcsec
            double delta_dec = (dec_deg - jpl.dec_deg) * 3600.0;  // arcsec
            double delta_dist = coords.dist - jpl.distance_au;
            double sep = angular_distance(ra_deg, dec_deg, jpl.ra_deg, jpl.dec_deg) * 3600.0;  // arcsec
            
            max_sep = std::max(max_sep, std::abs(sep));
            max_dist_err = std::max(max_dist_err, std::abs(delta_dist));
            max_ra_err = std::max(max_ra_err, std::abs(delta_ra));
            max_dec_err = std::max(max_dec_err, std::abs(delta_dec));
            
            // Determina giorno
            int day;
            if (jpl.jd < DateToJD(2025, 11, 27)) day = 26;
            else if (jpl.jd < DateToJD(2025, 11, 28)) day = 27;
            else if (jpl.jd < DateToJD(2025, 11, 29)) day = 28;
            else if (jpl.jd < DateToJD(2025, 11, 30)) day = 29;
            else day = 30;
            
            // Output
            std::cout << "2025-11-" << std::setfill('0') << std::setw(2) << day;
            std::cout << " " << std::setfill(' ');
            std::cout << std::fixed << std::setprecision(4);
            std::cout << " " << std::setw(9) << ra_deg << "° ";
            std::cout << " " << std::setw(9) << dec_deg << "° ";
            std::cout << " " << std::setw(9) << jpl.ra_deg << "° ";
            std::cout << " " << std::setw(9) << jpl.dec_deg << "°\n";
            
            std::cout << "             ";
            std::cout << std::setprecision(6) << std::setw(9) << coords.dist;
            std::cout << "  " << std::setprecision(2) << std::setw(10) << delta_ra;
            std::cout << "    " << std::setw(10) << delta_dec;
            std::cout << "    " << std::setprecision(0) << std::setw(10) << sep << "\n\n";
        }
        
        std::cout << "════════════════════════════════════════════════════════════\n\n";
        std::cout << "RISULTATI FINALI:\n";
        std::cout << std::fixed;
        std::cout << "  Max ΔRA:          " << std::setprecision(2) << max_ra_err << " arcsec\n";
        std::cout << "  Max ΔDec:         " << std::setprecision(2) << max_dec_err << " arcsec\n";
        std::cout << "  Max Separazione:  " << std::setprecision(2) << max_sep << " arcsec (" 
                  << std::setprecision(4) << (max_sep/3600.0) << "°)\n";
        std::cout << "  Max ΔDistanza:    " << std::setprecision(6) << max_dist_err << " AU\n\n";
        
        std::cout << "VALUTAZIONE:\n";
        if (max_sep < 0.1) {
            std::cout << "✅ ECCELLENTE: Errore < 0.1 arcsec (JPL-grade accuracy)\n";
        } else if (max_sep < 1.0) {
            std::cout << "✅ OTTIMO: Errore < 1 arcsec\n";
        } else if (max_sep < 10.0) {
            std::cout << "⚠️  BUONO: Errore < 10 arcsec\n";
        } else if (max_sep < 60.0) {
            std::cout << "⚠️  ACCETTABILE: Errore < 60 arcsec\n";
        } else {
            std::cout << "❌ INACCETTABILE: Errore > 60 arcsec\n";
        }
        
        std::cout << "\n════════════════════════════════════════════════════════════\n\n";
        
        std::cout << "NOTE:\n";
        std::cout << "- RKF78: 7th/8th order Runge-Kutta-Fehlberg (13 stages)\n";
        std::cout << "- Perturbazioni: Sole, 8 pianeti, Schwarzschild\n";
        std::cout << "- Propagazione: " << std::setprecision(1) << (jpl_data[4].jd - epoch_jd)/365.25 
                  << " anni (2018→2025)\n";
        std::cout << "- Frame: ICRF (equatoriale J2000)\n\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ ERRORE: " << e.what() << std::endl;
        return 1;
    }
}
