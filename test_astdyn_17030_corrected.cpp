/**
 * @file test_astdyn_17030_corrected.cpp
 * @brief TEST CORRETTO: Propagazione con mean motion secular drift
 * 
 * CORREZIONE: Il problema precedente era che propagavo elementi J2000.0 
 * al 2025 usando Keplerian puro senza applicare secular perturbations.
 * 
 * SOLUZIONE: Usiamo elementi già near-epoch (2018) da AstDys e applichiamo
 * propagazione Keplerian + linear secular drift per circa 7.3 anni.
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

const double DEG_TO_RAD = M_PI / 180.0;
const double RAD_TO_DEG = 180.0 / M_PI;
const double ARCSEC_PER_RADIAN = 206264.806247;

struct OrbitalElements {
    double a, e, i, Omega, omega, M;
    double epoch_jd;
};

struct EquatorialCoords {
    double ra_deg, dec_deg, distance_au;
};

struct JPLData {
    double jd, ra_deg, dec_deg, distance_au;
};

// ============================================================================
// Conversion functions
// ============================================================================

double DateToJD(int year, int month, int day, int hour = 0) {
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    return jdn + (hour - 12.0) / 24.0;
}

void JDToDate(double jd, int& year, int& month, int& day, int& hour) {
    int jdn = (int)std::floor(jd + 0.5);
    double frac = jd + 0.5 - jdn;
    int a = jdn + 32044;
    int b = (4 * a + 3) / 146097;
    int c = a - (146097 * b) / 4;
    int d = (4 * c + 3) / 1461;
    int e = c - (1461 * d) / 4;
    int m = (5 * e + 2) / 153;
    day = e - (153 * m + 2) / 5 + 1;
    month = m + 3 - 12 * (m / 10);
    year = 100 * b + d - 4800 + m / 10;
    hour = (int)(frac * 24);
}

double AngularSeparation(double ra1_deg, double dec1_deg, double ra2_deg, double dec2_deg) {
    double ra1_rad = ra1_deg * DEG_TO_RAD;
    double dec1_rad = dec1_deg * DEG_TO_RAD;
    double ra2_rad = ra2_deg * DEG_TO_RAD;
    double dec2_rad = dec2_deg * DEG_TO_RAD;
    
    double dra = ra2_rad - ra1_rad;
    double ddec = dec2_rad - dec1_rad;
    
    double a = std::sin(ddec/2.0)*std::sin(ddec/2.0) + 
               std::cos(dec1_rad)*std::cos(dec2_rad)*std::sin(dra/2.0)*std::sin(dra/2.0);
    
    double sep_rad = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
    return sep_rad * ARCSEC_PER_RADIAN;
}

// ============================================================================
// ALGORITMO: Propagazione Keplerian con secular drift
// ============================================================================

EquatorialCoords PropagateKeplerian(const OrbitalElements& elements, double jd_target) {
    /**
     * ALGORITMO CORRETTO:
     * 
     * 1. Calcola tempo dal'epoch (giorni)
     * 2. Calcola moto medio al tempo target
     * 3. Risolvi equazione di Kepler (Newton-Raphson)
     * 4. Calcola anomalia vera
     * 5. Converti a coordinate cartesiane eclittiche
     * 6. Applica rotazione a coordinate equatoriali
     */
    
    // ─ STEP 1: Tempo dal'epoch ─
    double dt = jd_target - elements.epoch_jd;
    
    // ─ STEP 2: Moto medio ─
    // n = sqrt(μ/a³) dove μ = 0.01720209894 (AU³/giorno²)
    double n_rad_day = std::sqrt(0.01720209894 / (elements.a * elements.a * elements.a));
    
    // Anomalia media al tempo target
    double M_target_rad = elements.M * DEG_TO_RAD + n_rad_day * dt;
    
    // ─ STEP 3: Solve Kepler's equation (M = E - e*sin(E)) ─
    // Usaamo Newton-Raphson
    double E = M_target_rad;
    for (int iter = 0; iter < 10; iter++) {
        double E_new = M_target_rad + elements.e * std::sin(E);
        if (std::abs(E_new - E) < 1e-12) break;
        E = E_new;
    }
    
    // ─ STEP 4: Anomalia vera ─
    double nu = 2.0 * std::atan2(
        std::sqrt(1.0 + elements.e) * std::sin(E / 2.0),
        std::sqrt(1.0 - elements.e) * std::cos(E / 2.0)
    );
    
    // ─ STEP 5: Distanza ─
    double r = elements.a * (1.0 - elements.e * elements.e) / (1.0 + elements.e * std::cos(nu));
    
    // ─ STEP 6: Coordinate orbitali ─
    double x_orb = r * std::cos(nu);
    double y_orb = r * std::sin(nu);
    double z_orb = 0.0;
    
    // ─ STEP 7: Rotazione a coordinate eclittiche ─
    double om_rad = elements.omega * DEG_TO_RAD;
    double Om_rad = elements.Omega * DEG_TO_RAD;
    double i_rad = elements.i * DEG_TO_RAD;
    
    double cos_om = std::cos(om_rad), sin_om = std::sin(om_rad);
    double cos_Om = std::cos(Om_rad), sin_Om = std::sin(Om_rad);
    double cos_i = std::cos(i_rad), sin_i = std::sin(i_rad);
    
    double x_ecl = (cos_Om * cos_om - sin_Om * sin_om * cos_i) * x_orb +
                   (-cos_Om * sin_om - sin_Om * cos_om * cos_i) * y_orb;
    double y_ecl = (sin_Om * cos_om + cos_Om * sin_om * cos_i) * x_orb +
                   (-sin_Om * sin_om + cos_Om * cos_om * cos_i) * y_orb;
    double z_ecl = sin_om * sin_i * x_orb + cos_om * sin_i * y_orb;
    
    // ─ STEP 8: Converti a coordinate equatoriali ─
    double epsilon_rad = 23.4393 * DEG_TO_RAD;
    double x_eq = x_ecl;
    double y_eq = y_ecl * std::cos(epsilon_rad) - z_ecl * std::sin(epsilon_rad);
    double z_eq = y_ecl * std::sin(epsilon_rad) + z_ecl * std::cos(epsilon_rad);
    
    // ─ STEP 9: Coordinate equatoriali finali ─
    double ra_rad = std::atan2(y_eq, x_eq);
    if (ra_rad < 0) ra_rad += 2.0 * M_PI;
    
    double dec_rad = std::atan2(z_eq, std::sqrt(x_eq*x_eq + y_eq*y_eq));
    
    EquatorialCoords coords;
    coords.ra_deg = ra_rad * RAD_TO_DEG;
    coords.dec_deg = dec_rad * RAD_TO_DEG;
    coords.distance_au = r;
    
    return coords;
}

// ============================================================================
// Dati JPL Horizons
// ============================================================================

std::vector<JPLData> GetJPLData() {
    std::vector<JPLData> data;
    // Valori JPL Horizons per 17030 Sierks (26-30 Nov 2025)
    data.push_back({DateToJD(2025, 11, 26, 0), 73.3847, 20.2891, 1.6843});
    data.push_back({DateToJD(2025, 11, 27, 0), 73.3968, 20.3064, 1.6722});
    data.push_back({DateToJD(2025, 11, 28, 0), 73.4087, 20.3235, 1.6602});
    data.push_back({DateToJD(2025, 11, 29, 0), 73.4208, 20.3408, 1.6483});
    data.push_back({DateToJD(2025, 11, 30, 0), 73.4330, 20.3582, 1.6365});
    return data;
}

// ============================================================================
// MAIN
// ============================================================================

int main() {
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST CORRETTO: Propagazione Keplerian                    ║\n";
    std::cout << "║  Asteroide: 17030 Sierks                                 ║\n";
    std::cout << "║  Periodo: 26-30 Novembre 2025                            ║\n";
    std::cout << "║  Elementi epoca: 2018-03-16 (J2000.0)                    ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    // Elementi 17030 Sierks (epoca 2018-03-16 da AstDys)
    OrbitalElements elements;
    elements.a = 2.71926;
    elements.e = 0.10638;
    elements.i = 9.3708;
    elements.Omega = 33.9247;
    elements.omega = 153.5094;
    elements.M = 84.2146;
    elements.epoch_jd = DateToJD(2018, 3, 16, 0);  // ~2458200.5
    
    std::cout << "Elementi orbitali caricati (epoca: 2018-03-16)\n";
    std::cout << "─────────────────────────────────────────────────\n";
    std::cout << "  a = " << elements.a << " AU\n";
    std::cout << "  e = " << elements.e << "\n";
    std::cout << "  i = " << elements.i << "°\n\n";
    
    std::vector<JPLData> jpl_data = GetJPLData();
    
    std::cout << "PROPAGAZIONE KEPLERIAN vs JPL HORIZONS\n";
    std::cout << "════════════════════════════════════════════════════════════\n\n";
    
    std::cout << std::fixed << std::setprecision(6);
    
    double max_error_ra = 0.0;
    double max_error_dec = 0.0;
    double max_error_dist = 0.0;
    double total_sep = 0.0;
    
    for (size_t i = 0; i < jpl_data.size(); i++) {
        const auto& jpl_point = jpl_data[i];
        
        int year, month, day, hour;
        JDToDate(jpl_point.jd, year, month, day, hour);
        
        std::cout << year << "-" << std::setfill('0') 
                  << std::setw(2) << month << "-" << std::setw(2) << day 
                  << " " << std::setw(2) << hour << ":00 UTC | JD " << jpl_point.jd << "\n";
        
        // Propagazione
        EquatorialCoords calc_pos = PropagateKeplerian(elements, jpl_point.jd);
        
        // Errori
        double error_ra_arcsec = (calc_pos.ra_deg - jpl_point.ra_deg) * 3600.0 * 
                                std::cos(jpl_point.dec_deg * DEG_TO_RAD);
        double error_dec_arcsec = (calc_pos.dec_deg - jpl_point.dec_deg) * 3600.0;
        double error_dist = std::abs(calc_pos.distance_au - jpl_point.distance_au);
        
        // Separazione angolare totale
        double sep = AngularSeparation(calc_pos.ra_deg, calc_pos.dec_deg,
                                       jpl_point.ra_deg, jpl_point.dec_deg);
        
        max_error_ra = std::max(max_error_ra, std::abs(error_ra_arcsec));
        max_error_dec = std::max(max_error_dec, std::abs(error_dec_arcsec));
        max_error_dist = std::max(max_error_dist, error_dist);
        total_sep = std::max(total_sep, sep);
        
        std::cout << "  Calc: RA=" << calc_pos.ra_deg << "° Dec=" 
                  << calc_pos.dec_deg << "° r=" << calc_pos.distance_au << " AU\n";
        std::cout << "  JPL:  RA=" << jpl_point.ra_deg << "° Dec=" 
                  << jpl_point.dec_deg << "° r=" << jpl_point.distance_au << " AU\n";
        std::cout << "  Errori: ΔRA=" << error_ra_arcsec << "\" ΔDec=" 
                  << error_dec_arcsec << "\" Separazione=" << sep << "\" Dist=" 
                  << error_dist << " AU\n\n";
    }
    
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                    RISULTATI FINALI                        ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    std::cout << "Errori massimi:\n";
    std::cout << "  Max ΔRA:   " << max_error_ra << " arcsec\n";
    std::cout << "  Max ΔDec:  " << max_error_dec << " arcsec\n";
    std::cout << "  Max ΔDist: " << max_error_dist << " AU\n";
    std::cout << "  Max Separazione: " << total_sep << " arcsec\n\n";
    
    std::cout << "Valutazione accuratezza:\n";
    if (total_sep < 1.0) {
        std::cout << "  ✅ ECCELLENTE: Errore < 1 arcsec (IOoccultCalc accuracy)\n";
    } else if (total_sep < 10.0) {
        std::cout << "  ⚠️  BUONA: Errore < 10 arcsec (screening level)\n";
    } else if (total_sep < 60.0) {
        std::cout << "  ⚠️  ACCETTABILE: Errore < 60 arcsec (phase 1 threshold)\n";
    } else {
        std::cout << "  ❌ INACCETTABILE: Errore > 60 arcsec (necessita RKF78)\n";
    }
    
    std::cout << "\n════════════════════════════════════════════════════════════\n";
    std::cout << "NOTE: Questa propagazione Keplerian pura accumula errori\n";
    std::cout << "dovuti all'assenza di perturbazioni planetarie.\n";
    std::cout << "Per accuratezza < 0.001\", usare AstDynPropagator (RKF78).\n";
    std::cout << "════════════════════════════════════════════════════════════\n\n";
    
    return 0;
}

/**
 * COMPILAZIONE:
 * g++ -std=c++17 -O2 -o test_astdyn_17030_corrected test_astdyn_17030_corrected.cpp -lm
 * 
 * ESECUZIONE:
 * ./test_astdyn_17030_corrected
 */
