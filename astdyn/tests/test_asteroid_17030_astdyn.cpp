/**
 * @file test_asteroid_17030_astdyn.cpp
 * @brief Test occultazione asteroide 17030 usando AstDynPropagator
 * 
 * Utilizza la classe AstDynPropagator con:
 * - Integratore RKF78 
 * - Perturbazioni planetarie (8 pianeti)
 * - Propagazione da elementi JPL Horizons (epoca 2018) al 28/11/2026
 * - Confronto con effemeridi JPL
 * - Calcolo distanza angolare da stella GAIA DR3 3411546266140512128
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <array>

// Definisci ASTDYN_LIB_MODE prima di includere astdyn_propagator
#define ASTDYN_LIB_MODE
#include "../tools/astdyn_propagator.cpp"

using namespace astdyn;

// Stella GAIA DR3 3411546266140512128 (epoca J2000.0)
const double STAR_RA_DEG = 73.4161003759929;
const double STAR_DEC_DEG = 20.3316626372542;
const double STAR_PMRA = 1.097;    // mas/yr
const double STAR_PMDEC = -0.155;  // mas/yr

struct JPLEphemeris {
    double mjd;
    double ra_deg;
    double dec_deg;
};

// Converti RA da HMS a gradi
double hms_to_deg(int h, int m, double s) {
    return (h + m/60.0 + s/3600.0) * 15.0;
}

// Converti Dec da DMS a gradi  
double dms_to_deg(int d, int m, double s) {
    double sign = (d >= 0) ? 1.0 : -1.0;
    return sign * (std::abs(d) + m/60.0 + s/3600.0);
}

// Effemeridi JPL per 17030 il 28/11/2026
std::vector<JPLEphemeris> get_jpl_ephemeris() {
    std::vector<JPLEphemeris> ephem;
    
    ephem.push_back({60642.0, hms_to_deg(10, 21, 18.80), dms_to_deg(11, 56, 16.8)});
    ephem.push_back({60642.0 + 5.0/1440.0, hms_to_deg(10, 21, 18.94), dms_to_deg(11, 56, 16.3)});
    ephem.push_back({60642.0 + 10.0/1440.0, hms_to_deg(10, 21, 19.08), dms_to_deg(11, 56, 15.7)});
    ephem.push_back({60642.0 + 15.0/1440.0, hms_to_deg(10, 21, 19.22), dms_to_deg(11, 56, 15.1)});
    ephem.push_back({60642.0 + 20.0/1440.0, hms_to_deg(10, 21, 19.35), dms_to_deg(11, 56, 14.6)});
    ephem.push_back({60642.0 + 25.0/1440.0, hms_to_deg(10, 21, 19.49), dms_to_deg(11, 56, 14.0)});
    ephem.push_back({60642.0 + 30.0/1440.0, hms_to_deg(10, 21, 19.63), dms_to_deg(11, 56, 13.4)});
    ephem.push_back({60642.0 + 35.0/1440.0, hms_to_deg(10, 21, 19.76), dms_to_deg(11, 56, 12.9)});
    ephem.push_back({60642.0 + 40.0/1440.0, hms_to_deg(10, 21, 19.90), dms_to_deg(11, 56, 12.3)});
    ephem.push_back({60642.0 + 45.0/1440.0, hms_to_deg(10, 21, 20.04), dms_to_deg(11, 56, 11.7)});
    ephem.push_back({60642.0 + 50.0/1440.0, hms_to_deg(10, 21, 20.17), dms_to_deg(11, 56, 11.2)});
    ephem.push_back({60642.0 + 55.0/1440.0, hms_to_deg(10, 21, 20.31), dms_to_deg(11, 56, 10.6)});
    ephem.push_back({60642.0 + 60.0/1440.0, hms_to_deg(10, 21, 20.45), dms_to_deg(11, 56, 10.1)});
    
    return ephem;
}

// Calcola distanza angolare tra due punti (gradi)
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

// Applica moto proprio alla stella
void apply_proper_motion(double ra_deg, double dec_deg, 
                        double pmra, double pmdec, 
                        double years, 
                        double& ra_new, double& dec_new) {
    const double MAS_TO_DEG = 1.0 / 3.6e6;
    const double DEG2RAD = M_PI / 180.0;
    ra_new = ra_deg + (pmra * MAS_TO_DEG * years) / std::cos(dec_deg * DEG2RAD);
    dec_new = dec_deg + pmdec * MAS_TO_DEG * years;
}

int main() {
    try {
        std::cout << "========================================\n";
        std::cout << " Test Asteroide 17030 con AstDynPropagator\n";
        std::cout << "========================================\n\n";
        
        // Elementi orbitali JPL Horizons (epoca 2018-Mar-16.00)
        // IAU76/J2000 helio. ecliptic osc. elements
        OrbitalElements elem;
        elem.name = "17030";
        elem.epoch = 2458193.5;  // JD
        elem.a = 3.173489964321051;  // AU
        elem.e = 0.04796607451625862;
        elem.i = 2.904309538190326;  // gradi
        elem.Omega = 104.1845838362649;  // gradi
        elem.omega = 102.1497438064497;  // gradi
        elem.M = 99.03517819281583;  // gradi (anomalia media all'epoca)
        
        std::cout << "Elementi orbitali JPL (epoca JD " << std::fixed << std::setprecision(1) 
                  << elem.epoch << " = 2018-Mar-16):\n";
        std::cout << std::setprecision(8);
        std::cout << "  a = " << elem.a << " AU\n";
        std::cout << "  e = " << elem.e << "\n";
        std::cout << "  i = " << elem.i << "°\n";
        std::cout << "  ω = " << elem.omega << "°\n";
        std::cout << "  Ω = " << elem.Omega << "°\n";
        std::cout << "  M = " << elem.M << "°\n\n";
        
        // Crea propagatore con tolleranza alta
        std::cout << "Creazione propagatore AstDynPropagator...\n";
        AstDynPropagator prop(1e-12);
        prop.usePlanets(true);     // Abilita perturbazioni planetarie
        prop.useAST17(false);      // Disabilita AST17 (troppo piccole)
        prop.useRelativity(true);  // Abilita relatività
        std::cout << "  Tolleranza: 1e-12\n";
        std::cout << "  Perturbazioni planetarie: ON\n";
        std::cout << "  Correzione relativistica: ON\n\n";
        
        // Effemeridi JPL di riferimento
        auto jpl_ephem = get_jpl_ephemeris();
        double target_jd = 2400000.5 + 60642.0;  // 28 Nov 2026 00:00 UTC
        
        // Stella GAIA con moto proprio
        double years_from_j2000 = (60642.0 - 51544.5) / 365.25;
        double star_ra, star_dec;
        apply_proper_motion(STAR_RA_DEG, STAR_DEC_DEG, STAR_PMRA, STAR_PMDEC, 
                           years_from_j2000, star_ra, star_dec);
        
        std::cout << "Stella GAIA DR3 3411546266140512128 (epoca 28/11/2026):\n";
        std::cout << std::setprecision(8);
        std::cout << "  RA  = " << star_ra << "° = " 
                  << (int)(star_ra/15) << "h " 
                  << (int)((star_ra/15 - (int)(star_ra/15))*60) << "m "
                  << std::setprecision(2)
                  << ((star_ra/15 - (int)(star_ra/15))*60 - (int)((star_ra/15 - (int)(star_ra/15))*60))*60 << "s\n";
        std::cout << std::setprecision(8);
        std::cout << "  Dec = " << star_dec << "° = "
                  << (int)star_dec << "° "
                  << (int)((star_dec - (int)star_dec)*60) << "' "
                  << std::setprecision(1)
                  << ((star_dec - (int)star_dec)*60 - (int)((star_dec - (int)star_dec)*60))*60 << "\"\n\n";
        
        std::cout << "========================================\n";
        std::cout << " Propagazione: " << elem.epoch << " → " << target_jd << "\n";
        std::cout << " Intervallo: " << std::fixed << std::setprecision(2) 
                  << (target_jd - elem.epoch) / 365.25 << " anni\n";
        std::cout << "========================================\n\n";
        
        std::cout << "Tempo UTC        RA (AstDyn)    Dec (AstDyn)    RA (JPL)        Dec (JPL)       Diff RA  Diff Dec  Dist Stella\n";
        std::cout << "---------------- --------------- --------------- --------------- --------------- -------- --------- ------------\n";
        
        double min_distance = 1e10;
        double closest_time = 0;
        
        // Propaga per ogni istante
        for (size_t idx = 0; idx < jpl_ephem.size(); idx++) {
            int minute = idx * 5;
            double jd = 2400000.5 + jpl_ephem[idx].mjd;
            
            // Propaga con AstDynPropagator
            EquatorialCoords coords = prop.propagateElements(elem, jd);
            
            double ast_ra = coords.ra;
            double ast_dec = coords.dec;
            
            // Differenze con JPL
            double diff_ra = (ast_ra - jpl_ephem[idx].ra_deg) * 3600.0;  // arcsec
            double diff_dec = (ast_dec - jpl_ephem[idx].dec_deg) * 3600.0;  // arcsec
            
            // Distanza angolare da stella
            double distance = angular_distance(ast_ra, ast_dec, star_ra, star_dec);
            double distance_arcsec = distance * 3600.0;
            
            if (distance < min_distance) {
                min_distance = distance;
                closest_time = minute;
            }
            
            // Output
            std::cout << std::setfill('0');
            std::cout << "28/11/2026 " 
                      << std::setw(2) << minute/60 << ":" 
                      << std::setw(2) << minute%60 << ":00   ";
            std::cout << std::setfill(' ');
            std::cout << std::fixed << std::setprecision(6);
            std::cout << std::setw(12) << ast_ra << "°  ";
            std::cout << std::setw(12) << ast_dec << "°  ";
            std::cout << std::setw(12) << jpl_ephem[idx].ra_deg << "°  ";
            std::cout << std::setw(12) << jpl_ephem[idx].dec_deg << "°  ";
            std::cout << std::setprecision(2);
            std::cout << std::setw(7) << diff_ra << "\"  ";
            std::cout << std::setw(8) << diff_dec << "\"  ";
            std::cout << std::setw(10) << distance_arcsec << "\"";
            
            if (distance_arcsec < 60.0) {
                std::cout << "  *** CLOSE ***";
            }
            std::cout << "\n";
        }
        
        std::cout << "\n========================================\n";
        std::cout << "RISULTATI:\n";
        std::cout << "Distanza minima da stella: " << std::fixed << std::setprecision(2) 
                  << min_distance * 3600.0 << " arcsec = " 
                  << std::setprecision(4) << min_distance << "°\n";
        std::cout << "Tempo minima distanza: ";
        std::cout << std::setfill('0') << std::setw(2) << (int)(closest_time/60) << ":" 
                  << std::setw(2) << ((int)closest_time)%60 << " UTC\n";
        std::cout << "========================================\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << std::endl;
        return 1;
    }
}
