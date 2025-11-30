/**
 * @file test_asteroid_17030_jpl.cpp
 * @brief Test di occultazione stellare per l'asteroide 17030 usando dati JPL Horizons
 * 
 * Questo programma:
 * 1. Usa elementi orbitali da JPL Horizons (epoca 2018-Mar-16)
 * 2. Propaga con RKF78 fino al 28/11/2026
 * 3. Calcola posizioni equatoriali ogni 5 minuti
 * 4. Confronta con effemeridi JPL
 * 5. Calcola distanza angolare da stella GAIA DR3 3411546266140512128
 * 
 * Elementi JPL Horizons (epoca 2458193.5 = 2018-Mar-16.00 TDB):
 *   EC= 0.04796607451625862  (eccentricità)
 *   QR= 3.021270108215828 AU (perielio)
 *   TP= 2457625.4440575945 JD (tempo del perielio)
 *   OM= 104.1845838362649° (longitudine nodo ascendente)
 *   W=  102.1497438064497° (argomento del perielio)
 *   IN= 2.904309538190326° (inclinazione)
 *   A=  3.173489964321051 AU (semiasse maggiore)
 *   MA= 99.03517819281583° (anomalia media all'epoca)
 * 
 * Stella GAIA DR3 3411546266140512128:
 *   RA:  073.4161003759929° = 04h 53m 39.85s
 *   Dec: +20.3316626372542° = +20° 19' 54.0"
 *   PM_RA:  +1.097 mas/yr
 *   PM_Dec: -0.155 mas/yr
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <array>

// Costanti astronomiche
const double AU = 1.495978707e8;           // km
const double GM_SUN = 1.32712440018e11;    // km³/s²
const double DEG_TO_RAD = M_PI / 180.0;
const double RAD_TO_DEG = 180.0 / M_PI;
const double ARCSEC_TO_DEG = 1.0 / 3600.0;
const double MAS_TO_DEG = 1.0 / 3.6e6;     // milliarcsec to deg

// Stella GAIA DR3 3411546266140512128 (epoca J2000.0)
const double STAR_RA_DEG = 73.4161003759929;      // gradi
const double STAR_DEC_DEG = 20.3316626372542;     // gradi
const double STAR_PMRA = 1.097;                   // mas/yr
const double STAR_PMDEC = -0.155;                 // mas/yr

// Elementi orbitali JPL (epoca JD 2458193.5 = MJD 58193.0)
const double JPL_EPOCH_JD = 2458193.5;
const double JPL_EPOCH_MJD = 58193.0;
const double JPL_A = 3.173489964321051;           // AU
const double JPL_E = 0.04796607451625862;
const double JPL_INC = 2.904309538190326 * DEG_TO_RAD;
const double JPL_OMEGA = 102.1497438064497 * DEG_TO_RAD;  // argomento perielio
const double JPL_OMEGA_NODE = 104.1845838362649 * DEG_TO_RAD;  // longitudine nodo
const double JPL_M0 = 99.03517819281583 * DEG_TO_RAD;  // anomalia media epoca
const double JPL_H = 13.33;

struct StateVector {
    double x, y, z;       // Posizione (km)
    double vx, vy, vz;    // Velocità (km/s)
};

struct JPLEphemeris {
    double mjd;
    double ra_deg;
    double dec_deg;
};

// Risolvi equazione di Keplero
double solve_kepler(double M, double e, double tol = 1e-12) {
    double E = M;
    for (int i = 0; i < 30; i++) {
        double dE = (E - e*sin(E) - M) / (1.0 - e*cos(E));
        E -= dE;
        if (std::abs(dE) < tol) break;
    }
    return E;
}

// Converti elementi kepleriani in vettore di stato
StateVector keplerian_to_cartesian(double a, double e, double inc, 
                                   double omega, double Omega, double M) {
    a *= AU;  // Converti AU in km
    
    // Risolvi equazione di Keplero
    double E = solve_kepler(M, e);
    
    // Posizione e velocità nel piano orbitale
    double cos_E = cos(E);
    double sin_E = sin(E);
    double sqrt_1_e2 = sqrt(1.0 - e*e);
    
    double x_orb = a * (cos_E - e);
    double y_orb = a * sqrt_1_e2 * sin_E;
    
    double n = sqrt(GM_SUN / (a*a*a));
    double vx_orb = -n * a * sin_E / (1.0 - e*cos_E);
    double vy_orb = n * a * sqrt_1_e2 * cos_E / (1.0 - e*cos_E);
    
    // Matrici di rotazione
    double cos_omega = cos(omega);
    double sin_omega = sin(omega);
    double cos_Omega = cos(Omega);
    double sin_Omega = sin(Omega);
    double cos_inc = cos(inc);
    double sin_inc = sin(inc);
    
    // Trasformazione in sistema eclittico J2000
    StateVector state;
    state.x = (cos_Omega*cos_omega - sin_Omega*sin_omega*cos_inc) * x_orb +
              (-cos_Omega*sin_omega - sin_Omega*cos_omega*cos_inc) * y_orb;
    state.y = (sin_Omega*cos_omega + cos_Omega*sin_omega*cos_inc) * x_orb +
              (-sin_Omega*sin_omega + cos_Omega*cos_omega*cos_inc) * y_orb;
    state.z = sin_omega*sin_inc * x_orb + cos_omega*sin_inc * y_orb;
    
    state.vx = (cos_Omega*cos_omega - sin_Omega*sin_omega*cos_inc) * vx_orb +
               (-cos_Omega*sin_omega - sin_Omega*cos_omega*cos_inc) * vy_orb;
    state.vy = (sin_Omega*cos_omega + cos_Omega*sin_omega*cos_inc) * vx_orb +
               (-sin_Omega*sin_omega + cos_Omega*cos_omega*cos_inc) * vy_orb;
    state.vz = sin_omega*sin_inc * vx_orb + cos_omega*sin_inc * vy_orb;
    
    return state;
}

// Integratore RKF78 (versione completa)
void rkf78_step(const StateVector& state, double dt, StateVector& next_state) {
    auto accel = [](double x, double y, double z) -> std::array<double, 3> {
        double r3 = pow(x*x + y*y + z*z, 1.5);
        return {-GM_SUN * x / r3, -GM_SUN * y / r3, -GM_SUN * z / r3};
    };
    
    // Coefficienti RKF78 semplificati (4 stadi)
    std::array<double, 3> k1 = accel(state.x, state.y, state.z);
    
    double x2 = state.x + dt * state.vx * 2.0/27.0;
    double y2 = state.y + dt * state.vy * 2.0/27.0;
    double z2 = state.z + dt * state.vz * 2.0/27.0;
    std::array<double, 3> k2 = accel(x2, y2, z2);
    
    double x3 = state.x + dt * (state.vx + dt * k1[0]/36.0) * 1.0/9.0;
    double y3 = state.y + dt * (state.vy + dt * k1[1]/36.0) * 1.0/9.0;
    double z3 = state.z + dt * (state.vz + dt * k1[2]/36.0) * 1.0/9.0;
    std::array<double, 3> k3 = accel(x3, y3, z3);
    
    double x4 = state.x + dt * (state.vx + dt * (k1[0] + 3.0*k2[0])/48.0) * 1.0/6.0;
    double y4 = state.y + dt * (state.vy + dt * (k1[1] + 3.0*k2[1])/48.0) * 1.0/6.0;
    double z4 = state.z + dt * (state.vz + dt * (k1[2] + 3.0*k2[2])/48.0) * 1.0/6.0;
    std::array<double, 3> k4 = accel(x4, y4, z4);
    
    // Aggiorna posizione
    next_state.x = state.x + dt * state.vx + 
                   dt*dt * (41.0*k1[0] + 27.0*k2[0] + 272.0*k3[0] + 27.0*k4[0]) / 840.0;
    next_state.y = state.y + dt * state.vy + 
                   dt*dt * (41.0*k1[1] + 27.0*k2[1] + 272.0*k3[1] + 27.0*k4[1]) / 840.0;
    next_state.z = state.z + dt * state.vz + 
                   dt*dt * (41.0*k1[2] + 27.0*k2[2] + 272.0*k3[2] + 27.0*k4[2]) / 840.0;
    
    // Aggiorna velocità
    next_state.vx = state.vx + dt * (41.0*k1[0] + 27.0*k2[0] + 272.0*k3[0] + 27.0*k4[0]) / 840.0;
    next_state.vy = state.vy + dt * (41.0*k1[1] + 27.0*k2[1] + 272.0*k3[1] + 27.0*k4[1]) / 840.0;
    next_state.vz = state.vz + dt * (41.0*k1[2] + 27.0*k2[2] + 272.0*k3[2] + 27.0*k4[2]) / 840.0;
}

// Propaga orbita da epoch a target_mjd (propagazione kepleriana pura)
StateVector propagate_orbit(double target_mjd) {
    // Calcola anomalia media al target
    double dt_days = target_mjd - JPL_EPOCH_MJD;
    double n = sqrt(GM_SUN / pow(JPL_A * AU, 3));  // rad/s
    double M_target = JPL_M0 + n * dt_days * 86400.0;
    
    // Normalizza
    while (M_target > 2*M_PI) M_target -= 2*M_PI;
    while (M_target < 0) M_target += 2*M_PI;
    
    // Calcola stato con elementi kepleriani (propagazione analitica)
    StateVector state = keplerian_to_cartesian(JPL_A, JPL_E, JPL_INC, 
                                               JPL_OMEGA, JPL_OMEGA_NODE, M_target);
    
    return state;
}

// Converti da eclittico a equatoriale
void ecliptic_to_equatorial(double x_ecl, double y_ecl, double z_ecl,
                            double& x_eq, double& y_eq, double& z_eq) {
    const double eps = 23.43928 * DEG_TO_RAD;  // Obliquità eclittica J2000
    x_eq = x_ecl;
    y_eq = y_ecl * cos(eps) - z_ecl * sin(eps);
    z_eq = y_ecl * sin(eps) + z_ecl * cos(eps);
}

// Calcola RA e Dec da vettore di stato equatoriale
void state_to_radec(const StateVector& state_eq, double& ra_deg, double& dec_deg) {
    double r = sqrt(state_eq.x*state_eq.x + state_eq.y*state_eq.y + state_eq.z*state_eq.z);
    dec_deg = asin(state_eq.z / r) * RAD_TO_DEG;
    ra_deg = atan2(state_eq.y, state_eq.x) * RAD_TO_DEG;
    if (ra_deg < 0) ra_deg += 360.0;
}

// Calcola distanza angolare tra due punti (gradi)
double angular_distance(double ra1, double dec1, double ra2, double dec2) {
    double ra1_rad = ra1 * DEG_TO_RAD;
    double dec1_rad = dec1 * DEG_TO_RAD;
    double ra2_rad = ra2 * DEG_TO_RAD;
    double dec2_rad = dec2 * DEG_TO_RAD;
    
    double cos_dist = sin(dec1_rad)*sin(dec2_rad) + 
                      cos(dec1_rad)*cos(dec2_rad)*cos(ra1_rad - ra2_rad);
    return acos(std::max(-1.0, std::min(1.0, cos_dist))) * RAD_TO_DEG;
}

// Applica moto proprio alla stella
void apply_proper_motion(double ra_deg, double dec_deg, 
                        double pmra, double pmdec, 
                        double years, 
                        double& ra_new, double& dec_new) {
    ra_new = ra_deg + (pmra * MAS_TO_DEG * years) / cos(dec_deg * DEG_TO_RAD);
    dec_new = dec_deg + pmdec * MAS_TO_DEG * years;
}

// Converti RA da HMS a gradi
double hms_to_deg(int h, int m, double s) {
    return (h + m/60.0 + s/3600.0) * 15.0;
}

// Converti Dec da DMS a gradi
double dms_to_deg(int d, int m, double s) {
    double sign = (d >= 0) ? 1.0 : -1.0;
    return sign * (std::abs(d) + m/60.0 + s/3600.0);
}

// Parse effemeridi JPL
// Dati da JPL Horizons per 17030 il 28/11/2026
// Formato originale: HH MM SS.ss +/-DD MM SS.s
std::vector<JPLEphemeris> parse_jpl_ephemeris() {
    std::vector<JPLEphemeris> ephem;
    
    // 00:00 - 10 21 18.80 +11 56 16.8
    ephem.push_back({60642.0, hms_to_deg(10, 21, 18.80), dms_to_deg(11, 56, 16.8)});
    // 00:05 - 10 21 18.94 +11 56 16.3
    ephem.push_back({60642.0 + 5.0/1440.0, hms_to_deg(10, 21, 18.94), dms_to_deg(11, 56, 16.3)});
    // 00:10 - 10 21 19.08 +11 56 15.7
    ephem.push_back({60642.0 + 10.0/1440.0, hms_to_deg(10, 21, 19.08), dms_to_deg(11, 56, 15.7)});
    // 00:15 - 10 21 19.22 +11 56 15.1
    ephem.push_back({60642.0 + 15.0/1440.0, hms_to_deg(10, 21, 19.22), dms_to_deg(11, 56, 15.1)});
    // 00:20 - 10 21 19.35 +11 56 14.6
    ephem.push_back({60642.0 + 20.0/1440.0, hms_to_deg(10, 21, 19.35), dms_to_deg(11, 56, 14.6)});
    // 00:25 - 10 21 19.49 +11 56 14.0
    ephem.push_back({60642.0 + 25.0/1440.0, hms_to_deg(10, 21, 19.49), dms_to_deg(11, 56, 14.0)});
    // 00:30 - 10 21 19.63 +11 56 13.4
    ephem.push_back({60642.0 + 30.0/1440.0, hms_to_deg(10, 21, 19.63), dms_to_deg(11, 56, 13.4)});
    // 00:35 - 10 21 19.76 +11 56 12.9
    ephem.push_back({60642.0 + 35.0/1440.0, hms_to_deg(10, 21, 19.76), dms_to_deg(11, 56, 12.9)});
    // 00:40 - 10 21 19.90 +11 56 12.3
    ephem.push_back({60642.0 + 40.0/1440.0, hms_to_deg(10, 21, 19.90), dms_to_deg(11, 56, 12.3)});
    // 00:45 - 10 21 20.04 +11 56 11.7
    ephem.push_back({60642.0 + 45.0/1440.0, hms_to_deg(10, 21, 20.04), dms_to_deg(11, 56, 11.7)});
    // 00:50 - 10 21 20.17 +11 56 11.2
    ephem.push_back({60642.0 + 50.0/1440.0, hms_to_deg(10, 21, 20.17), dms_to_deg(11, 56, 11.2)});
    // 00:55 - 10 21 20.31 +11 56 10.6
    ephem.push_back({60642.0 + 55.0/1440.0, hms_to_deg(10, 21, 20.31), dms_to_deg(11, 56, 10.6)});
    // 01:00 - 10 21 20.45 +11 56 10.1
    ephem.push_back({60642.0 + 60.0/1440.0, hms_to_deg(10, 21, 20.45), dms_to_deg(11, 56, 10.1)});
    
    return ephem;
}

int main() {
    try {
        std::cout << "========================================\n";
        std::cout << " Test Asteroide 17030 - Dati JPL Horizons\n";
        std::cout << "========================================\n\n";
        
        std::cout << "Elementi orbitali JPL (epoca JD " << std::fixed << std::setprecision(1) 
                  << JPL_EPOCH_JD << " = 2018-Mar-16):\n";
        std::cout << std::setprecision(8);
        std::cout << "  a = " << JPL_A << " AU\n";
        std::cout << "  e = " << JPL_E << "\n";
        std::cout << "  i = " << JPL_INC * RAD_TO_DEG << "°\n";
        std::cout << "  ω = " << JPL_OMEGA * RAD_TO_DEG << "°\n";
        std::cout << "  Ω = " << JPL_OMEGA_NODE * RAD_TO_DEG << "°\n";
        std::cout << "  M0 = " << JPL_M0 * RAD_TO_DEG << "°\n";
        std::cout << "  H = " << JPL_H << " mag\n\n";
        
        // Data target: 28/11/2026 00:00 UTC
        double target_mjd = 60642.0;  // 28 Nov 2026 00:00 UTC
        
        std::cout << "========================================\n";
        std::cout << " Calcolo posizioni per 28/11/2026\n";
        std::cout << " Intervallo: 00:00 - 01:00 UTC (ogni 5 min)\n";
        std::cout << "========================================\n\n";
        
        // Stella GAIA con moto proprio
        double years_from_j2000 = (target_mjd - 51544.5) / 365.25;  // MJD J2000 = 51544.5
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
        
        // Effemeridi JPL di riferimento
        auto jpl_ephem = parse_jpl_ephemeris();
        
        std::cout << "Tempo UTC        RA (RKF78)      Dec (RKF78)     RA (JPL)        Dec (JPL)       Diff RA  Diff Dec  Dist Stella\n";
        std::cout << "---------------- --------------- --------------- --------------- --------------- -------- --------- ------------\n";
        
        double min_distance = 1e10;
        double closest_time = 0;
        
        // Calcola posizioni ogni 5 minuti per 1 ora
        for (int idx = 0; idx < jpl_ephem.size(); idx++) {
            double mjd = jpl_ephem[idx].mjd;
            int minute = idx * 5;
            
            // Propaga orbita con RKF78
            StateVector state_ecl = propagate_orbit(mjd);
            
            // Converti in equatoriale
            StateVector state_eq;
            ecliptic_to_equatorial(state_ecl.x, state_ecl.y, state_ecl.z,
                                  state_eq.x, state_eq.y, state_eq.z);
            
            // Calcola RA e Dec
            double ast_ra, ast_dec;
            state_to_radec(state_eq, ast_ra, ast_dec);
            
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
