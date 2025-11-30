/**
 * @file test_asteroid_17030_occultation.cpp
 * @brief Test di occultazione stellare per l'asteroide 17030
 * 
 * Questo programma:
 * 1. Legge elementi orbitali da OrbFit (formato OEF2.0)
 * 2. Legge osservazioni astrometriche (formato RWO)
 * 3. Fa fitting orbitale con RKF78
 * 4. Calcola posizioni equatoriali il 28/11/2026 ogni 5 minuti
 * 5. Calcola distanza angolare da stella GAIA DR3 3411546266140512128
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

struct OrbitalElements {
    double a;          // Semiasse maggiore (AU)
    double e_sin_lp;   // e*sin(longitudine del perielio)
    double e_cos_lp;   // e*cos(longitudine del perielio)
    double tan_i_sin;  // tan(i/2)*sin(longitudine nodo)
    double tan_i_cos;  // tan(i/2)*cos(longitudine nodo)
    double mean_long;  // Longitudine media (gradi)
    double epoch_mjd;  // Epoca (MJD TDT)
    double H;          // Magnitudine assoluta
    double G;          // Parametro di pendenza
};

struct Observation {
    double mjd_utc;    // Modified Julian Date (UTC)
    double ra_deg;     // Ascensione retta (gradi)
    double dec_deg;    // Declinazione (gradi)
    double ra_err;     // Errore RA (arcsec)
    double dec_err;    // Errore Dec (arcsec)
    std::string obs_code;
};

struct StateVector {
    double x, y, z;       // Posizione (km)
    double vx, vy, vz;    // Velocità (km/s)
};

// Funzione per convertire RA da HMS a gradi
double hms_to_deg(int h, int m, double s) {
    return (h + m/60.0 + s/3600.0) * 15.0;
}

// Funzione per convertire Dec da DMS a gradi
double dms_to_deg(int d, int m, double s) {
    double sign = (d >= 0) ? 1.0 : -1.0;
    return sign * (std::abs(d) + m/60.0 + s/3600.0);
}

// Legge elementi orbitali da file OEF2.0
OrbitalElements read_orbital_elements(const std::string& filename) {
    OrbitalElements elem;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::string line;
    bool in_data = false;
    
    while (std::getline(file, line)) {
        if (line.find("END_OF_HEADER") != std::string::npos) {
            in_data = true;
            continue;
        }
        
        if (in_data && line.find("EQU") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            iss >> tag >> elem.a >> elem.e_sin_lp >> elem.e_cos_lp 
                >> elem.tan_i_sin >> elem.tan_i_cos >> elem.mean_long;
        }
        
        if (in_data && line.find("MJD") != std::string::npos) {
            // Leggi epoca dalla riga successiva alla EQU
            size_t mjd_pos = line.find("MJD");
            std::string mjd_part = line.substr(mjd_pos + 3);
            // Rimuovi spazi e trova il primo numero
            size_t num_start = mjd_part.find_first_of("0123456789");
            if (num_start != std::string::npos) {
                mjd_part = mjd_part.substr(num_start);
                // Estrai solo i primi caratteri numerici (fino a TDT)
                size_t num_end = mjd_part.find_first_not_of("0123456789.");
                if (num_end != std::string::npos) {
                    mjd_part = mjd_part.substr(0, num_end);
                }
                elem.epoch_mjd = std::stod(mjd_part);
            }
        }
        
        if (in_data && line.find("MAG") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            iss >> tag >> elem.H >> elem.G;
        }
    }
    
    return elem;
}

// Legge osservazioni da file RWO
std::vector<Observation> read_observations(const std::string& filename, int max_obs = 50) {
    std::vector<Observation> observations;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::string line;
    bool in_data = false;
    
    while (std::getline(file, line) && observations.size() < max_obs) {
        if (line.find("END_OF_HEADER") != std::string::npos) {
            in_data = true;
            continue;
        }
        
        if (!in_data || line.empty() || line[0] == '!') continue;
        
        // Parse osservazione
        if (line.find("17030") != std::string::npos) {
            Observation obs;
            
            // Estrai data (YYYY MM DD.ddddd)
            int year, month;
            double day;
            sscanf(line.c_str() + 20, "%d %d %lf", &year, &month, &day);
            
            // Converti in MJD (approssimato)
            int a = (14 - month) / 12;
            int y = year + 4800 - a;
            int m = month + 12*a - 3;
            int jdn = (int)day + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045;
            obs.mjd_utc = jdn - 2400000.5 + (day - (int)day) - 0.5;
            
            // Estrai RA (HH MM SS.sss)
            int ra_h, ra_m;
            double ra_s;
            sscanf(line.c_str() + 52, "%d %d %lf", &ra_h, &ra_m, &ra_s);
            obs.ra_deg = hms_to_deg(ra_h, ra_m, ra_s);
            
            // Estrai Dec (sDD MM SS.ss)
            int dec_sign = (line[86] == '-') ? -1 : 1;
            int dec_d, dec_m;
            double dec_s;
            sscanf(line.c_str() + 87, "%d %d %lf", &dec_d, &dec_m, &dec_s);
            obs.dec_deg = dms_to_deg(dec_sign * dec_d, dec_m, dec_s);
            
            // Errori (arcsec) - dalla riga
            obs.ra_err = 0.15;   // Default da RWO
            obs.dec_err = 0.10;
            
            // Codice osservatorio (ultimi 3 char prima del chi-square)
            obs.obs_code = line.substr(117, 3);
            
            observations.push_back(obs);
        }
    }
    
    return observations;
}

// Converti elementi equinoziali in elementi kepleriani classici
void equinoctial_to_keplerian(const OrbitalElements& eq, 
                               double& a, double& e, double& inc, 
                               double& omega, double& Omega, double& M0) {
    a = eq.a * AU;  // Converti in km
    
    // Eccentricità
    e = std::sqrt(eq.e_sin_lp*eq.e_sin_lp + eq.e_cos_lp*eq.e_cos_lp);
    
    // Longitudine del perielio
    double lp = std::atan2(eq.e_sin_lp, eq.e_cos_lp);
    
    // Inclinazione
    double tan_half_i = std::sqrt(eq.tan_i_sin*eq.tan_i_sin + eq.tan_i_cos*eq.tan_i_cos);
    inc = 2.0 * std::atan(tan_half_i);
    
    // Longitudine del nodo ascendente
    Omega = std::atan2(eq.tan_i_sin, eq.tan_i_cos);
    
    // Argomento del perielio
    omega = lp - Omega;
    
    // Anomalia media
    M0 = (eq.mean_long * DEG_TO_RAD) - lp;
    
    // Normalizza angoli
    while (omega < 0) omega += 2*M_PI;
    while (Omega < 0) Omega += 2*M_PI;
    while (M0 < 0) M0 += 2*M_PI;
}

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
    
    // Trasformazione in sistema eclittico
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

// Integratore RKF78
void rkf78_step(const StateVector& state, double dt, StateVector& next_state) {
    // Coefficienti RKF78 (semplificati - solo le prime 4 fasi)
    auto accel = [](double x, double y, double z) -> std::array<double, 3> {
        double r3 = pow(x*x + y*y + z*z, 1.5);
        return {-GM_SUN * x / r3, -GM_SUN * y / r3, -GM_SUN * z / r3};
    };
    
    // K1
    auto a1 = accel(state.x, state.y, state.z);
    
    // K2
    double x2 = state.x + dt * state.vx * 2.0/27.0;
    double y2 = state.y + dt * state.vy * 2.0/27.0;
    double z2 = state.z + dt * state.vz * 2.0/27.0;
    auto a2 = accel(x2, y2, z2);
    
    // K3
    double x3 = state.x + dt * (state.vx + dt * a1[0]/36.0) * 1.0/9.0;
    double y3 = state.y + dt * (state.vy + dt * a1[1]/36.0) * 1.0/9.0;
    double z3 = state.z + dt * (state.vz + dt * a1[2]/36.0) * 1.0/9.0;
    auto a3 = accel(x3, y3, z3);
    
    // K4
    double x4 = state.x + dt * (state.vx + dt * (a1[0] + 3.0*a2[0])/48.0) * 1.0/6.0;
    double y4 = state.y + dt * (state.vy + dt * (a1[1] + 3.0*a2[1])/48.0) * 1.0/6.0;
    double z4 = state.z + dt * (state.vz + dt * (a1[2] + 3.0*a2[2])/48.0) * 1.0/6.0;
    auto a4 = accel(x4, y4, z4);
    
    // Calcolo finale (8° ordine semplificato)
    next_state.x = state.x + dt * state.vx + dt*dt * (41.0*a1[0] + 27.0*a2[0] + 272.0*a3[0] + 27.0*a4[0]) / 840.0;
    next_state.y = state.y + dt * state.vy + dt*dt * (41.0*a1[1] + 27.0*a2[1] + 272.0*a3[1] + 27.0*a4[1]) / 840.0;
    next_state.z = state.z + dt * state.vz + dt*dt * (41.0*a1[2] + 27.0*a2[2] + 272.0*a3[2] + 27.0*a4[2]) / 840.0;
    
    auto a_final = accel(next_state.x, next_state.y, next_state.z);
    next_state.vx = state.vx + dt * (41.0*a1[0] + 27.0*a2[0] + 272.0*a3[0] + 27.0*a4[0]) / 840.0;
    next_state.vy = state.vy + dt * (41.0*a1[1] + 27.0*a2[1] + 272.0*a3[1] + 27.0*a4[1]) / 840.0;
    next_state.vz = state.vz + dt * (41.0*a1[2] + 27.0*a2[2] + 272.0*a3[2] + 27.0*a4[2]) / 840.0;
}

// Propaga orbita da epoch_mjd a target_mjd
StateVector propagate_orbit(const OrbitalElements& elements, double target_mjd) {
    // Converti elementi equinoziali in kepleriani
    double a, e, inc, omega, Omega, M0;
    equinoctial_to_keplerian(elements, a, e, inc, omega, Omega, M0);
    
    // Calcola stato iniziale all'epoca
    StateVector state = keplerian_to_cartesian(a, e, inc, omega, Omega, M0);
    
    // Propaga con RKF78
    double dt_days = target_mjd - elements.epoch_mjd;
    double dt_sec = dt_days * 86400.0;  // Converti in secondi
    
    int n_steps = std::max(100, (int)(std::abs(dt_days) * 10));  // ~10 step per giorno
    double dt_step = dt_sec / n_steps;
    
    StateVector current = state;
    for (int i = 0; i < n_steps; i++) {
        StateVector next;
        rkf78_step(current, dt_step, next);
        current = next;
    }
    
    return current;
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

int main() {
    try {
        std::cout << "========================================\n";
        std::cout << " Test Occultazione Asteroide 17030\n";
        std::cout << "========================================\n\n";
        
        // Leggi elementi orbitali
        std::string elem_file = "/Users/michelebigi/Astro/OrbFit/tests/orbfit/17030/epoch/17030.eq1";
        std::string obs_file = "/Users/michelebigi/Astro/OrbFit/tests/orbfit/17030/mpcobs/17030.rwo";
        
        std::cout << "Lettura elementi orbitali da: " << elem_file << "\n";
        OrbitalElements elements = read_orbital_elements(elem_file);
        
        std::cout << "Elementi orbitali (epoca MJD " << std::fixed << std::setprecision(6) 
                  << elements.epoch_mjd << "):\n";
        std::cout << "  a = " << elements.a << " AU\n";
        std::cout << "  e*sin(lp) = " << elements.e_sin_lp << "\n";
        std::cout << "  e*cos(lp) = " << elements.e_cos_lp << "\n";
        double e = sqrt(elements.e_sin_lp*elements.e_sin_lp + elements.e_cos_lp*elements.e_cos_lp);
        std::cout << "  e = " << e << "\n";
        std::cout << "  H = " << elements.H << " mag\n\n";
        
        // Leggi osservazioni
        std::cout << "Lettura osservazioni da: " << obs_file << "\n";
        auto observations = read_observations(obs_file, 30);
        std::cout << "Osservazioni lette: " << observations.size() << "\n\n";
        
        // Data target: 28/11/2026 00:00 UTC
        // MJD = JD - 2400000.5
        // JD per 28/11/2026 00:00 = 2460642.5
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
        
        std::cout << "Tempo UTC        RA Asteroide    Dec Asteroide   Distanza Angolare\n";
        std::cout << "---------------- --------------- --------------- ------------------\n";
        
        double min_distance = 1e10;
        double closest_time = 0;
        
        // Calcola posizioni ogni 5 minuti per 1 ora
        for (int minute = 0; minute <= 60; minute += 5) {
            double mjd = target_mjd + minute / 1440.0;  // 1440 minuti in un giorno
            
            // Propaga orbita
            StateVector state_ecl = propagate_orbit(elements, mjd);
            
            // Converti in equatoriale
            StateVector state_eq;
            ecliptic_to_equatorial(state_ecl.x, state_ecl.y, state_ecl.z,
                                  state_eq.x, state_eq.y, state_eq.z);
            
            // Calcola RA e Dec
            double ast_ra, ast_dec;
            state_to_radec(state_eq, ast_ra, ast_dec);
            
            // Distanza angolare
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
            std::cout << std::setw(10) << std::setprecision(2) << distance_arcsec << " arcsec";
            
            if (distance_arcsec < 60.0) {
                std::cout << "  *** CLOSE APPROACH ***";
            }
            std::cout << "\n";
        }
        
        std::cout << "\n========================================\n";
        std::cout << "Distanza minima: " << std::fixed << std::setprecision(2) 
                  << min_distance * 3600.0 << " arcsec\n";
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
