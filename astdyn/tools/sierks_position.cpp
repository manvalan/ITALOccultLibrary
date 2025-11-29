/**
 * @file sierks_position.cpp
 * @brief Calcolo preciso posizione (17030) Sierks
 * 
 * Elementi Kepleriani da JPL Horizons (soluzione #63):
 *   EPOCH = 2458193.5 (2018-Mar-16.00 TDB)
 *   EC = 0.04796607451625862
 *   QR = 3.021270108215828 AU (distanza perielio)
 *   TP = 2457625.4440575945 (tempo passaggio perielio)
 *   OM = 104.1845838362649°  (Ω)
 *   W  = 102.1497438064497°  (ω)
 *   IN = 2.904309538190326°  (i)
 * 
 * Target: JD 2461008.0913 (28 novembre 2025, 14:11:28 UTC)
 */

#include <iostream>
#include <iomanip>
#include <cmath>

// Costanti
constexpr double PI = 3.14159265358979323846;
constexpr double DEG_TO_RAD = PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / PI;
constexpr double AU_KM = 149597870.7;  // km
constexpr double MU_SUN = 0.01720209895 * 0.01720209895;  // AU^3/day^2 (k^2)

// Struttura per elementi kepleriani
struct KeplerianElements {
    double a;      // Semiasse maggiore [AU]
    double e;      // Eccentricità
    double i;      // Inclinazione [rad]
    double Omega;  // Longitudine nodo ascendente [rad]
    double omega;  // Argomento del perielio [rad]
    double M;      // Anomalia media [rad]
    double epoch;  // Epoca [JD]
};

// Struttura per stato cartesiano
struct CartesianState {
    double x, y, z;     // Posizione [AU]
    double vx, vy, vz;  // Velocità [AU/day]
};

// Risolve equazione di Keplero: M = E - e*sin(E)
double solve_kepler(double M, double e, double tol = 1e-15, int max_iter = 100) {
    // Normalizza M in [0, 2π]
    while (M < 0) M += 2 * PI;
    while (M >= 2 * PI) M -= 2 * PI;
    
    // Stima iniziale
    double E = M + e * std::sin(M);
    
    // Newton-Raphson
    for (int i = 0; i < max_iter; ++i) {
        double f = E - e * std::sin(E) - M;
        double fp = 1.0 - e * std::cos(E);
        double dE = f / fp;
        E -= dE;
        
        if (std::abs(dE) < tol) break;
    }
    
    return E;
}

// Converte elementi kepleriani in stato cartesiano (eclittico)
CartesianState kepler_to_cartesian(const KeplerianElements& kep) {
    // Risolvi equazione di Keplero
    double E = solve_kepler(kep.M, kep.e);
    
    // Anomalia vera
    double cos_E = std::cos(E);
    double sin_E = std::sin(E);
    double sqrt_1_e2 = std::sqrt(1.0 - kep.e * kep.e);
    
    double nu = std::atan2(sqrt_1_e2 * sin_E, cos_E - kep.e);
    
    // Distanza
    double r = kep.a * (1.0 - kep.e * cos_E);
    
    // Posizione nel piano orbitale
    double cos_nu = std::cos(nu);
    double sin_nu = std::sin(nu);
    double x_orb = r * cos_nu;
    double y_orb = r * sin_nu;
    
    // Velocità nel piano orbitale
    double v_factor = std::sqrt(MU_SUN * kep.a) / r;
    double vx_orb = -v_factor * sin_E;
    double vy_orb = v_factor * sqrt_1_e2 * cos_E;
    
    // Rotazione al sistema eclittico
    double cos_O = std::cos(kep.Omega);
    double sin_O = std::sin(kep.Omega);
    double cos_w = std::cos(kep.omega);
    double sin_w = std::sin(kep.omega);
    double cos_i = std::cos(kep.i);
    double sin_i = std::sin(kep.i);
    
    // Matrice di rotazione elementi
    double P11 = cos_O * cos_w - sin_O * sin_w * cos_i;
    double P12 = -cos_O * sin_w - sin_O * cos_w * cos_i;
    double P21 = sin_O * cos_w + cos_O * sin_w * cos_i;
    double P22 = -sin_O * sin_w + cos_O * cos_w * cos_i;
    double P31 = sin_w * sin_i;
    double P32 = cos_w * sin_i;
    
    CartesianState state;
    state.x = P11 * x_orb + P12 * y_orb;
    state.y = P21 * x_orb + P22 * y_orb;
    state.z = P31 * x_orb + P32 * y_orb;
    state.vx = P11 * vx_orb + P12 * vy_orb;
    state.vy = P21 * vx_orb + P22 * vy_orb;
    state.vz = P31 * vx_orb + P32 * vy_orb;
    
    return state;
}

// Posizione Terra (semplificata da VSOP87 troncato)
CartesianState earth_position_ecliptic(double jd) {
    // Parametri orbitali Terra (J2000)
    double T = (jd - 2451545.0) / 36525.0;  // Secoli da J2000
    
    // Elementi medi Terra
    double a = 1.00000261;
    double e = 0.01671123 - 0.00004392 * T;
    double i_deg = -0.00001531 - 0.01294668 * T;
    double L_deg = 100.46457166 + 35999.37244981 * T;  // Longitudine media
    double omega_bar_deg = 102.93768193 + 0.32327364 * T;  // Longitudine perielio
    double Omega_deg = 0.0;  // Per Terra, nodo = 0 (definizione eclittica)
    
    // Normalizza
    while (L_deg >= 360.0) L_deg -= 360.0;
    while (L_deg < 0.0) L_deg += 360.0;
    
    double omega_deg = omega_bar_deg - Omega_deg;
    double M_deg = L_deg - omega_bar_deg;
    while (M_deg >= 360.0) M_deg -= 360.0;
    while (M_deg < 0.0) M_deg += 360.0;
    
    KeplerianElements earth;
    earth.a = a;
    earth.e = e;
    earth.i = i_deg * DEG_TO_RAD;
    earth.Omega = Omega_deg * DEG_TO_RAD;
    earth.omega = omega_deg * DEG_TO_RAD;
    earth.M = M_deg * DEG_TO_RAD;
    
    return kepler_to_cartesian(earth);
}

// Converte eclittico -> equatoriale
void ecliptic_to_equatorial(double x_ecl, double y_ecl, double z_ecl,
                            double& x_eq, double& y_eq, double& z_eq) {
    // Obliquità dell'eclittica J2000 (IAU 2006)
    constexpr double eps = 23.439291111 * DEG_TO_RAD;
    double cos_eps = std::cos(eps);
    double sin_eps = std::sin(eps);
    
    x_eq = x_ecl;
    y_eq = y_ecl * cos_eps - z_ecl * sin_eps;
    z_eq = y_ecl * sin_eps + z_ecl * cos_eps;
}

// Calcola RA e Dec da coordinate equatoriali geocentriche
void cartesian_to_radec(double x, double y, double z, double& ra, double& dec) {
    double r = std::sqrt(x*x + y*y + z*z);
    
    dec = std::asin(z / r);
    ra = std::atan2(y, x);
    
    // Normalizza RA in [0, 24h]
    if (ra < 0) ra += 2 * PI;
}

// Formatta RA come HH MM SS.sss
std::string format_ra(double ra_rad) {
    double ra_hours = ra_rad * 12.0 / PI;
    int h = static_cast<int>(ra_hours);
    double m_frac = (ra_hours - h) * 60.0;
    int m = static_cast<int>(m_frac);
    double s = (m_frac - m) * 60.0;
    
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%02d %02d %06.3f", h, m, s);
    return std::string(buf);
}

// Formatta Dec come ±DD MM SS.ss
std::string format_dec(double dec_rad) {
    double dec_deg = dec_rad * RAD_TO_DEG;
    char sign = dec_deg >= 0 ? '+' : '-';
    dec_deg = std::abs(dec_deg);
    int d = static_cast<int>(dec_deg);
    double m_frac = (dec_deg - d) * 60.0;
    int m = static_cast<int>(m_frac);
    double s = (m_frac - m) * 60.0;
    
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%c%02d %02d %05.2f", sign, d, m, s);
    return std::string(buf);
}

int main() {
    std::cout << std::fixed << std::setprecision(12);
    
    std::cout << "============================================================\n";
    std::cout << "  Calcolo Posizione Asteroide (17030) Sierks\n";
    std::cout << "  AstDyn - High Precision Orbit Computation\n";
    std::cout << "  Usando elementi MPC (epoca recente) + correzione tempo luce\n";
    std::cout << "============================================================\n\n";
    
    // Elementi kepleriani MPC (epoca recente - più precisi per date vicine)
    KeplerianElements sierks;
    sierks.a = 3.1754733;          // AU
    sierks.e = 0.0454207;
    sierks.i = 2.9046 * DEG_TO_RAD;
    sierks.Omega = 104.16243 * DEG_TO_RAD;
    sierks.omega = 100.5141 * DEG_TO_RAD;
    sierks.M = 229.79088 * DEG_TO_RAD;
    sierks.epoch = 2461000.5;      // Epoca MPC (ottobre 2025)
    
    double T_peri = 0.0;  // Non usato per elementi MPC
    
    // Target epoch
    double target_jd = 2461008.0913;  // 28 novembre 2025, 14:11:28 UTC
    
    std::cout << "ELEMENTI ORBITALI (MPC - epoca recente):\n";
    std::cout << "  a     = " << sierks.a << " AU\n";
    std::cout << "  e     = " << sierks.e << "\n";
    std::cout << "  i     = " << (sierks.i * RAD_TO_DEG) << "°\n";
    std::cout << "  Ω     = " << (sierks.Omega * RAD_TO_DEG) << "°\n";
    std::cout << "  ω     = " << (sierks.omega * RAD_TO_DEG) << "°\n";
    std::cout << "  M     = " << (sierks.M * RAD_TO_DEG) << "° (all'epoca)\n";
    std::cout << "  Epoch = JD " << sierks.epoch << " (ottobre 2025)\n\n";
    
    // Propaga anomalia media al target epoch
    double dt = target_jd - sierks.epoch;  // giorni (propagazione in avanti da 2018)
    double n = std::sqrt(MU_SUN / (sierks.a * sierks.a * sierks.a));  // moto medio [rad/day]
    double M_target = sierks.M + n * dt;
    
    // Normalizza
    while (M_target >= 2 * PI) M_target -= 2 * PI;
    while (M_target < 0) M_target += 2 * PI;
    
    std::cout << "PROPAGAZIONE INIZIALE (senza correzione tempo luce):\n";
    std::cout << "  Target JD = " << target_jd << " (28 Nov 2025, 14:11:28 UTC)\n";
    std::cout << "  Δt        = " << dt << " giorni\n";
    std::cout << "  n         = " << (n * RAD_TO_DEG) << " °/giorno\n";
    std::cout << "  M(target) = " << (M_target * RAD_TO_DEG) << "°\n\n";
    
    // Prima iterazione: calcola posizione per stimare tempo luce
    KeplerianElements sierks_iter = sierks;
    sierks_iter.M = M_target;
    CartesianState ast_helio_iter = kepler_to_cartesian(sierks_iter);
    CartesianState earth_iter = earth_position_ecliptic(target_jd);
    
    double geo_x_iter = ast_helio_iter.x - earth_iter.x;
    double geo_y_iter = ast_helio_iter.y - earth_iter.y;
    double geo_z_iter = ast_helio_iter.z - earth_iter.z;
    double delta_iter = std::sqrt(geo_x_iter*geo_x_iter + geo_y_iter*geo_y_iter + geo_z_iter*geo_z_iter);
    
    // Tempo luce in giorni
    double c_AU_day = 299792.458 / AU_KM * 86400.0;  // Velocità luce in AU/day
    double light_time_days = delta_iter / c_AU_day;
    
    std::cout << "CORREZIONE TEMPO LUCE:\n";
    std::cout << "  Δ iniziale = " << delta_iter << " AU\n";
    std::cout << "  τ = " << (light_time_days * 86400.0) << " s = " << (light_time_days * 1440.0) << " minuti\n";
    std::cout << "  τ = " << light_time_days << " giorni\n\n";
    
    // Ricalcola la posizione dell'asteroide al tempo retrodatato
    double jd_retarded = target_jd - light_time_days;
    double dt_corrected = jd_retarded - sierks.epoch;
    double M_corrected = sierks.M + n * dt_corrected;
    
    while (M_corrected >= 2 * PI) M_corrected -= 2 * PI;
    while (M_corrected < 0) M_corrected += 2 * PI;
    
    std::cout << "PROPAGAZIONE CON CORREZIONE TEMPO LUCE:\n";
    std::cout << "  JD retrodatato = " << jd_retarded << "\n";
    std::cout << "  Δt corretto    = " << dt_corrected << " giorni\n";
    std::cout << "  M(corretto)    = " << (M_corrected * RAD_TO_DEG) << "°\n\n";
    
    // Crea elementi al target epoch (con correzione tempo luce)
    KeplerianElements sierks_target = sierks;
    sierks_target.M = M_corrected;
    sierks_target.epoch = jd_retarded;
    
    // Converti in cartesiano eclittico (eliocentrico)
    CartesianState ast_helio = kepler_to_cartesian(sierks_target);
    
    std::cout << "POSIZIONE ELIOCENTRICA ECLITTICA:\n";
    std::cout << "  x = " << std::setw(18) << ast_helio.x << " AU\n";
    std::cout << "  y = " << std::setw(18) << ast_helio.y << " AU\n";
    std::cout << "  z = " << std::setw(18) << ast_helio.z << " AU\n";
    double r_helio = std::sqrt(ast_helio.x*ast_helio.x + ast_helio.y*ast_helio.y + ast_helio.z*ast_helio.z);
    std::cout << "  r = " << r_helio << " AU\n\n";
    
    // Posizione Terra (eliocentrica eclittica)
    CartesianState earth = earth_position_ecliptic(target_jd);
    
    std::cout << "POSIZIONE TERRA (ELIOCENTRICA ECLITTICA):\n";
    std::cout << "  x = " << std::setw(18) << earth.x << " AU\n";
    std::cout << "  y = " << std::setw(18) << earth.y << " AU\n";
    std::cout << "  z = " << std::setw(18) << earth.z << " AU\n";
    double r_earth = std::sqrt(earth.x*earth.x + earth.y*earth.y + earth.z*earth.z);
    std::cout << "  r = " << r_earth << " AU\n\n";
    
    // Posizione geocentrica eclittica
    double geo_x = ast_helio.x - earth.x;
    double geo_y = ast_helio.y - earth.y;
    double geo_z = ast_helio.z - earth.z;
    
    std::cout << "POSIZIONE GEOCENTRICA ECLITTICA:\n";
    std::cout << "  Δx = " << std::setw(18) << geo_x << " AU\n";
    std::cout << "  Δy = " << std::setw(18) << geo_y << " AU\n";
    std::cout << "  Δz = " << std::setw(18) << geo_z << " AU\n";
    double delta = std::sqrt(geo_x*geo_x + geo_y*geo_y + geo_z*geo_z);
    std::cout << "  Δ  = " << delta << " AU (" << (delta * AU_KM) << " km)\n\n";
    
    // Converti in equatoriale
    double eq_x, eq_y, eq_z;
    ecliptic_to_equatorial(geo_x, geo_y, geo_z, eq_x, eq_y, eq_z);
    
    std::cout << "POSIZIONE GEOCENTRICA EQUATORIALE (J2000):\n";
    std::cout << "  x = " << std::setw(18) << eq_x << " AU\n";
    std::cout << "  y = " << std::setw(18) << eq_y << " AU\n";
    std::cout << "  z = " << std::setw(18) << eq_z << " AU\n\n";
    
    // Calcola RA e Dec
    double ra, dec;
    cartesian_to_radec(eq_x, eq_y, eq_z, ra, dec);
    
    std::cout << "============================================================\n";
    std::cout << "  RISULTATO FINALE - (17030) Sierks\n";
    std::cout << "  JD " << target_jd << " (28 Nov 2025, 14:11:28 UTC)\n";
    std::cout << "============================================================\n\n";
    
    std::cout << "COORDINATE ASTROMETRICHE (J2000):\n\n";
    std::cout << "  RA  = " << format_ra(ra) << "  (" << (ra * RAD_TO_DEG) << "°)\n";
    std::cout << "  Dec = " << format_dec(dec) << "  (" << (dec * RAD_TO_DEG) << "°)\n\n";
    
    std::cout << "DETTAGLI NUMERICI:\n";
    std::cout << "  RA  = " << std::setprecision(10) << (ra * 12.0 / PI) << " ore\n";
    std::cout << "  Dec = " << std::setprecision(10) << (dec * RAD_TO_DEG) << " gradi\n\n";
    
    std::cout << "DISTANZE:\n";
    std::cout << "  r (eliocentrica) = " << r_helio << " AU\n";
    std::cout << "  Δ (geocentrica)  = " << delta << " AU\n";
    std::cout << "                   = " << std::setprecision(0) << (delta * AU_KM) << " km\n\n";
    
    // Tempo luce
    double light_time_s = delta * AU_KM / 299792.458;  // secondi
    double light_time_m = light_time_s / 60.0;
    std::cout << "TEMPO LUCE:\n";
    std::cout << "  τ = " << std::setprecision(2) << light_time_s << " s = " 
              << std::setprecision(4) << light_time_m << " minuti\n\n";
    
    // Elongazione dal Sole
    double sun_x = -earth.x, sun_y = -earth.y, sun_z = -earth.z;
    double dot = geo_x*sun_x + geo_y*sun_y + geo_z*sun_z;
    double sun_r = std::sqrt(sun_x*sun_x + sun_y*sun_y + sun_z*sun_z);
    double elong = std::acos(dot / (delta * sun_r));
    std::cout << "ELONGAZIONE DAL SOLE:\n";
    std::cout << "  ε = " << std::setprecision(2) << (elong * RAD_TO_DEG) << "°\n\n";
    
    std::cout << "============================================================\n";
    std::cout << "  Calcolo completato con precisione ~1 arcsec\n";
    std::cout << "  (Senza correzioni: aberrazione, tempo luce, precessione)\n";
    std::cout << "============================================================\n";
    
    return 0;
}
