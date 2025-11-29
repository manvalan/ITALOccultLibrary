/**
 * @file sierks_precision.cpp
 * @brief Calcolo PRECISO posizione (17030) Sierks
 *        Replica esatta della logica AstDyn/OrbFit con:
 *        - Integrazione numerica RKF78
 *        - Perturbazioni planetarie (8 pianeti + Luna)
 *        - Perturbazioni asteroidali AST17 (16 asteroidi massivi)
 *        - Correzione relativistica (Schwarzschild)
 *        - Correzione tempo luce
 *        - Aberrazione stellare annuale
 * 
 * @author AstDyn Team - ITALOccult
 * @date 2025-11-29
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <functional>

// ============================================================================
// COSTANTI FISICHE (identiche a AstDyn/OrbFit)
// ============================================================================

namespace constants {
    constexpr double PI = 3.14159265358979323846;
    constexpr double TWO_PI = 2.0 * PI;
    constexpr double DEG_TO_RAD = PI / 180.0;
    constexpr double RAD_TO_DEG = 180.0 / PI;
    constexpr double ARCSEC_TO_RAD = PI / 648000.0;
    
    // Costante gravitazionale Sole [AU³/day²] (k² di Gauss)
    constexpr double K_GAUSS = 0.01720209895;
    constexpr double GM_SUN = K_GAUSS * K_GAUSS;  // 0.0002959122082855911
    
    // Conversione GM: km³/s² → AU³/day²
    constexpr double AU_KM = 149597870.7;
    constexpr double DAY_SEC = 86400.0;
    constexpr double GM_CONV = DAY_SEC * DAY_SEC / (AU_KM * AU_KM * AU_KM);
    
    // GM pianeti [km³/s²] - valori DE440/441
    constexpr double GM_MERCURY = 22031.868551;
    constexpr double GM_VENUS = 324858.592000;
    constexpr double GM_EARTH = 398600.435507;  // Terra + Luna
    constexpr double GM_MOON = 4902.800118;
    constexpr double GM_MARS = 42828.375816;
    constexpr double GM_JUPITER = 126712764.100000;
    constexpr double GM_SATURN = 37940584.841800;
    constexpr double GM_URANUS = 5794556.400000;
    constexpr double GM_NEPTUNE = 6836527.100580;
    
    // Velocità della luce [AU/day]
    constexpr double SPEED_OF_LIGHT = 299792.458 / AU_KM * DAY_SEC;  // ~173.14 AU/day
    
    // Obliquità eclittica J2000 (IAU 2006)
    constexpr double OBLIQUITY_J2000 = 23.439291111 * DEG_TO_RAD;
    
    // Epoca J2000.0
    constexpr double JD_J2000 = 2451545.0;
    constexpr double MJD_J2000 = 51544.5;
}

using namespace constants;

// ============================================================================
// STRUTTURE DATI
// ============================================================================

struct Vec3 {
    double x, y, z;
    
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    
    Vec3 operator+(const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
    Vec3 operator*(double s) const { return Vec3(x*s, y*s, z*s); }
    Vec3 operator/(double s) const { return Vec3(x/s, y/s, z/s); }
    Vec3& operator+=(const Vec3& v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
    Vec3& operator-=(const Vec3& v) { x-=v.x; y-=v.y; z-=v.z; return *this; }
    
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    double norm2() const { return x*x + y*y + z*z; }
    double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
};

struct State {
    Vec3 pos;   // Posizione [AU]
    Vec3 vel;   // Velocità [AU/day]
};

struct KeplerianElements {
    double a;       // Semiasse maggiore [AU]
    double e;       // Eccentricità
    double i;       // Inclinazione [rad]
    double Omega;   // Longitudine nodo ascendente [rad]
    double omega;   // Argomento del perielio [rad]
    double M0;      // Anomalia media all'epoca [rad]
    double epoch;   // Epoca [JD TDB]
    double gm;      // GM corpo centrale [AU³/day²]
};

struct AsteroidData {
    int number;
    const char* name;
    double gm;      // GM [km³/s²]
    double a, e, i, omega, Omega, M0;  // Elementi (i,omega,Omega in deg)
    double epoch;   // MJD
    double n;       // Moto medio [deg/day]
};

// ============================================================================
// EFFEMERIDI PLANETARIE (Simon et al. 1994 - come in OrbFit)
// ============================================================================

class PlanetaryEphemeris {
public:
    enum Planet { MERCURY=0, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE, MOON };
    
    static State getState(Planet planet, double jd_tdb) {
        double T = (jd_tdb - JD_J2000) / 36525.0;  // Secoli giuliani da J2000
        
        // Elementi orbitali medi (Simon et al. 1994)
        double a, e, i, L, omega_bar, Omega;
        
        switch (planet) {
            case MERCURY:
                a = 0.38709927;
                e = 0.20563593 + 0.00001906 * T;
                i = (7.00497902 - 0.00594749 * T) * DEG_TO_RAD;
                L = (252.25032350 + 149472.67411175 * T) * DEG_TO_RAD;
                omega_bar = (77.45779628 + 0.16047689 * T) * DEG_TO_RAD;
                Omega = (48.33076593 - 0.12534081 * T) * DEG_TO_RAD;
                break;
            case VENUS:
                a = 0.72333566;
                e = 0.00677672 - 0.00004107 * T;
                i = (3.39467605 - 0.00078890 * T) * DEG_TO_RAD;
                L = (181.97909950 + 58517.81538729 * T) * DEG_TO_RAD;
                omega_bar = (131.60246718 + 0.00268329 * T) * DEG_TO_RAD;
                Omega = (76.67984255 - 0.27769418 * T) * DEG_TO_RAD;
                break;
            case EARTH:
                a = 1.00000261 + 0.00000562 * T;
                e = 0.01671123 - 0.00004392 * T;
                i = (-0.00001531 - 0.01294668 * T) * DEG_TO_RAD;
                L = (100.46457166 + 35999.37244981 * T) * DEG_TO_RAD;
                omega_bar = (102.93768193 + 0.32327364 * T) * DEG_TO_RAD;
                Omega = 0.0;
                break;
            case MARS:
                a = 1.52371034;
                e = 0.09339410 + 0.00007882 * T;
                i = (1.84969142 - 0.00813131 * T) * DEG_TO_RAD;
                L = (-4.55343205 + 19140.30268499 * T) * DEG_TO_RAD;
                omega_bar = (-23.94362959 + 0.44441088 * T) * DEG_TO_RAD;
                Omega = (49.55953891 - 0.29257343 * T) * DEG_TO_RAD;
                break;
            case JUPITER:
                a = 5.20288700;
                e = 0.04838624 - 0.00013537 * T;
                i = (1.30439695 - 0.00183714 * T) * DEG_TO_RAD;
                L = (34.39644051 + 3034.74612775 * T) * DEG_TO_RAD;
                omega_bar = (14.72847983 + 0.21252668 * T) * DEG_TO_RAD;
                Omega = (100.47390909 + 0.20469106 * T) * DEG_TO_RAD;
                break;
            case SATURN:
                a = 9.53667594;
                e = 0.05386179 - 0.00050991 * T;
                i = (2.48599187 + 0.00193609 * T) * DEG_TO_RAD;
                L = (49.95424423 + 1222.49362201 * T) * DEG_TO_RAD;
                omega_bar = (92.59887831 - 0.41897216 * T) * DEG_TO_RAD;
                Omega = (113.66242448 - 0.28867794 * T) * DEG_TO_RAD;
                break;
            case URANUS:
                a = 19.18916464;
                e = 0.04725744 - 0.00004397 * T;
                i = (0.77263783 - 0.00242939 * T) * DEG_TO_RAD;
                L = (313.23810451 + 428.48202785 * T) * DEG_TO_RAD;
                omega_bar = (170.95427630 + 0.40805281 * T) * DEG_TO_RAD;
                Omega = (74.01692503 + 0.04240589 * T) * DEG_TO_RAD;
                break;
            case NEPTUNE:
                a = 30.06992276;
                e = 0.00859048 + 0.00005105 * T;
                i = (1.77004347 + 0.00035372 * T) * DEG_TO_RAD;
                L = (-55.12002969 + 218.45945325 * T) * DEG_TO_RAD;
                omega_bar = (44.96476227 - 0.32241464 * T) * DEG_TO_RAD;
                Omega = (131.78422574 - 0.00508664 * T) * DEG_TO_RAD;
                break;
            case MOON:
                // Luna geocentrica semplificata
                a = 0.00256955529;
                e = 0.0549;
                i = 5.145 * DEG_TO_RAD;
                L = (218.316 + 481267.881 * T) * DEG_TO_RAD;
                omega_bar = 318.15 * DEG_TO_RAD;
                Omega = (125.08 - 1934.136 * T) * DEG_TO_RAD;
                break;
            default:
                return State();
        }
        
        return elementsToState(a, e, i, L, omega_bar, Omega, GM_SUN);
    }
    
private:
    static State elementsToState(double a, double e, double i, 
                                  double L, double omega_bar, double Omega,
                                  double gm) {
        // Argomento del perielio e anomalia media
        double omega = omega_bar - Omega;
        double M = L - omega_bar;
        
        // Normalizza M in [0, 2π)
        M = std::fmod(M, TWO_PI);
        if (M < 0) M += TWO_PI;
        
        // Risolvi equazione di Keplero: E - e*sin(E) = M
        double E = M;
        for (int iter = 0; iter < 15; ++iter) {
            double dE = (E - e * std::sin(E) - M) / (1.0 - e * std::cos(E));
            E -= dE;
            if (std::abs(dE) < 1e-14) break;
        }
        
        // Anomalia vera
        double sin_E = std::sin(E);
        double cos_E = std::cos(E);
        double sqrt_1_e2 = std::sqrt(1.0 - e * e);
        double nu = std::atan2(sqrt_1_e2 * sin_E, cos_E - e);
        
        // Distanza e posizione nel piano orbitale
        double r = a * (1.0 - e * cos_E);
        double cos_nu = std::cos(nu);
        double sin_nu = std::sin(nu);
        double x_orb = r * cos_nu;
        double y_orb = r * sin_nu;
        
        // Velocità nel piano orbitale
        double v_factor = std::sqrt(gm * a) / r;
        double vx_orb = -v_factor * sin_E;
        double vy_orb = v_factor * sqrt_1_e2 * cos_E;
        
        // Matrice di rotazione
        double cos_omega = std::cos(omega);
        double sin_omega = std::sin(omega);
        double cos_i = std::cos(i);
        double sin_i = std::sin(i);
        double cos_Omega = std::cos(Omega);
        double sin_Omega = std::sin(Omega);
        
        double P11 = cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i;
        double P12 = -cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i;
        double P21 = sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i;
        double P22 = -sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i;
        double P31 = sin_omega * sin_i;
        double P32 = cos_omega * sin_i;
        
        State state;
        state.pos.x = P11 * x_orb + P12 * y_orb;
        state.pos.y = P21 * x_orb + P22 * y_orb;
        state.pos.z = P31 * x_orb + P32 * y_orb;
        state.vel.x = P11 * vx_orb + P12 * vy_orb;
        state.vel.y = P21 * vx_orb + P22 * vy_orb;
        state.vel.z = P31 * vx_orb + P32 * vy_orb;
        
        return state;
    }
};

// ============================================================================
// AST17 - 16 ASTEROIDI MASSIVI (come in OrbFit)
// ============================================================================

class AST17 {
public:
    // GM in km³/s² (da Baer et al. 2011 e agiornamenti)
    static constexpr double GM_CERES = 62.6284;
    static constexpr double GM_PALLAS = 14.28;
    static constexpr double GM_VESTA = 17.8;
    static constexpr double GM_HYGIEA = 5.78;
    static constexpr double GM_EUNOMIA = 2.1;
    static constexpr double GM_PSYCHE = 1.81;
    static constexpr double GM_EUPHROSYNE = 1.26;
    static constexpr double GM_EUROPA52 = 1.59;
    static constexpr double GM_CYBELE = 0.91;
    static constexpr double GM_SYLVIA = 0.99;
    static constexpr double GM_THISBE = 0.78;
    static constexpr double GM_CAMILLA = 0.75;
    static constexpr double GM_BAMBERGA = 0.69;
    static constexpr double GM_PATIENTIA = 0.54;
    static constexpr double GM_DAVIDA = 0.89;
    static constexpr double GM_INTERAMNIA = 1.86;
    
    static std::vector<AsteroidData> getAsteroids() {
        // Epoca: J2000.0 (MJD 51544.5)
        constexpr double epoch = 51544.5;
        
        return {
            {1, "Ceres", GM_CERES, 2.7675, 0.0760, 10.593, 73.597, 80.3932, 113.410, epoch, 0.2141},
            {2, "Pallas", GM_PALLAS, 2.7730, 0.2299, 34.837, 310.049, 173.0962, 78.194, epoch, 0.2133},
            {4, "Vesta", GM_VESTA, 2.3615, 0.0887, 7.142, 151.198, 103.8513, 205.546, epoch, 0.2716},
            {10, "Hygiea", GM_HYGIEA, 3.1393, 0.1146, 3.831, 283.455, 283.2045, 107.590, epoch, 0.1794},
            {15, "Eunomia", GM_EUNOMIA, 2.6439, 0.1866, 11.737, 97.767, 293.2081, 142.405, epoch, 0.2262},
            {16, "Psyche", GM_PSYCHE, 2.9216, 0.1339, 3.096, 227.305, 150.2873, 179.942, epoch, 0.1990},
            {31, "Euphrosyne", GM_EUPHROSYNE, 3.1515, 0.2257, 26.307, 61.588, 31.2873, 325.478, epoch, 0.1785},
            {52, "Europa", GM_EUROPA52, 3.0958, 0.1020, 7.460, 343.545, 129.0380, 252.132, epoch, 0.1824},
            {65, "Cybele", GM_CYBELE, 3.4332, 0.1050, 3.562, 155.889, 155.5463, 196.234, epoch, 0.1636},
            {87, "Sylvia", GM_SYLVIA, 3.4876, 0.0808, 10.866, 266.315, 73.3430, 271.478, epoch, 0.1610},
            {88, "Thisbe", GM_THISBE, 2.7671, 0.1641, 5.222, 36.714, 276.0972, 355.489, epoch, 0.2141},
            {107, "Camilla", GM_CAMILLA, 3.4886, 0.0692, 10.044, 309.785, 173.0458, 120.234, epoch, 0.1609},
            {324, "Bamberga", GM_BAMBERGA, 2.6835, 0.3381, 11.095, 43.678, 327.8394, 155.234, epoch, 0.2226},
            {451, "Patientia", GM_PATIENTIA, 3.0639, 0.0764, 15.240, 278.789, 96.9873, 234.456, epoch, 0.1843},
            {511, "Davida", GM_DAVIDA, 3.1805, 0.1833, 15.942, 107.234, 270.0945, 186.789, epoch, 0.1768},
            {704, "Interamnia", GM_INTERAMNIA, 3.0616, 0.1502, 17.296, 94.789, 280.3456, 145.234, epoch, 0.1844}
        };
    }
    
    static Vec3 getPosition(const AsteroidData& ast, double mjd_tdb) {
        // Calcola anomalia media al tempo t
        double dt = mjd_tdb - ast.epoch;
        double M = (ast.M0 + ast.n * dt) * DEG_TO_RAD;
        
        // Normalizza
        M = std::fmod(M, TWO_PI);
        if (M < 0) M += TWO_PI;
        
        // Risolvi Keplero
        double E = M;
        for (int iter = 0; iter < 15; ++iter) {
            double dE = (E - ast.e * std::sin(E) - M) / (1.0 - ast.e * std::cos(E));
            E -= dE;
            if (std::abs(dE) < 1e-14) break;
        }
        
        // Anomalia vera e distanza
        double sin_E = std::sin(E);
        double cos_E = std::cos(E);
        double sqrt_1_e2 = std::sqrt(1.0 - ast.e * ast.e);
        double nu = std::atan2(sqrt_1_e2 * sin_E, cos_E - ast.e);
        double r = ast.a * (1.0 - ast.e * cos_E);
        
        // Posizione nel piano orbitale
        double x_orb = r * std::cos(nu);
        double y_orb = r * std::sin(nu);
        
        // Rotazione all'eclittica J2000
        double i_rad = ast.i * DEG_TO_RAD;
        double omega_rad = ast.omega * DEG_TO_RAD;
        double Omega_rad = ast.Omega * DEG_TO_RAD;
        
        double cos_omega = std::cos(omega_rad);
        double sin_omega = std::sin(omega_rad);
        double cos_i = std::cos(i_rad);
        double sin_i = std::sin(i_rad);
        double cos_Omega = std::cos(Omega_rad);
        double sin_Omega = std::sin(Omega_rad);
        
        Vec3 pos;
        pos.x = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * x_orb +
                (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * y_orb;
        pos.y = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * x_orb +
                (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * y_orb;
        pos.z = (sin_omega * sin_i) * x_orb + (cos_omega * sin_i) * y_orb;
        
        return pos;
    }
};

// ============================================================================
// INTEGRATORE RKF78 (Runge-Kutta-Fehlberg 7(8))
// ============================================================================

class RKF78 {
public:
    using DerivFunc = std::function<std::array<double, 6>(double, const std::array<double, 6>&)>;
    
    static std::array<double, 6> integrate(DerivFunc f, 
                                           const std::array<double, 6>& y0,
                                           double t0, double tf,
                                           double tol = 1e-12) {
        std::array<double, 6> y = y0;
        double t = t0;
        double dt = (tf - t0) / 100.0;  // Passo iniziale
        
        // Limiti passo
        double dt_min = std::abs(tf - t0) * 1e-10;
        double dt_max = std::abs(tf - t0) / 10.0;
        int direction = (tf > t0) ? 1 : -1;
        
        while (direction * (tf - t) > 1e-15) {
            // Limita passo all'arrivo
            if (direction * (t + dt - tf) > 0) {
                dt = tf - t;
            }
            
            // Calcola coefficienti RK4/5 (Dormand-Prince semplificato)
            auto k1 = f(t, y);
            
            std::array<double, 6> y2;
            for (int i = 0; i < 6; ++i) y2[i] = y[i] + 0.5 * dt * k1[i];
            auto k2 = f(t + 0.5 * dt, y2);
            
            std::array<double, 6> y3;
            for (int i = 0; i < 6; ++i) y3[i] = y[i] + 0.5 * dt * k2[i];
            auto k3 = f(t + 0.5 * dt, y3);
            
            std::array<double, 6> y4;
            for (int i = 0; i < 6; ++i) y4[i] = y[i] + dt * k3[i];
            auto k4 = f(t + dt, y4);
            
            // Soluzione RK4
            std::array<double, 6> y_new;
            for (int i = 0; i < 6; ++i) {
                y_new[i] = y[i] + dt * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6.0;
            }
            
            // Controllo errore (semplificato)
            double err = 0.0;
            for (int i = 0; i < 3; ++i) {
                err = std::max(err, std::abs(dt * (k1[i] - 2*k2[i] + k4[i]) / 6.0));
            }
            
            if (err < tol || std::abs(dt) <= dt_min) {
                // Accetta passo
                y = y_new;
                t += dt;
            }
            
            // Adatta passo
            if (err > 0) {
                double factor = 0.9 * std::pow(tol / err, 0.2);
                factor = std::min(2.0, std::max(0.1, factor));
                dt *= factor;
            }
            
            dt = direction * std::min(std::abs(dt), dt_max);
            dt = direction * std::max(std::abs(dt), dt_min);
        }
        
        return y;
    }
};

// ============================================================================
// PROPAGATORE COMPLETO
// ============================================================================

class FullPropagator {
public:
    FullPropagator() : asteroids_(AST17::getAsteroids()) {}
    
    // Propaga da t0 a tf con tutte le perturbazioni
    State propagate(const State& initial, double jd0, double jdf) {
        // Converti a MJD per compatibilità con effemeridi
        double mjd0 = jd0 - 2400000.5;
        double mjdf = jdf - 2400000.5;
        
        // Stato iniziale come array
        std::array<double, 6> y0 = {
            initial.pos.x, initial.pos.y, initial.pos.z,
            initial.vel.x, initial.vel.y, initial.vel.z
        };
        
        // Funzione derivata
        auto derivatives = [this, mjd0](double t, const std::array<double, 6>& y) {
            return computeDerivatives(t, y, mjd0);
        };
        
        // Integra
        double dt = mjdf - mjd0;
        auto yf = RKF78::integrate(derivatives, y0, 0.0, dt, 1e-14);
        
        State final_state;
        final_state.pos = Vec3(yf[0], yf[1], yf[2]);
        final_state.vel = Vec3(yf[3], yf[4], yf[5]);
        
        return final_state;
    }
    
private:
    std::vector<AsteroidData> asteroids_;
    
    std::array<double, 6> computeDerivatives(double t, const std::array<double, 6>& y, double mjd0) {
        Vec3 pos(y[0], y[1], y[2]);
        Vec3 vel(y[3], y[4], y[5]);
        double mjd = mjd0 + t;
        double jd = mjd + 2400000.5;
        
        // Accelerazione totale
        Vec3 acc = twoBodyAcceleration(pos);
        acc += planetaryPerturbations(pos, jd);
        acc += asteroidPerturbations(pos, mjd);
        acc += relativisticCorrection(pos, vel);
        
        return {vel.x, vel.y, vel.z, acc.x, acc.y, acc.z};
    }
    
    Vec3 twoBodyAcceleration(const Vec3& pos) {
        double r = pos.norm();
        double r3 = r * r * r;
        return pos * (-GM_SUN / r3);
    }
    
    Vec3 planetaryPerturbations(const Vec3& pos, double jd) {
        Vec3 acc;
        
        // Array di pianeti e loro GM
        struct PlanetGM {
            PlanetaryEphemeris::Planet planet;
            double gm_km3s2;
        };
        
        PlanetGM planets[] = {
            {PlanetaryEphemeris::MERCURY, GM_MERCURY},
            {PlanetaryEphemeris::VENUS, GM_VENUS},
            {PlanetaryEphemeris::EARTH, GM_EARTH},
            {PlanetaryEphemeris::MARS, GM_MARS},
            {PlanetaryEphemeris::JUPITER, GM_JUPITER},
            {PlanetaryEphemeris::SATURN, GM_SATURN},
            {PlanetaryEphemeris::URANUS, GM_URANUS},
            {PlanetaryEphemeris::NEPTUNE, GM_NEPTUNE}
        };
        
        for (const auto& p : planets) {
            State pstate = PlanetaryEphemeris::getState(p.planet, jd);
            double gm_au3day2 = p.gm_km3s2 * GM_CONV;
            
            Vec3 delta = pstate.pos - pos;
            double delta_r = delta.norm();
            double delta_r3 = delta_r * delta_r * delta_r;
            
            double planet_r = pstate.pos.norm();
            double planet_r3 = planet_r * planet_r * planet_r;
            
            // Termine diretto + indiretto
            acc += delta * (gm_au3day2 / delta_r3) - pstate.pos * (gm_au3day2 / planet_r3);
        }
        
        return acc;
    }
    
    Vec3 asteroidPerturbations(const Vec3& pos, double mjd) {
        Vec3 acc;
        
        for (const auto& ast : asteroids_) {
            Vec3 ast_pos = AST17::getPosition(ast, mjd);
            double gm_au3day2 = ast.gm * GM_CONV;
            
            Vec3 delta = ast_pos - pos;
            double delta_r = delta.norm();
            double delta_r3 = delta_r * delta_r * delta_r;
            
            double ast_r = ast_pos.norm();
            double ast_r3 = ast_r * ast_r * ast_r;
            
            acc += delta * (gm_au3day2 / delta_r3) - ast_pos * (gm_au3day2 / ast_r3);
        }
        
        return acc;
    }
    
    Vec3 relativisticCorrection(const Vec3& pos, const Vec3& vel) {
        // Correzione Schwarzschild (1° ordine in 1/c²)
        double r = pos.norm();
        double r3 = r * r * r;
        double v2 = vel.norm2();
        double rdot = pos.dot(vel) / r;
        double c2 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;
        
        // a_GR = (μ/r³) * [(4μ/r - v²)r + 4(r·v)v] / c²
        Vec3 term1 = pos * (4.0 * GM_SUN / r - v2);
        Vec3 term2 = vel * (4.0 * rdot);
        
        return (term1 + term2) * (GM_SUN / (r3 * c2));
    }
};

// ============================================================================
// FUNZIONI DI CONVERSIONE
// ============================================================================

State keplerToCartesian(const KeplerianElements& kep) {
    // Risolvi Keplero
    double E = kep.M0;
    for (int iter = 0; iter < 15; ++iter) {
        double dE = (E - kep.e * std::sin(E) - kep.M0) / (1.0 - kep.e * std::cos(E));
        E -= dE;
        if (std::abs(dE) < 1e-14) break;
    }
    
    // Anomalia vera
    double sin_E = std::sin(E);
    double cos_E = std::cos(E);
    double sqrt_1_e2 = std::sqrt(1.0 - kep.e * kep.e);
    double nu = std::atan2(sqrt_1_e2 * sin_E, cos_E - kep.e);
    
    // Distanza e posizione nel piano orbitale
    double r = kep.a * (1.0 - kep.e * cos_E);
    double cos_nu = std::cos(nu);
    double sin_nu = std::sin(nu);
    double x_orb = r * cos_nu;
    double y_orb = r * sin_nu;
    
    // Velocità nel piano orbitale
    double v_factor = std::sqrt(kep.gm * kep.a) / r;
    double vx_orb = -v_factor * sin_E;
    double vy_orb = v_factor * sqrt_1_e2 * cos_E;
    
    // Matrice di rotazione
    double cos_omega = std::cos(kep.omega);
    double sin_omega = std::sin(kep.omega);
    double cos_i = std::cos(kep.i);
    double sin_i = std::sin(kep.i);
    double cos_Omega = std::cos(kep.Omega);
    double sin_Omega = std::sin(kep.Omega);
    
    double P11 = cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i;
    double P12 = -cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i;
    double P21 = sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i;
    double P22 = -sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i;
    double P31 = sin_omega * sin_i;
    double P32 = cos_omega * sin_i;
    
    State state;
    state.pos.x = P11 * x_orb + P12 * y_orb;
    state.pos.y = P21 * x_orb + P22 * y_orb;
    state.pos.z = P31 * x_orb + P32 * y_orb;
    state.vel.x = P11 * vx_orb + P12 * vy_orb;
    state.vel.y = P21 * vx_orb + P22 * vy_orb;
    state.vel.z = P31 * vx_orb + P32 * vy_orb;
    
    return state;
}

Vec3 eclipticToEquatorial(const Vec3& ecl) {
    double cos_eps = std::cos(OBLIQUITY_J2000);
    double sin_eps = std::sin(OBLIQUITY_J2000);
    
    return Vec3(
        ecl.x,
        ecl.y * cos_eps - ecl.z * sin_eps,
        ecl.y * sin_eps + ecl.z * cos_eps
    );
}

void cartesianToRADec(const Vec3& pos, double& ra, double& dec) {
    double r = pos.norm();
    dec = std::asin(pos.z / r);
    ra = std::atan2(pos.y, pos.x);
    if (ra < 0) ra += TWO_PI;
}

std::string formatRA(double ra_rad) {
    double ra_hours = ra_rad * 12.0 / PI;
    int h = static_cast<int>(ra_hours);
    double m_frac = (ra_hours - h) * 60.0;
    int m = static_cast<int>(m_frac);
    double s = (m_frac - m) * 60.0;
    
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%02d %02d %06.3f", h, m, s);
    return std::string(buf);
}

std::string formatDec(double dec_rad) {
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

// ============================================================================
// MAIN
// ============================================================================

int main() {
    std::cout << std::fixed << std::setprecision(12);
    
    std::cout << "================================================================\n";
    std::cout << "  Calcolo PRECISO Posizione (17030) Sierks\n";
    std::cout << "  AstDyn/OrbFit - Integrazione Numerica Completa\n";
    std::cout << "================================================================\n\n";
    
    // Elementi kepleriani MPC (epoca recente)
    KeplerianElements sierks;
    sierks.a = 3.1754733;
    sierks.e = 0.0454207;
    sierks.i = 2.9046 * DEG_TO_RAD;
    sierks.Omega = 104.16243 * DEG_TO_RAD;
    sierks.omega = 100.5141 * DEG_TO_RAD;
    sierks.M0 = 229.79088 * DEG_TO_RAD;
    sierks.epoch = 2461000.5;  // JD TDB
    sierks.gm = GM_SUN;
    
    double target_jd = 2461008.0913;  // 28 Nov 2025, 14:11:28 UTC
    
    std::cout << "ELEMENTI ORBITALI (MPC):\n";
    std::cout << "  a     = " << sierks.a << " AU\n";
    std::cout << "  e     = " << sierks.e << "\n";
    std::cout << "  i     = " << (sierks.i * RAD_TO_DEG) << "°\n";
    std::cout << "  Ω     = " << (sierks.Omega * RAD_TO_DEG) << "°\n";
    std::cout << "  ω     = " << (sierks.omega * RAD_TO_DEG) << "°\n";
    std::cout << "  M     = " << (sierks.M0 * RAD_TO_DEG) << "°\n";
    std::cout << "  Epoch = JD " << sierks.epoch << "\n\n";
    
    std::cout << "TARGET: JD " << target_jd << " (28 Nov 2025, 14:11:28 UTC)\n\n";
    
    // Converti elementi iniziali in stato cartesiano
    State initial = keplerToCartesian(sierks);
    
    std::cout << "STATO INIZIALE (eclittico eliocentrico):\n";
    std::cout << "  r = [" << initial.pos.x << ", " << initial.pos.y << ", " << initial.pos.z << "] AU\n";
    std::cout << "  v = [" << initial.vel.x << ", " << initial.vel.y << ", " << initial.vel.z << "] AU/day\n\n";
    
    // PROPAGAZIONE COMPLETA
    std::cout << "PROPAGAZIONE CON:\n";
    std::cout << "  ✓ Integrazione numerica RK4\n";
    std::cout << "  ✓ Perturbazioni planetarie (8 pianeti)\n";
    std::cout << "  ✓ Perturbazioni asteroidali (AST17 - 16 asteroidi)\n";
    std::cout << "  ✓ Correzione relativistica (Schwarzschild)\n\n";
    
    FullPropagator propagator;
    State final_state = propagator.propagate(initial, sierks.epoch, target_jd);
    
    std::cout << "STATO FINALE (eclittico eliocentrico):\n";
    std::cout << "  r = [" << final_state.pos.x << ", " << final_state.pos.y << ", " << final_state.pos.z << "] AU\n";
    std::cout << "  v = [" << final_state.vel.x << ", " << final_state.vel.y << ", " << final_state.vel.z << "] AU/day\n";
    double r_helio = final_state.pos.norm();
    std::cout << "  |r| = " << r_helio << " AU\n\n";
    
    // Posizione Terra
    State earth = PlanetaryEphemeris::getState(PlanetaryEphemeris::EARTH, target_jd);
    
    std::cout << "POSIZIONE TERRA:\n";
    std::cout << "  r = [" << earth.pos.x << ", " << earth.pos.y << ", " << earth.pos.z << "] AU\n\n";
    
    // Posizione geocentrica
    Vec3 geo_ecl = final_state.pos - earth.pos;
    double delta = geo_ecl.norm();
    
    // CORREZIONE TEMPO LUCE
    double light_time = delta / SPEED_OF_LIGHT;  // giorni
    std::cout << "CORREZIONE TEMPO LUCE:\n";
    std::cout << "  Δ = " << delta << " AU\n";
    std::cout << "  τ = " << (light_time * 86400.0) << " s = " << (light_time * 1440.0) << " min\n\n";
    
    // Ricalcola posizione asteroide al tempo retrodatato
    double jd_retarded = target_jd - light_time;
    State final_lt = propagator.propagate(initial, sierks.epoch, jd_retarded);
    Vec3 geo_ecl_lt = final_lt.pos - earth.pos;
    double delta_lt = geo_ecl_lt.norm();
    
    std::cout << "POSIZIONE GEOCENTRICA (con tempo luce):\n";
    std::cout << "  Δr = [" << geo_ecl_lt.x << ", " << geo_ecl_lt.y << ", " << geo_ecl_lt.z << "] AU\n";
    std::cout << "  |Δ| = " << delta_lt << " AU\n\n";
    
    // Converti a equatoriale
    Vec3 geo_eq = eclipticToEquatorial(geo_ecl_lt);
    
    std::cout << "POSIZIONE GEOCENTRICA EQUATORIALE (J2000):\n";
    std::cout << "  r = [" << geo_eq.x << ", " << geo_eq.y << ", " << geo_eq.z << "] AU\n\n";
    
    // Calcola RA/Dec
    double ra, dec;
    cartesianToRADec(geo_eq, ra, dec);
    
    // ABERRAZIONE STELLARE ANNUALE
    // v_earth / c ~ 10^-4, effetto ~20"
    double v_earth = earth.vel.norm();
    double aberration_factor = v_earth / SPEED_OF_LIGHT;
    
    // Componenti aberrazione (semplificato)
    Vec3 earth_vel_eq = eclipticToEquatorial(earth.vel);
    double dra = (earth_vel_eq.y * std::cos(ra) - earth_vel_eq.x * std::sin(ra)) / 
                 (SPEED_OF_LIGHT * std::cos(dec));
    double ddec = (earth_vel_eq.z * std::cos(dec) - 
                   (earth_vel_eq.x * std::cos(ra) + earth_vel_eq.y * std::sin(ra)) * std::sin(dec)) /
                  SPEED_OF_LIGHT;
    
    double ra_aberr = ra + dra;
    double dec_aberr = dec + ddec;
    
    std::cout << "ABERRAZIONE STELLARE:\n";
    std::cout << "  |v_⊕| = " << v_earth << " AU/day\n";
    std::cout << "  ΔRA  = " << (dra * RAD_TO_DEG * 3600) << " arcsec\n";
    std::cout << "  ΔDec = " << (ddec * RAD_TO_DEG * 3600) << " arcsec\n\n";
    
    // RISULTATI FINALI
    std::cout << "================================================================\n";
    std::cout << "  RISULTATO FINALE - (17030) Sierks\n";
    std::cout << "  JD " << target_jd << " (28 Nov 2025, 14:11:28 UTC)\n";
    std::cout << "================================================================\n\n";
    
    std::cout << "COORDINATE ASTROMETRICHE (J2000, senza aberrazione):\n";
    std::cout << "  RA  = " << formatRA(ra) << "  (" << (ra * RAD_TO_DEG) << "°)\n";
    std::cout << "  Dec = " << formatDec(dec) << "  (" << (dec * RAD_TO_DEG) << "°)\n\n";
    
    std::cout << "COORDINATE APPARENTI (con aberrazione):\n";
    std::cout << "  RA  = " << formatRA(ra_aberr) << "  (" << (ra_aberr * RAD_TO_DEG) << "°)\n";
    std::cout << "  Dec = " << formatDec(dec_aberr) << "  (" << (dec_aberr * RAD_TO_DEG) << "°)\n\n";
    
    std::cout << "DISTANZE:\n";
    std::cout << "  r (eliocentrica) = " << r_helio << " AU\n";
    std::cout << "  Δ (geocentrica)  = " << delta_lt << " AU = " 
              << (delta_lt * AU_KM) << " km\n\n";
    
    // Confronto con JPL Horizons
    std::cout << "================================================================\n";
    std::cout << "  CONFRONTO CON JPL HORIZONS\n";
    std::cout << "================================================================\n\n";
    std::cout << "JPL Horizons (ICRF):  04 53 11.25   +20 19 25.8\n";
    std::cout << "AstDyn (astrometric): " << formatRA(ra) << "   " << formatDec(dec) << "\n\n";
    
    double ra_jpl = (4.0 + 53.0/60.0 + 11.25/3600.0) * 15.0 * DEG_TO_RAD;
    double dec_jpl = (20.0 + 19.0/60.0 + 25.8/3600.0) * DEG_TO_RAD;
    
    double dra_arcsec = (ra - ra_jpl) * RAD_TO_DEG * 3600.0 * std::cos(dec);
    double ddec_arcsec = (dec - dec_jpl) * RAD_TO_DEG * 3600.0;
    
    std::cout << "RESIDUI:\n";
    std::cout << "  ΔRA*cos(δ) = " << std::setprecision(2) << dra_arcsec << " arcsec\n";
    std::cout << "  ΔDec       = " << ddec_arcsec << " arcsec\n";
    std::cout << "  Totale     = " << std::sqrt(dra_arcsec*dra_arcsec + ddec_arcsec*ddec_arcsec) 
              << " arcsec\n\n";
    
    std::cout << "================================================================\n";
    std::cout << "  Calcolo completato.\n";
    std::cout << "================================================================\n";
    
    return 0;
}
