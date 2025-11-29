/**
 * @file test_rkf78.cpp
 * @brief Integratore RKF78 (Runge-Kutta-Fehlberg 7(8)) COMPLETO
 *        Con coefficienti originali di Fehlberg (1968)
 * 
 * RKF78 è un metodo embedded di ordine 7 con stima dell'errore di ordine 8.
 * È lo stesso integratore usato in OrbFit e nei software di meccanica celeste.
 * 
 * Coefficienti da: Fehlberg, E. (1968) "Classical Fifth-, Sixth-, Seventh-, 
 *                  and Eighth-Order Runge-Kutta Formulas with Stepsize Control"
 *                  NASA TR R-287
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
#include <algorithm>

// ============================================================================
// COSTANTI
// ============================================================================

namespace constants {
    constexpr double PI = 3.14159265358979323846;
    constexpr double TWO_PI = 2.0 * PI;
    constexpr double DEG_TO_RAD = PI / 180.0;
    constexpr double RAD_TO_DEG = 180.0 / PI;
    
    constexpr double K_GAUSS = 0.01720209895;
    constexpr double GM_SUN = K_GAUSS * K_GAUSS;
    
    constexpr double AU_KM = 149597870.7;
    constexpr double DAY_SEC = 86400.0;
    constexpr double GM_CONV = DAY_SEC * DAY_SEC / (AU_KM * AU_KM * AU_KM);
    
    constexpr double GM_MERCURY = 22031.868551;
    constexpr double GM_VENUS = 324858.592000;
    constexpr double GM_EARTH = 398600.435507;
    constexpr double GM_MARS = 42828.375816;
    constexpr double GM_JUPITER = 126712764.100000;
    constexpr double GM_SATURN = 37940584.841800;
    constexpr double GM_URANUS = 5794556.400000;
    constexpr double GM_NEPTUNE = 6836527.100580;
    
    constexpr double SPEED_OF_LIGHT = 299792.458 / AU_KM * DAY_SEC;
    constexpr double OBLIQUITY_J2000 = 23.439291111 * DEG_TO_RAD;
    constexpr double JD_J2000 = 2451545.0;
}

using namespace constants;

// ============================================================================
// INTEGRATORE RKF78 COMPLETO (Fehlberg 1968)
// ============================================================================

/**
 * @class RKF78
 * @brief Integratore Runge-Kutta-Fehlberg 7(8) con controllo automatico del passo
 * 
 * Caratteristiche:
 * - 13 valutazioni della funzione per passo
 * - Ordine 7 per la soluzione, ordine 8 per la stima dell'errore
 * - Controllo adattivo del passo basato sull'errore locale
 * - Altamente efficiente per problemi di meccanica celeste
 */
template<size_t N>
class RKF78 {
public:
    using State = std::array<double, N>;
    using DerivFunc = std::function<State(double, const State&)>;
    
    struct IntegrationStats {
        int steps_accepted = 0;
        int steps_rejected = 0;
        int function_evals = 0;
        double min_step = 1e30;
        double max_step = 0.0;
    };
    
    /**
     * @brief Costruttore
     * @param rtol Tolleranza relativa (default 1e-12)
     * @param atol Tolleranza assoluta (default 1e-15)
     */
    RKF78(double rtol = 1e-12, double atol = 1e-15) 
        : rtol_(rtol), atol_(atol) {
        initializeCoefficients();
    }
    
    /**
     * @brief Integra da t0 a tf
     * @param f Funzione derivata f(t, y)
     * @param y0 Stato iniziale
     * @param t0 Tempo iniziale
     * @param tf Tempo finale
     * @param h0 Passo iniziale (opzionale)
     * @return Stato finale
     */
    State integrate(DerivFunc f, const State& y0, double t0, double tf, double h0 = 0.0) {
        stats_ = IntegrationStats();
        
        State y = y0;
        double t = t0;
        
        // Direzione di integrazione
        int direction = (tf > t0) ? 1 : -1;
        double dt = tf - t0;
        
        // Passo iniziale
        double h = (h0 != 0.0) ? std::abs(h0) : estimateInitialStep(f, y0, t0, dt);
        h = direction * h;
        
        // Limiti del passo
        double h_min = std::abs(dt) * 1e-12;
        double h_max = std::abs(dt) * 0.1;
        
        // Loop principale
        while (direction * (tf - t) > 1e-15 * std::abs(tf)) {
            // Limita il passo per arrivare esattamente a tf
            if (direction * (t + h - tf) > 0) {
                h = tf - t;
            }
            
            // Tenta un passo
            State y_new;
            double err;
            bool accepted = tryStep(f, t, y, h, y_new, err);
            
            if (accepted) {
                t += h;
                y = y_new;
                stats_.steps_accepted++;
            } else {
                stats_.steps_rejected++;
            }
            
            // Calcola nuovo passo
            double h_new = computeNewStep(h, err, accepted);
            
            // Applica limiti
            h_new = direction * std::clamp(std::abs(h_new), h_min, h_max);
            
            // Aggiorna statistiche
            stats_.min_step = std::min(stats_.min_step, std::abs(h));
            stats_.max_step = std::max(stats_.max_step, std::abs(h));
            
            h = h_new;
        }
        
        return y;
    }
    
    const IntegrationStats& getStats() const { return stats_; }
    
private:
    double rtol_, atol_;
    IntegrationStats stats_;
    
    // Coefficienti RKF78 (Fehlberg 1968)
    // 13 stadi
    static constexpr int STAGES = 13;
    
    // Nodi c_i
    double c_[STAGES];
    
    // Matrice A (lower triangular)
    double a_[STAGES][STAGES];
    
    // Pesi per soluzione di ordine 7
    double b7_[STAGES];
    
    // Pesi per soluzione di ordine 8 (per stima errore)
    double b8_[STAGES];
    
    void initializeCoefficients() {
        // Inizializza a zero
        for (int i = 0; i < STAGES; ++i) {
            c_[i] = 0.0;
            b7_[i] = 0.0;
            b8_[i] = 0.0;
            for (int j = 0; j < STAGES; ++j) {
                a_[i][j] = 0.0;
            }
        }
        
        // Nodi c_i (Fehlberg 1968, Table 5)
        c_[0]  = 0.0;
        c_[1]  = 2.0 / 27.0;
        c_[2]  = 1.0 / 9.0;
        c_[3]  = 1.0 / 6.0;
        c_[4]  = 5.0 / 12.0;
        c_[5]  = 1.0 / 2.0;
        c_[6]  = 5.0 / 6.0;
        c_[7]  = 1.0 / 6.0;
        c_[8]  = 2.0 / 3.0;
        c_[9]  = 1.0 / 3.0;
        c_[10] = 1.0;
        c_[11] = 0.0;
        c_[12] = 1.0;
        
        // Matrice A (coefficienti di Fehlberg)
        a_[1][0] = 2.0 / 27.0;
        
        a_[2][0] = 1.0 / 36.0;
        a_[2][1] = 1.0 / 12.0;
        
        a_[3][0] = 1.0 / 24.0;
        a_[3][2] = 1.0 / 8.0;
        
        a_[4][0] = 5.0 / 12.0;
        a_[4][2] = -25.0 / 16.0;
        a_[4][3] = 25.0 / 16.0;
        
        a_[5][0] = 1.0 / 20.0;
        a_[5][3] = 1.0 / 4.0;
        a_[5][4] = 1.0 / 5.0;
        
        a_[6][0] = -25.0 / 108.0;
        a_[6][3] = 125.0 / 108.0;
        a_[6][4] = -65.0 / 27.0;
        a_[6][5] = 125.0 / 54.0;
        
        a_[7][0] = 31.0 / 300.0;
        a_[7][4] = 61.0 / 225.0;
        a_[7][5] = -2.0 / 9.0;
        a_[7][6] = 13.0 / 900.0;
        
        a_[8][0] = 2.0;
        a_[8][3] = -53.0 / 6.0;
        a_[8][4] = 704.0 / 45.0;
        a_[8][5] = -107.0 / 9.0;
        a_[8][6] = 67.0 / 90.0;
        a_[8][7] = 3.0;
        
        a_[9][0] = -91.0 / 108.0;
        a_[9][3] = 23.0 / 108.0;
        a_[9][4] = -976.0 / 135.0;
        a_[9][5] = 311.0 / 54.0;
        a_[9][6] = -19.0 / 60.0;
        a_[9][7] = 17.0 / 6.0;
        a_[9][8] = -1.0 / 12.0;
        
        a_[10][0] = 2383.0 / 4100.0;
        a_[10][3] = -341.0 / 164.0;
        a_[10][4] = 4496.0 / 1025.0;
        a_[10][5] = -301.0 / 82.0;
        a_[10][6] = 2133.0 / 4100.0;
        a_[10][7] = 45.0 / 82.0;
        a_[10][8] = 45.0 / 164.0;
        a_[10][9] = 18.0 / 41.0;
        
        a_[11][0] = 3.0 / 205.0;
        a_[11][5] = -6.0 / 41.0;
        a_[11][6] = -3.0 / 205.0;
        a_[11][7] = -3.0 / 41.0;
        a_[11][8] = 3.0 / 41.0;
        a_[11][9] = 6.0 / 41.0;
        
        a_[12][0] = -1777.0 / 4100.0;
        a_[12][3] = -341.0 / 164.0;
        a_[12][4] = 4496.0 / 1025.0;
        a_[12][5] = -289.0 / 82.0;
        a_[12][6] = 2193.0 / 4100.0;
        a_[12][7] = 51.0 / 82.0;
        a_[12][8] = 33.0 / 164.0;
        a_[12][9] = 12.0 / 41.0;
        a_[12][11] = 1.0;
        
        // Pesi per soluzione di ordine 7
        b7_[0]  = 41.0 / 840.0;
        b7_[5]  = 34.0 / 105.0;
        b7_[6]  = 9.0 / 35.0;
        b7_[7]  = 9.0 / 35.0;
        b7_[8]  = 9.0 / 280.0;
        b7_[9]  = 9.0 / 280.0;
        b7_[10] = 41.0 / 840.0;
        
        // Pesi per soluzione di ordine 8 (usata per stima errore)
        b8_[5]  = 34.0 / 105.0;
        b8_[6]  = 9.0 / 35.0;
        b8_[7]  = 9.0 / 35.0;
        b8_[8]  = 9.0 / 280.0;
        b8_[9]  = 9.0 / 280.0;
        b8_[11] = 41.0 / 840.0;
        b8_[12] = 41.0 / 840.0;
    }
    
    double estimateInitialStep(DerivFunc& f, const State& y0, double t0, double dt) {
        // Stima del passo iniziale basata sulla scala del problema
        State f0 = f(t0, y0);
        stats_.function_evals++;
        
        double d0 = 0.0, d1 = 0.0;
        for (size_t i = 0; i < N; ++i) {
            double sc = atol_ + rtol_ * std::abs(y0[i]);
            d0 += (y0[i] / sc) * (y0[i] / sc);
            d1 += (f0[i] / sc) * (f0[i] / sc);
        }
        d0 = std::sqrt(d0 / N);
        d1 = std::sqrt(d1 / N);
        
        double h0 = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * d0 / d1;
        h0 = std::min(h0, 0.01 * std::abs(dt));
        
        return h0;
    }
    
    bool tryStep(DerivFunc& f, double t, const State& y, double h,
                 State& y_new, double& err) {
        // Array per i k_i
        std::array<State, STAGES> k;
        
        // Calcola k_0
        k[0] = f(t, y);
        stats_.function_evals++;
        
        // Calcola k_1 ... k_12
        for (int i = 1; i < STAGES; ++i) {
            State y_stage;
            for (size_t j = 0; j < N; ++j) {
                y_stage[j] = y[j];
                for (int l = 0; l < i; ++l) {
                    y_stage[j] += h * a_[i][l] * k[l][j];
                }
            }
            k[i] = f(t + c_[i] * h, y_stage);
            stats_.function_evals++;
        }
        
        // Calcola soluzione di ordine 7 e ordine 8
        State y7, y8;
        for (size_t j = 0; j < N; ++j) {
            y7[j] = y[j];
            y8[j] = y[j];
            for (int i = 0; i < STAGES; ++i) {
                y7[j] += h * b7_[i] * k[i][j];
                y8[j] += h * b8_[i] * k[i][j];
            }
        }
        
        // Stima dell'errore (differenza tra ordine 8 e 7)
        err = 0.0;
        for (size_t j = 0; j < N; ++j) {
            double sc = atol_ + rtol_ * std::max(std::abs(y[j]), std::abs(y8[j]));
            double e = (y8[j] - y7[j]) / sc;
            err += e * e;
        }
        err = std::sqrt(err / N);
        
        // Usa la soluzione di ordine 8 (local extrapolation)
        y_new = y8;
        
        return err <= 1.0;
    }
    
    double computeNewStep(double h, double err, bool accepted) {
        // Fattore di sicurezza
        constexpr double SAFETY = 0.9;
        
        // Limiti per il cambio di passo
        constexpr double FAC_MIN = 0.1;
        constexpr double FAC_MAX = 5.0;
        
        double fac;
        if (err == 0.0) {
            fac = FAC_MAX;
        } else {
            // Per RKF78: ordine = 8, quindi esponente = 1/8
            fac = SAFETY * std::pow(1.0 / err, 1.0 / 8.0);
            fac = std::clamp(fac, FAC_MIN, FAC_MAX);
        }
        
        // Se il passo è stato rifiutato, riduci di più
        if (!accepted) {
            fac = std::min(fac, 1.0);
        }
        
        return h * fac;
    }
};

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
    Vec3& operator+=(const Vec3& v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    double norm2() const { return x*x + y*y + z*z; }
    double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
};

struct OrbitalState {
    Vec3 pos;
    Vec3 vel;
};

struct AsteroidData {
    int number;
    const char* name;
    double gm, a, e, i, omega, Omega, M0, epoch, n;
};

// ============================================================================
// EFFEMERIDI PLANETARIE
// ============================================================================

class PlanetaryEphemeris {
public:
    enum Planet { MERCURY=0, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE };
    
    static OrbitalState getState(Planet planet, double jd_tdb) {
        double T = (jd_tdb - JD_J2000) / 36525.0;
        double a, e, i, L, omega_bar, Omega;
        
        switch (planet) {
            case MERCURY:
                a = 0.38709927; e = 0.20563593 + 0.00001906 * T;
                i = (7.00497902 - 0.00594749 * T) * DEG_TO_RAD;
                L = (252.25032350 + 149472.67411175 * T) * DEG_TO_RAD;
                omega_bar = (77.45779628 + 0.16047689 * T) * DEG_TO_RAD;
                Omega = (48.33076593 - 0.12534081 * T) * DEG_TO_RAD;
                break;
            case VENUS:
                a = 0.72333566; e = 0.00677672 - 0.00004107 * T;
                i = (3.39467605 - 0.00078890 * T) * DEG_TO_RAD;
                L = (181.97909950 + 58517.81538729 * T) * DEG_TO_RAD;
                omega_bar = (131.60246718 + 0.00268329 * T) * DEG_TO_RAD;
                Omega = (76.67984255 - 0.27769418 * T) * DEG_TO_RAD;
                break;
            case EARTH:
                a = 1.00000261 + 0.00000562 * T; e = 0.01671123 - 0.00004392 * T;
                i = (-0.00001531 - 0.01294668 * T) * DEG_TO_RAD;
                L = (100.46457166 + 35999.37244981 * T) * DEG_TO_RAD;
                omega_bar = (102.93768193 + 0.32327364 * T) * DEG_TO_RAD;
                Omega = 0.0;
                break;
            case MARS:
                a = 1.52371034; e = 0.09339410 + 0.00007882 * T;
                i = (1.84969142 - 0.00813131 * T) * DEG_TO_RAD;
                L = (-4.55343205 + 19140.30268499 * T) * DEG_TO_RAD;
                omega_bar = (-23.94362959 + 0.44441088 * T) * DEG_TO_RAD;
                Omega = (49.55953891 - 0.29257343 * T) * DEG_TO_RAD;
                break;
            case JUPITER:
                a = 5.20288700; e = 0.04838624 - 0.00013537 * T;
                i = (1.30439695 - 0.00183714 * T) * DEG_TO_RAD;
                L = (34.39644051 + 3034.74612775 * T) * DEG_TO_RAD;
                omega_bar = (14.72847983 + 0.21252668 * T) * DEG_TO_RAD;
                Omega = (100.47390909 + 0.20469106 * T) * DEG_TO_RAD;
                break;
            case SATURN:
                a = 9.53667594; e = 0.05386179 - 0.00050991 * T;
                i = (2.48599187 + 0.00193609 * T) * DEG_TO_RAD;
                L = (49.95424423 + 1222.49362201 * T) * DEG_TO_RAD;
                omega_bar = (92.59887831 - 0.41897216 * T) * DEG_TO_RAD;
                Omega = (113.66242448 - 0.28867794 * T) * DEG_TO_RAD;
                break;
            case URANUS:
                a = 19.18916464; e = 0.04725744 - 0.00004397 * T;
                i = (0.77263783 - 0.00242939 * T) * DEG_TO_RAD;
                L = (313.23810451 + 428.48202785 * T) * DEG_TO_RAD;
                omega_bar = (170.95427630 + 0.40805281 * T) * DEG_TO_RAD;
                Omega = (74.01692503 + 0.04240589 * T) * DEG_TO_RAD;
                break;
            case NEPTUNE:
                a = 30.06992276; e = 0.00859048 + 0.00005105 * T;
                i = (1.77004347 + 0.00035372 * T) * DEG_TO_RAD;
                L = (-55.12002969 + 218.45945325 * T) * DEG_TO_RAD;
                omega_bar = (44.96476227 - 0.32241464 * T) * DEG_TO_RAD;
                Omega = (131.78422574 - 0.00508664 * T) * DEG_TO_RAD;
                break;
            default: return OrbitalState();
        }
        return elementsToState(a, e, i, L, omega_bar, Omega, GM_SUN);
    }
    
private:
    static OrbitalState elementsToState(double a, double e, double i, double L, 
                                         double omega_bar, double Omega, double gm) {
        double omega = omega_bar - Omega;
        double M = std::fmod(L - omega_bar, TWO_PI);
        if (M < 0) M += TWO_PI;
        
        double E = M;
        for (int iter = 0; iter < 15; ++iter) {
            double dE = (E - e * std::sin(E) - M) / (1.0 - e * std::cos(E));
            E -= dE;
            if (std::abs(dE) < 1e-14) break;
        }
        
        double sin_E = std::sin(E), cos_E = std::cos(E);
        double sqrt_1_e2 = std::sqrt(1.0 - e * e);
        double nu = std::atan2(sqrt_1_e2 * sin_E, cos_E - e);
        double r = a * (1.0 - e * cos_E);
        
        double x_orb = r * std::cos(nu), y_orb = r * std::sin(nu);
        double v_factor = std::sqrt(gm * a) / r;
        double vx_orb = -v_factor * sin_E, vy_orb = v_factor * sqrt_1_e2 * cos_E;
        
        double cO = std::cos(Omega), sO = std::sin(Omega);
        double cw = std::cos(omega), sw = std::sin(omega);
        double ci = std::cos(i), si = std::sin(i);
        
        double P11 = cO*cw - sO*sw*ci, P12 = -cO*sw - sO*cw*ci;
        double P21 = sO*cw + cO*sw*ci, P22 = -sO*sw + cO*cw*ci;
        double P31 = sw*si, P32 = cw*si;
        
        OrbitalState s;
        s.pos = Vec3(P11*x_orb + P12*y_orb, P21*x_orb + P22*y_orb, P31*x_orb + P32*y_orb);
        s.vel = Vec3(P11*vx_orb + P12*vy_orb, P21*vx_orb + P22*vy_orb, P31*vx_orb + P32*vy_orb);
        return s;
    }
};

// ============================================================================
// AST17 - 16 ASTEROIDI MASSIVI
// ============================================================================

class AST17 {
public:
    static std::vector<AsteroidData> getAsteroids() {
        constexpr double epoch = 51544.5;
        return {
            {1, "Ceres", 62.6284, 2.7675, 0.0760, 10.593, 73.597, 80.3932, 113.410, epoch, 0.2141},
            {2, "Pallas", 14.28, 2.7730, 0.2299, 34.837, 310.049, 173.0962, 78.194, epoch, 0.2133},
            {4, "Vesta", 17.8, 2.3615, 0.0887, 7.142, 151.198, 103.8513, 205.546, epoch, 0.2716},
            {10, "Hygiea", 5.78, 3.1393, 0.1146, 3.831, 283.455, 283.2045, 107.590, epoch, 0.1794},
            {15, "Eunomia", 2.1, 2.6439, 0.1866, 11.737, 97.767, 293.2081, 142.405, epoch, 0.2262},
            {16, "Psyche", 1.81, 2.9216, 0.1339, 3.096, 227.305, 150.2873, 179.942, epoch, 0.1990},
            {31, "Euphrosyne", 1.26, 3.1515, 0.2257, 26.307, 61.588, 31.2873, 325.478, epoch, 0.1785},
            {52, "Europa", 1.59, 3.0958, 0.1020, 7.460, 343.545, 129.0380, 252.132, epoch, 0.1824},
            {65, "Cybele", 0.91, 3.4332, 0.1050, 3.562, 155.889, 155.5463, 196.234, epoch, 0.1636},
            {87, "Sylvia", 0.99, 3.4876, 0.0808, 10.866, 266.315, 73.3430, 271.478, epoch, 0.1610},
            {88, "Thisbe", 0.78, 2.7671, 0.1641, 5.222, 36.714, 276.0972, 355.489, epoch, 0.2141},
            {107, "Camilla", 0.75, 3.4886, 0.0692, 10.044, 309.785, 173.0458, 120.234, epoch, 0.1609},
            {324, "Bamberga", 0.69, 2.6835, 0.3381, 11.095, 43.678, 327.8394, 155.234, epoch, 0.2226},
            {451, "Patientia", 0.54, 3.0639, 0.0764, 15.240, 278.789, 96.9873, 234.456, epoch, 0.1843},
            {511, "Davida", 0.89, 3.1805, 0.1833, 15.942, 107.234, 270.0945, 186.789, epoch, 0.1768},
            {704, "Interamnia", 1.86, 3.0616, 0.1502, 17.296, 94.789, 280.3456, 145.234, epoch, 0.1844}
        };
    }
    
    static Vec3 getPosition(const AsteroidData& ast, double mjd_tdb) {
        double M = (ast.M0 + ast.n * (mjd_tdb - ast.epoch)) * DEG_TO_RAD;
        M = std::fmod(M, TWO_PI); if (M < 0) M += TWO_PI;
        
        double E = M;
        for (int iter = 0; iter < 15; ++iter) {
            double dE = (E - ast.e * std::sin(E) - M) / (1.0 - ast.e * std::cos(E));
            E -= dE;
            if (std::abs(dE) < 1e-14) break;
        }
        
        double sin_E = std::sin(E), cos_E = std::cos(E);
        double sqrt_1_e2 = std::sqrt(1.0 - ast.e * ast.e);
        double nu = std::atan2(sqrt_1_e2 * sin_E, cos_E - ast.e);
        double r = ast.a * (1.0 - ast.e * cos_E);
        double x_orb = r * std::cos(nu), y_orb = r * std::sin(nu);
        
        double i_rad = ast.i * DEG_TO_RAD;
        double omega_rad = ast.omega * DEG_TO_RAD;
        double Omega_rad = ast.Omega * DEG_TO_RAD;
        
        double cO = std::cos(Omega_rad), sO = std::sin(Omega_rad);
        double cw = std::cos(omega_rad), sw = std::sin(omega_rad);
        double ci = std::cos(i_rad), si = std::sin(i_rad);
        
        return Vec3(
            (cO*cw - sO*sw*ci)*x_orb + (-cO*sw - sO*cw*ci)*y_orb,
            (sO*cw + cO*sw*ci)*x_orb + (-sO*sw + cO*cw*ci)*y_orb,
            sw*si*x_orb + cw*si*y_orb
        );
    }
};

// ============================================================================
// PROPAGATORE CON RKF78
// ============================================================================

class FullPropagator {
public:
    FullPropagator(double rtol = 1e-12, double atol = 1e-15) 
        : asteroids_(AST17::getAsteroids()), integrator_(rtol, atol) {}
    
    OrbitalState propagate(const OrbitalState& initial, double jd0, double jdf) {
        mjd0_ = jd0 - 2400000.5;
        
        std::array<double, 6> y0 = {
            initial.pos.x, initial.pos.y, initial.pos.z,
            initial.vel.x, initial.vel.y, initial.vel.z
        };
        
        auto derivatives = [this](double t, const std::array<double, 6>& y) {
            return computeDerivatives(t, y);
        };
        
        double dt = jdf - jd0;
        auto yf = integrator_.integrate(derivatives, y0, 0.0, dt);
        
        OrbitalState final_state;
        final_state.pos = Vec3(yf[0], yf[1], yf[2]);
        final_state.vel = Vec3(yf[3], yf[4], yf[5]);
        return final_state;
    }
    
    const RKF78<6>::IntegrationStats& getStats() const { 
        return integrator_.getStats(); 
    }
    
private:
    std::vector<AsteroidData> asteroids_;
    RKF78<6> integrator_;
    double mjd0_;
    
    std::array<double, 6> computeDerivatives(double t, const std::array<double, 6>& y) {
        Vec3 pos(y[0], y[1], y[2]), vel(y[3], y[4], y[5]);
        double mjd = mjd0_ + t;
        double jd = mjd + 2400000.5;
        
        Vec3 acc = twoBodyAcceleration(pos);
        acc += planetaryPerturbations(pos, jd);
        acc += asteroidPerturbations(pos, mjd);
        acc += relativisticCorrection(pos, vel);
        
        return {vel.x, vel.y, vel.z, acc.x, acc.y, acc.z};
    }
    
    Vec3 twoBodyAcceleration(const Vec3& pos) {
        double r3 = std::pow(pos.norm(), 3);
        return pos * (-GM_SUN / r3);
    }
    
    Vec3 planetaryPerturbations(const Vec3& pos, double jd) {
        Vec3 acc;
        struct PG { PlanetaryEphemeris::Planet p; double gm; };
        PG planets[] = {
            {PlanetaryEphemeris::MERCURY, GM_MERCURY}, {PlanetaryEphemeris::VENUS, GM_VENUS},
            {PlanetaryEphemeris::EARTH, GM_EARTH}, {PlanetaryEphemeris::MARS, GM_MARS},
            {PlanetaryEphemeris::JUPITER, GM_JUPITER}, {PlanetaryEphemeris::SATURN, GM_SATURN},
            {PlanetaryEphemeris::URANUS, GM_URANUS}, {PlanetaryEphemeris::NEPTUNE, GM_NEPTUNE}
        };
        for (const auto& p : planets) {
            OrbitalState ps = PlanetaryEphemeris::getState(p.p, jd);
            double gm = p.gm * GM_CONV;
            Vec3 delta = ps.pos - pos;
            double dr3 = std::pow(delta.norm(), 3);
            double pr3 = std::pow(ps.pos.norm(), 3);
            acc += delta * (gm / dr3) - ps.pos * (gm / pr3);
        }
        return acc;
    }
    
    Vec3 asteroidPerturbations(const Vec3& pos, double mjd) {
        Vec3 acc;
        for (const auto& ast : asteroids_) {
            Vec3 ap = AST17::getPosition(ast, mjd);
            double gm = ast.gm * GM_CONV;
            Vec3 delta = ap - pos;
            double dr3 = std::pow(delta.norm(), 3);
            double ar3 = std::pow(ap.norm(), 3);
            acc += delta * (gm / dr3) - ap * (gm / ar3);
        }
        return acc;
    }
    
    Vec3 relativisticCorrection(const Vec3& pos, const Vec3& vel) {
        double r = pos.norm(), r3 = r*r*r;
        double v2 = vel.norm2();
        double rdot = pos.dot(vel) / r;
        double c2 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;
        Vec3 t1 = pos * (4.0 * GM_SUN / r - v2);
        Vec3 t2 = vel * (4.0 * rdot);
        return (t1 + t2) * (GM_SUN / (r3 * c2));
    }
};

// ============================================================================
// UTILITY
// ============================================================================

OrbitalState keplerToCartesian(double a, double e, double i, double Omega, 
                                double omega, double M0, double gm) {
    double E = M0;
    for (int iter = 0; iter < 15; ++iter) {
        double dE = (E - e * std::sin(E) - M0) / (1.0 - e * std::cos(E));
        E -= dE;
        if (std::abs(dE) < 1e-14) break;
    }
    
    double sin_E = std::sin(E), cos_E = std::cos(E);
    double sqrt_1_e2 = std::sqrt(1.0 - e * e);
    double nu = std::atan2(sqrt_1_e2 * sin_E, cos_E - e);
    double r = a * (1.0 - e * cos_E);
    
    double x_orb = r * std::cos(nu), y_orb = r * std::sin(nu);
    double v_factor = std::sqrt(gm * a) / r;
    double vx_orb = -v_factor * sin_E, vy_orb = v_factor * sqrt_1_e2 * cos_E;
    
    double cO = std::cos(Omega), sO = std::sin(Omega);
    double cw = std::cos(omega), sw = std::sin(omega);
    double ci = std::cos(i), si = std::sin(i);
    
    double P11 = cO*cw - sO*sw*ci, P12 = -cO*sw - sO*cw*ci;
    double P21 = sO*cw + cO*sw*ci, P22 = -sO*sw + cO*cw*ci;
    double P31 = sw*si, P32 = cw*si;
    
    OrbitalState s;
    s.pos = Vec3(P11*x_orb + P12*y_orb, P21*x_orb + P22*y_orb, P31*x_orb + P32*y_orb);
    s.vel = Vec3(P11*vx_orb + P12*vy_orb, P21*vx_orb + P22*vy_orb, P31*vx_orb + P32*vy_orb);
    return s;
}

Vec3 eclipticToEquatorial(const Vec3& ecl) {
    double ce = std::cos(OBLIQUITY_J2000), se = std::sin(OBLIQUITY_J2000);
    return Vec3(ecl.x, ecl.y*ce - ecl.z*se, ecl.y*se + ecl.z*ce);
}

void cartesianToRADec(const Vec3& pos, double& ra, double& dec) {
    double r = pos.norm();
    dec = std::asin(pos.z / r);
    ra = std::atan2(pos.y, pos.x);
    if (ra < 0) ra += TWO_PI;
}

std::string formatRA(double ra_rad) {
    double h = ra_rad * 12.0 / PI;
    int hi = (int)h; double m = (h - hi) * 60; int mi = (int)m;
    double s = (m - mi) * 60;
    char buf[32]; std::snprintf(buf, 32, "%02d %02d %06.3f", hi, mi, s);
    return buf;
}

std::string formatDec(double dec_rad) {
    double d = dec_rad * RAD_TO_DEG;
    char sign = d >= 0 ? '+' : '-'; d = std::abs(d);
    int di = (int)d; double m = (d - di) * 60; int mi = (int)m;
    double s = (m - mi) * 60;
    char buf[32]; std::snprintf(buf, 32, "%c%02d %02d %05.2f", sign, di, mi, s);
    return buf;
}

// ============================================================================
// MAIN - TEST RKF78
// ============================================================================

int main() {
    std::cout << std::fixed;
    
    std::cout << "================================================================\n";
    std::cout << "  TEST INTEGRATORE RKF78 (Fehlberg 1968)\n";
    std::cout << "  13 stadi, ordine 7(8), controllo adattivo del passo\n";
    std::cout << "================================================================\n\n";
    
    // Elementi kepleriani (17030) Sierks
    double a = 3.1754733;
    double e = 0.0454207;
    double i = 2.9046 * DEG_TO_RAD;
    double Omega = 104.16243 * DEG_TO_RAD;
    double omega = 100.5141 * DEG_TO_RAD;
    double M0 = 229.79088 * DEG_TO_RAD;
    double jd_epoch = 2461000.5;
    double jd_target = 2461008.0913;
    
    // Stato iniziale
    OrbitalState s0 = keplerToCartesian(a, e, i, Omega, omega, M0, GM_SUN);
    
    std::cout << "ASTEROIDE: (17030) Sierks\n";
    std::cout << "Epoca:  JD " << std::setprecision(4) << jd_epoch << "\n";
    std::cout << "Target: JD " << jd_target << "\n";
    std::cout << "Δt = " << (jd_target - jd_epoch) << " giorni\n\n";
    
    std::cout << "STATO INIZIALE:\n";
    std::cout << std::setprecision(15);
    std::cout << "  r = [" << s0.pos.x << ", " << s0.pos.y << ", " << s0.pos.z << "] AU\n";
    std::cout << "  v = [" << s0.vel.x << ", " << s0.vel.y << ", " << s0.vel.z << "] AU/day\n\n";
    
    // Propagatore con RKF78
    FullPropagator prop(1e-14, 1e-17);  // Alta precisione
    
    // STEP 1: Forward
    std::cout << "STEP 1: Propagazione AVANTI (epoca → target)\n";
    std::cout << "================================================\n";
    OrbitalState s1 = prop.propagate(s0, jd_epoch, jd_target);
    
    auto stats1 = prop.getStats();
    std::cout << "Passi accettati: " << stats1.steps_accepted << "\n";
    std::cout << "Passi rifiutati: " << stats1.steps_rejected << "\n";
    std::cout << "Valutazioni f:   " << stats1.function_evals << "\n";
    std::cout << "Passo min:       " << std::scientific << std::setprecision(3) << stats1.min_step << " giorni\n";
    std::cout << "Passo max:       " << stats1.max_step << " giorni\n\n";
    
    std::cout << std::fixed << std::setprecision(15);
    std::cout << "Stato finale:\n";
    std::cout << "  r = [" << s1.pos.x << ", " << s1.pos.y << ", " << s1.pos.z << "] AU\n";
    std::cout << "  v = [" << s1.vel.x << ", " << s1.vel.y << ", " << s1.vel.z << "] AU/day\n\n";
    
    // Calcola RA/Dec
    OrbitalState earth = PlanetaryEphemeris::getState(PlanetaryEphemeris::EARTH, jd_target);
    Vec3 geo = s1.pos - earth.pos;
    double lt = geo.norm() / SPEED_OF_LIGHT;
    OrbitalState s1_lt = prop.propagate(s0, jd_epoch, jd_target - lt);
    Vec3 geo_lt = s1_lt.pos - earth.pos;
    Vec3 geo_eq = eclipticToEquatorial(geo_lt);
    double ra1, dec1;
    cartesianToRADec(geo_eq, ra1, dec1);
    
    std::cout << "POSIZIONE (con correzione tempo luce):\n";
    std::cout << "  RA  = " << formatRA(ra1) << "\n";
    std::cout << "  Dec = " << formatDec(dec1) << "\n\n";
    
    // STEP 2: Backward
    std::cout << "STEP 2: Propagazione INDIETRO (target → epoca)\n";
    std::cout << "================================================\n";
    OrbitalState s2 = prop.propagate(s1, jd_target, jd_epoch);
    
    auto stats2 = prop.getStats();
    std::cout << "Passi accettati: " << stats2.steps_accepted << "\n";
    std::cout << "Passi rifiutati: " << stats2.steps_rejected << "\n";
    std::cout << "Valutazioni f:   " << stats2.function_evals << "\n\n";
    
    std::cout << "Stato dopo round-trip:\n";
    std::cout << "  r = [" << s2.pos.x << ", " << s2.pos.y << ", " << s2.pos.z << "] AU\n";
    std::cout << "  v = [" << s2.vel.x << ", " << s2.vel.y << ", " << s2.vel.z << "] AU/day\n\n";
    
    // VERIFICA
    std::cout << "================================================================\n";
    std::cout << "  VERIFICA ROUND-TRIP RKF78\n";
    std::cout << "================================================================\n\n";
    
    Vec3 dr = s2.pos - s0.pos;
    Vec3 dv = s2.vel - s0.vel;
    
    std::cout << "ERRORE POSIZIONE:\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  Δx = " << dr.x << " AU = " << std::fixed << std::setprecision(6) << (dr.x * AU_KM * 1000) << " m\n";
    std::cout << "  Δy = " << std::scientific << dr.y << " AU = " << std::fixed << (dr.y * AU_KM * 1000) << " m\n";
    std::cout << "  Δz = " << std::scientific << dr.z << " AU = " << std::fixed << (dr.z * AU_KM * 1000) << " m\n";
    std::cout << "  |Δr| = " << std::scientific << dr.norm() << " AU = " 
              << std::fixed << std::setprecision(3) << (dr.norm() * AU_KM * 1000) << " m\n\n";
    
    std::cout << "ERRORE VELOCITÀ:\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  |Δv| = " << dv.norm() << " AU/day = " 
              << std::fixed << std::setprecision(9) << (dv.norm() * AU_KM / DAY_SEC) << " km/s\n\n";
    
    // Confronto con JPL
    std::cout << "================================================================\n";
    std::cout << "  CONFRONTO CON JPL HORIZONS\n";
    std::cout << "================================================================\n\n";
    
    std::cout << "AL TARGET (JD " << std::setprecision(4) << jd_target << "):\n";
    std::cout << "  JPL Horizons: 04 53 11.25   +20 19 25.8\n";
    std::cout << "  RKF78:        " << formatRA(ra1) << "   " << formatDec(dec1) << "\n";
    
    double ra_jpl = (4.0 + 53.0/60.0 + 11.25/3600.0) * 15.0 * DEG_TO_RAD;
    double dec_jpl = (20.0 + 19.0/60.0 + 25.8/3600.0) * DEG_TO_RAD;
    double dra = (ra1 - ra_jpl) * RAD_TO_DEG * 3600.0 * std::cos(dec1);
    double ddec = (dec1 - dec_jpl) * RAD_TO_DEG * 3600.0;
    std::cout << "  Errore: ΔRA*cos(δ)=" << std::fixed << std::setprecision(2) << dra 
              << "\", ΔDec=" << ddec << "\" → Tot=" 
              << std::sqrt(dra*dra + ddec*ddec) << "\"\n\n";
    
    double err_m = dr.norm() * AU_KM * 1000;
    if (err_m < 1.0) {
        std::cout << "✅ RKF78 ROUND-TRIP: errore < 1 metro!\n";
    } else if (err_m < 1000.0) {
        std::cout << "✅ RKF78 ROUND-TRIP: errore = " << err_m << " m\n";
    } else {
        std::cout << "⚠️  RKF78 ROUND-TRIP: errore = " << (err_m/1000) << " km\n";
    }
    
    std::cout << "\n================================================================\n";
    
    return 0;
}
