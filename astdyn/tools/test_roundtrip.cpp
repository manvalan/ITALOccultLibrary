/**
 * @file test_roundtrip.cpp
 * @brief Test Round-Trip: propaga avanti e poi indietro all'epoca iniziale
 *        Verifica la precisione dell'integratore confrontando con JPL Horizons
 * 
 * Test:
 * 1. Parti dall'epoca MPC (JD 2461000.5)
 * 2. Propaga al target (JD 2461008.0913) - confronta con Horizons
 * 3. Propaga indietro all'epoca iniziale
 * 4. Verifica che torni allo stato iniziale
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
// COSTANTI (identiche a sierks_precision.cpp)
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
// STRUTTURE
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
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    double norm2() const { return x*x + y*y + z*z; }
    double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
};

struct State {
    Vec3 pos;
    Vec3 vel;
};

struct KeplerianElements {
    double a, e, i, Omega, omega, M0, epoch, gm;
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
    
    static State getState(Planet planet, double jd_tdb) {
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
            default: return State();
        }
        return elementsToState(a, e, i, L, omega_bar, Omega, GM_SUN);
    }
    
private:
    static State elementsToState(double a, double e, double i, double L, 
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
        
        State s;
        s.pos = Vec3(P11*x_orb + P12*y_orb, P21*x_orb + P22*y_orb, P31*x_orb + P32*y_orb);
        s.vel = Vec3(P11*vx_orb + P12*vy_orb, P21*vx_orb + P22*vy_orb, P31*vx_orb + P32*vy_orb);
        return s;
    }
};

// ============================================================================
// AST17
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
// INTEGRATORE RK4
// ============================================================================

class RK4 {
public:
    using DerivFunc = std::function<std::array<double, 6>(double, const std::array<double, 6>&)>;
    
    static std::array<double, 6> integrate(DerivFunc f, const std::array<double, 6>& y0,
                                           double t0, double tf, double h = 0.01) {
        std::array<double, 6> y = y0;
        double t = t0;
        int direction = (tf > t0) ? 1 : -1;
        h = direction * std::abs(h);
        
        while (direction * (tf - t) > 1e-12) {
            if (direction * (t + h - tf) > 0) h = tf - t;
            
            auto k1 = f(t, y);
            std::array<double, 6> y2, y3, y4;
            for (int i = 0; i < 6; ++i) y2[i] = y[i] + 0.5*h*k1[i];
            auto k2 = f(t + 0.5*h, y2);
            for (int i = 0; i < 6; ++i) y3[i] = y[i] + 0.5*h*k2[i];
            auto k3 = f(t + 0.5*h, y3);
            for (int i = 0; i < 6; ++i) y4[i] = y[i] + h*k3[i];
            auto k4 = f(t + h, y4);
            
            for (int i = 0; i < 6; ++i)
                y[i] += h * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6.0;
            t += h;
        }
        return y;
    }
};

// ============================================================================
// PROPAGATORE
// ============================================================================

class FullPropagator {
public:
    FullPropagator() : asteroids_(AST17::getAsteroids()) {}
    
    State propagate(const State& initial, double jd0, double jdf) {
        double mjd0 = jd0 - 2400000.5;
        std::array<double, 6> y0 = {
            initial.pos.x, initial.pos.y, initial.pos.z,
            initial.vel.x, initial.vel.y, initial.vel.z
        };
        
        auto derivatives = [this, mjd0](double t, const std::array<double, 6>& y) {
            return computeDerivatives(t, y, mjd0);
        };
        
        double dt = (jdf - jd0);
        auto yf = RK4::integrate(derivatives, y0, 0.0, dt, 0.01);
        
        State final_state;
        final_state.pos = Vec3(yf[0], yf[1], yf[2]);
        final_state.vel = Vec3(yf[3], yf[4], yf[5]);
        return final_state;
    }
    
private:
    std::vector<AsteroidData> asteroids_;
    
    std::array<double, 6> computeDerivatives(double t, const std::array<double, 6>& y, double mjd0) {
        Vec3 pos(y[0], y[1], y[2]), vel(y[3], y[4], y[5]);
        double mjd = mjd0 + t;
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
            State ps = PlanetaryEphemeris::getState(p.p, jd);
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

State keplerToCartesian(const KeplerianElements& kep) {
    double E = kep.M0;
    for (int iter = 0; iter < 15; ++iter) {
        double dE = (E - kep.e * std::sin(E) - kep.M0) / (1.0 - kep.e * std::cos(E));
        E -= dE;
        if (std::abs(dE) < 1e-14) break;
    }
    
    double sin_E = std::sin(E), cos_E = std::cos(E);
    double sqrt_1_e2 = std::sqrt(1.0 - kep.e * kep.e);
    double nu = std::atan2(sqrt_1_e2 * sin_E, cos_E - kep.e);
    double r = kep.a * (1.0 - kep.e * cos_E);
    
    double x_orb = r * std::cos(nu), y_orb = r * std::sin(nu);
    double v_factor = std::sqrt(kep.gm * kep.a) / r;
    double vx_orb = -v_factor * sin_E, vy_orb = v_factor * sqrt_1_e2 * cos_E;
    
    double cO = std::cos(kep.Omega), sO = std::sin(kep.Omega);
    double cw = std::cos(kep.omega), sw = std::sin(kep.omega);
    double ci = std::cos(kep.i), si = std::sin(kep.i);
    
    double P11 = cO*cw - sO*sw*ci, P12 = -cO*sw - sO*cw*ci;
    double P21 = sO*cw + cO*sw*ci, P22 = -sO*sw + cO*cw*ci;
    double P31 = sw*si, P32 = cw*si;
    
    State s;
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
// MAIN - TEST ROUND TRIP
// ============================================================================

int main() {
    std::cout << std::fixed << std::setprecision(15);
    
    std::cout << "================================================================\n";
    std::cout << "  TEST ROUND-TRIP: (17030) Sierks\n";
    std::cout << "  Epoca → Target → Epoca (verifica consistenza integratore)\n";
    std::cout << "================================================================\n\n";
    
    // Elementi kepleriani
    KeplerianElements sierks;
    sierks.a = 3.1754733;
    sierks.e = 0.0454207;
    sierks.i = 2.9046 * DEG_TO_RAD;
    sierks.Omega = 104.16243 * DEG_TO_RAD;
    sierks.omega = 100.5141 * DEG_TO_RAD;
    sierks.M0 = 229.79088 * DEG_TO_RAD;
    sierks.epoch = 2461000.5;
    sierks.gm = GM_SUN;
    
    double jd_epoch = 2461000.5;
    double jd_target = 2461008.0913;
    
    // Stato iniziale
    State s0 = keplerToCartesian(sierks);
    
    std::cout << "EPOCA INIZIALE: JD " << jd_epoch << "\n";
    std::cout << "================================================\n";
    std::cout << "Posizione: [" << s0.pos.x << ", " << s0.pos.y << ", " << s0.pos.z << "] AU\n";
    std::cout << "Velocità:  [" << s0.vel.x << ", " << s0.vel.y << ", " << s0.vel.z << "] AU/day\n";
    std::cout << "|r| = " << s0.pos.norm() << " AU\n\n";
    
    // Propagatore
    FullPropagator prop;
    
    // STEP 1: Epoca → Target
    std::cout << "STEP 1: Propaga Epoca → Target (Δt = +" << (jd_target - jd_epoch) << " giorni)\n";
    std::cout << "================================================\n";
    State s1 = prop.propagate(s0, jd_epoch, jd_target);
    std::cout << "Posizione: [" << s1.pos.x << ", " << s1.pos.y << ", " << s1.pos.z << "] AU\n";
    std::cout << "Velocità:  [" << s1.vel.x << ", " << s1.vel.y << ", " << s1.vel.z << "] AU/day\n";
    std::cout << "|r| = " << s1.pos.norm() << " AU\n\n";
    
    // Calcola RA/Dec al target
    State earth_t = PlanetaryEphemeris::getState(PlanetaryEphemeris::EARTH, jd_target);
    Vec3 geo_ecl = s1.pos - earth_t.pos;
    double delta = geo_ecl.norm();
    
    // Correzione tempo luce
    double lt = delta / SPEED_OF_LIGHT;
    State s1_lt = prop.propagate(s0, jd_epoch, jd_target - lt);
    Vec3 geo_ecl_lt = s1_lt.pos - earth_t.pos;
    Vec3 geo_eq = eclipticToEquatorial(geo_ecl_lt);
    double ra1, dec1;
    cartesianToRADec(geo_eq, ra1, dec1);
    
    std::cout << "POSIZIONE AL TARGET (con light-time):\n";
    std::cout << "  RA  = " << formatRA(ra1) << "\n";
    std::cout << "  Dec = " << formatDec(dec1) << "\n";
    std::cout << "  Δ   = " << geo_ecl_lt.norm() << " AU\n\n";
    
    // STEP 2: Target → Epoca (indietro)
    std::cout << "STEP 2: Propaga Target → Epoca (Δt = " << (jd_epoch - jd_target) << " giorni)\n";
    std::cout << "================================================\n";
    State s2 = prop.propagate(s1, jd_target, jd_epoch);
    std::cout << "Posizione: [" << s2.pos.x << ", " << s2.pos.y << ", " << s2.pos.z << "] AU\n";
    std::cout << "Velocità:  [" << s2.vel.x << ", " << s2.vel.y << ", " << s2.vel.z << "] AU/day\n";
    std::cout << "|r| = " << s2.pos.norm() << " AU\n\n";
    
    // VERIFICA ROUND-TRIP
    std::cout << "================================================================\n";
    std::cout << "  VERIFICA ROUND-TRIP\n";
    std::cout << "================================================================\n\n";
    
    Vec3 dr = s2.pos - s0.pos;
    Vec3 dv = s2.vel - s0.vel;
    
    std::cout << "ERRORE POSIZIONE:\n";
    std::cout << "  Δx = " << std::scientific << std::setprecision(6) << dr.x << " AU = " 
              << std::fixed << std::setprecision(3) << (dr.x * AU_KM) << " km\n";
    std::cout << "  Δy = " << std::scientific << std::setprecision(6) << dr.y << " AU = " 
              << std::fixed << std::setprecision(3) << (dr.y * AU_KM) << " km\n";
    std::cout << "  Δz = " << std::scientific << std::setprecision(6) << dr.z << " AU = " 
              << std::fixed << std::setprecision(3) << (dr.z * AU_KM) << " km\n";
    std::cout << "  |Δr| = " << std::scientific << std::setprecision(6) << dr.norm() << " AU = " 
              << std::fixed << std::setprecision(3) << (dr.norm() * AU_KM) << " km\n\n";
    
    std::cout << "ERRORE VELOCITÀ:\n";
    std::cout << "  Δvx = " << std::scientific << std::setprecision(6) << dv.x << " AU/day\n";
    std::cout << "  Δvy = " << std::scientific << std::setprecision(6) << dv.y << " AU/day\n";
    std::cout << "  Δvz = " << std::scientific << std::setprecision(6) << dv.z << " AU/day\n";
    std::cout << "  |Δv| = " << std::scientific << std::setprecision(6) << dv.norm() << " AU/day = "
              << std::fixed << std::setprecision(6) << (dv.norm() * AU_KM / DAY_SEC) << " km/s\n\n";
    
    // Calcola RA/Dec all'epoca (per confronto con Horizons)
    State earth_0 = PlanetaryEphemeris::getState(PlanetaryEphemeris::EARTH, jd_epoch);
    Vec3 geo_ecl_0 = s0.pos - earth_0.pos;
    double delta_0 = geo_ecl_0.norm();
    double lt_0 = delta_0 / SPEED_OF_LIGHT;
    
    // Per l'epoca iniziale, gli elementi sono già "al posto giusto"
    Vec3 geo_eq_0 = eclipticToEquatorial(geo_ecl_0);
    double ra0, dec0;
    cartesianToRADec(geo_eq_0, ra0, dec0);
    
    std::cout << "================================================================\n";
    std::cout << "  CONFRONTO CON JPL HORIZONS\n";
    std::cout << "================================================================\n\n";
    
    std::cout << "AL TARGET (JD " << std::fixed << std::setprecision(4) << jd_target << "):\n";
    std::cout << "  JPL Horizons: 04 53 11.25   +20 19 25.8\n";
    std::cout << "  AstDyn:       " << formatRA(ra1) << "   " << formatDec(dec1) << "\n";
    
    double ra_jpl = (4.0 + 53.0/60.0 + 11.25/3600.0) * 15.0 * DEG_TO_RAD;
    double dec_jpl = (20.0 + 19.0/60.0 + 25.8/3600.0) * DEG_TO_RAD;
    double dra = (ra1 - ra_jpl) * RAD_TO_DEG * 3600.0 * std::cos(dec1);
    double ddec = (dec1 - dec_jpl) * RAD_TO_DEG * 3600.0;
    std::cout << "  Errore: ΔRA*cos(δ)=" << std::fixed << std::setprecision(2) << dra 
              << "\", ΔDec=" << ddec << "\" → Tot=" 
              << std::sqrt(dra*dra + ddec*ddec) << "\"\n\n";
    
    std::cout << "ALL'EPOCA (JD " << std::fixed << std::setprecision(4) << jd_epoch << "):\n";
    std::cout << "  AstDyn:  " << formatRA(ra0) << "   " << formatDec(dec0) << "\n";
    std::cout << "  (Query Horizons per verifica)\n\n";
    
    std::cout << "================================================================\n";
    std::cout << "  TEST COMPLETATO\n";
    std::cout << "================================================================\n";
    
    double err_pos_km = dr.norm() * AU_KM;
    double err_vel_ms = dv.norm() * AU_KM / DAY_SEC * 1000.0;
    
    if (err_pos_km < 1.0 && err_vel_ms < 1.0) {
        std::cout << "✅ ROUND-TRIP OK: errore < 1 km in posizione\n";
    } else {
        std::cout << "⚠️  ROUND-TRIP: errore = " << err_pos_km << " km\n";
    }
    
    return 0;
}
