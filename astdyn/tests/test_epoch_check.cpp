#include <iostream>
#include <cmath>

const double AU = 1.495978707e8;
const double GM_SUN = 1.32712440018e11;
const double DEG_TO_RAD = M_PI / 180.0;
const double RAD_TO_DEG = 180.0 / M_PI;

const double JPL_A = 3.173489964321051;
const double JPL_E = 0.04796607451625862;
const double JPL_INC = 2.904309538190326 * DEG_TO_RAD;
const double JPL_OMEGA = 102.1497438064497 * DEG_TO_RAD;
const double JPL_OMEGA_NODE = 104.1845838362649 * DEG_TO_RAD;
const double JPL_M0 = 99.03517819281583 * DEG_TO_RAD;

double solve_kepler(double M, double e) {
    double E = M;
    for (int i = 0; i < 30; i++) {
        double dE = (E - e*sin(E) - M) / (1.0 - e*cos(E));
        E -= dE;
        if (std::abs(dE) < 1e-12) break;
    }
    return E;
}

int main() {
    // Test all'epoca (M = M0)
    double E = solve_kepler(JPL_M0, JPL_E);
    double a_km = JPL_A * AU;
    double sqrt_1_e2 = sqrt(1.0 - JPL_E*JPL_E);
    
    double x_orb = a_km * (cos(E) - JPL_E);
    double y_orb = a_km * sqrt_1_e2 * sin(E);
    
    double cos_omega = cos(JPL_OMEGA);
    double sin_omega = sin(JPL_OMEGA);
    double cos_Omega = cos(JPL_OMEGA_NODE);
    double sin_Omega = sin(JPL_OMEGA_NODE);
    double cos_inc = cos(JPL_INC);
    double sin_inc = sin(JPL_INC);
    
    // Eclittico
    double x_ecl = (cos_Omega*cos_omega - sin_Omega*sin_omega*cos_inc) * x_orb +
                   (-cos_Omega*sin_omega - sin_Omega*cos_omega*cos_inc) * y_orb;
    double y_ecl = (sin_Omega*cos_omega + cos_Omega*sin_omega*cos_inc) * x_orb +
                   (-sin_Omega*sin_omega + cos_Omega*cos_omega*cos_inc) * y_orb;
    double z_ecl = sin_omega*sin_inc * x_orb + cos_omega*sin_inc * y_orb;
    
    std::cout << "Posizione eclittica (km):\n";
    std::cout << "  X = " << x_ecl << "\n";
    std::cout << "  Y = " << y_ecl << "\n";
    std::cout << "  Z = " << z_ecl << "\n\n";
    
    // Equatoriale
    const double eps = 23.43928 * DEG_TO_RAD;
    double x_eq = x_ecl;
    double y_eq = y_ecl * cos(eps) - z_ecl * sin(eps);
    double z_eq = y_ecl * sin(eps) + z_ecl * cos(eps);
    
    double r = sqrt(x_eq*x_eq + y_eq*y_eq + z_eq*z_eq);
    double dec = asin(z_eq / r) * RAD_TO_DEG;
    double ra = atan2(y_eq, x_eq) * RAD_TO_DEG;
    if (ra < 0) ra += 360.0;
    
    std::cout << "RA/Dec all'epoca 2018-03-16:\n";
    std::cout << "  RA  = " << ra << "° (JPL: 323.5°)\n";
    std::cout << "  Dec = " << dec << "° (JPL: -15.56°)\n";
    
    return 0;
}
