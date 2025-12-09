/**
 * @file ResidualCalculator.cpp
 * @brief Implementation of residual calculator
 */

#include "astdyn/orbit_determination/ResidualCalculator.hpp"
#include <cmath>

namespace astdyn::orbit_determination {

ResidualCalculator::ResidualCalculator() {
    // Default: use VSOP87 for Earth position
    ephemeris_func_ = [](int body_id, double mjd) {
        // Simplified: Earth at origin (geocentric observations)
        return Eigen::Vector3d::Zero();
    };
}

void ResidualCalculator::cartesian_to_radec(
    const Eigen::Vector3d& r,
    double& ra_deg,
    double& dec_deg
) const {
    double x = r(0);
    double y = r(1);
    double z = r(2);
    
    double r_norm = r.norm();
    
    // RA = atan2(y, x)
    ra_deg = std::atan2(y, x) * 180.0 / M_PI;
    if (ra_deg < 0) ra_deg += 360.0;
    
    // Dec = asin(z / r)
    dec_deg = std::asin(z / r_norm) * 180.0 / M_PI;
}

Eigen::Vector3d ResidualCalculator::apply_light_time(
    const Eigen::Vector3d& r_helio,
    const Eigen::Vector3d& r_obs,
    double& light_time_days
) const {
    // Iterative light-time correction
    constexpr double c_au_per_day = 173.1446;  // Speed of light [AU/day]
    constexpr int max_iter = 3;
    
    Eigen::Vector3d r_topo = r_helio - r_obs;
    light_time_days = r_topo.norm() / c_au_per_day;
    
    // For now, simple approximation (no iteration)
    // Full version would propagate back by light_time_days
    
    return r_topo;
}

Eigen::Vector3d ResidualCalculator::get_earth_position(double mjd) const {
    // Use ephemeris function (body_id 3 = Earth)
    return ephemeris_func_(3, mjd);
}

Eigen::Vector3d ResidualCalculator::get_observatory_position(
    const std::string& code,
    double mjd
) const {
    // Simplified: observatory at geocenter
    // Full version would use MPC observatory database
    return Eigen::Vector3d::Zero();
}

Residual ResidualCalculator::compute_residual(
    const Observation& obs,
    const Eigen::Vector<double, 6>& state0,
    double epoch0_mjd
) {
    Residual res;
    res.epoch_mjd = obs.epoch_mjd;
    res.rejected = obs.rejected;
    
    // For now: simplified version without propagation
    // Assume state0 is already at obs epoch
    Eigen::Vector3d r_helio = state0.head<3>();
    
    // Get observer position
    Eigen::Vector3d r_earth = get_earth_position(obs.epoch_mjd);
    Eigen::Vector3d r_obs = get_observatory_position(obs.observatory_code, obs.epoch_mjd);
    
    // Topocentric position
    double light_time;
    Eigen::Vector3d r_topo = apply_light_time(r_helio, r_earth + r_obs, light_time);
    
    // Convert to RA/Dec
    double ra_comp, dec_comp;
    cartesian_to_radec(r_topo, ra_comp, dec_comp);
    
    res.ra_computed_deg = ra_comp;
    res.dec_computed_deg = dec_comp;
    
    // Compute residuals (O-C) in arcsec
    double ra_obs_deg = obs.ra_deg;
    double dec_obs_deg = obs.dec_deg;
    
    res.ra_residual_arcsec = (ra_obs_deg - ra_comp) * 3600.0 * std::cos(dec_obs_deg * M_PI / 180.0);
    res.dec_residual_arcsec = (dec_obs_deg - dec_comp) * 3600.0;
    
    return res;
}

std::vector<Residual> ResidualCalculator::compute_residuals(
    const std::vector<Observation>& observations,
    const Eigen::Vector<double, 6>& state0,
    double epoch0_mjd
) {
    std::vector<Residual> residuals;
    residuals.reserve(observations.size());
    
    for (const auto& obs : observations) {
        residuals.push_back(compute_residual(obs, state0, epoch0_mjd));
    }
    
    return residuals;
}

Observatory ResidualCalculator::get_observatory(const std::string& code) const {
    auto it = observatories_.find(code);
    if (it != observatories_.end()) {
        return it->second;
    }
    
    // Default: geocenter
    Observatory obs;
    obs.code = code;
    obs.position = Eigen::Vector3d::Zero();
    obs.longitude_deg = 0.0;
    obs.latitude_deg = 0.0;
    obs.altitude_m = 0.0;
    return obs;
}

void ResidualCalculator::load_observatories(const std::string& filename) {
    // TODO: Load MPC observatory database
    // For now, empty
}

} // namespace astdyn::orbit_determination
