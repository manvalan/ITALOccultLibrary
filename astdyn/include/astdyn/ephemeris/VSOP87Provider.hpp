#ifndef ASTDYN_VSOP87_PROVIDER_HPP
#define ASTDYN_VSOP87_PROVIDER_HPP

#include "astdyn/ephemeris/EphemerisProvider.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"

namespace astdyn::ephemeris {

class VSOP87Provider : public astdyn::ephemeris::EphemerisProvider {
public:
    VSOP87Provider() = default;
    
    Eigen::Vector3d getPosition(astdyn::ephemeris::CelestialBody body, double jd_tdb) override {
        return PlanetaryEphemeris::getPosition(body, jd_tdb);
    }
    
    Eigen::Vector3d getVelocity(astdyn::ephemeris::CelestialBody body, double jd_tdb) override {
        return PlanetaryEphemeris::getVelocity(body, jd_tdb);
    }
    
    std::string getName() const override {
        return "VSOP87";
    }
    
    double getAccuracy() const override {
        return 20.0;
    }
    
    bool isAvailable() const override {
        return true;
    }
};

} // namespace astdyn::ephemeris

#endif // ASTDYN_VSOP87_PROVIDER_HPP
