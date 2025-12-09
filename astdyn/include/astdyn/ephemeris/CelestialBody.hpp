/**
 * @file CelestialBody.hpp
 * @brief Common celestial body identifiers
 * @date 2025-12-09
 */

#ifndef ASTDYN_CELESTIAL_BODY_HPP
#define ASTDYN_CELESTIAL_BODY_HPP

namespace astdyn::ephemeris {

/**
 * @brief Celestial body identifiers
 * Compatible with VSOP87 body IDs and can be mapped to NAIF IDs
 */
enum class CelestialBody {
    SUN = 10,
    MERCURY = 1,
    VENUS = 2,
    EARTH = 3,
    MARS = 4,
    JUPITER = 5,
    SATURN = 6,
    URANUS = 7,
    NEPTUNE = 8,
    PLUTO = 9,
    MOON = 301
};

} // namespace astdyn::ephemeris

#endif // ASTDYN_CELESTIAL_BODY_HPP
