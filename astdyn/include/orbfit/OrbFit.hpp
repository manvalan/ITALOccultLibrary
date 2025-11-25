/**
 * @file OrbFit.hpp
 * @brief Main include file for OrbFit library
 * @author OrbFit C++ Conversion Team
 * @date 2025-11-23
 * 
 * Include this file to access all OrbFit functionality.
 */

#ifndef ORBFIT_HPP
#define ORBFIT_HPP

// Version and configuration
#include "orbfit/Version.hpp"
#include "orbfit/Config.hpp"

// Core types and constants
#include "orbfit/core/Types.hpp"
#include "orbfit/core/Constants.hpp"

// Utilities (to be added in Phase 2)
// #include "orbfit/utils/Logger.hpp"
// #include "orbfit/utils/StringUtils.hpp"

// Time (to be added in Phase 2)
// #include "orbfit/time/TimeScale.hpp"
// #include "orbfit/time/TimeConversions.hpp"

// Math (to be added in Phase 2)
// #include "orbfit/math/MathUtils.hpp"
// #include "orbfit/math/LinearAlgebra.hpp"

// Ephemeris (to be added in Phase 3)
// #include "orbfit/ephemeris/JPLEphemeris.hpp"

// Observations (to be added in Phase 4)
// #include "orbfit/observations/Observation.hpp"
// #include "orbfit/observations/ObservationReader.hpp"

// Orbital elements (to be added in Phase 5)
// #include "orbfit/orbit/KeplerianElements.hpp"
// #include "orbfit/orbit/StateVector.hpp"

// Propagation (to be added in Phase 6)
// #include "orbfit/propagation/OrbitPropagator.hpp"
// #include "orbfit/propagation/ForceModel.hpp"

namespace orbfit {

/**
 * @brief Initialize OrbFit library
 * 
 * Call this function once at program startup to initialize
 * the library and set up any global state.
 * 
 * @return true if initialization was successful
 */
inline bool initialize() {
    // Future: Initialize logging, load configuration files, etc.
    return true;
}

/**
 * @brief Shutdown OrbFit library
 * 
 * Call this function at program exit to clean up resources.
 */
inline void shutdown() {
    // Future: Clean up resources, close files, etc.
}

} // namespace orbfit

#endif // ORBFIT_HPP
