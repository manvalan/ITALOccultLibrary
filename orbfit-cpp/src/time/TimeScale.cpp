/**
 * @file TimeScale.cpp
 * @brief Implementation of time scale conversions
 * @author OrbFit C++ Conversion Team
 * @date 2025-11-23
 */

#include <orbfit/time/TimeScale.hpp>
#include <cmath>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <map>

namespace orbfit {
namespace time {

// Static data for leap seconds
static std::map<double, int> leap_second_table;
static bool leap_seconds_loaded = false;

// ============================================================================
// Julian Date Conversions
// ============================================================================

std::tuple<int, int, int, double> mjd_to_calendar(double mjd) {
    double jd = mjd_to_jd(mjd);
    
    // Algorithm from Fliegel and van Flandern (1968)
    int l = static_cast<int>(jd + 68569);
    int n = (4 * l) / 146097;
    l = l - (146097 * n + 3) / 4;
    int i = (4000 * (l + 1)) / 1461001;
    l = l - (1461 * i) / 4 + 31;
    int j = (80 * l) / 2447;
    int day = l - (2447 * j) / 80;
    l = j / 11;
    int month = j + 2 - 12 * l;
    int year = 100 * (n - 49) + i + l;
    
    double fraction = jd - std::floor(jd) + 0.5;
    if (fraction >= 1.0) {
        fraction -= 1.0;
        day += 1;
    }
    
    return std::make_tuple(year, month, day, fraction);
}

double calendar_to_mjd(int year, int month, int day, double fraction) {
    // Algorithm from Fliegel and van Flandern (1968)
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    
    int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    double jd = jdn - 0.5 + fraction;
    
    return jd_to_mjd(jd);
}

// ============================================================================
// Time Scale Conversions
// ============================================================================

double utc_to_tai(double mjd_utc, int leap_seconds) {
    return mjd_utc + leap_seconds / constants::SECONDS_PER_DAY;
}

double tai_to_utc(double mjd_tai, int leap_seconds) {
    return mjd_tai - leap_seconds / constants::SECONDS_PER_DAY;
}

double tai_to_tt(double mjd_tai) {
    return mjd_tai + TT_MINUS_TAI / constants::SECONDS_PER_DAY;
}

double tt_to_tai(double mjd_tt) {
    return mjd_tt - TT_MINUS_TAI / constants::SECONDS_PER_DAY;
}

double utc_to_tt(double mjd_utc) {
    return tai_to_tt(utc_to_tai(mjd_utc));
}

double tt_to_utc(double mjd_tt) {
    return tai_to_utc(tt_to_tai(mjd_tt));
}

double tt_to_tdb(double mjd_tt) {
    // Fairhead & Bretagnon (1990) approximation
    // TDB - TT ≈ 0.001658 sin(g) + 0.000014 sin(2g) seconds
    // where g is the mean anomaly of Earth
    
    double jd_tt = mjd_to_jd(mjd_tt);
    double T = julian_centuries_j2000(jd_tt);
    
    // Mean anomaly of Earth (radians)
    double g = constants::DEG_TO_RAD * (357.528 + 35999.050 * T);
    
    // TDB - TT in seconds
    double dt_seconds = 0.001658 * std::sin(g) + 0.000014 * std::sin(2.0 * g);
    
    return mjd_tt + dt_seconds / constants::SECONDS_PER_DAY;
}

double tdb_to_tt(double mjd_tdb) {
    // Iterative solution (Newton-Raphson)
    double mjd_tt = mjd_tdb;  // Initial guess
    
    for (int iter = 0; iter < 5; ++iter) {
        double mjd_tdb_calc = tt_to_tdb(mjd_tt);
        double error = mjd_tdb - mjd_tdb_calc;
        
        if (std::abs(error) < 1.0e-12) {
            break;
        }
        
        mjd_tt += error;
    }
    
    return mjd_tt;
}

double tai_to_gps(double mjd_tai) {
    return mjd_tai + GPS_MINUS_TAI / constants::SECONDS_PER_DAY;
}

double gps_to_tai(double mjd_gps) {
    return mjd_gps - GPS_MINUS_TAI / constants::SECONDS_PER_DAY;
}

// ============================================================================
// UT1 Conversions
// ============================================================================

double get_dut1(double mjd_utc) {
    // Simplified: return approximate value
    // TODO: Load from IERS finals.data file
    return 0.0;  // Placeholder
}

double utc_to_ut1(double mjd_utc) {
    double dut1 = get_dut1(mjd_utc);
    return mjd_utc + dut1 / constants::SECONDS_PER_DAY;
}

double ut1_to_utc(double mjd_ut1) {
    double dut1 = get_dut1(mjd_ut1);
    return mjd_ut1 - dut1 / constants::SECONDS_PER_DAY;
}

// ============================================================================
// Time Utilities
// ============================================================================

std::optional<double> parse_time_string(
    const std::string& time_string,
    TimeScale time_scale
) {
    // Simple parser for ISO 8601 format: YYYY-MM-DD or YYYY-MM-DDTHH:MM:SS
    std::istringstream iss(time_string);
    int year, month, day, hour = 0, minute = 0;
    double second = 0.0;
    char sep;
    
    iss >> year >> sep >> month >> sep >> day;
    
    if (iss.fail()) {
        return std::nullopt;
    }
    
    // Try to parse time part
    if (iss >> sep >> hour >> sep >> minute >> sep >> second) {
        // Successfully parsed time
    }
    
    double fraction = (hour + minute / 60.0 + second / 3600.0) / 24.0;
    double mjd = calendar_to_mjd(year, month, day, fraction);
    
    // Convert to requested time scale
    switch (time_scale) {
        case TimeScale::UTC:
            return mjd;
        case TimeScale::TT:
            return utc_to_tt(mjd);
        case TimeScale::TDB:
            return utc_to_tdb(mjd);
        default:
            return mjd;
    }
}

std::string format_time(
    double mjd,
    TimeScale time_scale,
    const std::string& format
) {
    auto [year, month, day, fraction] = mjd_to_calendar(mjd);
    
    double hours_frac = fraction * 24.0;
    int hours = static_cast<int>(hours_frac);
    double minutes_frac = (hours_frac - hours) * 60.0;
    int minutes = static_cast<int>(minutes_frac);
    double seconds = (minutes_frac - minutes) * 60.0;
    
    std::ostringstream oss;
    oss << std::setfill('0')
        << std::setw(4) << year << "-"
        << std::setw(2) << month << "-"
        << std::setw(2) << day << " "
        << std::setw(2) << hours << ":"
        << std::setw(2) << minutes << ":"
        << std::fixed << std::setprecision(3) << std::setw(6) << seconds;
    
    return oss.str();
}

double now(TimeScale time_scale) {
    // Get current system time
    auto now = std::chrono::system_clock::now();
    auto duration = now.time_since_epoch();
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();
    
    // Convert Unix time to MJD
    // Unix epoch (1970-01-01 00:00:00) = MJD 40587.0
    double mjd_utc = 40587.0 + seconds / constants::SECONDS_PER_DAY;
    
    // Convert to requested time scale
    switch (time_scale) {
        case TimeScale::UTC:
            return mjd_utc;
        case TimeScale::TT:
            return utc_to_tt(mjd_utc);
        case TimeScale::TDB:
            return utc_to_tdb(mjd_utc);
        case TimeScale::UT1:
            return utc_to_ut1(mjd_utc);
        default:
            return mjd_utc;
    }
}

// ============================================================================
// Leap Second Management
// ============================================================================

bool load_leap_seconds(const std::string& filepath) {
    // TODO: Implement leap second file loading
    // For now, initialize with known values
    leap_second_table.clear();
    
    // Some historical leap seconds (MJD, cumulative leap seconds)
    leap_second_table[41317.0] = 10;  // 1972-01-01
    leap_second_table[57754.0] = 37;  // 2017-01-01
    
    leap_seconds_loaded = true;
    return true;
}

int get_leap_seconds(double mjd_utc) {
    if (!leap_seconds_loaded) {
        return LEAP_SECONDS_2025;  // Default
    }
    
    // Find the leap second count for this date
    int leap_seconds = 10;  // Minimum value
    
    for (const auto& entry : leap_second_table) {
        if (mjd_utc >= entry.first) {
            leap_seconds = entry.second;
        } else {
            break;
        }
    }
    
    return leap_seconds;
}

} // namespace time
} // namespace orbfit
