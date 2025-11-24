/**
 * @file test_types.cpp
 * @brief Unit tests for type definitions and utilities
 */

#include <gtest/gtest.h>
#include <orbfit/core/Types.hpp>
#include <cmath>

using namespace orbfit;

class TypesTest : public ::testing::Test {
protected:
    static constexpr double TOLERANCE = 1.0e-10;
};

// ============================================================================
// Eigen Vector Tests
// ============================================================================

TEST_F(TypesTest, Vector3dCreation) {
    Vector3d v(1.0, 2.0, 3.0);
    EXPECT_DOUBLE_EQ(v(0), 1.0);
    EXPECT_DOUBLE_EQ(v(1), 2.0);
    EXPECT_DOUBLE_EQ(v(2), 3.0);
}

TEST_F(TypesTest, Vector6dCreation) {
    Vector6d state;
    state << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    
    EXPECT_EQ(state.size(), 6);
    EXPECT_DOUBLE_EQ(state(0), 1.0);
    EXPECT_DOUBLE_EQ(state(5), 6.0);
}

TEST_F(TypesTest, Matrix3dCreation) {
    Matrix3d m = Matrix3d::Identity();
    
    EXPECT_DOUBLE_EQ(m(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(m(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(m(2, 2), 1.0);
    EXPECT_DOUBLE_EQ(m(0, 1), 0.0);
}

// ============================================================================
// Special Values Tests
// ============================================================================

TEST_F(TypesTest, NaNDetection) {
    double nan_value = NaN;
    double normal_value = 1.0;
    
    EXPECT_TRUE(isNaN(nan_value));
    EXPECT_FALSE(isNaN(normal_value));
}

TEST_F(TypesTest, InfinityDetection) {
    double inf_value = Infinity;
    double normal_value = 1.0;
    
    EXPECT_FALSE(isFinite(inf_value));
    EXPECT_TRUE(isFinite(normal_value));
}

TEST_F(TypesTest, FiniteDetection) {
    EXPECT_TRUE(isFinite(0.0));
    EXPECT_TRUE(isFinite(1.0));
    EXPECT_TRUE(isFinite(-1.0));
    EXPECT_FALSE(isFinite(Infinity));
    EXPECT_FALSE(isFinite(-Infinity));
    EXPECT_FALSE(isFinite(NaN));
}

// ============================================================================
// Result Type Tests
// ============================================================================

TEST_F(TypesTest, ResultSuccess) {
    auto result = Result<double>::Success(42.0);
    
    EXPECT_TRUE(result.success);
    EXPECT_TRUE(result); // Test implicit bool conversion
    EXPECT_DOUBLE_EQ(result.value, 42.0);
    EXPECT_TRUE(result.error_message.empty());
}

TEST_F(TypesTest, ResultFailure) {
    auto result = Result<double>::Failure("Something went wrong");
    
    EXPECT_FALSE(result.success);
    EXPECT_FALSE(result); // Test implicit bool conversion
    EXPECT_EQ(result.error_message, "Something went wrong");
}

TEST_F(TypesTest, ResultWithComplexType) {
    auto result = Result<Vector3d>::Success(Vector3d(1.0, 2.0, 3.0));
    
    EXPECT_TRUE(result);
    EXPECT_DOUBLE_EQ(result.value(0), 1.0);
    EXPECT_DOUBLE_EQ(result.value(1), 2.0);
    EXPECT_DOUBLE_EQ(result.value(2), 3.0);
}

// ============================================================================
// Enum Tests
// ============================================================================

TEST_F(TypesTest, CoordinateSystemEnum) {
    CoordinateSystem cs = CoordinateSystem::ECLIPTIC_J2000;
    EXPECT_EQ(cs, CoordinateSystem::ECLIPTIC_J2000);
    EXPECT_NE(cs, CoordinateSystem::EQUATORIAL_J2000);
}

TEST_F(TypesTest, ElementTypeEnum) {
    ElementType et = ElementType::KEPLERIAN;
    EXPECT_EQ(et, ElementType::KEPLERIAN);
    EXPECT_NE(et, ElementType::CARTESIAN);
}

TEST_F(TypesTest, TimeScaleEnum) {
    TimeScale ts = TimeScale::TDB;
    EXPECT_EQ(ts, TimeScale::TDB);
    EXPECT_NE(ts, TimeScale::UTC);
}

TEST_F(TypesTest, ObservationTypeEnum) {
    ObservationType ot = ObservationType::OPTICAL_RA_DEC;
    EXPECT_EQ(ot, ObservationType::OPTICAL_RA_DEC);
    EXPECT_NE(ot, ObservationType::RADAR_RANGE);
}

TEST_F(TypesTest, IntegratorTypeEnum) {
    IntegratorType it = IntegratorType::RADAU15;
    EXPECT_EQ(it, IntegratorType::RADAU15);
    EXPECT_NE(it, IntegratorType::RUNGE_KUTTA_GAUSS);
}

TEST_F(TypesTest, ForceComponentEnum) {
    ForceComponent fc = ForceComponent::POINT_MASS;
    EXPECT_EQ(fc, ForceComponent::POINT_MASS);
    EXPECT_NE(fc, ForceComponent::YARKOVSKY);
}

// ============================================================================
// Type Alias Tests
// ============================================================================

TEST_F(TypesTest, TimeTypeAliases) {
    MJD mjd = 51544.5;
    JD jd = 2451545.0;
    TimeSpan span = 365.25;
    
    EXPECT_DOUBLE_EQ(mjd, 51544.5);
    EXPECT_DOUBLE_EQ(jd, 2451545.0);
    EXPECT_DOUBLE_EQ(span, 365.25);
}

TEST_F(TypesTest, AngleTypeAliases) {
    Radians rad = 3.14159;
    Degrees deg = 180.0;
    
    EXPECT_NEAR(rad, 3.14159, TOLERANCE);
    EXPECT_DOUBLE_EQ(deg, 180.0);
}

TEST_F(TypesTest, DistanceTypeAliases) {
    AU_Distance au_dist = 1.0;
    KM_Distance km_dist = 149597870.7;
    
    EXPECT_DOUBLE_EQ(au_dist, 1.0);
    EXPECT_DOUBLE_EQ(km_dist, 149597870.7);
}

TEST_F(TypesTest, VelocityTypeAliases) {
    AU_Per_Day au_vel = 0.01;
    KM_Per_Second km_vel = 30.0;
    
    EXPECT_DOUBLE_EQ(au_vel, 0.01);
    EXPECT_DOUBLE_EQ(km_vel, 30.0);
}

TEST_F(TypesTest, StringTypeAliases) {
    ObservatoryCode obs_code = "500";
    ObjectName obj_name = "2024 AA";
    
    EXPECT_EQ(obs_code, "500");
    EXPECT_EQ(obj_name, "2024 AA");
}

// Run all tests
// main() is provided by gtest_main
