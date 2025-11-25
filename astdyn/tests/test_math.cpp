/**
 * @file test_math.cpp
 * @brief Unit tests for MathUtils and LinearAlgebra modules
 */

#include <gtest/gtest.h>
#include "astdyn/math/MathUtils.hpp"
#include "astdyn/math/LinearAlgebra.hpp"

using namespace astdyn::math;

// ========== MathUtils Tests ==========

TEST(MathUtilsTest, CholeskyDecomposition) {
    Eigen::Matrix3d A;
    A << 4, 12, -16,
         12, 37, -43,
        -16, -43, 98;
    
    auto result = cholesky_decompose(A);
    ASSERT_TRUE(result.has_value());
    
    Eigen::Matrix3d L = result.value();
    Eigen::Matrix3d reconstructed = L * L.transpose();
    
    EXPECT_TRUE(A.isApprox(reconstructed, 1e-10));
}

TEST(MathUtilsTest, CholeskyInversion) {
    Eigen::Matrix3d A;
    A << 4, 2, 1,
         2, 5, 3,
         1, 3, 6;
    
    auto result = cholesky_invert(A);
    ASSERT_TRUE(result.has_value());
    
    Eigen::Matrix3d inv = result.value();
    Eigen::Matrix3d identity = A * inv;
    
    EXPECT_TRUE(identity.isApprox(Eigen::Matrix3d::Identity(), 1e-10));
}

TEST(MathUtilsTest, MatrixInversion) {
    Eigen::Matrix3d A;
    A << 1, 2, 3,
         0, 1, 4,
         5, 6, 0;
    
    auto result = matrix_invert(A);
    ASSERT_TRUE(result.has_value());
    
    Eigen::Matrix3d inv = result.value();
    Eigen::Matrix3d identity = A * inv;
    
    EXPECT_TRUE(identity.isApprox(Eigen::Matrix3d::Identity(), 1e-10));
}

TEST(MathUtilsTest, WeightedRMSNorm) {
    Eigen::VectorXd v(3);
    v << 1.0, 2.0, 3.0;
    
    Eigen::MatrixXd weight_matrix = Eigen::MatrixXd::Identity(3, 3);
    
    double norm = weighted_rms_norm(v, weight_matrix);
    double expected = std::sqrt((1.0 + 4.0 + 9.0) / 3.0);
    
    EXPECT_NEAR(norm, expected, 1e-10);
}

TEST(MathUtilsTest, BilinearForm) {
    Eigen::Vector3d x(1.0, 2.0, 3.0);
    Eigen::Vector3d y(4.0, 5.0, 6.0);
    Eigen::Matrix3d A = Eigen::Matrix3d::Identity();
    
    double result = bilinear_form(x, A, y);
    double expected = x.dot(y); // Since A is identity
    
    EXPECT_NEAR(result, expected, 1e-10);
}

TEST(MathUtilsTest, IsSymmetric) {
    Eigen::Matrix3d symmetric;
    symmetric << 1, 2, 3,
                 2, 4, 5,
                 3, 5, 6;
    
    Eigen::Matrix3d asymmetric;
    asymmetric << 1, 2, 3,
                  4, 5, 6,
                  7, 8, 9;
    
    EXPECT_TRUE(is_symmetric(symmetric));
    EXPECT_FALSE(is_symmetric(asymmetric));
}

TEST(MathUtilsTest, IsPositiveDefinite) {
    Eigen::Matrix3d pd;
    pd << 4, 2, 1,
          2, 5, 3,
          1, 3, 6;
    
    Eigen::Matrix3d npd;
    npd << 1, 2, 3,
           2, 1, 4,
           3, 4, 1;
    
    EXPECT_TRUE(is_positive_definite(pd));
    EXPECT_FALSE(is_positive_definite(npd));
}

TEST(MathUtilsTest, StatisticalFunctions) {
    Eigen::VectorXd data(5);
    data << 1.0, 2.0, 3.0, 4.0, 5.0;
    
    EXPECT_NEAR(data.mean(), 3.0, 1e-10);
    EXPECT_NEAR(standard_deviation(data), std::sqrt(2.5), 1e-10);
    EXPECT_NEAR(median(data), 3.0, 1e-10);
}

TEST(MathUtilsTest, NormalizeAngle) {
    // Test with positive normalization [0, 2π)
    EXPECT_NEAR(normalize_angle_positive(3.0 * M_PI), M_PI, 1e-10);
    EXPECT_NEAR(normalize_angle_positive(-M_PI), M_PI, 1e-10);
    EXPECT_NEAR(normalize_angle_positive(0.5 * M_PI), 0.5 * M_PI, 1e-10);
}

// ========== LinearAlgebra Tests ==========

TEST(LinearAlgebraTest, QRDecomposition) {
    Eigen::MatrixXd A(3, 3);
    A << 12, -51, 4,
          6, 167, -68,
         -4, 24, -41;
    
    auto qr = qr_decompose(A);
    
    Eigen::MatrixXd reconstructed = qr.Q * qr.R;
    EXPECT_TRUE(A.isApprox(reconstructed, 1e-10));
    
    // Q should be orthogonal
    Eigen::MatrixXd QtQ = qr.Q.transpose() * qr.Q;
    EXPECT_TRUE(QtQ.isApprox(Eigen::MatrixXd::Identity(3, 3), 1e-10));
}

TEST(LinearAlgebraTest, SVDDecomposition) {
    Eigen::MatrixXd A(3, 2);
    A << 1, 2,
         3, 4,
         5, 6;
    
    auto svd = svd_decompose(A);
    
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(3, 2);
    S(0, 0) = svd.singular_values(0);
    S(1, 1) = svd.singular_values(1);
    
    Eigen::MatrixXd reconstructed = svd.U * S * svd.V.transpose();
    EXPECT_TRUE(A.isApprox(reconstructed, 1e-10));
}

TEST(LinearAlgebraTest, Pseudoinverse) {
    Eigen::MatrixXd A(3, 2);
    A << 1, 2,
         3, 4,
         5, 6;
    
    Eigen::MatrixXd pinv = pseudoinverse(A);
    
    // A * pinv * A should equal A
    Eigen::MatrixXd reconstructed = A * pinv * A;
    EXPECT_TRUE(A.isApprox(reconstructed, 1e-8));
}

TEST(LinearAlgebraTest, EigenvalueSymmetric) {
    Eigen::MatrixXd A(3, 3);
    A << 4, 1, 0,
         1, 4, 1,
         0, 1, 4;
    
    auto eigen = eigen_symmetric(A);
    
    // Check: A * v = lambda * v for first eigenpair
    Eigen::VectorXd v0 = eigen.eigenvectors.col(0);
    double lambda0 = eigen.eigenvalues(0);
    
    Eigen::VectorXd lhs = A * v0;
    Eigen::VectorXd rhs = lambda0 * v0;
    
    EXPECT_TRUE(lhs.isApprox(rhs, 1e-10));
}

TEST(LinearAlgebraTest, LeastSquares) {
    // y = 2x + 1
    Eigen::MatrixXd A(3, 2);
    A << 1, 1,
         1, 2,
         1, 3;
    
    Eigen::VectorXd b(3);
    b << 3, 5, 7;
    
    auto result = least_squares(A, b);
    ASSERT_TRUE(result.has_value());
    
    Eigen::VectorXd x = result.value();
    EXPECT_NEAR(x(0), 1.0, 1e-10); // intercept
    EXPECT_NEAR(x(1), 2.0, 1e-10); // slope
}

TEST(LinearAlgebraTest, LinearRegression) {
    Eigen::VectorXd x(4);
    x << 1, 2, 3, 4;
    
    Eigen::VectorXd y(4);
    y << 2, 4, 6, 8; // y = 2x
    
    auto result = linear_regression(x, y);
    ASSERT_TRUE(result.has_value());
    
    double slope = result.value().first;
    double intercept = result.value().second;
    
    EXPECT_NEAR(slope, 2.0, 1e-10);
    EXPECT_NEAR(intercept, 0.0, 1e-10);
}

TEST(LinearAlgebraTest, MatrixRank) {
    Eigen::MatrixXd full_rank(3, 3);
    full_rank << 1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;
    
    Eigen::MatrixXd rank_deficient(3, 3);
    rank_deficient << 1, 2, 3,
                      2, 4, 6,
                      3, 6, 9;
    
    EXPECT_EQ(matrix_rank(full_rank), 3);
    EXPECT_EQ(matrix_rank(rank_deficient), 1);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
