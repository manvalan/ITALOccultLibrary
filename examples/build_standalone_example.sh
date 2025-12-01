#!/bin/bash
# build_standalone_example.sh
# Script per compilare l'esempio standalone di integrazione AstDyn

set -e  # Exit on error

echo "========================================"
echo "  Building AstDyn Integration Example"
echo "========================================"
echo ""

# Directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
TEMPLATES_DIR="$ROOT_DIR/templates_ioccultcalc"
EXAMPLES_DIR="$ROOT_DIR/examples"
BUILD_DIR="$EXAMPLES_DIR/build"

# Check prerequisites
echo "[1/6] Checking prerequisites..."

# Check compiler
if ! command -v g++ &> /dev/null; then
    echo "ERROR: g++ not found. Please install C++ compiler."
    exit 1
fi
echo "  ✓ g++ found: $(g++ --version | head -n1)"

# Check Eigen3
if ! pkg-config --exists eigen3; then
    echo "ERROR: Eigen3 not found. Install with: brew install eigen"
    exit 1
fi
EIGEN_INCLUDE=$(pkg-config --cflags eigen3)
echo "  ✓ Eigen3 found: $(pkg-config --modversion eigen3)"

# Check AstDyn
ASTDYN_INCLUDE="/usr/local/include"
ASTDYN_LIB="/usr/local/lib"
if [ ! -f "$ASTDYN_INCLUDE/astdyn/AstDynPropagator.hpp" ]; then
    echo "ERROR: AstDyn headers not found in $ASTDYN_INCLUDE/astdyn/"
    echo "Please install AstDyn library first."
    exit 1
fi
echo "  ✓ AstDyn headers found"

if [ ! -f "$ASTDYN_LIB/libastdyn.a" ] && [ ! -f "$ASTDYN_LIB/libastdyn.dylib" ]; then
    echo "ERROR: AstDyn library not found in $ASTDYN_LIB/"
    echo "Please install AstDyn library first."
    exit 1
fi
echo "  ✓ AstDyn library found"

# Create build directory
echo ""
echo "[2/6] Creating build directory..."
mkdir -p "$BUILD_DIR"
echo "  ✓ Build directory: $BUILD_DIR"

# Compile eq1_parser.cpp
echo ""
echo "[3/6] Compiling eq1_parser.cpp..."
g++ -std=c++17 -c \
    -I"$TEMPLATES_DIR/include" \
    $EIGEN_INCLUDE \
    -O3 -Wall -Wextra \
    "$TEMPLATES_DIR/src/eq1_parser.cpp" \
    -o "$BUILD_DIR/eq1_parser.o"
echo "  ✓ eq1_parser.o created"

# Compile orbital_conversions.cpp
echo ""
echo "[4/6] Compiling orbital_conversions.cpp..."
g++ -std=c++17 -c \
    -I"$TEMPLATES_DIR/include" \
    $EIGEN_INCLUDE \
    -O3 -Wall -Wextra \
    "$TEMPLATES_DIR/src/orbital_conversions.cpp" \
    -o "$BUILD_DIR/orbital_conversions.o"
echo "  ✓ orbital_conversions.o created"

# Compile astdyn_wrapper.cpp
echo ""
echo "[5/6] Compiling astdyn_wrapper.cpp..."
g++ -std=c++17 -c \
    -I"$TEMPLATES_DIR/include" \
    -I"$ASTDYN_INCLUDE" \
    $EIGEN_INCLUDE \
    -O3 -Wall -Wextra \
    "$TEMPLATES_DIR/src/astdyn_wrapper.cpp" \
    -o "$BUILD_DIR/astdyn_wrapper.o"
echo "  ✓ astdyn_wrapper.o created"

# Link executable
echo ""
echo "[6/6] Linking test_astdyn_integration..."
g++ -std=c++17 \
    -I"$TEMPLATES_DIR/include" \
    -I"$ASTDYN_INCLUDE" \
    $EIGEN_INCLUDE \
    -O3 -Wall -Wextra \
    "$EXAMPLES_DIR/test_astdyn_integration_standalone.cpp" \
    "$BUILD_DIR/eq1_parser.o" \
    "$BUILD_DIR/orbital_conversions.o" \
    "$BUILD_DIR/astdyn_wrapper.o" \
    -L"$ASTDYN_LIB" \
    -lastdyn \
    -o "$BUILD_DIR/test_astdyn_integration"

echo "  ✓ Executable created: $BUILD_DIR/test_astdyn_integration"

# Summary
echo ""
echo "========================================"
echo "  BUILD SUCCESSFUL"
echo "========================================"
echo ""
echo "To run the test:"
echo "  cd $BUILD_DIR"
echo "  ./test_astdyn_integration ../../astdyn/data/17030.eq1 2460643.77083"
echo ""
echo "Expected output:"
echo "  - Parsed eq1 elements"
echo "  - Keplerian conversion"
echo "  - Cartesian states (Ecliptic + ICRF)"
echo "  - Propagation results"
echo "  - Performance metrics"
echo ""
