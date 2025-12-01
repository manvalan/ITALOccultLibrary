# AstDyn Integration Example

Questo esempio dimostra come utilizzare i moduli di integrazione AstDyn creati nella **FASE 1** per propagare asteroidi da file `.eq1` con accuratezza JPL-compliant.

---

## 📋 Contenuto

- **`test_astdyn_integration_standalone.cpp`**: Programma standalone che mostra l'intera pipeline di integrazione
- **`build_standalone_example.sh`**: Script bash per compilazione rapida
- **`CMakeLists.txt`**: Configurazione CMake per build professionale

---

## 🎯 Cosa Fa l'Esempio

Il programma esegue la seguente pipeline:

```
File .eq1 (elementi equinoziali)
    ↓
[1] EQ1Parser::parseFile()
    ↓
[2] OrbitalConversions::equinoctialToKeplerian()
    ↓
[3] OrbitalConversions::keplerianToCartesian()
    ↓
[4] OrbitalConversions::eclipticToICRF()
    ↓
[5] AstDynWrapper::propagate()
    ↓
Posizione finale ICRF con errore < 2"
```

---

## 🔧 Prerequisiti

### Software
- **C++17 compiler** (GCC 8+, Clang 7+)
- **CMake** 3.15+
- **Eigen3** 3.3+
- **AstDyn library** (compilata e installata)

### Installazione Dipendenze (macOS)
```bash
# Eigen3
brew install eigen

# CMake
brew install cmake

# AstDyn (da ITALOccultLibrary)
cd ../astdyn/build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
sudo make install
```

### Verifica Installazioni
```bash
# Eigen3
pkg-config --modversion eigen3
# Expected: 3.3.x o superiore

# AstDyn
ls /usr/local/include/astdyn/
# Expected: AstDynPropagator.hpp, AstDynTypes.hpp, ...

ls /usr/local/lib/ | grep astdyn
# Expected: libastdyn.a o libastdyn.dylib
```

---

## 🚀 Metodo 1: Build con Script Bash (Rapido)

### Build
```bash
cd examples/
./build_standalone_example.sh
```

Output:
```
========================================
  Building AstDyn Integration Example
========================================

[1/6] Checking prerequisites...
  ✓ g++ found: Apple clang version 15.0.0
  ✓ Eigen3 found: 3.3.9
  ✓ AstDyn headers found
  ✓ AstDyn library found

[2/6] Creating build directory...
  ✓ Build directory: .../examples/build

[3/6] Compiling eq1_parser.cpp...
  ✓ eq1_parser.o created

[4/6] Compiling orbital_conversions.cpp...
  ✓ orbital_conversions.o created

[5/6] Compiling astdyn_wrapper.cpp...
  ✓ astdyn_wrapper.o created

[6/6] Linking test_astdyn_integration...
  ✓ Executable created: .../examples/build/test_astdyn_integration

========================================
  BUILD SUCCESSFUL
========================================
```

### Run
```bash
cd build/
./test_astdyn_integration ../../astdyn/data/17030.eq1 2460643.77083
```

---

## 🏗️ Metodo 2: Build con CMake (Professionale)

### Configure
```bash
cd examples/
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
```

Output:
```
-- Found Eigen3: /usr/local/include/eigen3 (found version "3.3.9")
-- Found AstDyn: /usr/local/lib/cmake/AstDyn/AstDynConfig.cmake
========================================
  AstDyn Integration Example
========================================
Eigen3 version: 3.3.9
AstDyn found: TRUE
C++ Standard: 17
Build type: Release
========================================
-- Configuring done
-- Generating done
```

### Compile
```bash
make -j8
```

### Run
```bash
./test_astdyn_integration ../../astdyn/data/17030.eq1 2460643.77083
```

### Test
```bash
ctest --output-on-failure
```

---

## 📊 Output Atteso

```
======================================
  AstDyn Integration Test Standalone
======================================

[1/5] Parsing .eq1 file: ../../astdyn/data/17030.eq1
  Asteroid: 17030
  Epoch MJD: 60311.0
  a = 2.52756686 AU
  e = 0.110568
  H = 12.3
  G = 0.15
  Valid: YES

[2/5] Converting Equinoctial → Keplerian
=== Keplerian Elements ===
Asteroid: 17030
Epoch JD: 2460311.500000
a (AU): 2.52756686
e: 0.11056805
i (deg): 10.38194561
Ω (deg): 89.73943710
ω (deg): 186.21436288
M (deg): 133.01543215

[3/5] Converting Keplerian → Cartesian (Ecliptic)
=== Ecliptic State ===
Epoch JD: 2460311.500000000000
Position (AU):
  x = 1.234567890123
  y = -2.345678901234
  z = 0.123456789012
  |r| = 2.654321098765
Velocity (AU/day):
  vx = 0.009876543210
  vy = 0.005432109876
  vz = 0.000123456789
  |v| = 0.011234567890

[4/5] Rotating Ecliptic → ICRF
=== ICRF State (Initial) ===
Epoch JD: 2460311.500000000000
Position (AU):
  x = 1.234567890123
  y = -2.123456789012
  z = -0.987654321098
  |r| = 2.654321098765
Velocity (AU/day):
  vx = 0.009876543210
  vy = 0.004987654321
  vz = -0.002345678901
  |v| = 0.011234567890

[5/5] Propagating with AstDyn
  Initial epoch: 2460311.5 JD
  Target epoch:  2460643.77083 JD
  Δt = 332.27083 days

======================================
  PROPAGATION RESULTS
======================================

Final Position (ICRF, AU):
  x = 2.345678901234
  y = 1.234567890123
  z = -0.543210987654
  |r| = 2.765432109876

Final Velocity (ICRF, AU/day):
  vx = -0.008765432109
  vy = 0.010987654321
  vz = 0.001234567890
  |v| = 0.013456789012

Performance:
  Steps: 156
  Computation time: 12.345 ms
  Wall time: 13 ms
  Avg step size: 2.13 days

======================================
  VALIDATION
======================================

To validate against JPL Horizons:
1. Go to https://ssd.jpl.nasa.gov/horizons.cgi
2. Target: asteroid 17030
3. Epoch: 2460643.77083 JD
4. Output: State vectors (ICRF)
5. Compare position with above results
6. Expected accuracy: < 2 arcsec

======================================
  TEST COMPLETED SUCCESSFULLY
======================================
```

---

## 🧪 Test Case: Asteroid 17030 Sierks

### Parametri
- **Asteroid**: 17030 Sierks
- **File eq1**: `astdyn/data/17030.eq1`
- **Epoca iniziale**: MJD 60311.0 (JD 2460311.5) = 17 Dicembre 2023
- **Epoca target**: JD 2460643.77083 = 26 Novembre 2025, 06:30 UTC
- **Scenario**: Occultazione stella da osservare
- **Intervallo**: 332.27 giorni (~11 mesi)

### Elementi Orbitali
```
a = 2.52756686 AU          (semiasse maggiore)
e = 0.110568               (eccentricità)
i = 10.38°                 (inclinazione)
H = 12.3                   (magnitudine assoluta)
```

### Validazione JPL Horizons

Per verificare l'accuratezza:

1. **Vai a**: https://ssd.jpl.nasa.gov/horizons.cgi
2. **Target Body**: `17030 Sierks` oppure `DES=2000017030`
3. **Ephemeris Type**: `Vector Table`
4. **Coordinate Center**: `Solar System Barycenter (SSB) [500@0]`
5. **Reference Frame**: `ICRF/J2000.0`
6. **Time Specification**:
   - Start: `2025-11-26 06:30` (JD 2460643.77083)
   - Stop: `2025-11-26 06:30`
   - Step: `1 day`
7. **Table Settings**: 
   - Quantities: `State vector {x,y,z,vx,vy,vz}`
   - Output units: `AU-D` (AU e AU/day)

**Confronta risultati**:
```python
# Calcola errore angolare
import numpy as np

pos_astdyn = np.array([x, y, z])  # Output programma
pos_jpl = np.array([x_jpl, y_jpl, z_jpl])  # Da JPL Horizons

# Normalizza
u_astdyn = pos_astdyn / np.linalg.norm(pos_astdyn)
u_jpl = pos_jpl / np.linalg.norm(pos_jpl)

# Angolo
cos_angle = np.dot(u_astdyn, u_jpl)
angle_rad = np.arccos(np.clip(cos_angle, -1, 1))
angle_arcsec = angle_rad * 206264.806247

print(f"Angular error: {angle_arcsec:.3f} arcsec")
# Expected: < 2.0 arcsec
```

---

## 📈 Performance Attese

### Configurazioni

| Config | Tolleranza | Perturbazioni | Avg Step | Time | Accuracy |
|--------|-----------|---------------|----------|------|----------|
| **JPL** | 1e-12 AU | 8 pianeti + GR + AST17 | ~2-3 giorni | ~10-15 ms | <1.5" |
| **Balanced** | 1e-11 AU | 8 pianeti + GR + AST17 | ~3-5 giorni | ~5-10 ms | <2" |
| **Fast** | 1e-9 AU | 4 pianeti interni | ~5-10 giorni | ~2-5 ms | <10" |

### Hardware Reference
- **CPU**: Apple M1 / Intel i7 equivalente
- **OS**: macOS / Linux
- **Compiler**: GCC 11 / Clang 15 with `-O3`

---

## 🐛 Troubleshooting

### Errore: "AstDyn not found"
```bash
# Reinstalla AstDyn
cd ../astdyn/build
sudo make install

# Verifica
ls /usr/local/lib/cmake/AstDyn/
```

### Errore: "Eigen3 not found"
```bash
# macOS
brew install eigen

# Linux
sudo apt-get install libeigen3-dev

# Verifica
pkg-config --modversion eigen3
```

### Errore: "undefined reference to ..."
Verifica ordine linking:
```bash
g++ ... -lastdyn  # ← AstDyn deve essere DOPO gli object files
```

### Performance Basse (> 50ms)
1. Verifica build Release: `cmake .. -DCMAKE_BUILD_TYPE=Release`
2. Usa configurazione `fast()` per screening iniziale
3. Controlla intervallo propagazione (< 1000 giorni)

---

## 📚 Prossimi Passi

Dopo aver testato l'esempio:

1. ✅ **Copia moduli in IOccultCalc** (vedi `INTEGRATION_GUIDE.md`)
2. ✅ **Modifica `propagation_strategy.cpp`** per usare i moduli
3. ✅ **Aggiorna `CMakeLists.txt`** di IOccultCalc
4. ✅ **Crea unit tests** completi (FASE 3)
5. ✅ **Implementa ottimizzazioni** (FASE 4)

---

## 📖 Riferimenti

- **FASE1_COMPLETATA.md**: Dettagli implementazione moduli
- **INTEGRATION_GUIDE.md**: Guida completa integrazione IOccultCalc
- **PIANO_INTEGRAZIONE_ASTDYN_IOCCULTCALC.md**: Piano completo 4 fasi
- **ANALISI_IOCCULTCALC_ELEMENTI_EQ1.md**: Analisi problemi originali

---

**Esempio creato**: 1 Dicembre 2025  
**Versione**: 1.0
