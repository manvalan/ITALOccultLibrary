# PARTE 1 - Capitolo 2: Piani di Riferimento e Trasformazioni
# PART 1 - Chapter 2: Reference Planes and Transformations

---

## 🇮🇹 Italiano

### 2.1 Sistemi di Riferimento Inerziali

Un **sistema di riferimento inerziale** è un sistema di coordinate che non è soggetto ad accelerazioni. In astrodinamica, utilizziamo diversi sistemi di riferimento quasi-inerziali, poiché un sistema perfettamente inerziale è un'idealizzazione teorica.

#### 2.1.1 ICRS - International Celestial Reference System

L'**ICRS** (International Celestial Reference System) è il sistema di riferimento celeste standard definito dalla IAU. Le sue caratteristiche sono:

- **Origine**: Baricentro del Sistema Solare
- **Orientamento**: Definito da osservazioni VLBI di quasar extragalattici
- **Stabilità**: Non è soggetto a precessione o nutazione
- **Epoca di riferimento**: Nessuna (sistema non dinamico)

L'ICRS è materializzato dall'**ICRF** (International Celestial Reference Frame), un catalogo di sorgenti radio extragalattiche con posizioni note.

#### 2.1.2 J2000 - Sistema Equatoriale Medio

Il sistema **J2000** è definito come:

- **Origine**: Baricentro del Sistema Solare
- **Piano XY**: Piano equatoriale medio all'epoca J2000.0
- **Asse X**: Direzione del punto vernale medio (equinozio) al J2000.0
- **Asse Z**: Perpendicolare all'equatore medio, verso il polo nord celeste
- **Asse Y**: Completa la terna destrorsa

#### 2.1.3 Relazione ICRS ↔ J2000

La differenza tra ICRS e J2000 è molto piccola (circa 17.3 milliarcsec in orientamento). La trasformazione è data da un **frame bias**:

**Matrice di trasformazione J2000 → ICRS:**

```
       ⎡  0.9999999999999928  -0.0000000707827974  -0.0000000805621715 ⎤
B =    ⎢  0.0000000707827948   0.9999999999999969  -0.0000000330604145 ⎥
       ⎣  0.0000000805621738   0.0000000330604088   0.9999999999999962 ⎦
```

Per la maggior parte delle applicazioni, ICRS e J2000 possono essere considerati equivalenti.

### 2.2 Matrici di Rotazione Fondamentali

Le trasformazioni tra sistemi di coordinate si basano su **matrici di rotazione** attorno agli assi principali.

#### 2.2.1 Rotazione attorno all'asse X

```
         ⎡ 1     0         0     ⎤
Rx(θ) =  ⎢ 0   cos(θ)    sin(θ)  ⎥
         ⎣ 0  -sin(θ)    cos(θ)  ⎦
```

#### 2.2.2 Rotazione attorno all'asse Y

```
         ⎡ cos(θ)   0  -sin(θ) ⎤
Ry(θ) =  ⎢   0      1     0    ⎥
         ⎣ sin(θ)   0   cos(θ) ⎦
```

#### 2.2.3 Rotazione attorno all'asse Z

```
         ⎡  cos(θ)   sin(θ)   0 ⎤
Rz(θ) =  ⎢ -sin(θ)   cos(θ)   0 ⎥
         ⎣    0        0      1 ⎦
```

#### 2.2.4 Proprietà delle Matrici di Rotazione

Le matrici di rotazione sono **matrici ortogonali**:

- Determinante = 1 (rotazione propria)
- R⁻¹ = Rᵀ (l'inversa è uguale alla trasposta)
- R(-θ) = Rᵀ(θ)

### 2.3 Trasformazioni tra Sistemi di Riferimento

#### 2.3.1 Equatoriale J2000 ↔ Eclittica J2000

La trasformazione utilizza l'obliquità dell'eclittica ε₀ = 23.439291°:

**Implementazione in C++:**

```cpp
namespace orbfit::coordinates {
    // Matrice da J2000 a Eclittica
    static Matrix3d j2000_to_ecliptic() {
        constexpr double epsilon0 = 23.439291 * constants::DEG_TO_RAD;
        return rotation_x(epsilon0);
    }
    
    // Matrice da Eclittica a J2000
    static Matrix3d ecliptic_to_j2000() {
        constexpr double epsilon0 = 23.439291 * constants::DEG_TO_RAD;
        return rotation_x(-epsilon0);
    }
}
```

**Esempio numerico:**

Per ε₀ = 23.439291°:

```
cos(ε₀) = 0.9174821
sin(ε₀) = 0.3977772

           ⎡ 1.0000000   0.0000000   0.0000000 ⎤
Rx(ε₀) =   ⎢ 0.0000000   0.9174821   0.3977772 ⎥
           ⎣ 0.0000000  -0.3977772   0.9174821 ⎦
```

#### 2.3.2 Sistema Inerziale ↔ Sistema Terrestre (ITRF)

La trasformazione da sistema inerziale (J2000) a sistema terrestre (ITRF) richiede:

1. **Precessione**: Movimento secolare dell'asse di rotazione terrestre
2. **Nutazione**: Oscillazioni periodiche dell'asse terrestre
3. **Earth Rotation Angle (ERA)**: Rotazione giornaliera della Terra
4. **Polar Motion**: Movimento del polo rispetto alla crosta terrestre

**Trasformazione semplificata (solo ERA):**

```cpp
static Matrix3d j2000_to_itrf_simple(double mjd_ut1) {
    // Earth Rotation Angle
    double T = mjd_ut1 - constants::MJD2000;
    double era = 2.0 * constants::PI * 
                 (0.7790572732640 + 1.00273781191135448 * T);
    
    // Normalizzazione a [0, 2π)
    era = std::fmod(era, 2.0 * constants::PI);
    if (era < 0.0) era += 2.0 * constants::PI;
    
    // Rotazione attorno a Z
    return rotation_z(-era);
}
```

#### 2.3.3 Greenwich Mean Sidereal Time (GMST)

Il **GMST** rappresenta l'angolo orario del punto vernale medio al meridiano di Greenwich.

**Formula (IAU 2000):**

```
GMST = 24110.54841 + 8640184.812866·T + 0.093104·T² - 6.2×10⁻⁶·T³
```

dove T è in secoli giuliani dal J2000.0.

**Implementazione:**

```cpp
static double gmst(double mjd_ut1) {
    double T = (mjd_ut1 - constants::MJD2000) / 36525.0;
    
    // GMST a 0h UT1 (in secondi)
    double gmst0 = 24110.54841 + 8640184.812866 * T 
                 + 0.093104 * T * T - 6.2e-6 * T * T * T;
    
    // Converti in giorni e aggiungi la frazione del giorno
    gmst0 /= constants::SECONDS_PER_DAY;
    double frac_day = mjd_ut1 - std::floor(mjd_ut1);
    double gmst_days = gmst0 + frac_day * 1.00273790935;
    
    // Converti in radianti
    return 2.0 * constants::PI * (gmst_days - std::floor(gmst_days));
}
```

### 2.4 Trasformazione degli Stati Orbitali

Quando si trasformano vettori di **stato orbitale** (posizione + velocità) tra sistemi rotanti, è necessario considerare il termine di **Coriolis**.

#### 2.4.1 Trasformazione della Velocità

Per un sistema rotante con velocità angolare ω:

```
v_rot = R · v_inert - ω × (R · r_inert)
```

**Implementazione:**

```cpp
static Vector3d transform_velocity(
    const Vector3d& pos, const Vector3d& vel,
    FrameType from, FrameType to,
    double mjd_ut1
) {
    Matrix3d R = get_transformation(from, to, mjd_ut1);
    Vector3d vel_rotated = R * vel;
    
    bool from_rotating = (from == FrameType::ITRF);
    bool to_rotating = (to == FrameType::ITRF);
    
    if (from_rotating != to_rotating) {
        // Velocità angolare terrestre [rad/s]
        constexpr double omega_earth = 7.292115e-5;
        Vector3d omega(0.0, 0.0, omega_earth);
        
        Vector3d pos_rotated = R * pos;
        
        if (to_rotating) {
            vel_rotated -= omega.cross(pos_rotated);
        } else {
            vel_rotated += omega.cross(pos_rotated);
        }
    }
    
    return vel_rotated;
}
```

### 2.5 Precessione e Nutazione

#### 2.5.1 Precessione

La **precessione** è il lento movimento conico dell'asse di rotazione terrestre causato dalle attrazioni gravitazionali di Sole e Luna sul rigonfiamento equatoriale terrestre.

**Periodo di precessione**: ~25,772 anni (un ciclo completo)

**Variazione dell'obliquità**: ~0.47" per anno

Il modello di precessione IAU 2006 fornisce:

```
ε(T) = ε₀ - 46.836769"·T - 0.0001831"·T² + 0.00200340"·T³ - ...
```

dove T è in secoli giuliani dal J2000.0.

#### 2.5.2 Nutazione

La **nutazione** rappresenta le oscillazioni periodiche dell'asse terrestre sovrapposte al moto di precessione.

**Componenti principali:**
- Periodo di 18.6 anni (moto retrogrado dei nodi lunari)
- Ampiezza: ~9.2" in nutazione in longitudine
- Ampiezza: ~6.9" in nutazione in obliquità

### 2.6 Esempio Pratico con Dati Reali

#### Esempio 2.1: Trasformazione di coordinate di un satellite

**Dati:** Satellite in orbita geostazionaria
- Posizione (ITRF, MJD 60000.5): r = [42164.0, 0.0, 0.0] km
- Velocità (ITRF): v = [0.0, 3.0747, 0.0] km/s

**Obiettivo:** Trasformare in coordinate J2000

**Calcolo del GMST:**
```
MJD = 60000.5
T = (60000.5 - 51544.5) / 36525.0 = 0.2315
GMST ≈ 5.832 rad = 334.12°
```

**Matrice di rotazione (ITRF → J2000):**
```
           ⎡  0.4087  -0.9127   0 ⎤
R(GMST) =  ⎢  0.9127   0.4087   0 ⎥
           ⎣    0        0      1 ⎦
```

**Posizione in J2000:**
```
r_J2000 = R · r_ITRF
        = [0.4087×42164, 0.9127×42164, 0]
        = [17234.5, 38489.0, 0.0] km
```

**Velocità in J2000 (con Coriolis):**
```
ω = [0, 0, 7.292×10⁻⁵] rad/s
v_J2000 = R · v_ITRF + ω × r_J2000
        ≈ [-0.0257, 1.257, 0.0] + [2.807, -1.257, 0.0]
        = [2.781, 0.0, 0.0] km/s
```

### 2.7 Schema Riassuntivo delle Trasformazioni

```
                    ┌──────────────────┐
                    │      ICRS        │
                    │  (Baricentrico)  │
                    └────────┬─────────┘
                             │ Frame Bias (17 mas)
                             ↓
                    ┌──────────────────┐
                    │   J2000 (EME)    │
                    │  Equatoriale     │
                    └────────┬─────────┘
                             │
              ┌──────────────┼──────────────┐
              │              │              │
              ↓              ↓              ↓
     ┌────────────┐   ┌────────────┐  ┌────────────┐
     │  Eclittica │   │   TEME     │  │    MOD     │
     │   J2000    │   │  (True)    │  │  (Mean)    │
     └────────────┘   └──────┬─────┘  └────────────┘
                             │
                             ↓ ERA
                    ┌──────────────────┐
                    │      ITRF        │
                    │  (Terrestre)     │
                    └──────────────────┘
```

**Legenda:**
- ICRS: International Celestial Reference System
- EME: Earth Mean Equator
- MOD: Mean of Date
- TEME: True Equator Mean Equinox
- ITRF: International Terrestrial Reference Frame
- ERA: Earth Rotation Angle

---

## 🇬🇧 English

### 2.1 Inertial Reference Systems

An **inertial reference system** is a coordinate system that is not subject to accelerations. In astrodynamics, we use various quasi-inertial reference systems, since a perfectly inertial system is a theoretical idealization.

#### 2.1.1 ICRS - International Celestial Reference System

The **ICRS** (International Celestial Reference System) is the standard celestial reference system defined by the IAU. Its characteristics are:

- **Origin**: Barycenter of the Solar System
- **Orientation**: Defined by VLBI observations of extragalactic quasars
- **Stability**: Not subject to precession or nutation
- **Reference epoch**: None (non-dynamic system)

The ICRS is materialized by the **ICRF** (International Celestial Reference Frame), a catalog of extragalactic radio sources with known positions.

#### 2.1.2 J2000 - Mean Equatorial System

The **J2000** system is defined as:

- **Origin**: Barycenter of the Solar System
- **XY Plane**: Mean equatorial plane at epoch J2000.0
- **X Axis**: Direction of the mean vernal point (equinox) at J2000.0
- **Z Axis**: Perpendicular to the mean equator, toward the celestial north pole
- **Y Axis**: Completes the right-handed triad

#### 2.1.3 ICRS ↔ J2000 Relationship

The difference between ICRS and J2000 is very small (about 17.3 milliarcsec in orientation). The transformation is given by a **frame bias**:

**Transformation matrix J2000 → ICRS:**

```
       ⎡  0.9999999999999928  -0.0000000707827974  -0.0000000805621715 ⎤
B =    ⎢  0.0000000707827948   0.9999999999999969  -0.0000000330604145 ⎥
       ⎣  0.0000000805621738   0.0000000330604088   0.9999999999999962 ⎦
```

For most applications, ICRS and J2000 can be considered equivalent.

### 2.2 Fundamental Rotation Matrices

Transformations between coordinate systems are based on **rotation matrices** around the principal axes.

#### 2.2.1 Rotation about the X axis

```
         ⎡ 1     0         0     ⎤
Rx(θ) =  ⎢ 0   cos(θ)    sin(θ)  ⎥
         ⎣ 0  -sin(θ)    cos(θ)  ⎦
```

#### 2.2.2 Rotation about the Y axis

```
         ⎡ cos(θ)   0  -sin(θ) ⎤
Ry(θ) =  ⎢   0      1     0    ⎥
         ⎣ sin(θ)   0   cos(θ) ⎦
```

#### 2.2.3 Rotation about the Z axis

```
         ⎡  cos(θ)   sin(θ)   0 ⎤
Rz(θ) =  ⎢ -sin(θ)   cos(θ)   0 ⎥
         ⎣    0        0      1 ⎦
```

#### 2.2.4 Properties of Rotation Matrices

Rotation matrices are **orthogonal matrices**:

- Determinant = 1 (proper rotation)
- R⁻¹ = Rᵀ (the inverse equals the transpose)
- R(-θ) = Rᵀ(θ)

### 2.3 Transformations Between Reference Systems

#### 2.3.1 Equatorial J2000 ↔ Ecliptic J2000

The transformation uses the obliquity of the ecliptic ε₀ = 23.439291°:

**C++ Implementation:**

```cpp
namespace orbfit::coordinates {
    // Matrix from J2000 to Ecliptic
    static Matrix3d j2000_to_ecliptic() {
        constexpr double epsilon0 = 23.439291 * constants::DEG_TO_RAD;
        return rotation_x(epsilon0);
    }
    
    // Matrix from Ecliptic to J2000
    static Matrix3d ecliptic_to_j2000() {
        constexpr double epsilon0 = 23.439291 * constants::DEG_TO_RAD;
        return rotation_x(-epsilon0);
    }
}
```

**Numerical example:**

For ε₀ = 23.439291°:

```
cos(ε₀) = 0.9174821
sin(ε₀) = 0.3977772

           ⎡ 1.0000000   0.0000000   0.0000000 ⎤
Rx(ε₀) =   ⎢ 0.0000000   0.9174821   0.3977772 ⎥
           ⎣ 0.0000000  -0.3977772   0.9174821 ⎦
```

#### 2.3.2 Inertial System ↔ Terrestrial System (ITRF)

The transformation from inertial system (J2000) to terrestrial system (ITRF) requires:

1. **Precession**: Secular motion of Earth's rotation axis
2. **Nutation**: Periodic oscillations of Earth's axis
3. **Earth Rotation Angle (ERA)**: Daily rotation of Earth
4. **Polar Motion**: Movement of the pole relative to Earth's crust

**Simplified transformation (ERA only):**

```cpp
static Matrix3d j2000_to_itrf_simple(double mjd_ut1) {
    // Earth Rotation Angle
    double T = mjd_ut1 - constants::MJD2000;
    double era = 2.0 * constants::PI * 
                 (0.7790572732640 + 1.00273781191135448 * T);
    
    // Normalize to [0, 2π)
    era = std::fmod(era, 2.0 * constants::PI);
    if (era < 0.0) era += 2.0 * constants::PI;
    
    // Rotation about Z
    return rotation_z(-era);
}
```

#### 2.3.3 Greenwich Mean Sidereal Time (GMST)

**GMST** represents the hour angle of the mean vernal point at the Greenwich meridian.

**Formula (IAU 2000):**

```
GMST = 24110.54841 + 8640184.812866·T + 0.093104·T² - 6.2×10⁻⁶·T³
```

where T is in Julian centuries from J2000.0.

**Implementation:**

```cpp
static double gmst(double mjd_ut1) {
    double T = (mjd_ut1 - constants::MJD2000) / 36525.0;
    
    // GMST at 0h UT1 (in seconds)
    double gmst0 = 24110.54841 + 8640184.812866 * T 
                 + 0.093104 * T * T - 6.2e-6 * T * T * T;
    
    // Convert to days and add day fraction
    gmst0 /= constants::SECONDS_PER_DAY;
    double frac_day = mjd_ut1 - std::floor(mjd_ut1);
    double gmst_days = gmst0 + frac_day * 1.00273790935;
    
    // Convert to radians
    return 2.0 * constants::PI * (gmst_days - std::floor(gmst_days));
}
```

### 2.4 Orbital State Transformation

When transforming **orbital state vectors** (position + velocity) between rotating systems, the **Coriolis** term must be considered.

#### 2.4.1 Velocity Transformation

For a rotating system with angular velocity ω:

```
v_rot = R · v_inert - ω × (R · r_inert)
```

**Implementation:**

```cpp
static Vector3d transform_velocity(
    const Vector3d& pos, const Vector3d& vel,
    FrameType from, FrameType to,
    double mjd_ut1
) {
    Matrix3d R = get_transformation(from, to, mjd_ut1);
    Vector3d vel_rotated = R * vel;
    
    bool from_rotating = (from == FrameType::ITRF);
    bool to_rotating = (to == FrameType::ITRF);
    
    if (from_rotating != to_rotating) {
        // Earth angular velocity [rad/s]
        constexpr double omega_earth = 7.292115e-5;
        Vector3d omega(0.0, 0.0, omega_earth);
        
        Vector3d pos_rotated = R * pos;
        
        if (to_rotating) {
            vel_rotated -= omega.cross(pos_rotated);
        } else {
            vel_rotated += omega.cross(pos_rotated);
        }
    }
    
    return vel_rotated;
}
```

### 2.5 Precession and Nutation

#### 2.5.1 Precession

**Precession** is the slow conical motion of Earth's rotation axis caused by gravitational attractions of the Sun and Moon on Earth's equatorial bulge.

**Precession period**: ~25,772 years (one complete cycle)

**Obliquity variation**: ~0.47" per year

The IAU 2006 precession model provides:

```
ε(T) = ε₀ - 46.836769"·T - 0.0001831"·T² + 0.00200340"·T³ - ...
```

where T is in Julian centuries from J2000.0.

#### 2.5.2 Nutation

**Nutation** represents the periodic oscillations of Earth's axis superimposed on the precession motion.

**Main components:**
- 18.6-year period (retrograde motion of lunar nodes)
- Amplitude: ~9.2" in nutation in longitude
- Amplitude: ~6.9" in nutation in obliquity

### 2.6 Practical Example with Real Data

#### Example 2.1: Transformation of satellite coordinates

**Data:** Geostationary satellite
- Position (ITRF, MJD 60000.5): r = [42164.0, 0.0, 0.0] km
- Velocity (ITRF): v = [0.0, 3.0747, 0.0] km/s

**Objective:** Transform to J2000 coordinates

**GMST Calculation:**
```
MJD = 60000.5
T = (60000.5 - 51544.5) / 36525.0 = 0.2315
GMST ≈ 5.832 rad = 334.12°
```

**Rotation matrix (ITRF → J2000):**
```
           ⎡  0.4087  -0.9127   0 ⎤
R(GMST) =  ⎢  0.9127   0.4087   0 ⎥
           ⎣    0        0      1 ⎦
```

**Position in J2000:**
```
r_J2000 = R · r_ITRF
        = [0.4087×42164, 0.9127×42164, 0]
        = [17234.5, 38489.0, 0.0] km
```

**Velocity in J2000 (with Coriolis):**
```
ω = [0, 0, 7.292×10⁻⁵] rad/s
v_J2000 = R · v_ITRF + ω × r_J2000
        ≈ [-0.0257, 1.257, 0.0] + [2.807, -1.257, 0.0]
        = [2.781, 0.0, 0.0] km/s
```

### 2.7 Transformation Summary Diagram

```
                    ┌──────────────────┐
                    │      ICRS        │
                    │  (Barycentric)   │
                    └────────┬─────────┘
                             │ Frame Bias (17 mas)
                             ↓
                    ┌──────────────────┐
                    │   J2000 (EME)    │
                    │   Equatorial     │
                    └────────┬─────────┘
                             │
              ┌──────────────┼──────────────┐
              │              │              │
              ↓              ↓              ↓
     ┌────────────┐   ┌────────────┐  ┌────────────┐
     │  Ecliptic  │   │   TEME     │  │    MOD     │
     │   J2000    │   │  (True)    │  │  (Mean)    │
     └────────────┘   └──────┬─────┘  └────────────┘
                             │
                             ↓ ERA
                    ┌──────────────────┐
                    │      ITRF        │
                    │  (Terrestrial)   │
                    └──────────────────┘
```

**Legend:**
- ICRS: International Celestial Reference System
- EME: Earth Mean Equator
- MOD: Mean of Date
- TEME: True Equator Mean Equinox
- ITRF: International Terrestrial Reference Frame
- ERA: Earth Rotation Angle

---

## Riferimenti / References

1. IERS Conventions (2010). IERS Technical Note No. 36.

2. Kaplan, G.H. (2005). "The IAU Resolutions on Astronomical Reference Systems, Time Scales, and Earth Rotation Models". USNO Circular No. 179.

3. Capitaine, N. et al. (2003). "Expressions for IAU 2000 precession quantities". Astronomy & Astrophysics, 412, 567-586.

---

*Documento creato per il progetto ITALOccultLibrary*
*Document created for the ITALOccultLibrary project*

© 2025 Michele Bigi (mikbigi@gmail.com)
