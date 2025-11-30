# PARTE 1 - Capitolo 1: Principi di Base sui Sistemi di Coordinate
# PART 1 - Chapter 1: Basic Principles of Coordinate Systems

---

## 🇮🇹 Italiano

### 1.1 Introduzione ai Sistemi di Coordinate

Un **sistema di coordinate** è un metodo matematico per identificare univocamente la posizione di un punto nello spazio. In astrodinamica, utilizziamo diversi sistemi di coordinate a seconda del contesto e delle esigenze computazionali.

I sistemi di coordinate possono essere classificati in base a:

- **Dimensionalità**: 2D (sulla superficie di una sfera) o 3D (nello spazio)
- **Geometria**: Cartesiane, sferiche, cilindriche
- **Riferimento**: Geocentriche, eliocentriche, baricentriche
- **Rotazione**: Inerziali (fisse) o rotanti (solidali con un corpo)

### 1.2 Coordinate Cartesiane

Le **coordinate cartesiane** rappresentano la posizione di un punto tramite tre numeri reali (x, y, z) che indicano le distanze lungo tre assi ortogonali.

#### 1.2.1 Definizione Matematica

Dato un sistema di riferimento con origine O e tre assi ortogonali X, Y, Z, la posizione di un punto P è definita dal **vettore posizione**:

```
    ⎡ x ⎤
r = ⎢ y ⎥
    ⎣ z ⎦
```

La **distanza** dall'origine (raggio vettore) è:

```
r = |r| = √(x² + y² + z²)
```

#### 1.2.2 Implementazione nel Codice

Nel software OrbFit C++, i vettori cartesiani sono implementati usando la libreria Eigen:

```cpp
namespace orbfit {
    // Tipo vettore 3D
    using Vector3d = Eigen::Vector3d;
    
    // Esempio di utilizzo
    Vector3d position(1.0, 0.0, 0.0);  // 1 AU sull'asse X
    double radius = position.norm();    // Calcolo del modulo
}
```

### 1.3 Coordinate Sferiche

Le **coordinate sferiche** rappresentano un punto tramite:

- **r**: distanza radiale dall'origine
- **θ** (theta): angolo polare (colatitudine)
- **φ** (phi): angolo azimutale (longitudine)

#### 1.3.1 Conversione Cartesiane ↔ Sferiche

**Da Cartesiane a Sferiche:**

```
r = √(x² + y² + z²)

θ = arccos(z / r)

φ = atan2(y, x)
```

**Da Sferiche a Cartesiane:**

```
x = r · sin(θ) · cos(φ)
y = r · sin(θ) · sin(φ)
z = r · cos(θ)
```

#### 1.3.2 Diagramma

```
                Z
                ↑
                |      P(r, θ, φ)
                |     /
                |    /  θ
                |   /  ↙
                |  /
                | /__________ Y
               O/    φ
              /
             /
            X
```

### 1.4 Coordinate Equatoriali Celesti

Le **coordinate equatoriali** sono il sistema di riferimento fondamentale in astronomia. Si basano sulla proiezione dell'equatore terrestre sulla sfera celeste.

#### 1.4.1 Componenti

- **Ascensione Retta (α o RA)**: Angolo misurato verso est lungo l'equatore celeste, a partire dal punto vernale (γ). Intervallo: [0°, 360°) o [0h, 24h).

- **Declinazione (δ o Dec)**: Angolo misurato dal piano equatoriale verso i poli celesti. Intervallo: [-90°, +90°].

- **Distanza (r)**: Distanza dall'osservatore (opzionale, per coordinate 3D).

#### 1.4.2 Formule di Conversione

**Da Coordinate Equatoriali a Cartesiane (sistema equatoriale):**

```
x = r · cos(δ) · cos(α)
y = r · cos(δ) · sin(α)
z = r · sin(δ)
```

**Da Cartesiane a Coordinate Equatoriali:**

```
r = √(x² + y² + z²)

α = atan2(y, x)           [normalizzato a [0, 2π)]

δ = arcsin(z / r)         [oppure: arctan(z / √(x² + y²))]
```

### 1.5 Coordinate Eclittiche

Le **coordinate eclittiche** usano il piano dell'orbita terrestre (eclittica) come riferimento fondamentale.

#### 1.5.1 Componenti

- **Longitudine eclittica (λ)**: Angolo misurato lungo l'eclittica dal punto vernale. Intervallo: [0°, 360°).

- **Latitudine eclittica (β)**: Angolo misurato perpendicolarmente all'eclittica. Intervallo: [-90°, +90°].

#### 1.5.2 Obliquità dell'Eclittica

L'**obliquità dell'eclittica** (ε) è l'angolo tra il piano equatoriale e il piano eclittico. Al J2000.0:

```
ε₀ = 23.439291° = 23° 26' 21.448"
```

Questo valore varia lentamente nel tempo a causa della precessione e nutazione.

#### 1.5.3 Trasformazione Equatoriale ↔ Eclittica

La trasformazione tra sistemi equatoriali ed eclittici si effettua tramite una **rotazione attorno all'asse X** dell'angolo ε:

**Matrice di rotazione da Equatoriale a Eclittico:**

```
        ⎡ 1    0        0     ⎤
Rx(ε) = ⎢ 0   cos(ε)  sin(ε) ⎥
        ⎣ 0  -sin(ε)  cos(ε) ⎦
```

**Trasformazione:**

```
r_ecl = Rx(ε) · r_eq
```

**Trasformazione inversa:**

```
r_eq = Rx(-ε) · r_ecl = Rx(ε)ᵀ · r_ecl
```

### 1.6 Esempio Pratico con Dati Reali

#### Esempio 1.1: Conversione delle coordinate di Sirio

**Dati osservativi di Sirio (α Canis Majoris):**
- Ascensione Retta: α = 6h 45m 08.917s = 101.2872°
- Declinazione: δ = -16° 42' 58.02" = -16.7161°
- Distanza: d = 8.6 anni luce = 2.638 parsec

**Calcolo delle coordinate cartesiane equatoriali:**

Convertiamo la distanza in unità astronomiche:
```
d = 2.638 pc × 206265 AU/pc = 544,128 AU
```

Calcoliamo le componenti cartesiane:
```
x = d · cos(δ) · cos(α)
  = 544,128 · cos(-16.7161°) · cos(101.2872°)
  = 544,128 · 0.9576 · (-0.1956)
  = -101,888 AU

y = d · cos(δ) · sin(α)
  = 544,128 · 0.9576 · 0.9807
  = 510,957 AU

z = d · sin(δ)
  = 544,128 · (-0.2876)
  = -156,495 AU
```

**Verifica:**
```
r = √(x² + y² + z²)
  = √(101,888² + 510,957² + 156,495²)
  ≈ 544,128 AU ✓
```

### 1.7 Riepilogo delle Formule

| Conversione | Formula |
|-------------|---------|
| Cartesiane → Sferiche (r) | r = √(x² + y² + z²) |
| Cartesiane → Sferiche (θ) | θ = arccos(z/r) |
| Cartesiane → Sferiche (φ) | φ = atan2(y, x) |
| Sferiche → Cartesiane (x) | x = r·sin(θ)·cos(φ) |
| Sferiche → Cartesiane (y) | y = r·sin(θ)·sin(φ) |
| Sferiche → Cartesiane (z) | z = r·cos(θ) |
| Equatoriali → RA | α = atan2(y, x) |
| Equatoriali → Dec | δ = arcsin(z/r) |

---

## 🇬🇧 English

### 1.1 Introduction to Coordinate Systems

A **coordinate system** is a mathematical method for uniquely identifying the position of a point in space. In astrodynamics, we use different coordinate systems depending on context and computational requirements.

Coordinate systems can be classified based on:

- **Dimensionality**: 2D (on a sphere's surface) or 3D (in space)
- **Geometry**: Cartesian, spherical, cylindrical
- **Reference**: Geocentric, heliocentric, barycentric
- **Rotation**: Inertial (fixed) or rotating (body-fixed)

### 1.2 Cartesian Coordinates

**Cartesian coordinates** represent a point's position through three real numbers (x, y, z) indicating distances along three orthogonal axes.

#### 1.2.1 Mathematical Definition

Given a reference system with origin O and three orthogonal axes X, Y, Z, a point P's position is defined by the **position vector**:

```
    ⎡ x ⎤
r = ⎢ y ⎥
    ⎣ z ⎦
```

The **distance** from the origin (radius vector) is:

```
r = |r| = √(x² + y² + z²)
```

#### 1.2.2 Code Implementation

In OrbFit C++ software, Cartesian vectors are implemented using the Eigen library:

```cpp
namespace orbfit {
    // 3D vector type
    using Vector3d = Eigen::Vector3d;
    
    // Usage example
    Vector3d position(1.0, 0.0, 0.0);  // 1 AU on X axis
    double radius = position.norm();    // Magnitude calculation
}
```

### 1.3 Spherical Coordinates

**Spherical coordinates** represent a point through:

- **r**: radial distance from origin
- **θ** (theta): polar angle (colatitude)
- **φ** (phi): azimuthal angle (longitude)

#### 1.3.1 Cartesian ↔ Spherical Conversion

**From Cartesian to Spherical:**

```
r = √(x² + y² + z²)

θ = arccos(z / r)

φ = atan2(y, x)
```

**From Spherical to Cartesian:**

```
x = r · sin(θ) · cos(φ)
y = r · sin(θ) · sin(φ)
z = r · cos(θ)
```

#### 1.3.2 Diagram

```
                Z
                ↑
                |      P(r, θ, φ)
                |     /
                |    /  θ
                |   /  ↙
                |  /
                | /__________ Y
               O/    φ
              /
             /
            X
```

### 1.4 Celestial Equatorial Coordinates

**Equatorial coordinates** are the fundamental reference system in astronomy. They are based on the projection of Earth's equator onto the celestial sphere.

#### 1.4.1 Components

- **Right Ascension (α or RA)**: Angle measured eastward along the celestial equator, starting from the vernal point (γ). Range: [0°, 360°) or [0h, 24h).

- **Declination (δ or Dec)**: Angle measured from the equatorial plane toward the celestial poles. Range: [-90°, +90°].

- **Distance (r)**: Distance from observer (optional, for 3D coordinates).

#### 1.4.2 Conversion Formulas

**From Equatorial Coordinates to Cartesian (equatorial system):**

```
x = r · cos(δ) · cos(α)
y = r · cos(δ) · sin(α)
z = r · sin(δ)
```

**From Cartesian to Equatorial Coordinates:**

```
r = √(x² + y² + z²)

α = atan2(y, x)           [normalized to [0, 2π)]

δ = arcsin(z / r)         [or: arctan(z / √(x² + y²))]
```

### 1.5 Ecliptic Coordinates

**Ecliptic coordinates** use Earth's orbital plane (ecliptic) as the fundamental reference.

#### 1.5.1 Components

- **Ecliptic longitude (λ)**: Angle measured along the ecliptic from the vernal point. Range: [0°, 360°).

- **Ecliptic latitude (β)**: Angle measured perpendicular to the ecliptic. Range: [-90°, +90°].

#### 1.5.2 Obliquity of the Ecliptic

The **obliquity of the ecliptic** (ε) is the angle between the equatorial plane and the ecliptic plane. At J2000.0:

```
ε₀ = 23.439291° = 23° 26' 21.448"
```

This value varies slowly over time due to precession and nutation.

#### 1.5.3 Equatorial ↔ Ecliptic Transformation

The transformation between equatorial and ecliptic systems is performed through a **rotation about the X axis** by angle ε:

**Rotation matrix from Equatorial to Ecliptic:**

```
        ⎡ 1    0        0     ⎤
Rx(ε) = ⎢ 0   cos(ε)  sin(ε) ⎥
        ⎣ 0  -sin(ε)  cos(ε) ⎦
```

**Transformation:**

```
r_ecl = Rx(ε) · r_eq
```

**Inverse transformation:**

```
r_eq = Rx(-ε) · r_ecl = Rx(ε)ᵀ · r_ecl
```

### 1.6 Practical Example with Real Data

#### Example 1.1: Converting Sirius Coordinates

**Observational data for Sirius (α Canis Majoris):**
- Right Ascension: α = 6h 45m 08.917s = 101.2872°
- Declination: δ = -16° 42' 58.02" = -16.7161°
- Distance: d = 8.6 light-years = 2.638 parsec

**Calculating equatorial Cartesian coordinates:**

Convert distance to astronomical units:
```
d = 2.638 pc × 206265 AU/pc = 544,128 AU
```

Calculate Cartesian components:
```
x = d · cos(δ) · cos(α)
  = 544,128 · cos(-16.7161°) · cos(101.2872°)
  = 544,128 · 0.9576 · (-0.1956)
  = -101,888 AU

y = d · cos(δ) · sin(α)
  = 544,128 · 0.9576 · 0.9807
  = 510,957 AU

z = d · sin(δ)
  = 544,128 · (-0.2876)
  = -156,495 AU
```

**Verification:**
```
r = √(x² + y² + z²)
  = √(101,888² + 510,957² + 156,495²)
  ≈ 544,128 AU ✓
```

### 1.7 Formula Summary

| Conversion | Formula |
|------------|---------|
| Cartesian → Spherical (r) | r = √(x² + y² + z²) |
| Cartesian → Spherical (θ) | θ = arccos(z/r) |
| Cartesian → Spherical (φ) | φ = atan2(y, x) |
| Spherical → Cartesian (x) | x = r·sin(θ)·cos(φ) |
| Spherical → Cartesian (y) | y = r·sin(θ)·sin(φ) |
| Spherical → Cartesian (z) | z = r·cos(θ) |
| Equatorial → RA | α = atan2(y, x) |
| Equatorial → Dec | δ = arcsin(z/r) |

---

## Riferimenti / References

1. Smart, W.M. (1977). *Textbook on Spherical Astronomy*. Cambridge University Press.

2. Green, R.M. (1985). *Spherical Astronomy*. Cambridge University Press.

3. IAU SOFA Library - Standards of Fundamental Astronomy. http://www.iausofa.org/

---

*Documento creato per il progetto ITALOccultLibrary*
*Document created for the ITALOccultLibrary project*

© 2025 Michele Bigi (mikbigi@gmail.com)
