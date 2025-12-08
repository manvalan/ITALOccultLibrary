# PARTE 2 - Capitolo 2: UTC e Framework Temporali Astrodinamici
# PART 2 - Chapter 2: UTC and Astrodynamic Time Frameworks

---

## 🇮🇹 Italiano

### 2.1 Panoramica delle Scale Temporali

In astronomia e astrodinamica si utilizzano diverse **scale temporali**, ognuna ottimizzata per uno scopo specifico. La scelta della scala corretta è fondamentale per la precisione dei calcoli.

#### 2.1.1 Classificazione delle Scale Temporali

| Categoria | Scale | Uso Principale |
|-----------|-------|----------------|
| **Atomiche** | TAI, GPS | Standard di riferimento, navigazione |
| **Dinamiche** | TT, TDB, TCB | Effemeridi, meccanica celeste |
| **Civili** | UTC, UT1 | Uso quotidiano, navigazione |
| **Rotazione** | UT1, GMST | Orientamento Terra-spazio |

### 2.2 UTC - Coordinated Universal Time

L'**UTC** (Coordinated Universal Time) è la scala temporale civile internazionale. È mantenuta sincronizzata con la rotazione terrestre tramite l'inserimento di **secondi intercalari** (leap seconds).

#### 2.2.1 Caratteristiche

- **Base**: Orologi atomici (come TAI)
- **Sincronizzazione**: Mantenuta entro 0.9 secondi da UT1
- **Correzione**: Secondi intercalari aggiunti a fine giugno o dicembre
- **Uso**: Standard civile mondiale, telecomunicazioni

#### 2.2.2 Secondi Intercalari (Leap Seconds)

Dal 1972, i secondi intercalari sono stati inseriti per mantenere UTC entro ±0.9 s da UT1:

| Data | Leap Seconds Totali |
|------|---------------------|
| 1972-01-01 | 10 |
| 1980-01-01 | 19 |
| 1990-01-01 | 25 |
| 2000-01-01 | 32 |
| 2017-01-01 | 37 |
| 2025 (attuale) | 37 |

**Nota:** L'IAU ha approvato l'abolizione dei leap seconds dal 2035.

#### 2.2.3 Gestione Leap Seconds nel Codice

```cpp
namespace orbfit::time {

    // Costante: leap seconds al 2025
    constexpr int LEAP_SECONDS_2025 = 37;
    
    // Carica tabella leap seconds
    bool load_leap_seconds(const std::string& filepath);
    
    // Ottieni leap seconds per una data
    int get_leap_seconds(double mjd_utc);
    
}
```

### 2.3 TAI - International Atomic Time

Il **TAI** (Temps Atomique International) è la scala temporale atomica fondamentale, mantenuta dal BIPM (Bureau International des Poids et Mesures).

#### 2.3.1 Caratteristiche

- **Definizione**: Media ponderata di ~500 orologi atomici nel mondo
- **Unità**: Secondo SI (basato sulla transizione del cesio-133)
- **Continuità**: Nessun leap second, conteggio continuo
- **Origine**: 1 gennaio 1958, 0:00 UT

#### 2.3.2 Relazione TAI ↔ UTC

```
TAI = UTC + ΔAT
```

Dove ΔAT = numero di leap seconds (37 nel 2025).

**Implementazione:**

```cpp
// UTC → TAI
double utc_to_tai(double mjd_utc, int leap_seconds = LEAP_SECONDS_2025) {
    return mjd_utc + leap_seconds / constants::SECONDS_PER_DAY;
}

// TAI → UTC
double tai_to_utc(double mjd_tai, int leap_seconds = LEAP_SECONDS_2025) {
    return mjd_tai - leap_seconds / constants::SECONDS_PER_DAY;
}
```

### 2.4 TT - Terrestrial Time

Il **TT** (Terrestrial Time), precedentemente TDT (Terrestrial Dynamical Time), è la scala temporale usata per le effemeridi geocentriche.

#### 2.4.1 Caratteristiche

- **Uso**: Effemeridi geocentriche, meccanica celeste
- **Definizione**: TT = TAI + 32.184 s
- **Origine della costante**: Continuità storica con ET (Ephemeris Time)
- **Proprietà**: Scala uniforme, non soggetta a leap seconds

#### 2.4.2 Formule di Conversione

```
TT = TAI + 32.184 s
TT = UTC + ΔAT + 32.184 s
```

**Implementazione:**

```cpp
// Costante TT - TAI
constexpr double TT_MINUS_TAI = 32.184;  // secondi

// TAI → TT
double tai_to_tt(double mjd_tai) {
    return mjd_tai + TT_MINUS_TAI / constants::SECONDS_PER_DAY;
}

// TT → TAI
double tt_to_tai(double mjd_tt) {
    return mjd_tt - TT_MINUS_TAI / constants::SECONDS_PER_DAY;
}

// UTC → TT (conversione diretta)
double utc_to_tt(double mjd_utc) {
    return tai_to_tt(utc_to_tai(mjd_utc));
}

// TT → UTC (conversione diretta)
double tt_to_utc(double mjd_tt) {
    return tai_to_utc(tt_to_tai(mjd_tt));
}
```

### 2.5 TDB - Barycentric Dynamical Time

Il **TDB** (Barycentric Dynamical Time) è la scala temporale usata per le effemeridi baricentriche del Sistema Solare.

#### 2.5.1 Caratteristiche

- **Uso**: Effemeridi JPL, meccanica celeste baricentrica
- **Definizione**: Scala relativistica al baricentro del Sistema Solare
- **Relazione con TT**: Differenza periodica di ~1.6 ms (max)

#### 2.5.2 Formula di Conversione TT ↔ TDB

La differenza TDB - TT è data da una serie di termini periodici. L'approssimazione di Fairhead & Bretagnon (1990):

```
TDB - TT ≈ 0.001658 × sin(g) + 0.000014 × sin(2g)  [secondi]
```

Dove **g** è l'anomalia media della Terra:

```
g = 357.528° + 35999.050° × T  [in gradi]
T = secoli giuliani da J2000.0
```

**Implementazione:**

```cpp
double tt_to_tdb(double mjd_tt) {
    double jd_tt = mjd_to_jd(mjd_tt);
    double T = julian_centuries_j2000(jd_tt);
    
    // Anomalia media della Terra (radianti)
    double g = constants::DEG_TO_RAD * (357.528 + 35999.050 * T);
    
    // TDB - TT in secondi
    double dt_seconds = 0.001658 * std::sin(g) + 0.000014 * std::sin(2.0 * g);
    
    return mjd_tt + dt_seconds / constants::SECONDS_PER_DAY;
}
```

#### 2.5.3 Conversione Inversa TDB → TT

La conversione inversa richiede un metodo iterativo:

```cpp
double tdb_to_tt(double mjd_tdb) {
    double mjd_tt = mjd_tdb;  // Prima approssimazione
    
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
```

### 2.6 UT1 - Universal Time

L'**UT1** è la scala temporale legata alla rotazione effettiva della Terra.

#### 2.6.1 Caratteristiche

- **Definizione**: Angolo di rotazione della Terra rispetto alle stelle
- **Irregolarità**: Variazioni dovute a fenomeni geofisici
- **Uso**: Trasformazioni tra sistemi celesti e terrestri

#### 2.6.2 DUT1 = UT1 - UTC

La differenza tra UT1 e UTC è chiamata **DUT1**:

```
UT1 = UTC + DUT1
```

- **Range**: -0.9 s < DUT1 < +0.9 s (per definizione)
- **Pubblicazione**: IERS Bulletin A (settimanale)

**Implementazione semplificata:**

```cpp
// Ottieni DUT1 per una data (richiede dati IERS)
double get_dut1(double mjd_utc) {
    // TODO: Caricare da file finals.data IERS
    return 0.0;  // Placeholder
}

// UTC → UT1
double utc_to_ut1(double mjd_utc) {
    double dut1 = get_dut1(mjd_utc);
    return mjd_utc + dut1 / constants::SECONDS_PER_DAY;
}
```

### 2.7 GPS Time

Il **GPS Time** è la scala temporale del sistema di navigazione GPS.

#### 2.7.1 Caratteristiche

- **Origine**: 6 gennaio 1980, 0:00 UTC (quando UTC = TAI - 19 s)
- **Continuità**: Nessun leap second
- **Relazione**: GPS = TAI - 19 s

**Implementazione:**

```cpp
constexpr double GPS_MINUS_TAI = -19.0;  // secondi

// TAI → GPS
double tai_to_gps(double mjd_tai) {
    return mjd_tai + GPS_MINUS_TAI / constants::SECONDS_PER_DAY;
}

// GPS → TAI
double gps_to_tai(double mjd_gps) {
    return mjd_gps - GPS_MINUS_TAI / constants::SECONDS_PER_DAY;
}
```

### 2.8 Schema delle Relazioni tra Scale Temporali

```
                     ┌─────────────────────────────────────────┐
                     │              TAI                        │
                     │   (International Atomic Time)           │
                     └─────────────┬───────────────────────────┘
                                   │
         ┌─────────────────────────┼─────────────────────────┐
         │                         │                         │
         ↓ +32.184 s               ↓ -ΔAT                    ↓ -19 s
    ┌─────────┐              ┌─────────┐              ┌─────────┐
    │   TT    │              │   UTC   │              │   GPS   │
    │ (Terr.) │              │ (Civil) │              │ (Navig.)│
    └────┬────┘              └────┬────┘              └─────────┘
         │                        │
         ↓ +periodic              ↓ +DUT1
    ┌─────────┐              ┌─────────┐
    │   TDB   │              │   UT1   │
    │ (Baryc.)│              │ (Rotat.)│
    └─────────┘              └─────────┘
```

**Relazioni chiave:**

| Da | A | Formula |
|----|---|---------|
| UTC | TAI | TAI = UTC + ΔAT (37 s nel 2025) |
| TAI | TT | TT = TAI + 32.184 s |
| TT | TDB | TDB ≈ TT + 0.0017 sin(g) s |
| UTC | UT1 | UT1 = UTC + DUT1 |
| TAI | GPS | GPS = TAI - 19 s |

### 2.9 Esempi Pratici con Dati Reali

#### Esempio 2.1: Conversione completa di una data

**Problema:** Per il 1 gennaio 2025, 12:00:00 UTC, calcolare TT, TDB, TAI.

**Dati:**
- Data: 1 gennaio 2025, 12:00:00 UTC
- MJD(UTC) = 60676.5
- ΔAT = 37 s (leap seconds cumulativi)

**Calcoli:**

1. **UTC → TAI:**
   ```
   MJD(TAI) = MJD(UTC) + 37/86400
            = 60676.5 + 0.000428
            = 60676.500428
   ```

2. **TAI → TT:**
   ```
   MJD(TT) = MJD(TAI) + 32.184/86400
           = 60676.500428 + 0.000373
           = 60676.500801
   ```

3. **TT → TDB:**
   ```
   T = (60676.500801 - 51544.5) / 36525 = 0.250074 secoli
   g = 357.528 + 35999.050 × 0.250074 = 9357.0° ≈ 357.0° (mod 360)
   g = 6.231 rad
   
   Δt = 0.001658 × sin(6.231) + 0.000014 × sin(12.462)
      = 0.001658 × (-0.052) + 0.000014 × (-0.104)
      ≈ -0.000088 s
   
   MJD(TDB) = 60676.500801 - 0.000088/86400
            ≈ 60676.500800
   ```

**Risultati:**

| Scala | MJD | Differenza da UTC |
|-------|-----|-------------------|
| UTC | 60676.500000 | 0 |
| TAI | 60676.500428 | +37.0 s |
| TT | 60676.500801 | +69.184 s |
| TDB | 60676.500800 | ~+69.18 s |

#### Esempio 2.2: Calcolo dell'ora corrente

**Implementazione per ottenere l'ora corrente in diverse scale:**

```cpp
#include <chrono>

double now(TimeScale time_scale) {
    // Ottieni tempo sistema
    auto now = std::chrono::system_clock::now();
    auto duration = now.time_since_epoch();
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();
    
    // Unix epoch (1970-01-01 00:00:00) = MJD 40587.0
    double mjd_utc = 40587.0 + seconds / constants::SECONDS_PER_DAY;
    
    // Converti alla scala richiesta
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
```

### 2.10 Riepilogo delle Formule di Conversione

| Conversione | Formula |
|-------------|---------|
| UTC → TAI | TAI = UTC + ΔAT |
| TAI → UTC | UTC = TAI - ΔAT |
| TAI → TT | TT = TAI + 32.184 s |
| TT → TAI | TAI = TT - 32.184 s |
| UTC → TT | TT = UTC + ΔAT + 32.184 s |
| TT → TDB | TDB ≈ TT + 0.001658·sin(g) s |
| UTC → UT1 | UT1 = UTC + DUT1 |
| TAI → GPS | GPS = TAI - 19 s |

### 2.11 Costanti di Riferimento

```cpp
namespace orbfit::time {

    // Costanti di conversione tra scale temporali
    constexpr int LEAP_SECONDS_2025 = 37;     // ΔAT al 2025
    constexpr double TT_MINUS_TAI = 32.184;   // secondi
    constexpr double GPS_MINUS_TAI = -19.0;   // secondi
    
}

namespace orbfit::constants {
    
    // Costanti temporali generali
    constexpr double SECONDS_PER_DAY = 86400.0;
    constexpr double DAYS_PER_CENTURY = 36525.0;
    constexpr double JD2000 = 2451545.0;
    constexpr double MJD2000 = 51544.5;
    
}
```

---

## 🇬🇧 English

### 2.1 Overview of Time Scales

In astronomy and astrodynamics, several **time scales** are used, each optimized for a specific purpose. Choosing the correct scale is fundamental for calculation precision.

#### 2.1.1 Time Scale Classification

| Category | Scales | Primary Use |
|----------|--------|-------------|
| **Atomic** | TAI, GPS | Reference standard, navigation |
| **Dynamic** | TT, TDB, TCB | Ephemerides, celestial mechanics |
| **Civil** | UTC, UT1 | Daily use, navigation |
| **Rotation** | UT1, GMST | Earth-space orientation |

### 2.2 UTC - Coordinated Universal Time

**UTC** (Coordinated Universal Time) is the international civil time scale. It is kept synchronized with Earth's rotation through the insertion of **leap seconds**.

#### 2.2.1 Characteristics

- **Base**: Atomic clocks (like TAI)
- **Synchronization**: Maintained within 0.9 seconds of UT1
- **Correction**: Leap seconds added at end of June or December
- **Use**: World civil standard, telecommunications

#### 2.2.2 Leap Seconds

Since 1972, leap seconds have been inserted to keep UTC within ±0.9 s of UT1:

| Date | Total Leap Seconds |
|------|-------------------|
| 1972-01-01 | 10 |
| 1980-01-01 | 19 |
| 1990-01-01 | 25 |
| 2000-01-01 | 32 |
| 2017-01-01 | 37 |
| 2025 (current) | 37 |

**Note:** The IAU has approved the abolition of leap seconds from 2035.

#### 2.2.3 Leap Seconds Management in Code

```cpp
namespace orbfit::time {

    // Constant: leap seconds as of 2025
    constexpr int LEAP_SECONDS_2025 = 37;
    
    // Load leap seconds table
    bool load_leap_seconds(const std::string& filepath);
    
    // Get leap seconds for a date
    int get_leap_seconds(double mjd_utc);
    
}
```

### 2.3 TAI - International Atomic Time

**TAI** (Temps Atomique International) is the fundamental atomic time scale, maintained by BIPM (Bureau International des Poids et Mesures).

#### 2.3.1 Characteristics

- **Definition**: Weighted average of ~500 atomic clocks worldwide
- **Unit**: SI second (based on cesium-133 transition)
- **Continuity**: No leap seconds, continuous counting
- **Origin**: January 1, 1958, 0:00 UT

#### 2.3.2 TAI ↔ UTC Relationship

```
TAI = UTC + ΔAT
```

Where ΔAT = number of leap seconds (37 in 2025).

**Implementation:**

```cpp
// UTC → TAI
double utc_to_tai(double mjd_utc, int leap_seconds = LEAP_SECONDS_2025) {
    return mjd_utc + leap_seconds / constants::SECONDS_PER_DAY;
}

// TAI → UTC
double tai_to_utc(double mjd_tai, int leap_seconds = LEAP_SECONDS_2025) {
    return mjd_tai - leap_seconds / constants::SECONDS_PER_DAY;
}
```

### 2.4 TT - Terrestrial Time

**TT** (Terrestrial Time), formerly TDT (Terrestrial Dynamical Time), is the time scale used for geocentric ephemerides.

#### 2.4.1 Characteristics

- **Use**: Geocentric ephemerides, celestial mechanics
- **Definition**: TT = TAI + 32.184 s
- **Constant origin**: Historical continuity with ET (Ephemeris Time)
- **Property**: Uniform scale, not subject to leap seconds

#### 2.4.2 Conversion Formulas

```
TT = TAI + 32.184 s
TT = UTC + ΔAT + 32.184 s
```

**Implementation:**

```cpp
// TT - TAI constant
constexpr double TT_MINUS_TAI = 32.184;  // seconds

// TAI → TT
double tai_to_tt(double mjd_tai) {
    return mjd_tai + TT_MINUS_TAI / constants::SECONDS_PER_DAY;
}

// TT → TAI
double tt_to_tai(double mjd_tt) {
    return mjd_tt - TT_MINUS_TAI / constants::SECONDS_PER_DAY;
}

// UTC → TT (direct conversion)
double utc_to_tt(double mjd_utc) {
    return tai_to_tt(utc_to_tai(mjd_utc));
}

// TT → UTC (direct conversion)
double tt_to_utc(double mjd_tt) {
    return tai_to_utc(tt_to_tai(mjd_tt));
}
```

### 2.5 TDB - Barycentric Dynamical Time

**TDB** (Barycentric Dynamical Time) is the time scale used for barycentric Solar System ephemerides.

#### 2.5.1 Characteristics

- **Use**: JPL ephemerides, barycentric celestial mechanics
- **Definition**: Relativistic scale at the Solar System barycenter
- **Relationship with TT**: Periodic difference of ~1.6 ms (max)

#### 2.5.2 TT ↔ TDB Conversion Formula

The TDB - TT difference is given by a series of periodic terms. The Fairhead & Bretagnon (1990) approximation:

```
TDB - TT ≈ 0.001658 × sin(g) + 0.000014 × sin(2g)  [seconds]
```

Where **g** is Earth's mean anomaly:

```
g = 357.528° + 35999.050° × T  [in degrees]
T = Julian centuries from J2000.0
```

**Implementation:**

```cpp
double tt_to_tdb(double mjd_tt) {
    double jd_tt = mjd_to_jd(mjd_tt);
    double T = julian_centuries_j2000(jd_tt);
    
    // Earth mean anomaly (radians)
    double g = constants::DEG_TO_RAD * (357.528 + 35999.050 * T);
    
    // TDB - TT in seconds
    double dt_seconds = 0.001658 * std::sin(g) + 0.000014 * std::sin(2.0 * g);
    
    return mjd_tt + dt_seconds / constants::SECONDS_PER_DAY;
}
```

#### 2.5.3 Inverse Conversion TDB → TT

The inverse conversion requires an iterative method:

```cpp
double tdb_to_tt(double mjd_tdb) {
    double mjd_tt = mjd_tdb;  // First approximation
    
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
```

### 2.6 UT1 - Universal Time

**UT1** is the time scale linked to Earth's actual rotation.

#### 2.6.1 Characteristics

- **Definition**: Earth's rotation angle relative to the stars
- **Irregularity**: Variations due to geophysical phenomena
- **Use**: Transformations between celestial and terrestrial systems

#### 2.6.2 DUT1 = UT1 - UTC

The difference between UT1 and UTC is called **DUT1**:

```
UT1 = UTC + DUT1
```

- **Range**: -0.9 s < DUT1 < +0.9 s (by definition)
- **Publication**: IERS Bulletin A (weekly)

**Simplified implementation:**

```cpp
// Get DUT1 for a date (requires IERS data)
double get_dut1(double mjd_utc) {
    // TODO: Load from IERS finals.data file
    return 0.0;  // Placeholder
}

// UTC → UT1
double utc_to_ut1(double mjd_utc) {
    double dut1 = get_dut1(mjd_utc);
    return mjd_utc + dut1 / constants::SECONDS_PER_DAY;
}
```

### 2.7 GPS Time

**GPS Time** is the time scale of the GPS navigation system.

#### 2.7.1 Characteristics

- **Origin**: January 6, 1980, 0:00 UTC (when UTC = TAI - 19 s)
- **Continuity**: No leap seconds
- **Relationship**: GPS = TAI - 19 s

**Implementation:**

```cpp
constexpr double GPS_MINUS_TAI = -19.0;  // seconds

// TAI → GPS
double tai_to_gps(double mjd_tai) {
    return mjd_tai + GPS_MINUS_TAI / constants::SECONDS_PER_DAY;
}

// GPS → TAI
double gps_to_tai(double mjd_gps) {
    return mjd_gps - GPS_MINUS_TAI / constants::SECONDS_PER_DAY;
}
```

### 2.8 Time Scale Relationship Diagram

```
                     ┌─────────────────────────────────────────┐
                     │              TAI                        │
                     │   (International Atomic Time)           │
                     └─────────────┬───────────────────────────┘
                                   │
         ┌─────────────────────────┼─────────────────────────┐
         │                         │                         │
         ↓ +32.184 s               ↓ -ΔAT                    ↓ -19 s
    ┌─────────┐              ┌─────────┐              ┌─────────┐
    │   TT    │              │   UTC   │              │   GPS   │
    │ (Terr.) │              │ (Civil) │              │ (Navig.)│
    └────┬────┘              └────┬────┘              └─────────┘
         │                        │
         ↓ +periodic              ↓ +DUT1
    ┌─────────┐              ┌─────────┐
    │   TDB   │              │   UT1   │
    │ (Baryc.)│              │ (Rotat.)│
    └─────────┘              └─────────┘
```

**Key relationships:**

| From | To | Formula |
|------|-----|---------|
| UTC | TAI | TAI = UTC + ΔAT (37 s in 2025) |
| TAI | TT | TT = TAI + 32.184 s |
| TT | TDB | TDB ≈ TT + 0.0017 sin(g) s |
| UTC | UT1 | UT1 = UTC + DUT1 |
| TAI | GPS | GPS = TAI - 19 s |

### 2.9 Practical Examples with Real Data

#### Example 2.1: Complete Date Conversion

**Problem:** For January 1, 2025, 12:00:00 UTC, calculate TT, TDB, TAI.

**Data:**
- Date: January 1, 2025, 12:00:00 UTC
- MJD(UTC) = 60676.5
- ΔAT = 37 s (cumulative leap seconds)

**Calculations:**

1. **UTC → TAI:**
   ```
   MJD(TAI) = MJD(UTC) + 37/86400
            = 60676.5 + 0.000428
            = 60676.500428
   ```

2. **TAI → TT:**
   ```
   MJD(TT) = MJD(TAI) + 32.184/86400
           = 60676.500428 + 0.000373
           = 60676.500801
   ```

3. **TT → TDB:**
   ```
   T = (60676.500801 - 51544.5) / 36525 = 0.250074 centuries
   g = 357.528 + 35999.050 × 0.250074 = 9357.0° ≈ 357.0° (mod 360)
   g = 6.231 rad
   
   Δt = 0.001658 × sin(6.231) + 0.000014 × sin(12.462)
      = 0.001658 × (-0.052) + 0.000014 × (-0.104)
      ≈ -0.000088 s
   
   MJD(TDB) = 60676.500801 - 0.000088/86400
            ≈ 60676.500800
   ```

**Results:**

| Scale | MJD | Difference from UTC |
|-------|-----|---------------------|
| UTC | 60676.500000 | 0 |
| TAI | 60676.500428 | +37.0 s |
| TT | 60676.500801 | +69.184 s |
| TDB | 60676.500800 | ~+69.18 s |

#### Example 2.2: Current Time Calculation

**Implementation to get current time in different scales:**

```cpp
#include <chrono>

double now(TimeScale time_scale) {
    // Get system time
    auto now = std::chrono::system_clock::now();
    auto duration = now.time_since_epoch();
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();
    
    // Unix epoch (1970-01-01 00:00:00) = MJD 40587.0
    double mjd_utc = 40587.0 + seconds / constants::SECONDS_PER_DAY;
    
    // Convert to requested scale
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
```

### 2.10 Conversion Formula Summary

| Conversion | Formula |
|------------|---------|
| UTC → TAI | TAI = UTC + ΔAT |
| TAI → UTC | UTC = TAI - ΔAT |
| TAI → TT | TT = TAI + 32.184 s |
| TT → TAI | TAI = TT - 32.184 s |
| UTC → TT | TT = UTC + ΔAT + 32.184 s |
| TT → TDB | TDB ≈ TT + 0.001658·sin(g) s |
| UTC → UT1 | UT1 = UTC + DUT1 |
| TAI → GPS | GPS = TAI - 19 s |

### 2.11 Reference Constants

```cpp
namespace orbfit::time {

    // Time scale conversion constants
    constexpr int LEAP_SECONDS_2025 = 37;     // ΔAT as of 2025
    constexpr double TT_MINUS_TAI = 32.184;   // seconds
    constexpr double GPS_MINUS_TAI = -19.0;   // seconds
    
}

namespace orbfit::constants {
    
    // General time constants
    constexpr double SECONDS_PER_DAY = 86400.0;
    constexpr double DAYS_PER_CENTURY = 36525.0;
    constexpr double JD2000 = 2451545.0;
    constexpr double MJD2000 = 51544.5;
    
}
```

---

## Riferimenti / References

1. McCarthy, D.D., & Seidelmann, P.K. (2009). *TIME: From Earth Rotation to Atomic Physics*. Wiley-VCH.

2. IERS Conventions (2010). IERS Technical Note No. 36, Chapter 10.

3. Fairhead, L., & Bretagnon, P. (1990). "An Analytical Formula for the Time Transformation TB-TT". Astronomy & Astrophysics, 229, 240-247.

4. USNO: ftp://maia.usno.navy.mil/ser7/finals.all (UT1-UTC data)

5. BIPM: https://www.bipm.org/en/time-metrology (TAI reference)

---

*Documento creato per il progetto ITALOccultLibrary*
*Document created for the ITALOccultLibrary project*

© 2025 Michele Bigi (mikbigi@gmail.com)
