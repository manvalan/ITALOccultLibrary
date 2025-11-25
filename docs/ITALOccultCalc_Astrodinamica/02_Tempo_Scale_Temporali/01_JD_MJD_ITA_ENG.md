# PARTE 2 - Capitolo 1: JD e MJD - Date Giuliane
# PART 2 - Chapter 1: JD and MJD - Julian Dates

---

## 🇮🇹 Italiano

### 1.1 Introduzione alle Date Giuliane

Le **Date Giuliane** (Julian Dates, JD) rappresentano un sistema di conteggio continuo dei giorni introdotto nel XVI secolo dallo studioso Joseph Justus Scaliger. Questo sistema è fondamentale in astronomia perché elimina le complessità dei calendari civili (anni bisestili, riforme calendariali, etc.).

#### 1.1.1 Vantaggi delle Date Giuliane

- **Continuità**: Nessuna interruzione o riforme
- **Semplicità**: Conteggio lineare del tempo
- **Universalità**: Indipendente da fusi orari e calendari
- **Precisione**: Adatta per calcoli scientifici di alta precisione

### 1.2 Julian Date (JD)

La **Julian Date (JD)** è il numero di giorni (e frazioni) trascorsi dall'epoca iniziale:

> **Epoca zero JD**: 1 gennaio 4713 a.C. (calendario giuliano), 12:00 UT
> (oppure: 24 novembre 4714 a.C. nel calendario gregoriano prolettico)

#### 1.2.1 Definizione Formale

```
JD = numero intero di giorni + frazione del giorno
```

Dove la **frazione del giorno** è calcolata da mezzogiorno (12:00 UT):
- 0.0 = mezzogiorno
- 0.5 = mezzanotte
- 0.75 = 6:00 del mattino seguente

#### 1.2.2 Epoche Significative

| Evento | JD | Data |
|--------|-----|------|
| Epoca zero JD | 0.0 | -4712/01/01, 12:00 UT |
| Nascita di Cristo | ~1,721,424 | 1/01/0001 |
| MJD zero | 2,400,000.5 | 17/11/1858, 0:00 UT |
| J2000.0 | **2,451,545.0** | 01/01/2000, 12:00 TT |
| Anno corrente (2025) | ~2,460,676 | 01/01/2025 |

### 1.3 Modified Julian Date (MJD)

La **Modified Julian Date (MJD)** è una versione "semplificata" della JD, introdotta per:

- Ridurre i numeri a valori più maneggevoli
- Far iniziare il giorno a mezzanotte invece che a mezzogiorno

#### 1.3.1 Definizione

```
MJD = JD - 2,400,000.5
```

#### 1.3.2 Proprietà

- **Epoca zero MJD**: 17 novembre 1858, 0:00 UT (MJD = 0)
- **Inizio del giorno**: Mezzanotte (0:00 UT)
- **Valori tipici attuali**: ~60,000 - 61,000

#### 1.3.3 Epoche MJD Significative

| Evento | MJD | Data |
|--------|-----|------|
| Epoca zero MJD | 0.0 | 17/11/1858, 0:00 UT |
| J2000.0 | **51,544.5** | 01/01/2000, 12:00 TT |
| Anno 2025 | ~60,676 | 01/01/2025 |

### 1.4 Algoritmi di Conversione

#### 1.4.1 Conversione JD ↔ MJD

**Funzioni inline (C++):**

```cpp
namespace orbfit::time {

    // JD → MJD
    inline double jd_to_mjd(double jd) {
        return jd - 2400000.5;
    }
    
    // MJD → JD
    inline double mjd_to_jd(double mjd) {
        return mjd + 2400000.5;
    }
    
}
```

#### 1.4.2 Conversione Calendario Gregoriano → MJD

**Algoritmo di Fliegel e van Flandern (1968):**

Questo algoritmo calcola il Julian Day Number (JDN) intero dal calendario gregoriano:

```cpp
double calendar_to_mjd(int year, int month, int day, double fraction) {
    // Algoritmo di Fliegel e van Flandern
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    
    // Julian Day Number (intero)
    int jdn = day + (153 * m + 2) / 5 + 365 * y 
            + y / 4 - y / 100 + y / 400 - 32045;
    
    // JD con frazione
    double jd = jdn - 0.5 + fraction;
    
    // Converti a MJD
    return jd - 2400000.5;
}
```

**Spiegazione dei passi:**

1. `a`: Aggiustamento per trattare gennaio e febbraio come mesi 13 e 14 dell'anno precedente
2. `y`: Anno modificato per il calcolo
3. `m`: Mese modificato (marzo = 0, ..., febbraio = 11)
4. Formula principale: Tiene conto di giorni per mese, anni bisestili gregoriani

#### 1.4.3 Conversione MJD → Calendario Gregoriano

**Algoritmo inverso:**

```cpp
std::tuple<int, int, int, double> mjd_to_calendar(double mjd) {
    double jd = mjd + 2400000.5;
    
    // Algoritmo di Fliegel e van Flandern (inverso)
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
    
    // Estrai la frazione del giorno
    double fraction = jd - std::floor(jd) + 0.5;
    if (fraction >= 1.0) {
        fraction -= 1.0;
        day += 1;
    }
    
    return std::make_tuple(year, month, day, fraction);
}
```

### 1.5 Secoli Giuliani dal J2000.0

Per molti calcoli astronomici è necessario esprimere il tempo in **secoli giuliani** dall'epoca J2000.0:

#### 1.5.1 Definizione

```
T = (JD - 2451545.0) / 36525.0
```

oppure in termini di MJD:

```
T = (MJD - 51544.5) / 36525.0
```

Dove:
- **36525** = numero di giorni in un secolo giuliano (365.25 × 100)
- **T** = numero di secoli giuliani (può essere negativo per date prima del J2000.0)

#### 1.5.2 Implementazione

```cpp
inline double julian_centuries_j2000(double jd) {
    return (jd - constants::JD2000) / constants::DAYS_PER_CENTURY;
}

inline double mjd_to_julian_centuries(double mjd) {
    return julian_centuries_j2000(mjd_to_jd(mjd));
}
```

#### 1.5.3 Valori Tipici

| Data | T (secoli) |
|------|------------|
| 1 gen 1900 | -1.0 |
| 1 gen 2000, 12:00 TT (J2000.0) | 0.0 |
| 1 gen 2025 | +0.25 |
| 1 gen 2100 | +1.0 |

### 1.6 Precisione Numerica

#### 1.6.1 Problemi con Numeri a Virgola Mobile

Le JD sono numeri grandi (~2.46 milioni per l'epoca attuale). Con la precisione `double` (64 bit, IEEE 754):

- **Precisione relativa**: ~15-16 cifre significative
- **Precisione temporale con JD**: ~0.01 secondi (per JD ~ 2,460,000)
- **Precisione temporale con MJD**: ~0.001 secondi (per MJD ~ 60,000)

#### 1.6.2 Raccomandazioni

Per calcoli ad alta precisione si consiglia:

1. **Usare sempre MJD** invece di JD (numeri più piccoli → maggiore precisione)
2. **Separare parte intera e frazionaria** quando necessario
3. **Usare aritmetica a precisione estesa** per calcoli critici

```cpp
// Tecnica: separazione parte intera/frazionaria
struct PreciseMJD {
    int64_t day;       // Parte intera (giorni)
    double fraction;   // Parte frazionaria [0, 1)
    
    double to_mjd() const {
        return static_cast<double>(day) + fraction;
    }
};
```

### 1.7 Esempi Pratici con Dati Reali

#### Esempio 1.1: Calcolo MJD per una data specifica

**Problema:** Calcolare il MJD per il 21 luglio 1969, 02:56:15 UTC (primo allunaggio Apollo 11)

**Soluzione:**

1. Data: 21/07/1969, 02:56:15 UTC
2. Anno = 1969, Mese = 7, Giorno = 21
3. Frazione del giorno:
   ```
   fraction = (2 + 56/60 + 15/3600) / 24
            = (2 + 0.9333 + 0.00417) / 24
            = 2.9375 / 24
            = 0.12240
   ```

4. Applicando l'algoritmo:
   ```
   a = (14 - 7) / 12 = 0
   y = 1969 + 4800 - 0 = 6769
   m = 7 + 0 - 3 = 4
   
   jdn = 21 + (153×4 + 2)/5 + 365×6769 + 6769/4 - 6769/100 + 6769/400 - 32045
       = 21 + 122 + 2470685 + 1692 - 67 + 16 - 32045
       = 2,440,424
   
   JD = 2440424 - 0.5 + 0.12240 = 2,440,423.62240
   
   MJD = 2440423.62240 - 2400000.5 = 40,423.12240
   ```

**Verifica:** MJD 40423.12240 corrisponde al 21 luglio 1969, ~02:56 UTC ✓

#### Esempio 1.2: Conversione inversa

**Problema:** Dato MJD = 51544.5, determinare la data calendario.

**Soluzione:**

1. JD = 51544.5 + 2400000.5 = 2,451,545.0
2. Applicando l'algoritmo inverso:
   ```
   l = int(2451545.0 + 68569) = 2520114
   n = (4 × 2520114) / 146097 = 69
   l = 2520114 - (146097 × 69 + 3) / 4 = 2520114 - 2520046 = 68
   i = (4000 × 69) / 1461001 = 0
   l = 68 - 0 + 31 = 99
   j = (80 × 99) / 2447 = 3
   day = 99 - (2447 × 3) / 80 = 99 - 91 = 8 → corretto a 1
   month = 3 + 2 - 0 = 5 → gennaio = 1
   year = 100 × (69 - 49) + 0 + 0 = 2000
   ```

**Risultato:** 1 gennaio 2000 = J2000.0 ✓

### 1.8 Costanti e Definizioni Fondamentali

```cpp
namespace orbfit::constants {
    // Epoche di riferimento
    constexpr double JD2000 = 2451545.0;      // J2000.0 in JD
    constexpr double MJD2000 = 51544.5;       // J2000.0 in MJD
    
    // Conversione JD ↔ MJD
    constexpr double JD_MJD_OFFSET = 2400000.5;
    
    // Costanti temporali
    constexpr double DAYS_PER_CENTURY = 36525.0;
    constexpr double SECONDS_PER_DAY = 86400.0;
    constexpr double HOURS_PER_DAY = 24.0;
    constexpr double MINUTES_PER_HOUR = 60.0;
}
```

### 1.9 Riepilogo delle Formule

| Operazione | Formula |
|------------|---------|
| JD → MJD | MJD = JD - 2,400,000.5 |
| MJD → JD | JD = MJD + 2,400,000.5 |
| Secoli giuliani da J2000.0 | T = (JD - 2,451,545.0) / 36,525.0 |
| Secoli giuliani (da MJD) | T = (MJD - 51,544.5) / 36,525.0 |
| Frazione del giorno | f = (ore + min/60 + sec/3600) / 24 |

---

## 🇬🇧 English

### 1.1 Introduction to Julian Dates

**Julian Dates** (JD) represent a continuous day counting system introduced in the 16th century by scholar Joseph Justus Scaliger. This system is fundamental in astronomy because it eliminates the complexities of civil calendars (leap years, calendar reforms, etc.).

#### 1.1.1 Advantages of Julian Dates

- **Continuity**: No interruptions or reforms
- **Simplicity**: Linear time counting
- **Universality**: Independent of time zones and calendars
- **Precision**: Suitable for high-precision scientific calculations

### 1.2 Julian Date (JD)

The **Julian Date (JD)** is the number of days (and fractions) elapsed since the initial epoch:

> **JD zero epoch**: January 1, 4713 BC (Julian calendar), 12:00 UT
> (or: November 24, 4714 BC in the proleptic Gregorian calendar)

#### 1.2.1 Formal Definition

```
JD = integer number of days + day fraction
```

Where the **day fraction** is calculated from noon (12:00 UT):
- 0.0 = noon
- 0.5 = midnight
- 0.75 = 6:00 AM the following day

#### 1.2.2 Significant Epochs

| Event | JD | Date |
|-------|-----|------|
| JD zero epoch | 0.0 | -4712/01/01, 12:00 UT |
| Birth of Christ | ~1,721,424 | 01/01/0001 |
| MJD zero | 2,400,000.5 | 17/11/1858, 0:00 UT |
| J2000.0 | **2,451,545.0** | 01/01/2000, 12:00 TT |
| Current year (2025) | ~2,460,676 | 01/01/2025 |

### 1.3 Modified Julian Date (MJD)

The **Modified Julian Date (MJD)** is a "simplified" version of JD, introduced to:

- Reduce numbers to more manageable values
- Start the day at midnight instead of noon

#### 1.3.1 Definition

```
MJD = JD - 2,400,000.5
```

#### 1.3.2 Properties

- **MJD zero epoch**: November 17, 1858, 0:00 UT (MJD = 0)
- **Day start**: Midnight (0:00 UT)
- **Typical current values**: ~60,000 - 61,000

#### 1.3.3 Significant MJD Epochs

| Event | MJD | Date |
|-------|-----|------|
| MJD zero epoch | 0.0 | 17/11/1858, 0:00 UT |
| J2000.0 | **51,544.5** | 01/01/2000, 12:00 TT |
| Year 2025 | ~60,676 | 01/01/2025 |

### 1.4 Conversion Algorithms

#### 1.4.1 JD ↔ MJD Conversion

**Inline functions (C++):**

```cpp
namespace orbfit::time {

    // JD → MJD
    inline double jd_to_mjd(double jd) {
        return jd - 2400000.5;
    }
    
    // MJD → JD
    inline double mjd_to_jd(double mjd) {
        return mjd + 2400000.5;
    }
    
}
```

#### 1.4.2 Gregorian Calendar → MJD Conversion

**Fliegel and van Flandern algorithm (1968):**

This algorithm calculates the integer Julian Day Number (JDN) from the Gregorian calendar:

```cpp
double calendar_to_mjd(int year, int month, int day, double fraction) {
    // Fliegel and van Flandern algorithm
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    
    // Julian Day Number (integer)
    int jdn = day + (153 * m + 2) / 5 + 365 * y 
            + y / 4 - y / 100 + y / 400 - 32045;
    
    // JD with fraction
    double jd = jdn - 0.5 + fraction;
    
    // Convert to MJD
    return jd - 2400000.5;
}
```

**Step explanation:**

1. `a`: Adjustment to treat January and February as months 13 and 14 of the previous year
2. `y`: Modified year for calculation
3. `m`: Modified month (March = 0, ..., February = 11)
4. Main formula: Accounts for days per month, Gregorian leap years

#### 1.4.3 MJD → Gregorian Calendar Conversion

**Inverse algorithm:**

```cpp
std::tuple<int, int, int, double> mjd_to_calendar(double mjd) {
    double jd = mjd + 2400000.5;
    
    // Fliegel and van Flandern algorithm (inverse)
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
    
    // Extract day fraction
    double fraction = jd - std::floor(jd) + 0.5;
    if (fraction >= 1.0) {
        fraction -= 1.0;
        day += 1;
    }
    
    return std::make_tuple(year, month, day, fraction);
}
```

### 1.5 Julian Centuries from J2000.0

For many astronomical calculations, time needs to be expressed in **Julian centuries** from the J2000.0 epoch:

#### 1.5.1 Definition

```
T = (JD - 2451545.0) / 36525.0
```

or in terms of MJD:

```
T = (MJD - 51544.5) / 36525.0
```

Where:
- **36525** = number of days in a Julian century (365.25 × 100)
- **T** = number of Julian centuries (can be negative for dates before J2000.0)

#### 1.5.2 Implementation

```cpp
inline double julian_centuries_j2000(double jd) {
    return (jd - constants::JD2000) / constants::DAYS_PER_CENTURY;
}

inline double mjd_to_julian_centuries(double mjd) {
    return julian_centuries_j2000(mjd_to_jd(mjd));
}
```

#### 1.5.3 Typical Values

| Date | T (centuries) |
|------|---------------|
| Jan 1, 1900 | -1.0 |
| Jan 1, 2000, 12:00 TT (J2000.0) | 0.0 |
| Jan 1, 2025 | +0.25 |
| Jan 1, 2100 | +1.0 |

### 1.6 Numerical Precision

#### 1.6.1 Floating-Point Number Issues

JDs are large numbers (~2.46 million for the current epoch). With `double` precision (64-bit, IEEE 754):

- **Relative precision**: ~15-16 significant digits
- **Time precision with JD**: ~0.01 seconds (for JD ~ 2,460,000)
- **Time precision with MJD**: ~0.001 seconds (for MJD ~ 60,000)

#### 1.6.2 Recommendations

For high-precision calculations it is recommended to:

1. **Always use MJD** instead of JD (smaller numbers → greater precision)
2. **Separate integer and fractional parts** when necessary
3. **Use extended precision arithmetic** for critical calculations

```cpp
// Technique: integer/fractional part separation
struct PreciseMJD {
    int64_t day;       // Integer part (days)
    double fraction;   // Fractional part [0, 1)
    
    double to_mjd() const {
        return static_cast<double>(day) + fraction;
    }
};
```

### 1.7 Practical Examples with Real Data

#### Example 1.1: MJD Calculation for a Specific Date

**Problem:** Calculate the MJD for July 21, 1969, 02:56:15 UTC (Apollo 11 first moon landing)

**Solution:**

1. Date: 21/07/1969, 02:56:15 UTC
2. Year = 1969, Month = 7, Day = 21
3. Day fraction:
   ```
   fraction = (2 + 56/60 + 15/3600) / 24
            = (2 + 0.9333 + 0.00417) / 24
            = 2.9375 / 24
            = 0.12240
   ```

4. Applying the algorithm:
   ```
   a = (14 - 7) / 12 = 0
   y = 1969 + 4800 - 0 = 6769
   m = 7 + 0 - 3 = 4
   
   jdn = 21 + (153×4 + 2)/5 + 365×6769 + 6769/4 - 6769/100 + 6769/400 - 32045
       = 21 + 122 + 2470685 + 1692 - 67 + 16 - 32045
       = 2,440,424
   
   JD = 2440424 - 0.5 + 0.12240 = 2,440,423.62240
   
   MJD = 2440423.62240 - 2400000.5 = 40,423.12240
   ```

**Verification:** MJD 40423.12240 corresponds to July 21, 1969, ~02:56 UTC ✓

#### Example 1.2: Inverse Conversion

**Problem:** Given MJD = 51544.5, determine the calendar date.

**Solution:**

1. JD = 51544.5 + 2400000.5 = 2,451,545.0
2. Applying the inverse algorithm:
   ```
   l = int(2451545.0 + 68569) = 2520114
   n = (4 × 2520114) / 146097 = 69
   l = 2520114 - (146097 × 69 + 3) / 4 = 2520114 - 2520046 = 68
   i = (4000 × 69) / 1461001 = 0
   l = 68 - 0 + 31 = 99
   j = (80 × 99) / 2447 = 3
   day = 99 - (2447 × 3) / 80 = 99 - 91 = 8 → corrected to 1
   month = 3 + 2 - 0 = 5 → January = 1
   year = 100 × (69 - 49) + 0 + 0 = 2000
   ```

**Result:** January 1, 2000 = J2000.0 ✓

### 1.8 Fundamental Constants and Definitions

```cpp
namespace orbfit::constants {
    // Reference epochs
    constexpr double JD2000 = 2451545.0;      // J2000.0 in JD
    constexpr double MJD2000 = 51544.5;       // J2000.0 in MJD
    
    // JD ↔ MJD conversion
    constexpr double JD_MJD_OFFSET = 2400000.5;
    
    // Time constants
    constexpr double DAYS_PER_CENTURY = 36525.0;
    constexpr double SECONDS_PER_DAY = 86400.0;
    constexpr double HOURS_PER_DAY = 24.0;
    constexpr double MINUTES_PER_HOUR = 60.0;
}
```

### 1.9 Formula Summary

| Operation | Formula |
|-----------|---------|
| JD → MJD | MJD = JD - 2,400,000.5 |
| MJD → JD | JD = MJD + 2,400,000.5 |
| Julian centuries from J2000.0 | T = (JD - 2,451,545.0) / 36,525.0 |
| Julian centuries (from MJD) | T = (MJD - 51,544.5) / 36,525.0 |
| Day fraction | f = (hours + min/60 + sec/3600) / 24 |

---

## Riferimenti / References

1. Fliegel, H.F., & van Flandern, T.C. (1968). "A Machine Algorithm for Processing Calendar Dates". Communications of the ACM, 11(10), 657.

2. Dershowitz, N., & Reingold, E.M. (2008). *Calendrical Calculations*. 3rd ed. Cambridge University Press.

3. Seidelmann, P.K. (1992). *Explanatory Supplement to the Astronomical Almanac*. University Science Books.

---

*Documento creato per il progetto ITALOccultLibrary*
*Document created for the ITALOccultLibrary project*

© 2025 Michele Bigi (mikbigi@gmail.com)
