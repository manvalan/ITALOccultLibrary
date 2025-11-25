# PARTE 0: Introduzione Generale
# PART 0: General Introduction

---

## 🇮🇹 Italiano

### 0.1 Prefazione

Benvenuti nel manuale scientifico **ITALOccultCalc Astrodinamica**, una guida completa ai principi matematici e fisici che governano il calcolo delle effemeridi astronomiche e la determinazione orbitale.

Questo manuale accompagna la libreria software **ITALOccultLibrary** (OrbFit C++), un porting moderno in C++17 del classico software di determinazione orbitale OrbFit, originariamente sviluppato in Fortran dal team del Prof. Andrea Milani presso l'Università di Pisa.

### 0.2 Obiettivi del Manuale

Gli obiettivi principali di questo manuale sono:

1. **Fornire basi teoriche solide**: Spiegare i concetti fondamentali di astrodinamica in modo chiaro e rigoroso
2. **Documentare gli algoritmi**: Presentare le formule matematiche implementate nella libreria
3. **Offrire esempi pratici**: Dimostrare l'applicazione dei concetti con dati reali e verificabili
4. **Servire come riferimento**: Essere una risorsa utile per studenti, ricercatori e sviluppatori

### 0.3 Struttura del Manuale

Il manuale è organizzato nelle seguenti **PARTI** principali:

| Parte | Titolo | Contenuto |
|-------|--------|-----------|
| **0** | Introduzione Generale | Panoramica e obiettivi |
| **1** | Sistemi di Coordinate | Principi base, piani di riferimento, trasformazioni |
| **2** | Tempo e Scale Temporali | JD, MJD, UTC, TAI, TT, TDB |

Ogni parte è suddivisa in **Capitoli**, e ogni capitolo può contenere più **Paragrafi** o **Sezioni**.

### 0.4 Convenzioni Utilizzate

#### 0.4.1 Notazione Matematica

Nel manuale adottiamo le seguenti convenzioni:

- **Vettori**: Indicati in grassetto (es. **r**, **v**) o con freccia (es. r⃗)
- **Matrici**: Lettere maiuscole in grassetto (es. **R**, **M**)
- **Scalari**: Lettere normali (es. *t*, *θ*, *a*)
- **Unità di misura**: Sistema SI, con eccezioni astronomiche standard (AU, gradi, arcsec)

#### 0.4.2 Costanti Astronomiche

Le costanti astronomiche utilizzate seguono gli standard IAU (International Astronomical Union):

| Costante | Simbolo | Valore | Unità |
|----------|---------|--------|-------|
| Unità Astronomica | AU | 149,597,870.700 | km |
| Velocità della Luce | c | 299,792.458 | km/s |
| Costante Gravitazionale Solare | GM☉ | 1.32712440018 × 10¹¹ | km³/s² |
| Obliquità dell'eclittica (J2000) | ε₀ | 23.439291 | gradi |
| Raggio Equatoriale Terrestre | R⊕ | 6,378.137 | km |

#### 0.4.3 Epoche di Riferimento

- **J2000.0**: 1 gennaio 2000, 12:00:00 TT (JD 2451545.0)
- **MJD**: Modified Julian Date = JD − 2,400,000.5
- **Epoca standard**: J2000.0 (equinozio medio)

### 0.5 Software di Riferimento

Questo manuale documenta la libreria **OrbFit C++**, le cui caratteristiche principali sono:

- **Linguaggio**: C++17 moderno
- **Librerie**: Eigen3 (algebra lineare), Boost (utilities)
- **Build System**: CMake 3.15+
- **Test Framework**: Google Test
- **Documentazione API**: Doxygen

La struttura del codice segue una architettura modulare:

```
orbfit-cpp/
├── include/orbfit/    # Header pubblici
│   ├── core/          # Tipi e costanti
│   ├── time/          # Scale temporali
│   ├── coordinates/   # Sistemi di riferimento
│   ├── ephemeris/     # Effemeridi
│   └── propagation/   # Propagazione orbitale
└── src/               # Implementazioni
```

### 0.6 Prerequisiti

Per comprendere appieno questo manuale, si consiglia familiarità con:

- **Matematica**: Algebra lineare, trigonometria sferica, calcolo differenziale
- **Fisica**: Meccanica classica, leggi di Keplero, gravitazione
- **Astronomia**: Coordinate celesti, sistemi di tempo, meccanica celeste di base
- **Programmazione**: C++ moderno (per gli esempi di codice)

---

## 🇬🇧 English

### 0.1 Preface

Welcome to the **ITALOccultCalc Astrodynamics** scientific manual, a comprehensive guide to the mathematical and physical principles governing astronomical ephemeris calculation and orbit determination.

This manual accompanies the **ITALOccultLibrary** software library (OrbFit C++), a modern C++17 port of the classic OrbFit orbit determination software, originally developed in Fortran by Prof. Andrea Milani's team at the University of Pisa.

### 0.2 Manual Objectives

The main objectives of this manual are:

1. **Provide solid theoretical foundations**: Explain fundamental astrodynamics concepts clearly and rigorously
2. **Document algorithms**: Present the mathematical formulas implemented in the library
3. **Offer practical examples**: Demonstrate concept application with real, verifiable data
4. **Serve as reference**: Be a useful resource for students, researchers, and developers

### 0.3 Manual Structure

The manual is organized into the following main **PARTS**:

| Part | Title | Content |
|------|-------|---------|
| **0** | General Introduction | Overview and objectives |
| **1** | Coordinate Systems | Basic principles, reference planes, transformations |
| **2** | Time and Time Scales | JD, MJD, UTC, TAI, TT, TDB |

Each part is divided into **Chapters**, and each chapter may contain multiple **Paragraphs** or **Sections**.

### 0.4 Conventions Used

#### 0.4.1 Mathematical Notation

In this manual we adopt the following conventions:

- **Vectors**: Indicated in bold (e.g., **r**, **v**) or with arrow (e.g., r⃗)
- **Matrices**: Uppercase bold letters (e.g., **R**, **M**)
- **Scalars**: Normal letters (e.g., *t*, *θ*, *a*)
- **Units**: SI system, with standard astronomical exceptions (AU, degrees, arcsec)

#### 0.4.2 Astronomical Constants

The astronomical constants used follow IAU (International Astronomical Union) standards:

| Constant | Symbol | Value | Unit |
|----------|--------|-------|------|
| Astronomical Unit | AU | 149,597,870.700 | km |
| Speed of Light | c | 299,792.458 | km/s |
| Solar Gravitational Constant | GM☉ | 1.32712440018 × 10¹¹ | km³/s² |
| Obliquity of the ecliptic (J2000) | ε₀ | 23.439291 | degrees |
| Earth Equatorial Radius | R⊕ | 6,378.137 | km |

#### 0.4.3 Reference Epochs

- **J2000.0**: January 1, 2000, 12:00:00 TT (JD 2451545.0)
- **MJD**: Modified Julian Date = JD − 2,400,000.5
- **Standard epoch**: J2000.0 (mean equinox)

### 0.5 Reference Software

This manual documents the **OrbFit C++** library, whose main characteristics are:

- **Language**: Modern C++17
- **Libraries**: Eigen3 (linear algebra), Boost (utilities)
- **Build System**: CMake 3.15+
- **Test Framework**: Google Test
- **API Documentation**: Doxygen

The code structure follows a modular architecture:

```
orbfit-cpp/
├── include/orbfit/    # Public headers
│   ├── core/          # Types and constants
│   ├── time/          # Time scales
│   ├── coordinates/   # Reference systems
│   ├── ephemeris/     # Ephemerides
│   └── propagation/   # Orbit propagation
└── src/               # Implementations
```

### 0.6 Prerequisites

To fully understand this manual, familiarity with the following is recommended:

- **Mathematics**: Linear algebra, spherical trigonometry, differential calculus
- **Physics**: Classical mechanics, Kepler's laws, gravitation
- **Astronomy**: Celestial coordinates, time systems, basic celestial mechanics
- **Programming**: Modern C++ (for code examples)

---

## Bibliografia / References

1. Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed. Microcosm Press.

2. Montenbruck, O., & Gill, E. (2012). *Satellite Orbits: Models, Methods and Applications*. Springer.

3. Urban, S. E., & Seidelmann, P. K. (2013). *Explanatory Supplement to the Astronomical Almanac*. 3rd ed. University Science Books.

4. Milani, A., & Gronchi, G. F. (2010). *Theory of Orbit Determination*. Cambridge University Press.

5. IERS Conventions (2010). IERS Technical Note No. 36.

---

*Documento creato per il progetto ITALOccultLibrary*
*Document created for the ITALOccultLibrary project*

© 2025 Michele Bigi (mikbigi@gmail.com)
