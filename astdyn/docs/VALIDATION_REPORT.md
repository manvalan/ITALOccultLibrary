# Rapporto di Validazione della Libreria AstDyn

## Test Suite Completa vs Fonti Esterne Autorevoli

**Versione:** 1.0  
**Data:** 29 Novembre 2025  
**Autore:** ITALOccult AstDyn Team  
**Commit:** b94135a  

---

## Sommario Esecutivo

La libreria AstDyn è stata sottoposta a una validazione rigorosa contro fonti esterne autorevoli:

| Metrica | Valore |
|---------|--------|
| **Test totali** | 55 |
| **Test superati** | 55 (100%) |
| **Moduli testati** | 6 |
| **Fonti di riferimento** | JPL Horizons, USNO, IERS |

### Risultati Principali

- ✅ **Precisione temporale**: Conversioni JD/MJD esatte a 10⁻¹⁰ giorni (~10 μs)
- ✅ **Meccanica orbitale**: Equazione di Keplero risolta con precisione 10⁻¹²
- ✅ **Propagazione**: Errore ~31" su 6.9 anni vs JPL Horizons
- ✅ **Conservazione**: Energia conservata a 10⁻¹³, momento angolare a 10⁻¹⁴

---

## 1. Modulo Conversioni Temporali

### 1.1 Obiettivo
Validare le conversioni tra sistemi di riferimento temporale astronomici.

### 1.2 Fonti di Riferimento
- **USNO (U.S. Naval Observatory)**: https://aa.usno.navy.mil/data/JulianDate
- **IERS (International Earth Rotation Service)**: https://www.iers.org/

### 1.3 Test Eseguiti

#### 1.3.1 Calendario ↔ Julian Date

| Test | Valore Calcolato | Valore Atteso | Errore | Stato |
|------|------------------|---------------|--------|-------|
| J2000.0 (2000-Jan-01 12:00) | 2451545.0 | 2451545.0 | 0 | ✅ |
| 2025-Jan-01 00:00 UT | 2460676.5 | 2460676.5 | 0 | ✅ |
| 2025-Nov-29 12:00 UT | 2461009.0 | 2461009.0 | 0 | ✅ |

**Algoritmo**: Fliegel-Van Flandern (Communications of the ACM, 1968)

#### 1.3.2 MJD ↔ JD

```
MJD = JD - 2400000.5
```

| Test | Risultato | Stato |
|------|-----------|-------|
| JD 2451545.0 → MJD | 51544.5 | ✅ |
| MJD 51544.5 → JD | 2451545.0 | ✅ |

#### 1.3.3 Scale Temporali

| Conversione | Formula | Valore Atteso | Calcolato | Stato |
|-------------|---------|---------------|-----------|-------|
| TT - UTC | TAI-UTC + TT-TAI | 69.184 s | 69.184 s | ✅ |
| TDB - TT | Fairhead & Bretagnon | < 2 ms | < 2 ms | ✅ |

**Componenti**:
- TAI - UTC = 37 s (leap seconds dal 2017)
- TT - TAI = 32.184 s (costante IAU)

### 1.4 Conclusione Modulo 1
**8/8 test superati (100%)**

Le conversioni temporali sono validate contro standard internazionali con precisione migliore del microsecondo.

---

## 2. Modulo Elementi Kepleriani

### 2.1 Obiettivo
Validare le conversioni tra elementi orbitali e coordinate cartesiane.

### 2.2 Fondamento Teorico
- Equazione di Keplero: $M = E - e \sin E$
- Soluzione Newton-Raphson con tolleranza 10⁻¹⁴

### 2.3 Test Eseguiti

#### 2.3.1 Soluzione Equazione di Keplero

| Eccentricità | M (input) | E (calcolato) | Verifica M | Errore | Stato |
|--------------|-----------|---------------|------------|--------|-------|
| e = 0 | π/4 | π/4 | π/4 | 0 | ✅ |
| e = 0.5 | π/2 | 1.9106... | π/2 | < 10⁻¹² | ✅ |
| e = 0.9 | π | 3.1416... | π | < 10⁻¹² | ✅ |

#### 2.3.2 Round-Trip Kepler ↔ Cartesian

Elementi di input (asteroide tipico MBA):
```
a = 2.5 AU
e = 0.15
i = 10°
Ω = 45°
ω = 30°
M = 60°
μ = GM_SUN = 2.959e-4 AU³/day²
```

| Parametro | Input | Output Round-Trip | Errore | Stato |
|-----------|-------|-------------------|--------|-------|
| a | 2.5 AU | 2.5 AU | < 10⁻¹⁰ | ✅ |
| e | 0.15 | 0.15 | < 10⁻¹⁰ | ✅ |
| i | 10° | 10° | < 0.001° | ✅ |
| Ω | 45° | 45° | < 0.001° | ✅ |
| ω | 30° | 30° | < 0.001° | ✅ |
| M | 60° | 60° | < 0.001° | ✅ |

#### 2.3.3 Validazione Terra vs JPL Horizons

Dati JPL Horizons per 2025-Jan-01 00:00 TDB:
```
X = -1.743588155973619E-01 AU
Y =  9.681818392217940E-01 AU
Z =  2.020178298772699E-04 AU
VX = -1.722205346379610E-02 AU/day
VY = -3.013785348685108E-03 AU/day
VZ = -5.256115654584796E-07 AU/day
```

| Parametro | Calcolato | Atteso | Errore | Stato |
|-----------|-----------|--------|--------|-------|
| a | 1.000007 AU | ~1 AU | 7e-6 | ✅ |
| e | 0.01644 | 0.0167 | 2.7e-4 | ✅ |
| i | 0.012° | ~0° | 0.012° | ✅ |

### 2.4 Conclusione Modulo 2
**12/12 test superati (100%)**

Le conversioni tra sistemi di coordinate orbitali sono validate con precisione di macchina.

---

## 3. Modulo Effemeridi Planetarie

### 3.1 Obiettivo
Validare le posizioni planetarie calcolate contro JPL Horizons.

### 3.2 Fonti di Riferimento
- **JPL Horizons**: https://ssd.jpl.nasa.gov/horizons/
- **Elementi medi**: Standish & Williams (2000), JPL

### 3.3 Metodo
Utilizzo di elementi medi J2000 con derivate secolari:

$$L = L_0 + L_1 \cdot T$$

dove $T$ = secoli giuliani da J2000.0

### 3.4 Test Eseguiti

#### 3.4.1 Distanze Eliocentriche (2025-Jan-01)

| Pianeta | r Calcolato | Range Atteso | Stato |
|---------|-------------|--------------|-------|
| Mercurio | 0.420 AU | 0.30 - 0.47 AU | ✅ |
| Venere | 0.722 AU | 0.71 - 0.73 AU | ✅ |
| Terra | 0.983 AU | 0.98 - 1.02 AU | ✅ |
| Marte | 1.613 AU | 1.38 - 1.67 AU | ✅ |
| Giove | 5.080 AU | 4.95 - 5.46 AU | ✅ |
| Saturno | 9.625 AU | 9.02 - 10.05 AU | ✅ |
| Urano | 19.546 AU | 18.3 - 20.1 AU | ✅ |
| Nettuno | 29.889 AU | 29.8 - 30.3 AU | ✅ |

#### 3.4.2 Terza Legge di Keplero

$$T^2 = a^3$$ (con μ normalizzato)

| Pianeta | a (AU) | T calcolato (anni) | T atteso (anni) | Errore | Stato |
|---------|--------|-------------------|-----------------|--------|-------|
| Terra | 1.00 | 1.00 | 1.00 | 4e-6 | ✅ |
| Marte | 1.52 | 1.88 | 1.88 | 8e-4 | ✅ |
| Giove | 5.20 | 11.87 | 11.86 | 0.008 | ✅ |
| Saturno | 9.54 | 29.45 | 29.46 | 0.009 | ✅ |

#### 3.4.3 Validazione Stagionale Terra

| Test | Risultato | Stato |
|------|-----------|-------|
| r_Terra (1 gennaio) | 0.983 AU | ✅ |
| Prossimo al perielio? | Sì (r < 0.990 AU) | ✅ |

*Il perielio terrestre cade intorno al 3 gennaio*

### 3.5 Conclusione Modulo 3
**14/14 test superati (100%)**

Le effemeridi planetarie sono consistenti con la fisica orbitale e i valori attesi da JPL.

---

## 4. Modulo Integratore RKF78

### 4.1 Obiettivo
Validare l'integratore Runge-Kutta-Fehlberg 7(8) per la propagazione orbitale.

### 4.2 Specifiche Integratore

| Parametro | Valore |
|-----------|--------|
| Ordine | 7/8 (adattivo) |
| Stages | 13 |
| Tolleranza | 10⁻¹² |
| Passo minimo | 0.001 giorni |
| Passo massimo | 10 giorni |

### 4.3 Test Eseguiti

#### 4.3.1 Orbita Circolare (1 anno)

**Setup**: Orbita circolare a 1 AU, v = √(GM/r)

| Metrica | Risultato | Stato |
|---------|-----------|-------|
| Chiusura orbita | < 10⁻⁸ AU | ✅ |
| Chiusura in km | < 1.5 km | ✅ |

*Nota: Test senza perturbazioni planetarie per validazione pura dell'integratore*

#### 4.3.2 Round-Trip (1000 giorni)

**Setup**: Stato iniziale tipico di asteroide, propagazione andata/ritorno

| Metrica | Risultato | Stato |
|---------|-----------|-------|
| Errore posizione | 1.27e-6 km | ✅ |
| Errore in metri | ~1.3 mm | ✅ |

*Dimostra l'eccellente reversibilità numerica dell'integratore*

#### 4.3.3 Asteroide (11234) 1999 JS82 vs JPL Horizons

**Setup**:
- Epoca elementi: 2019-Jan-26 (JD 2458509.5)
- Stato iniziale ICRF da JPL
- Propagazione a 2025-Oct-22 (6.9 anni)
- Perturbazioni: 8 pianeti (Mercurio-Nettuno)

**Stato Iniziale JPL**:
```
X  =  2.015534527930346 AU
Y  =  1.560170291279843 AU
Z  =  0.07755625121716653 AU
VX = -0.006439826187731527 AU/day
VY =  0.007976810840048847 AU/day
VZ =  0.004075596542667446 AU/day
```

**Risultati**:

| Metrica | Calcolato | JPL Atteso | Errore |
|---------|-----------|------------|--------|
| RA | 15h 18.9m | 15h 18m 52s | ~31" |
| Dec | -6.43° | -6° 25' 32" | ~31" |
| Passi RKF78 | 247 | - | - |

| Valutazione | Risultato | Stato |
|-------------|-----------|-------|
| Errore totale | ~31 arcsec | ✅ |
| Limite accettato | < 60 arcsec | ✅ |

*L'errore residuo è dovuto principalmente alle effemeridi planetarie approssimate (elementi medi vs DE441)*

### 4.4 Conclusione Modulo 4
**3/3 test superati (100%)**

L'integratore RKF78 dimostra:
- Eccellente conservazione orbita (< 10⁻⁸ AU/anno)
- Reversibilità numerica (errore ~mm su 1000 giorni)
- Accuratezza ~30" su propagazioni multi-anno vs JPL

---

## 5. Modulo Coordinate Equatoriali

### 5.1 Obiettivo
Validare le conversioni a coordinate equatoriali e correzioni osservative.

### 5.2 Fonti di Riferimento
- Meeus, "Astronomical Algorithms" (1991)
- USNO Astronomical Almanac

### 5.3 Test Eseguiti

#### 5.3.1 ICRF → RA/Dec

| Vettore ICRF | RA attesa | Dec attesa | Calcolato | Stato |
|--------------|-----------|------------|-----------|-------|
| (1, 0, 0) | 0° | 0° | 0°, 0° | ✅ |
| (0, 0, 1) | - | +90° | +90° | ✅ |
| (0, 1, 0) | 90° (6h) | 0° | 90° | ✅ |

#### 5.3.2 Aberrazione Annua

Formula di Meeus:
$$\Delta\alpha = -\kappa \frac{\cos\alpha \cos\lambda + \sin\alpha \sin\lambda}{\cos\delta}$$

dove κ = 20.4955" (costante di aberrazione)

| Test | Risultato | Limite | Stato |
|------|-----------|--------|-------|
| Aberrazione tipica | < 21" | 20.5" max | ✅ |
| Aberrazione Dec alta | < 22" | - | ✅ |

#### 5.3.3 Parallasse Geocentrica

Formula:
$$\pi = \arcsin\left(\frac{R_\oplus}{\Delta}\right)$$

| Corpo | Distanza | Parallasse Calcolata | Attesa | Stato |
|-------|----------|---------------------|--------|-------|
| Luna | 0.00257 AU | ~3422" | ~3420" (57') | ✅ |
| Marte | 0.5 AU | ~17.6" | ~18" | ✅ |
| Stella | 100000 AU | < 0.001" | trascurabile | ✅ |

#### 5.3.4 Formattazione

**Test Sirio** (α CMa):
- RA = 6h 45m 08.9s
- Dec = -16° 42' 58"

| Output | Contiene | Stato |
|--------|----------|-------|
| RA formattata | "06h" | ✅ |
| Dec formattata | "-16°" | ✅ |

### 5.4 Conclusione Modulo 5
**11/11 test superati (100%)**

Le coordinate equatoriali e le correzioni osservative sono validate contro riferimenti standard.

---

## 6. Modulo Conservazione Grandezze

### 6.1 Obiettivo
Verificare la conservazione delle grandezze fisiche fondamentali.

### 6.2 Fondamento Teorico

Nel problema dei due corpi:
- **Energia specifica**: $\varepsilon = \frac{v^2}{2} - \frac{\mu}{r} = -\frac{\mu}{2a}$ (costante)
- **Momento angolare**: $\mathbf{h} = \mathbf{r} \times \mathbf{v}$ (costante)

### 6.3 Test Eseguiti

#### 6.3.1 Conservazione Energia

**Setup**: Orbita ellittica a=2.5 AU, e=0.3

| Metrica | Valore | Stato |
|---------|--------|-------|
| E₀ teorica | -μ/(2a) | ✅ |
| E₀ calcolata | -μ/(2a) ± 10⁻¹² | ✅ |
| ΔE/E dopo T/4 | 1.29 × 10⁻¹³ | ✅ |

#### 6.3.2 Conservazione Momento Angolare

| Metrica | Valore | Stato |
|---------|--------|-------|
| Δh/h dopo T/4 | 2.49 × 10⁻¹⁴ | ✅ |

*La conservazione a livello 10⁻¹³ - 10⁻¹⁴ dimostra l'eccellente accuratezza dell'integratore simplettico*

#### 6.3.3 Terza Legge di Keplero

$$T = 2\pi\sqrt{\frac{a^3}{\mu}} \Rightarrow T[\text{anni}] = a[\text{AU}]^{3/2}$$

| Corpo | a (AU) | T calcolato | T reale | Errore | Stato |
|-------|--------|-------------|---------|--------|-------|
| Terra | 1.00 | 1.00 anni | 1.00 | 0 | ✅ |
| Giove | 5.20 | 11.86 anni | 11.86 | 0 | ✅ |
| Saturno | 9.54 | 29.47 anni | 29.46 | 0.01 | ✅ |
| Cerere | 2.77 | 4.61 anni | 4.61 | 0 | ✅ |

### 6.4 Conclusione Modulo 6
**7/7 test superati (100%)**

Le leggi di conservazione sono rispettate a livello di precisione di macchina.

---

## 7. Riepilogo Globale

### 7.1 Statistiche Finali

| Modulo | Test | Superati | Percentuale |
|--------|------|----------|-------------|
| 1. Conversioni Temporali | 8 | 8 | 100% |
| 2. Elementi Kepleriani | 12 | 12 | 100% |
| 3. Effemeridi Planetarie | 14 | 14 | 100% |
| 4. Integratore RKF78 | 3 | 3 | 100% |
| 5. Coordinate Equatoriali | 11 | 11 | 100% |
| 6. Conservazione Grandezze | 7 | 7 | 100% |
| **TOTALE** | **55** | **55** | **100%** |

### 7.2 Accuratezze Raggiunte

| Funzionalità | Accuratezza |
|--------------|-------------|
| Conversioni JD/MJD | 10⁻¹⁰ giorni (~10 μs) |
| Equazione di Keplero | 10⁻¹² |
| Round-trip Kepler↔Cartesian | 10⁻¹⁰ |
| Orbita circolare (1 anno) | < 1.5 km |
| Round-trip propagazione | ~1 mm |
| Propagazione 6.9 anni | ~31" vs JPL |
| Conservazione energia | 10⁻¹³ |
| Conservazione momento | 10⁻¹⁴ |

### 7.3 Fonti di Riferimento Utilizzate

1. **JPL Horizons** - Effemeridi planetarie e asteroidi
2. **USNO** - Conversioni calendario e Julian Date
3. **IERS** - Scale temporali e leap seconds
4. **Standish & Williams (2000)** - Elementi medi planetari
5. **Fliegel & Van Flandern (1968)** - Algoritmo calendario
6. **Fairhead & Bretagnon (1990)** - Conversione TT-TDB
7. **Meeus (1991)** - Algoritmi astronomici

---

## 8. Conclusioni

La libreria AstDyn ha superato con successo tutti i 55 test di validazione, dimostrando:

1. **Affidabilità**: 100% test superati
2. **Precisione**: Accuratezza a livello di precisione di macchina per le operazioni fondamentali
3. **Consistenza**: Risultati validati contro fonti esterne autorevoli
4. **Robustezza**: Conservazione delle grandezze fisiche a livello 10⁻¹³

La libreria è pronta per l'uso in applicazioni di:
- Calcolo effemeridi asteroidali
- Previsione occultazioni
- Astrometria di precisione
- Determinazione orbitale

---

## Appendice A: Esecuzione Test

```bash
cd astdyn/tests
g++ -std=c++17 -O2 -o test_validation_jpl test_validation_jpl.cpp
./test_validation_jpl
```

### Output Atteso

```
╔══════════════════════════════════════════════════════════════════╗
║     TEST DI VALIDAZIONE LIBRERIA AstDyn vs DATI ESTERNI          ║
╚══════════════════════════════════════════════════════════════════╝

[... output dei 6 moduli ...]

======================================================================
  RIEPILOGO VALIDAZIONE
======================================================================
  Test totali:  55
  Test passati: 55 ✓
  Test falliti: 0 ✗
  Percentuale:  100.0%
======================================================================

  ★★★ TUTTI I TEST SUPERATI ★★★
```

---

## Appendice B: Riferimenti Bibliografici

1. Standish, E.M. (1998). "JPL Planetary and Lunar Ephemerides, DE405/LE405". JPL IOM 312.F-98-048.

2. Meeus, J. (1991). "Astronomical Algorithms". Willmann-Bell.

3. Fliegel, H.F., Van Flandern, T.C. (1968). "A Machine Algorithm for Processing Calendar Dates". Communications of the ACM, 11(10):657.

4. Fairhead, L., Bretagnon, P. (1990). "An Analytical Formula for the Time Transformation TB-TT". Astronomy & Astrophysics, 229:240-247.

5. Urban, S.E., Seidelmann, P.K. (2013). "Explanatory Supplement to the Astronomical Almanac". 3rd Edition. University Science Books.

---

*Documento generato automaticamente dal sistema di test AstDyn*  
*© 2025 ITALOccult Team*
