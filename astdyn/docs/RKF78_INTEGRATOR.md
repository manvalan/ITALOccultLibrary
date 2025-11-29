# RKF78 Integrator - Documentazione Completa

## Runge-Kutta-Fehlberg 7(8) per Propagazione Orbitale

**Autore**: AstDyS Team  
**Data**: 29 Novembre 2025  
**Versione**: 1.0

---

## 1. Panoramica

L'integratore **RKF78** (Runge-Kutta-Fehlberg di ordine 7 con stima dell'errore di ordine 8) è un metodo numerico ad alta precisione per la propagazione orbitale di corpi celesti.

### Caratteristiche Principali

| Proprietà | Valore |
|-----------|--------|
| **Ordine di propagazione** | 7 |
| **Ordine stima errore** | 8 |
| **Numero di stadi** | 13 |
| **Tipo** | Embedded, adattivo |
| **Valutazioni per passo** | 13 |
| **Riferimento** | Fehlberg (1968) NASA TR R-287 |

---

## 2. Formulazione Matematica

### 2.1 Schema Generale

Per un sistema $\dot{y} = f(t, y)$, il metodo calcola:

$$y_{n+1} = y_n + h \sum_{i=1}^{13} b_i k_i$$

dove gli stadi sono:

$$k_i = f\left(t_n + c_i h, \, y_n + h \sum_{j=1}^{i-1} a_{ij} k_j\right)$$

### 2.2 Stima dell'Errore

L'errore locale è stimato come differenza tra soluzioni di ordine 7 e 8:

$$\epsilon = h \sum_{i=1}^{13} (b_i^{(7)} - b_i^{(8)}) k_i$$

### 2.3 Controllo Adattivo del Passo

Il nuovo passo è calcolato come:

$$h_{new} = h \cdot \min\left(5, \max\left(0.1, 0.9 \cdot \left(\frac{\text{tol}}{|\epsilon|}\right)^{1/8}\right)\right)$$

---

## 3. Coefficienti di Fehlberg

### 3.1 Nodi $c_i$

```cpp
static constexpr double c[13] = {
    0.0,
    2.0/27.0,
    1.0/9.0,
    1.0/6.0,
    5.0/12.0,
    1.0/2.0,
    5.0/6.0,
    1.0/6.0,
    2.0/3.0,
    1.0/3.0,
    1.0,
    0.0,
    1.0
};
```

### 3.2 Pesi Ordine 7 ($b^{(7)}$)

```cpp
static constexpr double b7[13] = {
    41.0/840.0,
    0.0,
    0.0,
    0.0,
    0.0,
    34.0/105.0,
    9.0/35.0,
    9.0/35.0,
    9.0/280.0,
    9.0/280.0,
    41.0/840.0,
    0.0,
    0.0
};
```

### 3.3 Pesi Ordine 8 ($b^{(8)}$)

```cpp
static constexpr double b8[13] = {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    34.0/105.0,
    9.0/35.0,
    9.0/35.0,
    9.0/280.0,
    9.0/280.0,
    0.0,
    41.0/840.0,
    41.0/840.0
};
```

---

## 4. API dell'Integratore

### 4.1 Classe `RKF78Integrator`

```cpp
class RKF78Integrator {
public:
    // Costruttore con tolleranza (default 1e-12)
    explicit RKF78Integrator(double tol = 1e-12);
    
    // Imposta tolleranza
    void setTolerance(double tol);
    
    // Imposta limiti passo
    void setStepLimits(double h_min, double h_max);
    
    // Integrazione principale
    State integrate(
        const State& y0,           // Stato iniziale [r, v]
        double t0,                 // Tempo iniziale (JD)
        double t1,                 // Tempo finale (JD)
        AccelFunction accel,       // Funzione accelerazione
        IntegrationStats* stats    // Statistiche (opzionale)
    );
    
    // Singolo passo adattivo
    State step_adaptive(
        const State& y,
        double t,
        double& h,                 // Passo (in/out)
        AccelFunction accel,
        bool& accepted             // Passo accettato?
    );
};
```

### 4.2 Struttura `State`

```cpp
struct State {
    std::array<double, 3> r;  // Posizione [AU]
    std::array<double, 3> v;  // Velocità [AU/day]
    
    // Operatori aritmetici
    State operator+(const State& other) const;
    State operator*(double scalar) const;
    
    // Norma posizione
    double norm_r() const;
};
```

### 4.3 Struttura `IntegrationStats`

```cpp
struct IntegrationStats {
    int steps_accepted;     // Passi accettati
    int steps_rejected;     // Passi rifiutati
    int func_evaluations;   // Valutazioni f(t,y)
    double h_min;           // Passo minimo usato
    double h_max;           // Passo massimo usato
};
```

### 4.4 Tipo `AccelFunction`

```cpp
using AccelFunction = std::function<
    std::array<double, 3>(double t, const std::array<double, 3>& r, 
                          const std::array<double, 3>& v)
>;
```

---

## 5. Modello Dinamico

### 5.1 Equazione del Moto

L'accelerazione totale è:

$$\ddot{\mathbf{r}} = \ddot{\mathbf{r}}_\odot + \ddot{\mathbf{r}}_{\text{pianeti}} + \ddot{\mathbf{r}}_{\text{AST17}} + \ddot{\mathbf{r}}_{\text{rel}}$$

### 5.2 Termine Kepleriano (Sole)

$$\ddot{\mathbf{r}}_\odot = -\frac{GM_\odot}{r^3} \mathbf{r}$$

con $GM_\odot = 0.01720209895^2$ AU³/day² (costante di Gauss).

### 5.3 Perturbazioni Planetarie

$$\ddot{\mathbf{r}}_{\text{pianeti}} = \sum_{i=1}^{8} GM_i \left( \frac{\mathbf{r}_i - \mathbf{r}}{|\mathbf{r}_i - \mathbf{r}|^3} - \frac{\mathbf{r}_i}{r_i^3} \right)$$

**Pianeti inclusi** (effemeridi Simon et al. 1994):

| Pianeta | GM [AU³/day²] | Precisione |
|---------|---------------|------------|
| Mercurio | 4.9125e-11 | ~1" |
| Venere | 7.2435e-10 | ~1" |
| Terra-Luna | 8.9970e-10 | ~1" |
| Marte | 9.5495e-11 | ~1" |
| Giove | 2.8253e-07 | ~1" |
| Saturno | 8.4597e-08 | ~1" |
| Urano | 1.2920e-08 | ~5" |
| Nettuno | 1.5244e-08 | ~10" |

### 5.4 Perturbazioni AST17

$$\ddot{\mathbf{r}}_{\text{AST17}} = \sum_{j=1}^{16} GM_j \left( \frac{\mathbf{r}_j - \mathbf{r}}{|\mathbf{r}_j - \mathbf{r}|^3} - \frac{\mathbf{r}_j}{r_j^3} \right)$$

**Asteroidi AST17**:

| # | Nome | GM [AU³/day²] | Massa relativa |
|---|------|---------------|----------------|
| 1 | Ceres | 1.392e-13 | 1.000 |
| 2 | Pallas | 3.036e-14 | 0.218 |
| 4 | Vesta | 3.803e-14 | 0.273 |
| 10 | Hygiea | 1.231e-14 | 0.088 |
| 704 | Interamnia | 5.467e-15 | 0.039 |
| 511 | Davida | 5.467e-15 | 0.039 |
| 52 | Europa | 3.893e-15 | 0.028 |
| 15 | Eunomia | 4.565e-15 | 0.033 |
| 16 | Psyche | 3.389e-15 | 0.024 |
| 3 | Juno | 3.893e-15 | 0.028 |
| 87 | Sylvia | 2.155e-15 | 0.015 |
| 88 | Thisbe | 2.491e-15 | 0.018 |
| 31 | Euphrosyne | 2.323e-15 | 0.017 |
| 324 | Bamberga | 1.515e-15 | 0.011 |
| 451 | Patientia | 1.851e-15 | 0.013 |
| 65 | Cybele | 1.683e-15 | 0.012 |

### 5.5 Correzione Relativistica (Schwarzschild)

$$\ddot{\mathbf{r}}_{\text{rel}} = \frac{GM_\odot}{c^2 r^3} \left[ \left(4\frac{GM_\odot}{r} - v^2\right)\mathbf{r} + 4(\mathbf{r}\cdot\mathbf{v})\mathbf{v} \right]$$

con $c = 173.1446$ AU/day.

---

## 6. Esempio d'Uso

### 6.1 Propagazione Base

```cpp
#include "rkf78_integrator.hpp"

int main() {
    // Stato iniziale (asteroide Sierks)
    State y0;
    y0.r = {1.082716234178542, 3.086569667750047, -0.091582857678774};
    y0.v = {-0.008923488396027, 0.002807706967566, 0.000404145908519};
    
    // Epoche
    double jd0 = 2461000.5;    // Epoca elementi
    double jd1 = 2461008.0913; // Target
    
    // Crea integratore
    RKF78Integrator rk78(1e-12);  // Tolleranza 1e-12
    
    // Funzione accelerazione (include pianeti, AST17, relatività)
    auto accel = [](double t, const auto& r, const auto& v) {
        return compute_acceleration(t, r, v);
    };
    
    // Integra
    IntegrationStats stats;
    State y1 = rk78.integrate(y0, jd0, jd1, accel, &stats);
    
    // Risultati
    std::cout << "Posizione finale: " << y1.r[0] << ", " 
              << y1.r[1] << ", " << y1.r[2] << " AU\n";
    std::cout << "Passi: " << stats.steps_accepted << "\n";
    std::cout << "Valutazioni: " << stats.func_evaluations << "\n";
    
    return 0;
}
```

### 6.2 Output Tipico

```
================================================================
  TEST INTEGRATORE RKF78 (Fehlberg 1968)
================================================================

ASTEROIDE: (17030) Sierks
Epoca:  JD 2461000.5000
Target: JD 2461008.0913
Δt = 7.5913 giorni

Passi accettati: 12
Passi rifiutati: 0
Valutazioni f:   157
Passo min:       7.591e-02 giorni
Passo max:       7.591e-01 giorni

POSIZIONE FINALE:
  RA  = 04 53 11.597
  Dec = +20 19 26.23

CONFRONTO CON JPL HORIZONS:
  Errore: ΔRA*cos(δ)=4.88", ΔDec=0.43" → Tot=4.89"

ROUND-TRIP TEST:
  Errore posizione: 4.4e-16 AU = 0.07 mm ✅
  Errore velocità:  7.2e-18 AU/day = 0.01 nm/s ✅
================================================================
```

---

## 7. Confronto con Altri Integratori

| Integratore | Ordine | Stadi | Adattivo | Precisione | Efficienza |
|-------------|--------|-------|----------|------------|------------|
| RK4 | 4 | 4 | No | ★★★ | ★★★★ |
| RKF45 | 4(5) | 6 | Sì | ★★★★ | ★★★★ |
| **RKF78** | **7(8)** | **13** | **Sì** | **★★★★★** | **★★★★** |
| Bulirsch-Stoer | Variabile | Variabile | Sì | ★★★★★ | ★★★ |

### Vantaggi RKF78

1. **Alta precisione**: Ordine 7 con stima errore ordine 8
2. **Efficienza**: Solo 12 passi per propagazioni di settimane
3. **Robustezza**: Controllo adattivo del passo
4. **Reversibilità**: Errore round-trip a livello di macchina

---

## 8. Limitazioni e Avvertenze

### 8.1 Effemeridi Planetarie

Le posizioni planetarie usano le formule analitiche di **Simon et al. (1994)**:
- Precisione: 1-20" a seconda del pianeta
- Per massima precisione, usare effemeridi JPL DE441

### 8.2 Elementi Orbitali

L'accuratezza finale dipende dalla qualità degli elementi orbitali:
- Elementi MPC: precisione tipica ~5-10"
- Per sub-arcsecond, usare elementi JPL o fit differenziale

### 8.3 Perturbazioni Aggiuntive

Non incluse in questa versione:
- Yarkovsky effect
- Radiazione solare
- Maree terrestri
- Oblateness planetaria (J2)

---

## 9. Riferimenti

1. **Fehlberg, E.** (1968). "Classical fifth-, sixth-, seventh-, and eighth-order Runge-Kutta formulas with stepsize control." NASA TR R-287.

2. **Simon, J.L., et al.** (1994). "Numerical expressions for precession formulae and mean elements for the Moon and the planets." Astronomy & Astrophysics, 282, 663-683.

3. **Standish, E.M.** (1998). "JPL Planetary and Lunar Ephemerides, DE405/LE405." JPL IOM 312.F-98-048.

4. **Farnocchia, D., et al.** (2015). "Asteroid masses from the analysis of their close encounters with other asteroids." Icarus, 245, 94-111.

---

## 10. File Sorgente

| File | Descrizione |
|------|-------------|
| `tools/test_rkf78.cpp` | Implementazione completa RKF78 |
| `tools/sierks_precision.cpp` | Propagatore con RK4 |
| `tools/test_roundtrip.cpp` | Test consistenza round-trip |

---

**© 2025 AstDyS Team - Dipartimento di Matematica, Università di Pisa**
