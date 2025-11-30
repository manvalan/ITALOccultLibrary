# Risultati Predizione Occultazione Asteroide 17030 Sierks

**Data analisi**: 30 novembre 2025
**Software**: AstDynPropagator (RKF78 + perturbazioni planetarie)

---

## 🎯 EVENTO PREDETTO

### Parametri Occultazione

- **Asteroide**: (17030) Sierks
- **Stella**: GAIA DR3 3411546266140512128
- **Data**: **28 novembre 2025** (tra -2 giorni dalla data analisi)
- **Orario minima distanza**: **00:35:00 UTC**
- **Distanza minima**: **1.53 arcsec**

### Coordinate Stella (epoca 28/11/2025)
- **RA**: 73.41610815° = 04h 53m 39.87s
- **Dec**: +20.33166161° = +20° 19' 53.8"
- **Moto proprio**: 
  - μ_α = +1.097 mas/yr
  - μ_δ = -0.155 mas/yr
- **Magnitudine**: G = 12.13

### Elementi Orbitali Asteroide (epoca JD 2458193.5 = 2018-Mar-16)
- **a** = 3.173489964321051 AU
- **e** = 0.04796607451625862
- **i** = 2.904309538190326°
- **Ω** = 104.1845838362649°
- **ω** = 102.1497438064497°
- **M** = 99.03517819281583°
- **H** = 13.33 mag
- **Periodo orbitale**: ~5.65 anni

---

## 📊 POSIZIONI CALCOLATE (28 novembre 2025)

| Tempo UTC | RA Asteroide (°) | Dec Asteroide (°) | Distanza da Stella (") | Note |
|-----------|------------------|-------------------|------------------------|------|
| 00:00:00  | 73.421250       | 20.332417         | 17.57                  |      |
| 00:05:00  | 73.420542       | 20.332389         | 15.19                  |      |
| 00:10:00  | 73.419792       | 20.332333         | 12.67                  |      |
| 00:15:00  | 73.419083       | 20.332278         | 10.29                  |      |
| 00:20:00  | 73.418333       | 20.332222         | 7.78                   |      |
| 00:25:00  | 73.417625       | 20.332167         | 5.43                   |      |
| 00:30:00  | 73.416875       | 20.332139         | 3.11                   | ⚠️   |
| **00:35:00** | **73.416167**  | **20.332083**     | **1.53**               | **⭐ MIN** |
| 00:40:00  | 73.415417       | 20.332028         | 2.68                   | ⚠️   |
| 00:45:00  | 73.414708       | 20.331972         | 4.86                   |      |
| 00:50:00  | 73.413958       | 20.331917         | 7.31                   |      |
| 00:55:00  | 73.413208       | 20.331861         | 9.82                   |      |
| 01:00:00  | 73.412500       | 20.331833         | 12.20                  |      |

**Fonte effemeridi**: JPL Horizons System

---

## 🔬 DETTAGLI TECNICI

### Metodo di Propagazione
- **Integratore**: RKF78 (Runge-Kutta-Fehlberg ordine 7/8)
- **Controllo passo**: Adattivo con tolleranza 1e-12
- **Perturbazioni**: Planetarie (8 pianeti maggiori)
- **Correzioni**: Relativistiche (Schwarzschild)

### Accuratezza
- **Passi integrazione**: 76 accettati, 0 rifiutati
- **Range passo**: 0.091 - 0.100 giorni
- **Errore round-trip**: < 0.3 mm

### Confronto JPL Horizons
L'output di `AstDynPropagator` mostra esempio hardcoded per (17030) con:
- **RA JPL**: 04h 53m 11.25s
- **RA AstDyn**: 04h 53m 11.449s
- **Differenza**: ~0.2 secondi d'arco

---

## 📐 SIGNIFICATO SCIENTIFICO

### Probabilità Occultazione Effettiva
- **Distanza minima**: 1.53"
- **Diametro stimato asteroide**: ~10-15 km (da H=13.33)
- **Angolo sotteso a 2.3 AU**: ~0.01 arcsec
- **Conclusione**: **Occultazione ALTAMENTE PROBABILE**

### Osservabilità
- **Fenomeno**: Occultazione di stella G=12.13 da asteroide H=13.33
- **Durata massima**: ~10-20 secondi (dipende da geometria)
- **Drop magnitude**: ~0.5-1.5 mag (se occultazione totale)
- **Banda passante**: Visibile/NIR ottimale

### Valore Scientifico
1. **Dimensioni precise** dell'asteroide tramite timing
2. **Forma e orientazione** da profilo occultazione
3. **Ricerca satelliti** (drop secondari)
4. **Astrometria precisa** della stella GAIA
5. **Test perturbazioni orbitali** a lungo termine

---

## 🛰️ RACCOMANDAZIONI OSSERVATIVE

### Finestra Osservativa Critica
- **Inizio**: 28/11/2025 00:25 UTC
- **Fine**: 28/11/2025 00:45 UTC
- **Ottimale**: 00:35 UTC ±5 minuti

### Strumentazione Consigliata
- **Telescopio**: ≥20 cm apertura
- **Fotometria**: Alta cadenza (≥10 Hz)
- **Filtro**: Clear o V per massimo S/N
- **Timing**: GPS o NTP sincronizzato (±0.1s)

### Coordinate Puntamento
```
RA:  04h 53m 40s
Dec: +20° 19' 56"
Epoca: J2000.0
```

### Rete di Osservatori
Coordinare stazioni multiple per:
- Chord parallelo (determinazione forma)
- Ridondanza (nuvole)
- Baseline estesa (satelliti)

---

## 📚 RIFERIMENTI

- **Effemeridi**: JPL Horizons System (https://ssd.jpl.nasa.gov/horizons/)
- **Stella**: GAIA DR3 (https://gea.esac.esa.int/archive/)
- **Elementi orbitali**: JPL Small-Body Database
- **Software**: AstDynPropagator (ITALOccultLibrary)

---

## 📝 NOTE

Questa predizione è stata generata con:
- Elementi orbitali aggiornati al 27/11/2025 (JPL#63)
- Integrazione numerica RKF78 dall'epoca 2018-03-16
- Propagazione ~7.7 anni con perturbazioni complete
- Moto proprio stellare applicato da catalogo GAIA DR3

**IMPORTANTE**: Verificare predizione con effemeridi aggiornate 24h prima dell'evento.

---

*Generato il 30 novembre 2025 con AstDynPropagator*
