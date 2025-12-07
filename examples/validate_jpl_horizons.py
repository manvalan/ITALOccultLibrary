#!/usr/bin/env python3
"""
Query JPL Horizons per asteroide 17030 Sierks
Confronta risultati con output AstDyn
"""

import sys

# Risultati da AstDyn (UFFICIALI - elementi da AstDyS database)
astdyn_results = {
    'mjd': 61007.0,  # 28 Nov 2025
    'x': 1.020031376556,
    'y': 3.105582988568,
    'z': -0.088735066765,
    'vx': -0.008985735657,
    'vy': 0.002624448538,
    'vz': 0.000409483640
}

print("=" * 70)
print("  VALIDAZIONE ASTDYN - Confronto con JPL Horizons")
print("=" * 70)
print()

print("RISULTATI ASTDYN:")
print(f"  Epoca: MJD {astdyn_results['mjd']}")
print(f"  Posizione (AU):")
print(f"    X  = {astdyn_results['x']:18.12f}")
print(f"    Y  = {astdyn_results['y']:18.12f}")
print(f"    Z  = {astdyn_results['z']:18.12f}")
print(f"  Velocità (AU/day):")
print(f"    VX = {astdyn_results['vx']:18.12f}")
print(f"    VY = {astdyn_results['vy']:18.12f}")
print(f"    VZ = {astdyn_results['vz']:18.12f}")
print()

print("=" * 70)
print("  ISTRUZIONI PER CONFRONTO MANUALE CON JPL HORIZONS")
print("=" * 70)
print()
print("1. Apri JPL Horizons Web Interface:")
print("   https://ssd.jpl.nasa.gov/horizons/app.html")
print()
print("2. Configura query:")
print("   - Target Body: 17030 Sierks")
print("   - Observer Location: @sun (Heliocentric)")
print("   - Time Specification:")
print(f"     Start: 2025-Nov-28 00:00 UTC (MJD {astdyn_results['mjd']})")
print("     Stop:  2025-Nov-28 00:00 UTC")
print("     Step:  1 day")
print("   - Table Settings:")
print("     Type: Vector Table")
print("     Reference Frame: ICRF")
print("     Coordinate Type: Cartesian")
print("     Units: AU and AU/day")
print()
print("3. Confronta valori ottenuti")
print()
print("4. Calcola errore:")
print("   error = sqrt((x_jpl - x_astdyn)² + (y_jpl - y_astdyn)² + (z_jpl - z_astdyn)²)")
print()

print("=" * 70)
print("  QUERY HORIZONS API (se disponibile)")
print("=" * 70)
print()

try:
    from astroquery.jplhorizons import Horizons
    
    print("Interrogando JPL Horizons API...")
    print()
    
    # Query JPL Horizons
    obj = Horizons(id='17030', 
                   location='@sun',  # Heliocentric
                   epochs=astdyn_results['mjd'])
    
    # Ottieni vettori
    vectors = obj.vectors(refplane='earth')  # ICRF frame
    
    print("RISULTATI JPL HORIZONS:")
    print(f"  Epoca: {vectors['datetime_jd'][0]} JD")
    print(f"  Posizione (AU):")
    print(f"    X  = {vectors['x'][0]:18.12f}")
    print(f"    Y  = {vectors['y'][0]:18.12f}")
    print(f"    Z  = {vectors['z'][0]:18.12f}")
    print(f"  Velocità (AU/day):")
    print(f"    VX = {vectors['vx'][0]:18.12f}")
    print(f"    VY = {vectors['vy'][0]:18.12f}")
    print(f"    VZ = {vectors['vz'][0]:18.12f}")
    print()
    
    # Calcola differenze
    dx = vectors['x'][0] - astdyn_results['x']
    dy = vectors['y'][0] - astdyn_results['y']
    dz = vectors['z'][0] - astdyn_results['z']
    dvx = vectors['vx'][0] - astdyn_results['vx']
    dvy = vectors['vy'][0] - astdyn_results['vy']
    dvz = vectors['vz'][0] - astdyn_results['vz']
    
    pos_error = (dx**2 + dy**2 + dz**2)**0.5
    vel_error = (dvx**2 + dvy**2 + dvz**2)**0.5
    
    print("=" * 70)
    print("  DIFFERENZE (JPL - AstDyn)")
    print("=" * 70)
    print(f"  ΔX  = {dx:18.12f} AU")
    print(f"  ΔY  = {dy:18.12f} AU")
    print(f"  ΔZ  = {dz:18.12f} AU")
    print(f"  ΔVX = {dvx:18.12f} AU/day")
    print(f"  ΔVY = {dvy:18.12f} AU/day")
    print(f"  ΔVZ = {dvz:18.12f} AU/day")
    print()
    print(f"  Errore posizione: {pos_error:.12f} AU = {pos_error * 149597870.7:.3f} km")
    print(f"  Errore velocità:  {vel_error:.12f} AU/day")
    print()
    
    # Converti in arcosecondi (approssimazione)
    # Per un oggetto a ~3 AU, 1 AU ≈ 688" (dipende dalla geometria)
    arcsec_approx = pos_error * 688.0  # Approssimazione molto grezza
    print(f"  Errore angolare (stima): ~{arcsec_approx:.2f} arcsec")
    print()
    
    # Valutazione
    if pos_error < 1e-6:
        status = "✅ ECCELLENTE"
        comment = "Precisione sub-km, adatta per occultazioni"
    elif pos_error < 1e-4:
        status = "✅ OTTIMO"
        comment = "Precisione sufficiente per la maggior parte degli usi"
    elif pos_error < 1e-2:
        status = "⚠️  ACCETTABILE"
        comment = "Precisione adeguata per propagazioni brevi"
    else:
        status = "❌ INSUFFICIENTE"
        comment = "Errore troppo grande, verificare configurazione"
    
    print("=" * 70)
    print(f"  VALIDAZIONE: {status}")
    print("=" * 70)
    print(f"  {comment}")
    print()

except ImportError:
    print("❌ Modulo astroquery non disponibile.")
    print()
    print("Per installare:")
    print("  pip install astroquery")
    print()
    print("In alternativa, usa il confronto manuale sopra indicato.")
    print()

except Exception as e:
    print(f"❌ Errore durante query JPL Horizons: {e}")
    print()
    print("Usa il confronto manuale sopra indicato.")
    print()

print("=" * 70)
