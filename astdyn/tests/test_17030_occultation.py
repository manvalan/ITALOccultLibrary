#!/usr/bin/env python3
"""
Test occultazione asteroide 17030 usando astdyn_propagator
Calcola posizioni per il 28/11/2026 00:00-01:00 ogni 5 minuti
Confronta con effemeridi JPL e calcola distanza da stella GAIA
"""

import subprocess
import math
from datetime import datetime, timedelta

# Elementi orbitali JPL Horizons (epoca 2018-Mar-16.00 = JD 2458193.5)
ELEM_17030 = {
    'name': '17030',
    'epoch': 2458193.5,
    'a': 3.173489964321051,
    'e': 0.04796607451625862,
    'i': 2.904309538190326,
    'Omega': 104.1845838362649,
    'omega': 102.1497438064497,
    'M': 99.03517819281583
}

# Stella GAIA DR3 3411546266140512128 (J2000.0)
STAR = {
    'ra': 73.4161003759929,     # gradi
    'dec': 20.3316626372542,    # gradi
    'pmra': 1.097,              # mas/yr
    'pmdec': -0.155             # mas/yr
}

# Effemeridi JPL per 28/11/2025 (!!!)
JPL_EPHEM = [
    (0,  '04 53 41.10', '+20 19 56.7'),
    (5,  '04 53 40.93', '+20 19 56.6'),
    (10, '04 53 40.75', '+20 19 56.4'),
    (15, '04 53 40.58', '+20 19 56.2'),
    (20, '04 53 40.40', '+20 19 56.0'),
    (25, '04 53 40.23', '+20 19 55.8'),
    (30, '04 53 40.05', '+20 19 55.7'),
    (35, '04 53 39.88', '+20 19 55.5'),
    (40, '04 53 39.70', '+20 19 55.3'),
    (45, '04 53 39.53', '+20 19 55.1'),
    (50, '04 53 39.35', '+20 19 54.9'),
    (55, '04 53 39.17', '+20 19 54.7'),
    (60, '04 53 39.00', '+20 19 54.6'),
]

def hms_to_deg(hms_str):
    """Converti HH MM SS.ss in gradi"""
    parts = hms_str.strip().split()
    h, m, s = float(parts[0]), float(parts[1]), float(parts[2])
    return (h + m/60.0 + s/3600.0) * 15.0

def dms_to_deg(dms_str):
    """Converti +/-DD MM SS.s in gradi"""
    parts = dms_str.strip().split()
    sign = 1 if parts[0][0] != '-' else -1
    d, m, s = abs(float(parts[0])), float(parts[1]), float(parts[2])
    return sign * (d + m/60.0 + s/3600.0)

def angular_distance(ra1, dec1, ra2, dec2):
    """Calcola distanza angolare in gradi"""
    ra1_rad = math.radians(ra1)
    dec1_rad = math.radians(dec1)
    ra2_rad = math.radians(ra2)
    dec2_rad = math.radians(dec2)
    
    cos_dist = (math.sin(dec1_rad) * math.sin(dec2_rad) +
                math.cos(dec1_rad) * math.cos(dec2_rad) * 
                math.cos(ra1_rad - ra2_rad))
    return math.degrees(math.acos(max(-1, min(1, cos_dist))))

def apply_proper_motion(ra, dec, pmra, pmdec, years):
    """Applica moto proprio"""
    ra_new = ra + (pmra / 3.6e6 * years) / math.cos(math.radians(dec))
    dec_new = dec + pmdec / 3.6e6 * years
    return ra_new, dec_new

def main():
    print("=" * 70)
    print(" Test Occultazione Asteroide 17030")
    print(" Usando astdyn_propagator (RKF78 + perturbazioni)")
    print("=" * 70)
    print()
    
    # Data target: 28/11/2025 00:00 UTC
    target_date = datetime(2025, 11, 28, 0, 0, 0)
    target_jd = 2460277.5  # JD per 28/11/2025 00:00
    
    # Stella con moto proprio (28/11/2025 = MJD 60277.0)
    years_from_j2000 = (60277.0 - 51544.5) / 365.25
    star_ra, star_dec = apply_proper_motion(
        STAR['ra'], STAR['dec'], 
        STAR['pmra'], STAR['pmdec'], 
        years_from_j2000
    )
    
    print(f"Stella GAIA DR3 3411546266140512128 (epoca 28/11/2025):")
    print(f"  RA  = {star_ra:.8f}°")
    print(f"  Dec = {star_dec:.8f}°")
    print()
    
    print("Nota: astdyn_propagator ha esempio hardcoded per asteroide 17030")
    print("      Qui mostriamo solo il calcolo della distanza angolare")
    print()
    
    # Compila ed esegui astdyn_propagator
    print("Compilazione astdyn_propagator...")
    compile_cmd = [
        'g++', '-std=c++17', '-O2', 
        '../tools/astdyn_propagator.cpp',
        '-o', 'astdyn_propagator_tmp',
        '-lm'
    ]
    result = subprocess.run(compile_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERRORE compilazione: {result.stderr}")
        return
    
    print("Esecuzione astdyn_propagator...\n")
    run_cmd = ['./astdyn_propagator_tmp']
    result = subprocess.run(run_cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        # Mostra output
        print(result.stdout)
    else:
        print(f"ERRORE: {result.stderr}")
    
    # Calcola distanze con effemeridi JPL
    print("\n" + "=" * 70)
    print(" Confronto con JPL e distanza da stella")
    print("=" * 70)
    print()
    print("Tempo UTC        RA (JPL)        Dec (JPL)       Dist da Stella")
    print("-" * 70)
    
    min_dist = float('inf')
    closest_time = 0
    
    for minute, ra_str, dec_str in JPL_EPHEM:
        jpl_ra = hms_to_deg(ra_str)
        jpl_dec = dms_to_deg(dec_str)
        
        # Distanza da stella
        dist = angular_distance(jpl_ra, jpl_dec, star_ra, star_dec)
        dist_arcsec = dist * 3600
        
        if dist < min_dist:
            min_dist = dist
            closest_time = minute
        
        h = minute // 60
        m = minute % 60
        close_marker = "  *** POSSIBILE OCCULTAZIONE ***" if dist_arcsec < 5.0 else ""
        print(f"28/11/2025 {h:02d}:{m:02d}:00   {jpl_ra:12.6f}°  {jpl_dec:12.6f}°  {dist_arcsec:10.2f}\"{close_marker}")
    
    print()
    print("=" * 70)
    print(f"Distanza minima da stella: {min_dist * 3600:.2f} arcsec = {min_dist:.4f}°")
    print(f"Tempo minima distanza: {closest_time // 60:02d}:{closest_time % 60:02d} UTC")
    print("=" * 70)
    
    # Cleanup
    subprocess.run(['rm', '-f', 'astdyn_propagator_tmp'], capture_output=True)

if __name__ == '__main__':
    main()
