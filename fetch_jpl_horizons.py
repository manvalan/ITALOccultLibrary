#!/usr/bin/env python3
"""
Scarica e formatta dati JPL Horizons per asteroide 17030 (Sierks)
Genera tabella di confronto AstDyn vs Chebyshev vs JPL
"""

import urllib.request
import urllib.parse
import json
from datetime import datetime, timedelta
import sys

def get_jpl_horizons_data():
    """
    Scarica dati effemeridi da JPL Horizons via API
    Asteroide 17030 (Sierks)
    Data: 2025-11-25 to 2025-11-30, step: 6 hours
    """
    
    # Parametri Horizons
    params = {
        'COMMAND': '17030',           # Asteroid 17030 Sierks
        'EPHEM_TYPE': 'VECTORS',      # Coordinate vector output
        'OUT_UNITS': 'AU-D',          # AU and AU/day
        'COORD_TYPE': 'GEODETIC',     # Geodetic/ICRF
        'REF_PLANE': 'ECLIPTIC',      # Ecliptic
        'CENTER': '@sun',             # Barycentric Solar System
        'START_TIME': '2025-11-25 00:00:00',
        'STOP_TIME': '2025-11-30 23:59:00',
        'STEP_SIZE': '6h',            # 6 hour intervals
        'VEC_LABELS': 'NO',
        'VEC_CORR': 'NONE',
        'CAL_FORMAT': 'BOTH',
        'TIME_DIGITS': 'MINUTES',
        'CSV_FORMAT': 'YES',
    }
    
    # URL Horizons API
    url = 'https://ssd.jpl.nasa.gov/api/horizons.api'
    
    try:
        print("📡 Scaricamento dati JPL Horizons...")
        query_string = urllib.parse.urlencode(params)
        full_url = f"{url}?{query_string}"
        
        # Scarica dati
        with urllib.request.urlopen(full_url, timeout=30) as response:
            data = response.read().decode('utf-8')
            
        print("✓ Dati scaricati")
        return data
        
    except Exception as e:
        print(f"✗ Errore download Horizons: {e}")
        print("  Usando dati locali...")
        return None

def parse_horizons_csv(data):
    """Parse CSV output from Horizons"""
    lines = data.split('\n')
    
    results = []
    in_data = False
    
    for line in lines:
        # Cerca inizio dati
        if '$$SOE' in line:
            in_data = True
            continue
        if '$$EOE' in line:
            break
            
        if in_data and line.strip():
            # Parse line: 2025-Nov-25 00:00:00.0000   9.009851844833686E-01   3.165000000000000E+00...
            parts = line.split()
            if len(parts) >= 8:
                try:
                    date_str = f"{parts[0]} {parts[1]}"
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    vx = float(parts[5])
                    vy = float(parts[6])
                    vz = float(parts[7])
                    
                    results.append({
                        'date': date_str,
                        'x': x, 'y': y, 'z': z,
                        'vx': vx, 'vy': vy, 'vz': vz
                    })
                except ValueError:
                    continue
    
    return results

def gregorian_to_mjd(date_str):
    """Convert 'YYYY-Mon-DD HH:MM' to MJD"""
    try:
        dt = datetime.strptime(date_str, '%Y-%b-%d %H:%M:%S.%f')
    except ValueError:
        dt = datetime.strptime(date_str, '%Y-%b-%d %H:%M')
    
    # JD from datetime
    a = (14 - dt.month) // 12
    y = dt.year + 4800 - a
    m = dt.month + 12 * a - 3
    
    jd = dt.day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
    
    # Add time of day
    fraction = dt.hour / 24.0 + dt.minute / 1440.0 + dt.second / 86400.0
    
    # MJD = JD - 2400000.5
    mjd = jd - 2400000.5 + fraction
    
    return mjd

def generate_cpp_data(data):
    """Generate C++ vector initialization"""
    
    print("\n📋 Dati JPL Horizons (convertiti in C++):\n")
    print("std::vector<JPLData> data = {")
    
    for entry in data:
        mjd = gregorian_to_mjd(entry['date'])
        print(f"    // {entry['date']} = MJD {mjd:.2f}")
        print(f"    {{{mjd:.2f}, {entry['x']:.5f}, {entry['y']:.5f}, {entry['z']:.5f}, "
              f"{entry['vx']:.5f}, {entry['vy']:.5f}, {entry['vz']:.5f}}},")
    
    print("};")

def generate_markdown_table(data):
    """Generate Markdown table"""
    print("\n\n## 📊 Dati JPL Horizons - Asteroide 17030 Sierks (25-30 Nov 2025)\n")
    print("| Data/Ora | MJD TDB | X (AU) | Y (AU) | Z (AU) | VX (AU/d) | VY (AU/d) | VZ (AU/d) |")
    print("|----------|---------|--------|--------|--------|-----------|-----------|-----------|")
    
    for entry in data:
        mjd = gregorian_to_mjd(entry['date'])
        print(f"| {entry['date']} | {mjd:.4f} | {entry['x']:.6f} | {entry['y']:.6f} | {entry['z']:.6f} | "
              f"{entry['vx']:.6f} | {entry['vy']:.6f} | {entry['vz']:.6f} |")

if __name__ == '__main__':
    # Scarica dati
    csv_data = get_jpl_horizons_data()
    
    if csv_data:
        # Parse
        data = parse_horizons_csv(csv_data)
        
        if data:
            print(f"✓ {len(data)} effemeridi scaricate")
            
            # Genera output
            generate_cpp_data(data)
            generate_markdown_table(data)
            
            # Salva CSV
            with open('jpl_horizons_17030.csv', 'w') as f:
                f.write("Date,X_AU,Y_AU,Z_AU,VX_AU_d,VY_AU_d,VZ_AU_d,MJD_TDB\n")
                for entry in data:
                    mjd = gregorian_to_mjd(entry['date'])
                    f.write(f"{entry['date']},{entry['x']},{entry['y']},{entry['z']},"
                           f"{entry['vx']},{entry['vy']},{entry['vz']},{mjd:.4f}\n")
            
            print("\n✓ Dati salvati in jpl_horizons_17030.csv")
