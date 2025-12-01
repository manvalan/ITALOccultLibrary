/**
 * @file eq1_parser.cpp
 * @brief Implementazione parser file .eq1
 * @author IOccultCalc Integration Team
 * @date 1 Dicembre 2025
 */

#include "eq1_parser.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cmath>

namespace ioccultcalc {

// ============================================================================
// Funzioni di utility
// ============================================================================

void EQ1Parser::trim(std::string& str) {
    // Rimuovi whitespace a sinistra
    str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
    
    // Rimuovi whitespace a destra
    str.erase(std::find_if(str.rbegin(), str.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), str.end());
}

double EQ1Parser::extractMJD(const std::string& mjd_line) {
    // Cerca "MJD" nella linea
    size_t mjd_pos = mjd_line.find("MJD");
    if (mjd_pos == std::string::npos) {
        throw EQ1ParseException("MJD keyword not found in line: " + mjd_line);
    }
    
    // Estrai parte dopo "MJD"
    std::string mjd_part = mjd_line.substr(mjd_pos + 3);
    trim(mjd_part);
    
    // Trova primo carattere numerico
    size_t num_start = mjd_part.find_first_of("0123456789");
    if (num_start == std::string::npos) {
        throw EQ1ParseException("No numeric value found after MJD");
    }
    
    mjd_part = mjd_part.substr(num_start);
    
    // Estrai solo la parte numerica (fino a TDT/UTC o spazio)
    size_t num_end = mjd_part.find_first_not_of("0123456789.");
    if (num_end != std::string::npos) {
        mjd_part = mjd_part.substr(0, num_end);
    }
    
    try {
        return std::stod(mjd_part);
    } catch (const std::exception& e) {
        throw EQ1ParseException("Invalid MJD value: " + mjd_part);
    }
}

// ============================================================================
// Parse functions
// ============================================================================

bool EQ1Parser::parseLine(const std::string& line, 
                          EquinoctialElements& elem) {
    // Linea vuota o commento
    if (line.empty() || line[0] == '!' || line[0] == '#') {
        return false;
    }
    
    // Elementi equinoziali: EQU a h k p q lambda
    if (line.find("EQU") != std::string::npos) {
        std::istringstream iss(line);
        std::string tag;
        iss >> tag >> elem.a >> elem.h >> elem.k 
            >> elem.p >> elem.q >> elem.lambda;
        
        if (iss.fail()) {
            throw EQ1ParseException("Invalid EQU line format: " + line);
        }
        return true;
    }
    
    // Epoca: MJD xxxxx.xxx TDT
    if (line.find("MJD") != std::string::npos) {
        elem.epoch_mjd = extractMJD(line);
        return true;
    }
    
    // Magnitudine: MAG H G
    if (line.find("MAG") != std::string::npos) {
        std::istringstream iss(line);
        std::string tag;
        iss >> tag >> elem.H >> elem.G;
        
        if (iss.fail()) {
            // MAG è opzionale, ignora errori
            elem.H = 15.0;  // Default
            elem.G = 0.15;
        }
        return true;
    }
    
    // Numero/nome asteroide (linea numerica senza keyword)
    if (!line.empty() && std::isdigit(line[0])) {
        // Controlla che non sia una linea EQU/MJD/MAG
        if (line.find("EQU") == std::string::npos &&
            line.find("MJD") == std::string::npos &&
            line.find("MAG") == std::string::npos) {
            elem.name = line;
            return true;
        }
    }
    
    return false;
}

EquinoctialElements EQ1Parser::parseFile(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw EQ1ParseException("Cannot open file: " + filepath);
    }
    
    // Leggi tutto il contenuto
    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();
    
    return parseString(buffer.str());
}

EquinoctialElements EQ1Parser::parseString(const std::string& content) {
    EquinoctialElements elem;
    std::istringstream stream(content);
    std::string line;
    bool in_data_section = false;
    bool found_equ = false;
    bool found_mjd = false;
    
    while (std::getline(stream, line)) {
        // Trim whitespace
        trim(line);
        
        // Skip linee vuote
        if (line.empty()) continue;
        
        // Skip commenti
        if (line[0] == '!' || line[0] == '#') continue;
        
        // Header end marker
        if (line.find("END_OF_HEADER") != std::string::npos) {
            in_data_section = true;
            continue;
        }
        
        // Se non siamo ancora nella sezione dati, skip
        if (!in_data_section) continue;
        
        // Fine dati
        if (line.find("END") != std::string::npos && 
            line.find("END_OF_HEADER") == std::string::npos) {
            break;
        }
        
        // Parse linea
        try {
            if (parseLine(line, elem)) {
                // Traccia quali campi critici abbiamo trovato
                if (line.find("EQU") != std::string::npos) found_equ = true;
                if (line.find("MJD") != std::string::npos) found_mjd = true;
            }
        } catch (const EQ1ParseException& e) {
            // Propaga eccezioni di parsing
            throw;
        }
    }
    
    // Validazione: campi critici devono essere presenti
    if (!found_equ) {
        throw EQ1ParseException("EQU line not found in file");
    }
    
    if (!found_mjd) {
        throw EQ1ParseException("MJD line not found in file");
    }
    
    // Validazione elementi
    if (!elem.isValid()) {
        std::ostringstream oss;
        oss << "Invalid orbital elements: a=" << elem.a 
            << " e=" << elem.getEccentricity() 
            << " epoch_mjd=" << elem.epoch_mjd;
        throw EQ1ParseException(oss.str());
    }
    
    // Se nome non specificato, usa "Unknown"
    if (elem.name.empty()) {
        elem.name = "Unknown";
    }
    
    return elem;
}

} // namespace ioccultcalc
