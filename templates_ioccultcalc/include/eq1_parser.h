/**
 * @file eq1_parser.h
 * @brief Parser per file .eq1 formato OEF2.0 (AstDyS/OrbFit)
 * @author IOccultCalc Integration Team
 * @date 1 Dicembre 2025
 * 
 * Questo parser legge file di elementi orbitali nel formato OEF2.0
 * (Orbital Elements Format) usato da AstDyS e OrbFit.
 * 
 * Formato file .eq1:
 *   format  = 'OEF2.0'
 *   rectype = 'ML'
 *   refsys  = ECLM J2000
 *   END_OF_HEADER
 *   17030
 *   EQU   a   h   k   p   q   lambda
 *   MJD   epoch_mjd TDT
 *   MAG   H   G
 *   END
 */

#ifndef IOCCULTCALC_EQ1_PARSER_H
#define IOCCULTCALC_EQ1_PARSER_H

#include <string>
#include <stdexcept>

namespace ioccultcalc {

/**
 * @struct EquinoctialElements
 * @brief Elementi orbitali equinoziali (formato AstDyS/OrbFit)
 * 
 * Sistema equinoziale (Broucke & Cefola):
 * - a: semiasse maggiore [AU]
 * - h = e·sin(ϖ) dove ϖ = Ω + ω (longitudine perielio)
 * - k = e·cos(ϖ)
 * - p = tan(i/2)·sin(Ω) dove Ω = longitudine nodo ascendente
 * - q = tan(i/2)·cos(Ω)
 * - λ = Ω + ω + M (longitudine media) [gradi]
 * 
 * Frame di riferimento: ECLM J2000 (eclittica media J2000)
 */
struct EquinoctialElements {
    double a;           ///< Semiasse maggiore [AU]
    double h;           ///< e*sin(longitudine perielio)
    double k;           ///< e*cos(longitudine perielio)
    double p;           ///< tan(i/2)*sin(longitudine nodo)
    double q;           ///< tan(i/2)*cos(longitudine nodo)
    double lambda;      ///< Longitudine media [gradi]
    double epoch_mjd;   ///< Epoca [Modified Julian Date, TDT]
    double H;           ///< Magnitudine assoluta
    double G;           ///< Parametro di pendenza (Bowell)
    std::string name;   ///< Nome/numero asteroide
    
    /// Costruttore di default
    EquinoctialElements() 
        : a(0.0), h(0.0), k(0.0), p(0.0), q(0.0), lambda(0.0),
          epoch_mjd(0.0), H(0.0), G(0.15), name("") {}
    
    /// Conversione MJD → JD
    double getEpochJD() const { return epoch_mjd + 2400000.5; }
    
    /// Calcola eccentricità
    double getEccentricity() const { return std::sqrt(h*h + k*k); }
    
    /// Validazione elementi
    bool isValid() const {
        return a > 0.0 && 
               a < 1000.0 &&  // Limite ragionevole
               getEccentricity() < 1.0 &&  // Orbita ellittica
               epoch_mjd > 0.0;
    }
};

/**
 * @class EQ1Parser
 * @brief Parser per file .eq1 formato OEF2.0
 * 
 * Esempio d'uso:
 * @code
 *   try {
 *       auto elements = EQ1Parser::parseFile("17030.eq1");
 *       std::cout << "Asteroide: " << elements.name << "\n";
 *       std::cout << "a = " << elements.a << " AU\n";
 *       std::cout << "e = " << elements.getEccentricity() << "\n";
 *   } catch (const std::exception& e) {
 *       std::cerr << "Errore: " << e.what() << "\n";
 *   }
 * @endcode
 */
class EQ1Parser {
public:
    /**
     * @brief Legge elementi da file .eq1
     * @param filepath Percorso file .eq1
     * @return Elementi equinoziali parsed
     * @throws std::runtime_error Se file non trovato o formato invalido
     */
    static EquinoctialElements parseFile(const std::string& filepath);
    
    /**
     * @brief Legge elementi da stringa
     * @param content Contenuto file .eq1 come stringa
     * @return Elementi equinoziali parsed
     * @throws std::runtime_error Se formato invalido
     */
    static EquinoctialElements parseString(const std::string& content);
    
private:
    /**
     * @brief Parse singola linea del file
     * @param line Linea da parsare
     * @param elem Struct elementi da riempire (in/out)
     * @return true se linea parsed con successo
     */
    static bool parseLine(const std::string& line, 
                         EquinoctialElements& elem);
    
    /**
     * @brief Estrae valore numerico da stringa MJD
     * @param mjd_line Linea contenente "MJD xxxxx.xxx TDT"
     * @return Valore MJD
     */
    static double extractMJD(const std::string& mjd_line);
    
    /**
     * @brief Trim whitespace da stringa
     * @param str Stringa da pulire (in/out)
     */
    static void trim(std::string& str);
};

/**
 * @class EQ1ParseException
 * @brief Eccezione specifica per errori di parsing .eq1
 */
class EQ1ParseException : public std::runtime_error {
public:
    explicit EQ1ParseException(const std::string& msg) 
        : std::runtime_error("EQ1 Parse Error: " + msg) {}
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_EQ1_PARSER_H
