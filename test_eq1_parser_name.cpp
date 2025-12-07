/**
 * Test per verificare cosa legge il parser dal file .eq1
 */

#include <iostream>
#include <string>
#include <astdyn/io/parsers/OrbFitEQ1Parser.hpp>

int main() {
    try {
        astdyn::io::parsers::OrbFitEQ1Parser parser;
        
        std::cout << "Parsing astdyn/data/17030.eq1..." << std::endl;
        auto elements = parser.parse("astdyn/data/17030.eq1");
        
        std::cout << "\nRisultati parsing:" << std::endl;
        std::cout << "  object_name: '" << elements.object_name << "'" << std::endl;
        std::cout << "  lunghezza: " << elements.object_name.length() << " caratteri" << std::endl;
        std::cout << "  vuoto? " << (elements.object_name.empty() ? "SI" : "NO") << std::endl;
        std::cout << "  epoch_mjd_tdb: " << elements.epoch_mjd_tdb << std::endl;
        std::cout << "  semi_major_axis: " << elements.semi_major_axis << " AU" << std::endl;
        
        if (elements.object_name.empty()) {
            std::cout << "\n⚠️  PROBLEMA: object_name è vuoto!" << std::endl;
            std::cout << "Il parser non sta estraendo il nome dall'intestazione del file .eq1" << std::endl;
        } else {
            std::cout << "\n✓ Nome estratto correttamente" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
