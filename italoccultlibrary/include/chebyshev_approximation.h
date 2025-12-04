/**
 * @file chebyshev_approximation.h
 * @brief Approssimazione mediante polinomi di Chebyshev per traiettorie asteroidali
 * @author ITALOccultLibrary Development Team
 * @date 4 Dicembre 2025
 * 
 * Questa classe fornisce un'approssimazione mediante polinomi di Chebyshev
 * per le coordinate cartesiane della traiettoria asteroidale su un intervallo
 * temporale. I polinomi di Chebyshev sono ottimi per l'approssimazione di
 * funzioni lisce e consentono valutazione veloce della posizione.
 * 
 * Utilizzo tipico:
 * @code
 * // Crea approssimazione con 10 coefficienti per componente
 * ChebyshevApproximation approx(10);
 * 
 * // Propaga asteroide e raccoglie stati
 * std::vector<AsteroidState> states = integrator.propagateToEpochs(epochs);
 * 
 * // Fitta i polinomi di Chebyshev
 * approx.fit(states, start_epoch, end_epoch);
 * 
 * // Valuta posizione in un'epoca arbitraria
 * Eigen::Vector3d position = approx.evaluatePosition(query_epoch);
 * @endcode
 */

#ifndef ITALOCCULTLIB_CHEBYSHEV_APPROXIMATION_H
#define ITALOCCULTLIB_CHEBYSHEV_APPROXIMATION_H

#include <vector>
#include <array>
#include <Eigen/Core>
#include <memory>

namespace ioccultcalc {

// Forward declaration - definito in italoccult_integration.h
struct AsteroidState;

/**
 * @struct ChebyshevCoefficients
 * @brief Contenitore per i coefficienti di Chebyshev di una componente
 */
struct ChebyshevCoefficients {
    std::vector<double> coefficients;  ///< Coefficienti a0, a1, a2, ...
    double t_min;                      ///< Valore minimo del parametro t
    double t_max;                      ///< Valore massimo del parametro t
    
    /**
     * @brief Normalizza un tempo nell'intervallo [-1, 1] per la valutazione
     * @param t_epoch Epoca in MJD
     * @return Parametro normalizzato nel range [-1, 1]
     */
    double normalize(double t_epoch) const {
        return 2.0 * (t_epoch - t_min) / (t_max - t_min) - 1.0;
    }
    
    /**
     * @brief Denormalizza un parametro dall'intervallo [-1, 1] a tempo MJD
     * @param t_norm Parametro normalizzato
     * @return Epoca in MJD
     */
    double denormalize(double t_norm) const {
        return t_min + (t_norm + 1.0) * (t_max - t_min) / 2.0;
    }
};

/**
 * @class ChebyshevApproximation
 * @brief Approssimazione di traiettorie asteroidali mediante polinomi di Chebyshev
 * 
 * Questa classe calcola i coefficienti dei polinomi di Chebyshev per
 * approssimare la posizione e velocità di un asteroide su un intervallo
 * temporale. Fornisce metodi efficienti per:
 * 
 * - Fitting dei coefficienti da un set di stati
 * - Valutazione della posizione in un'epoca arbitraria
 * - Valutazione della velocità mediante derivata
 * - Analisi dell'errore di approssimazione
 * - Calcolo dell'energia orbitale
 */
class ChebyshevApproximation {
public:
    /**
     * @brief Costruttore
     * @param num_coefficients Numero di coefficienti di Chebyshev per componente (default: 10)
     */
    explicit ChebyshevApproximation(size_t num_coefficients = 10);
    
    /**
     * @brief Destruttore
     */
    ~ChebyshevApproximation() = default;
    
    /**
     * @brief Fitta i polinomi di Chebyshev su un set di stati asteroidali
     * @param states Vettore di stati (posizione + velocità)
     * @param start_epoch Epoca iniziale (MJD)
     * @param end_epoch Epoca finale (MJD)
     * @return true se il fitting è riuscito, false altrimenti
     * 
     * Il fitting utilizza un metodo ai minimi quadrati per trovare i
     * coefficienti di Chebyshev che meglio approssimano i dati.
     */
    bool fit(const std::vector<Eigen::Vector3d>& positions,
             double start_epoch,
             double end_epoch);
    
    /**
     * @brief Valuta la posizione in un'epoca arbitraria
     * @param epoch Epoca in MJD dove valutare
     * @return Posizione (X, Y, Z) in AU nel frame ICRF
     * @throw std::runtime_error Se i polinomi non sono stati fittati
     */
    Eigen::Vector3d evaluatePosition(double epoch) const;
    
    /**
     * @brief Valuta la velocità mediante derivata numerica dei polinomi
     * @param epoch Epoca in MJD dove valutare
     * @return Velocità (VX, VY, VZ) in AU/day nel frame ICRF
     * @throw std::runtime_error Se i polinomi non sono stati fittati
     */
    Eigen::Vector3d evaluateVelocity(double epoch) const;
    
    /**
     * @brief Valuta sia posizione che velocità
     * @param epoch Epoca in MJD
     * @return Coppia (posizione, velocità)
     */
    std::pair<Eigen::Vector3d, Eigen::Vector3d> evaluate(double epoch) const {
        return {evaluatePosition(epoch), evaluateVelocity(epoch)};
    }
    
    /**
     * @brief Calcola l'errore di approssimazione rispetto ai dati di training
     * @return RMS dell'errore di posizione (AU) per le tre componenti
     */
    Eigen::Vector3d getApproximationError() const;
    
    /**
     * @brief Valuta l'energia orbitale specifica
     * @param epoch Epoca in MJD
     * @return Energia specifica (AU²/day²)
     */
    double evaluateEnergy(double epoch) const;
    
    /**
     * @brief Valuta il momento angolare orbitale
     * @param epoch Epoca in MJD
     * @return Momento angolare (AU²/day)
     */
    Eigen::Vector3d evaluateAngularMomentum(double epoch) const;
    
    /**
     * @brief Restituisce i coefficienti di Chebyshev per una componente
     * @param component 0=X, 1=Y, 2=Z
     * @return Struct ChebyshevCoefficients
     * @throw std::out_of_range Se component non è 0, 1 o 2
     */
    const ChebyshevCoefficients& getCoefficients(size_t component) const;
    
    /**
     * @brief Numero di coefficienti per componente
     * @return Numero di coefficienti
     */
    size_t getNumCoefficients() const { return num_coefficients_; }
    
    /**
     * @brief Intervallo temporale di validità dell'approssimazione
     * @return Coppia (start_epoch, end_epoch) in MJD
     */
    std::pair<double, double> getTimeInterval() const {
        return {start_epoch_, end_epoch_};
    }
    
    /**
     * @brief Verifica se i polinomi sono stati fittati
     * @return true se i polinomi sono disponibili
     */
    bool isFitted() const { return fitted_; }
    
    /**
     * @brief Ripristina lo stato (cancella i coefficienti)
     */
    void reset();
    
    /**
     * @brief Salva i coefficienti su file
     * @param filename Percorso del file di output
     * @return true se il salvataggio è riuscito
     */
    bool saveToFile(const std::string& filename) const;
    
    /**
     * @brief Carica i coefficienti da file
     * @param filename Percorso del file da caricare
     * @return true se il caricamento è riuscito
     */
    bool loadFromFile(const std::string& filename);
    
    /**
     * @brief Statistica sulla qualità dell'approssimazione
     * @return Stringa con informazioni dettagliate
     */
    std::string getStatistics() const;

private:
    size_t num_coefficients_;                           ///< Numero di coefficienti
    std::array<ChebyshevCoefficients, 3> coefficients_; ///< Coefficienti X, Y, Z
    std::vector<double> training_epochs_;               ///< Epoche di training
    std::array<Eigen::VectorXd, 3> training_data_;      ///< Dati di training
    double start_epoch_;                                ///< Epoca iniziale dell'intervallo
    double end_epoch_;                                  ///< Epoca finale dell'intervallo
    bool fitted_;                                       ///< Flag di fitting completato
    
    /**
     * @brief Calcola il polinomio di Chebyshev di ordine n
     * @param n Ordine del polinomio
     * @param t Parametro nel range [-1, 1]
     * @return T_n(t)
     */
    static double chebyshevT(size_t n, double t);
    
    /**
     * @brief Calcola la derivata del polinomio di Chebyshev di ordine n
     * @param n Ordine del polinomio
     * @param t Parametro nel range [-1, 1]
     * @return dT_n(t)/dt
     */
    static double chebyshevTDerivative(size_t n, double t);
    
    /**
     * @brief Effettua il fitting ai minimi quadrati
     * @param data Dati da fittare (per una componente)
     * @param epochs Epoche corrispondenti ai dati
     * @return Vettore dei coefficienti di Chebyshev
     */
    std::vector<double> fitComponent(
        const std::vector<double>& data,
        const std::vector<double>& epochs);
    
    /**
     * @brief Valuta un polinomio di Chebyshev ai coefficienti dati
     * @param coeffs Vettore di coefficienti
     * @param t Parametro nel range [-1, 1]
     * @return Valore del polinomio
     */
    static double evaluateChebyshev(
        const std::vector<double>& coeffs,
        double t);
    
    /**
     * @brief Valuta la derivata di un polinomio di Chebyshev
     * @param coeffs Vettore di coefficienti
     * @param t Parametro nel range [-1, 1]
     * @return Valore della derivata rispetto a t
     */
    static double evaluateChebyshevDerivative(
        const std::vector<double>& coeffs,
        double t);
};

} // namespace ioccultcalc

#endif // ITALOCCULTLIB_CHEBYSHEV_APPROXIMATION_H
