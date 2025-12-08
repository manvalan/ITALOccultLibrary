/**
 * @file chebyshev_approximation.cpp
 * @brief Implementazione dell'approssimazione mediante polinomi di Chebyshev
 * @author ITALOccultLibrary Development Team
 * @date 4 Dicembre 2025
 */

#include "chebyshev_approximation.h"
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <Eigen/Dense>

namespace ioccultcalc {

ChebyshevApproximation::ChebyshevApproximation(size_t num_coefficients)
    : num_coefficients_(num_coefficients),
      start_epoch_(0.0),
      end_epoch_(0.0),
      fitted_(false) {
    
    // Inizializza i coefficienti
    for (size_t i = 0; i < 3; ++i) {
        coefficients_[i].coefficients.resize(num_coefficients, 0.0);
        coefficients_[i].t_min = -1.0;
        coefficients_[i].t_max = 1.0;
    }
}

double ChebyshevApproximation::chebyshevT(size_t n, double t) {
    // Calcola T_n(t) ricorsivamente usando la relazione:
    // T_0(t) = 1
    // T_1(t) = t
    // T_n(t) = 2*t*T_{n-1}(t) - T_{n-2}(t)
    
    if (n == 0) return 1.0;
    if (n == 1) return t;
    
    double t0 = 1.0;
    double t1 = t;
    double tn;
    
    for (size_t i = 2; i <= n; ++i) {
        tn = 2.0 * t * t1 - t0;
        t0 = t1;
        t1 = tn;
    }
    
    return t1;
}

double ChebyshevApproximation::chebyshevTDerivative(size_t n, double t) {
    // Calcola dT_n(t)/dt usando la relazione:
    // dT_0(t)/dt = 0
    // dT_1(t)/dt = 1
    // dT_n(t)/dt = n*U_{n-1}(t) dove U è il polinomio di Chebyshev di seconda specie
    // 
    // Alternativamente: dT_n(t)/dt = n * T_{n-1}(t) / (1 - t²) per |t| != 1
    // Ma usiamo una formula più stabile:
    // dT_n/dt = 2*n * sum_{k=0}^{n-1, n-k dispari} k*T_{n-2k-1}(t)
    
    if (n == 0) return 0.0;
    if (n == 1) return 1.0;
    
    // Formula ricorsiva: dT_n/dt = 2*t * dT_{n-1}/dt - dT_{n-2}/dt + 2*T_{n-1}(t)
    double dt0 = 0.0;
    double dt1 = 1.0;
    double dtn;
    double t0 = 1.0;
    double t1 = t;
    
    for (size_t i = 2; i <= n; ++i) {
        dtn = 2.0 * t * dt1 - dt0 + 2.0 * t1;
        dt0 = dt1;
        dt1 = dtn;
        
        double tn = 2.0 * t * t1 - t0;
        t0 = t1;
        t1 = tn;
    }
    
    return dt1;
}

std::vector<double> ChebyshevApproximation::fitComponent(
    const std::vector<double>& data,
    const std::vector<double>& epochs) {
    
    if (data.empty() || data.size() != epochs.size()) {
        throw std::runtime_error("Data e epochs devono avere la stessa lunghezza");
    }
    
    size_t n_data = data.size();
    
    // Costruisci matrice di Vandermonde dei polinomi di Chebyshev
    Eigen::MatrixXd A(n_data, num_coefficients_);
    Eigen::VectorXd b(n_data);
    
    for (size_t i = 0; i < n_data; ++i) {
        // Normalizza l'epoca
        double t = coefficients_[0].normalize(epochs[i]);
        
        // Calcola polinomi di Chebyshev
        for (size_t j = 0; j < num_coefficients_; ++j) {
            A(i, j) = chebyshevT(j, t);
        }
        
        b(i) = data[i];
    }
    
    // Risolvi il sistema ai minimi quadrati
    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
    
    // Estrai i coefficienti
    std::vector<double> coeffs(num_coefficients_);
    for (size_t i = 0; i < num_coefficients_; ++i) {
        coeffs[i] = x(i);
    }
    
    return coeffs;
}

bool ChebyshevApproximation::fit(const std::vector<Eigen::Vector3d>& positions,
                                  double start_epoch,
                                  double end_epoch) {
    
    if (positions.empty()) {
        throw std::runtime_error("Posizioni non possono essere vuote");
        return false;
    }
    
    if (start_epoch >= end_epoch) {
        throw std::runtime_error("start_epoch deve essere minore di end_epoch");
        return false;
    }
    
    start_epoch_ = start_epoch;
    end_epoch_ = end_epoch;
    
    // ============================================================================
    // CORREZIONE IMPORTANTE (4 Dec 2025):
    // I dati in input DEVONO provenire da propagazione con:
    // 1. RKF78 integrator (7th-8th order, tolerance 1e-12)
    // 2. TUTTE le perturbazioni attive (8 pianeti + asteroids + relativity)
    // 3. Frame conversion ECLM→ICRF applicata
    // 4. Coordinate barycentriche ICRF J2000.0
    // ============================================================================
    
    // Prepara i dati per il fitting
    std::vector<double> epochs;
    std::vector<double> x_data, y_data, z_data;
    
    for (size_t i = 0; i < positions.size(); ++i) {
        double t = start_epoch + i * (end_epoch - start_epoch) / (positions.size() - 1);
        epochs.push_back(t);
        x_data.push_back(positions[i].x());
        y_data.push_back(positions[i].y());
        z_data.push_back(positions[i].z());
    }
    
    training_epochs_ = epochs;
    
    // Aggiorna gli intervalli di normalizzazione
    for (size_t i = 0; i < 3; ++i) {
        coefficients_[i].t_min = start_epoch;
        coefficients_[i].t_max = end_epoch;
    }
    
    try {
        // Fitta i tre componenti
        coefficients_[0].coefficients = fitComponent(x_data, epochs);
        coefficients_[1].coefficients = fitComponent(y_data, epochs);
        coefficients_[2].coefficients = fitComponent(z_data, epochs);
        
        // Memorizza i dati di training
        training_data_[0] = Eigen::Map<Eigen::VectorXd>(x_data.data(), x_data.size());
        training_data_[1] = Eigen::Map<Eigen::VectorXd>(y_data.data(), y_data.size());
        training_data_[2] = Eigen::Map<Eigen::VectorXd>(z_data.data(), z_data.size());
        
        fitted_ = true;
        return true;
        
    } catch (const std::exception&) {
        return false;
    }
}

// Implementazione rinviata - vedi dichiarazione nel .h
// Questa viene definita dove AsteroidState è completo

double ChebyshevApproximation::evaluateChebyshev(
    const std::vector<double>& coeffs,
    double t) {
    
    if (coeffs.empty()) return 0.0;
    
    double result = 0.0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * chebyshevT(i, t);
    }
    
    return result;
}

double ChebyshevApproximation::evaluateChebyshevDerivative(
    const std::vector<double>& coeffs,
    double t) {
    
    if (coeffs.empty()) return 0.0;
    
    double result = 0.0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * chebyshevTDerivative(i, t);
    }
    
    return result;
}

Eigen::Vector3d ChebyshevApproximation::evaluatePosition(double epoch) const {
    if (!fitted_) {
        throw std::runtime_error("ChebyshevApproximation: I polinomi non sono stati fittati");
    }
    
    if (epoch < start_epoch_ || epoch > end_epoch_) {
        throw std::out_of_range("Epoca fuori dall'intervallo di approssimazione");
    }
    
    double t = coefficients_[0].normalize(epoch);
    
    return Eigen::Vector3d(
        evaluateChebyshev(coefficients_[0].coefficients, t),
        evaluateChebyshev(coefficients_[1].coefficients, t),
        evaluateChebyshev(coefficients_[2].coefficients, t)
    );
}

Eigen::Vector3d ChebyshevApproximation::evaluateVelocity(double epoch) const {
    if (!fitted_) {
        throw std::runtime_error("ChebyshevApproximation: I polinomi non sono stati fittati");
    }
    
    if (epoch < start_epoch_ || epoch > end_epoch_) {
        throw std::out_of_range("Epoca fuori dall'intervallo di approssimazione");
    }
    
    double t = coefficients_[0].normalize(epoch);
    
    // Derivata rispetto a t (normalizzato)
    double dt_dt_norm = 2.0 / (end_epoch_ - start_epoch_);
    
    // Velocità = d(posizione)/dt = d(posizione)/d(t_norm) * d(t_norm)/dt
    Eigen::Vector3d vel_dt_norm(
        evaluateChebyshevDerivative(coefficients_[0].coefficients, t),
        evaluateChebyshevDerivative(coefficients_[1].coefficients, t),
        evaluateChebyshevDerivative(coefficients_[2].coefficients, t)
    );
    
    return vel_dt_norm * dt_dt_norm;
}

Eigen::Vector3d ChebyshevApproximation::getApproximationError() const {
    if (!fitted_ || training_epochs_.empty()) {
        return Eigen::Vector3d::Zero();
    }
    
    Eigen::Vector3d rms_error = Eigen::Vector3d::Zero();
    
    for (size_t i = 0; i < training_epochs_.size(); ++i) {
        Eigen::Vector3d predicted = evaluatePosition(training_epochs_[i]);
        Eigen::Vector3d actual(
            training_data_[0](i),
            training_data_[1](i),
            training_data_[2](i)
        );
        Eigen::Vector3d error = predicted - actual;
        
        rms_error(0) += error(0) * error(0);
        rms_error(1) += error(1) * error(1);
        rms_error(2) += error(2) * error(2);
    }
    
    size_t n = training_epochs_.size();
    rms_error(0) = std::sqrt(rms_error(0) / n);
    rms_error(1) = std::sqrt(rms_error(1) / n);
    rms_error(2) = std::sqrt(rms_error(2) / n);
    
    return rms_error;
}

double ChebyshevApproximation::evaluateEnergy(double epoch) const {
    Eigen::Vector3d pos = evaluatePosition(epoch);
    Eigen::Vector3d vel = evaluateVelocity(epoch);
    
    double r = pos.norm();
    double v2 = vel.squaredNorm();
    
    // Energia specifica = v²/2 - μ/r
    // Approssimiamo con μ = 1 (unità astronomiche)
    double energy = 0.5 * v2 - 1.0 / r;
    
    return energy;
}

Eigen::Vector3d ChebyshevApproximation::evaluateAngularMomentum(double epoch) const {
    Eigen::Vector3d pos = evaluatePosition(epoch);
    Eigen::Vector3d vel = evaluateVelocity(epoch);
    
    return pos.cross(vel);
}

const ChebyshevCoefficients& ChebyshevApproximation::getCoefficients(size_t component) const {
    if (component > 2) {
        throw std::out_of_range("Component deve essere 0 (X), 1 (Y) o 2 (Z)");
    }
    return coefficients_[component];
}

void ChebyshevApproximation::reset() {
    fitted_ = false;
    training_epochs_.clear();
    for (size_t i = 0; i < 3; ++i) {
        coefficients_[i].coefficients.assign(num_coefficients_, 0.0);
        training_data_[i].resize(0);
    }
    start_epoch_ = 0.0;
    end_epoch_ = 0.0;
}

bool ChebyshevApproximation::saveToFile(const std::string& filename) const {
    if (!fitted_) return false;
    
    std::ofstream file(filename);
    if (!file.is_open()) return false;
    
    file << std::scientific << std::setprecision(15);
    file << "# Coefficienti di Chebyshev per traiettoria asteroidale\n";
    file << "# Epoche di validità (MJD):\n";
    file << "# start_epoch = " << start_epoch_ << "\n";
    file << "# end_epoch = " << end_epoch_ << "\n";
    file << "# num_coefficients = " << num_coefficients_ << "\n";
    file << "\n";
    
    // Salva coefficienti per ogni componente
    for (size_t comp = 0; comp < 3; ++comp) {
        char axis = (comp == 0) ? 'X' : (comp == 1) ? 'Y' : 'Z';
        file << "# Componente " << axis << "\n";
        
        for (size_t i = 0; i < coefficients_[comp].coefficients.size(); ++i) {
            file << coefficients_[comp].coefficients[i] << "\n";
        }
        file << "\n";
    }
    
    return file.good();
}

bool ChebyshevApproximation::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    std::string line;
    double start_epoch, end_epoch;
    size_t num_coeff;
    
    // Leggi header
    while (std::getline(file, line)) {
        if (line.find("start_epoch") != std::string::npos) {
            std::istringstream iss(line);
            std::string hash, key, equals;
            iss >> hash >> key >> equals >> start_epoch;
        }
        if (line.find("end_epoch") != std::string::npos) {
            std::istringstream iss(line);
            std::string hash, key, equals;
            iss >> hash >> key >> equals >> end_epoch;
        }
        if (line.find("num_coefficients") != std::string::npos) {
            std::istringstream iss(line);
            std::string hash, key, equals;
            iss >> hash >> key >> equals >> num_coeff;
        }
        if (line.empty() || line[0] != '#') break;
    }
    
    // Leggi coefficienti
    try {
        for (size_t comp = 0; comp < 3; ++comp) {
            size_t i = 0;
            while (i < num_coefficients_ && std::getline(file, line)) {
                if (line.empty() || line[0] == '#') continue;
                
                try {
                    coefficients_[comp].coefficients[i] = std::stod(line);
                    i++;
                } catch (...) {
                    continue;
                }
            }
        }
        
        start_epoch_ = start_epoch;
        end_epoch_ = end_epoch;
        
        // Update normalization intervals for all components
        for (size_t i = 0; i < 3; ++i) {
            coefficients_[i].t_min = start_epoch;
            coefficients_[i].t_max = end_epoch;
        }

        fitted_ = true;
        return true;
        
    } catch (const std::exception&) {
        return false;
    }
}

std::string ChebyshevApproximation::getStatistics() const {
    std::ostringstream oss;
    
    oss << "=== Statistiche Approssimazione Chebyshev ===\n\n";
    
    if (!fitted_) {
        oss << "Status: NON FITTATO\n";
        return oss.str();
    }
    
    oss << "Status: FITTATO\n";
    oss << "Num coefficienti: " << num_coefficients_ << "\n";
    oss << "Intervallo temporale: MJD " << std::fixed << std::setprecision(2)
        << start_epoch_ << " - " << end_epoch_ << "\n";
    oss << "Durata: " << (end_epoch_ - start_epoch_) << " giorni\n";
    
    if (!training_epochs_.empty()) {
        oss << "Num punti di training: " << training_epochs_.size() << "\n";
        
        Eigen::Vector3d error = getApproximationError();
        oss << "\nErrore RMS (AU):\n";
        oss << "  X: " << std::scientific << std::setprecision(3) << error.x() << "\n";
        oss << "  Y: " << std::scientific << error.y() << "\n";
        oss << "  Z: " << std::scientific << error.z() << "\n";
        oss << "  Media: " << std::scientific << error.norm() / 3.0 << "\n";
    }
    
    return oss.str();
}

} // namespace ioccultcalc
