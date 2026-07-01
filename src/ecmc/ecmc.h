#pragma once

#include <array>
#include <cmath>
#include <random>

#include "../field/field.h"
#include "../geometry/geometry.h"

/**
 * @brief Structure to store ECMC params.
 */
struct ECMCParams {
    double theta_sample = 100;
    double theta_refresh = 50;
    int algo = 0;       // 0: classic ecmc 1: topological lifting
    double eta = 1e-6;  // Only used if algo=1
};

/**
 * @brief Structure to store ECMC chain state for persistence.
 */
struct LocalChainState {
    int site;
    int mu;
    int epsilon;
    bool initialized = false;

    // Budget de refresh persistant
    double theta_parcouru_refresh = 0.0;

    // Compteur de lifts
    size_t event_counter = 0;
    size_t lift_counter = 0;
};

/**
 * @brief Structure to store ECMC distributions.
 */
struct Distributions {
    std::uniform_int_distribution<int> random_dir;
    std::uniform_int_distribution<int> random_eps;
    // Constructeur
    Distributions() : random_dir(0, 1), random_eps(0, 1) {};
};

/**
 * @brief Optimized function to solve the reject equation. (algo 0, 1)
 */
void solve_reject_fast(double A, double B, double& gamma, double& reject, int epsilon);

/**
 * @brief Computes the values of the 2 plaquettes attached to a link. (algo 0, 1)
 */
void compute_plaquettes(const GaugeField& field, const Geometry& geo, int site, int mu,
                        std::array<double, 2>& list_plaquettes);

/**
 * @brief Computes the values of the reject angles for each plaquette attached to a link. (algo 0,
 * 1)
 */
void compute_reject_angles_fast(const std::array<double, 2>& list_plaquettes, int epsilon,
                                const double& beta, std::array<double, 2>& reject_angles,
                                std::mt19937_64& rng);

/**
 * @brief Uses the tower of probability method to choose an index in a length 3 array. (algo 0, 1)
 */
size_t selectVariable_norev(const std::array<double, 3>& probas, std::mt19937_64& rng);

/**
 * @brief Updates a link. (algo 0, 1)
 */
void update(GaugeField& field, int site, int mu, double theta, int epsilon);

// =========== Algo 0 : Forward lifts ===========

/**
 * @brief Lift without backtracking from a link. (algo 0)
 */
std::pair<std::pair<int, int>, int> lift_improved_fast_norev(const GaugeField& field,
                                                             const Geometry& geo, int site, int mu,
                                                             int j, std::mt19937_64& rng);

// =========== Algo 1 : Topology-driven lifts ===========

/**
 * @brief Returns the bistochastic version of a 4x4 matrix. (algo 1)
 */
void sinkhorn_knopp(double W[4][4], int iterations = 15);

/**
 * @brief Returns the sum of the absolute value of local topological charge of the 2 plaquettes
 * attached to a link normalized to [0,1]. (algo 1)
 */
double get_topological_variation(const GaugeField& field, const Geometry& geo, int site, int mu);

/**
 * @brief Performs a topology-driven lift from a link. (algo 1)
 */
std::pair<std::pair<int, int>, int> lift_topological(const GaugeField& field, const Geometry& geo,
                                                     int site, int mu, int j,
                                                     const ECMCParams& params,
                                                     std::mt19937_64& rng);

// =========== Algo 2 : Plaquette ECMC ===========

namespace algo2 {

/**
 * @brief computes the 4 plaquettes around the plaquette at a site. Order : right, top, left,
 * bottom.
 */
void compute_plaquettes(const GaugeField& field, const Geometry& geo, int site,
                        std::array<double, 4>& list_plaquettes);

void compute_reject_angles_fast(const std::array<double, 4>& list_plaquettes, int epsilon,
                                const double& beta, std::array<double, 4>& reject_angles,
                                std::mt19937_64& rng);
}  // namespace algo2

// =========== All Algos ===========

/**
 * @brief Selects a random site in the lattice. (all algos)
 */
int random_site(const Geometry& geo, std::mt19937_64& rng);

/**
 * @brief Performs the ECMC algorithm until theta_sample is reached. (all algos)
 */
void ecmc_sample(LocalChainState& state, GaugeField& field, double beta, Distributions& d,
                 const Geometry& geo, const ECMCParams& params, std::mt19937_64& rng);
