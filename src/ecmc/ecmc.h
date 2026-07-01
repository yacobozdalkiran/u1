#pragma once

#include <array>
#include <cmath>
#include <random>

#include "../field/field.h"
#include "../geometry/geometry.h"

struct ECMCParams {
    double theta_sample = 100;
    double theta_refresh = 50;
    bool use_topological_lifting = false;
    double eta = 1e-6;
};

// Saves the state of a local Event-Chain
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

struct Distributions {
    std::uniform_int_distribution<int> random_dir;
    std::uniform_int_distribution<int> random_eps;
    // Constructeur 
    Distributions(const ECMCParams& p) : random_dir(0, 1), random_eps(0, 1) {};
};

#pragma omp declare simd
inline void solve_reject_fast(double A, double B, double& gamma, double& reject, int epsilon) {
    // Change A depending on epsilon
    A = (epsilon == -1) ? -A : A;

    // std::hypot est souvent mieux vectorisé par SVML
    double R = std::hypot(A, B);
    double invR = 1.0 / R;
    double period = 2.0 * R;

    double discarded_number = std::floor(gamma / period);
    gamma -= discarded_number * period;

    double phi = std::atan2(-A, B);
    phi += (phi < 0.0) ? (2.0 * M_PI) : 0.0;

    double alpha;
    double p1 = R - A;

    // Le compilateur Intel transforme ce bloc en "masking" SIMD
    if (phi < (M_PI * 0.5) || phi > (M_PI * 1.5)) {
        alpha = (gamma > p1) ? (gamma - p1) * invR - 1.0 : (gamma + A) * invR;
    } else {
        alpha = gamma * invR - 1.0;
    }

    // Clamp (std::clamp est vectorisable en C++17, ou version manuelle)
    alpha = (alpha > 1.0) ? 1.0 : ((alpha < -1.0) ? -1.0 : alpha);

    double theta = phi + std::asin(alpha);

    // Normalisation 2*PI sans if/else
    theta += (theta < 0.0) ? (2.0 * M_PI) : 0.0;
    theta -= (theta >= 2.0 * M_PI) ? (2.0 * M_PI) : 0.0;

    reject = theta + 2.0 * M_PI * discarded_number;
}

void compute_plaquettes(const GaugeField& field, const Geometry& geo, int site, int mu,
                        std::array<double, 2>& list_plaquettes);
void compute_reject_angles_fast(const GaugeField& field, int site, int mu,
                                const std::array<double, 2>& list_plaquettes, int epsilon,
                                const double& beta, std::array<double, 2>& reject_angles,
                                std::mt19937_64& rng);
size_t selectVariable_norev(const std::array<double, 3>& probas, std::mt19937_64& rng);

std::pair<std::pair<int, int>, int> lift_improved_fast_norev(const GaugeField& field,
                                                                const Geometry& geo, int site,
                                                                int mu, int j,
                                                                std::mt19937_64& rng);

// Nouvelles fonctions pour le lifting topologique
void sinkhorn_knopp(double W[4][4], int iterations = 15);
double get_topological_variation(const GaugeField& field, const Geometry& geo, int site, int mu);
std::pair<std::pair<int, int>, int> lift_topological(const GaugeField& field, const Geometry& geo,
                                                   int site, int mu, int j,
                                                   const ECMCParams& params, std::mt19937_64& rng);

void update(GaugeField& field, int site, int mu, double theta, int epsilon);
int random_site(const Geometry& geo, std::mt19937_64& rng);
void ecmc_sample(LocalChainState& state, GaugeField& field, double beta, Distributions& d, const Geometry& geo,
                 const ECMCParams& params, std::mt19937_64& rng);
