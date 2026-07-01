#include "ecmc.h"

#include <algorithm>

#include "../obs/observables.h"

void compute_plaquettes(const GaugeField& field, const Geometry& geo, int site, int mu,
                        std::array<double, 2>& list_plaquettes) {
    int x = site;                            // x
    int xmu = geo.get_neighbor(x, mu, pos);  // x+mu
    int nu = 1 - mu;
    // Plaquette forward
    int xnu = geo.get_neighbor(x, nu, pos);  // x+nu
    const auto& U0 = field.get_link(xmu, nu);
    const auto& U1 = field.get_link(xnu, mu);
    const auto& U2 = field.get_link(x, nu);
    list_plaquettes[0] = field.get_link(site, mu) + U0 - U2 - U1;

    // Plaquette backward
    int xmunu = geo.get_neighbor(xmu, nu, neg);  // x+mu-nu
    int xmnu = geo.get_neighbor(x, nu, neg);     // x-nu
    auto V0 = field.get_link(xmunu, nu);
    auto V1 = field.get_link(xmnu, mu);
    auto V2 = field.get_link(xmnu, nu);
    list_plaquettes[1] = field.get_link(site, mu) - V1 - V0 + V2;
}

void compute_reject_angles_fast(const std::array<double, 2>& list_plaquettes, int epsilon,
                                const double& beta, std::array<double, 2>& reject_angles,
                                std::mt19937_64& rng) {
    static std::uniform_real_distribution<double> unif01_g(0.0, 1.0);
    for (int i = 0; i < 2; i++) {
        double gamma = -std::log(1.0 - unif01_g(rng));
        double A = -epsilon * beta * std::cos(list_plaquettes[i]);
        double B = epsilon * beta * std::sin(list_plaquettes[i]);
        solve_reject_fast(A, B, gamma, reject_angles[i], epsilon);
    }
}

size_t selectVariable_norev(const std::array<double, 3>& probas, std::mt19937_64& rng) {
    static std::uniform_real_distribution<double> unif01(0.0, 1.0);
    double r = unif01(rng);
    if (r < probas[0]) return 0;
    if (r < probas[0] + probas[1]) return 1;
    return 2;
}

size_t selectVariable4(const std::array<double, 4>& probas, std::mt19937_64& rng) {
    static std::uniform_real_distribution<double> unif01(0.0, 1.0);
    double r = unif01(rng);
    double cumulative = 0;
    for (size_t i = 0; i < 4; ++i) {
        cumulative += probas[i];
        if (r < cumulative) return i;
    }
    return 3;
}

std::pair<std::pair<int, int>, int> lift_improved_fast_norev(const Geometry& geo, int site, int mu,
                                                             int j, std::mt19937_64& rng) {
    std::array<double, 3> probas = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    size_t index_lift = selectVariable_norev(probas, rng);
    int epsilon = 1;
    if (j == 0) {
        if (index_lift == 0) epsilon = -1;
    } else {
        if (index_lift == 2) epsilon = -1;
    }
    return std::make_pair(geo.get_link_staple(site, mu, j, index_lift), epsilon);
}

void sinkhorn_knopp(double W[4][4], int iterations) {
    for (int it = 0; it < iterations; ++it) {
        for (int i = 0; i < 4; ++i) {
            double sum = 0;
            for (int j = 0; j < 4; ++j) sum += W[i][j];
            if (sum > 1e-15) {
                double invSum = 1.0 / sum;
                for (int j = 0; j < 4; ++j) W[i][j] *= invSum;
            }
        }
        for (int j = 0; j < 4; ++j) {
            double sum = 0;
            for (int i = 0; i < 4; ++i) sum += W[i][j];
            if (sum > 1e-15) {
                double invSum = 1.0 / sum;
                for (int i = 0; i < 4; ++i) W[i][j] *= invSum;
            }
        }
    }
}

double get_topological_variation(const GaugeField& field, const Geometry& geo, int site, int mu) {
    double phi_fwd = std::abs(plaquette_phase(field, geo, site));
    int nu = 1 - mu;
    int site_prev_nu = geo.get_neighbor(site, nu, neg);
    double phi_bwd = std::abs(plaquette_phase(field, geo, site_prev_nu));
    return (phi_fwd + phi_bwd) / (2.0 * M_PI);
}

std::pair<std::pair<int, int>, int> lift_topological(const GaugeField& field, const Geometry& geo,
                                                     int site, int mu, int j,
                                                     const ECMCParams& params,
                                                     std::mt19937_64& rng) {
    std::pair<int, int> links[4];
    links[0] = {site, mu};
    for (int k = 0; k < 3; ++k) links[k + 1] = geo.get_link_staple(site, mu, j, k);

    double a[4];
    for (int k = 0; k < 4; ++k)
        a[k] = get_topological_variation(field, geo, links[k].first, links[k].second);

    double W[4][4];
    for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 4; ++k) {
            if (i == k)
                W[i][k] = 0.0;
            else
                W[i][k] = a[i] + a[k] + params.eta;
        }
    }
    sinkhorn_knopp(W, 15);

    std::array<double, 4> probas;
    for (int k = 0; k < 4; ++k) probas[k] = W[0][k];
    size_t next_idx = selectVariable4(probas, rng);

    int signs[2][4] = {{1, 1, -1, -1}, {1, -1, -1, 1}};
    int epsilon_factor = 1;
    if (signs[j][0] == signs[j][next_idx]) epsilon_factor = -1;

    return std::make_pair(links[next_idx], epsilon_factor);
}

void update(GaugeField& field, int site, int mu, double theta, int epsilon) {
    double new_value = field.get_link(site, mu) + epsilon * theta;
    field.set_link(site, mu, new_value);
}

int random_site(const Geometry& geo, std::mt19937_64& rng) {
    static std::uniform_int_distribution<int> random_coord(0, geo.V - 1);
    return random_coord(rng);
}

void ecmc_sample(LocalChainState& state, GaugeField& field, double beta, Distributions& d,
                 const Geometry& geo, const ECMCParams& params, std::mt19937_64& rng) {
    if (!state.initialized) {
        state.site = random_site(geo, rng);
        state.mu = d.random_dir(rng);
        state.epsilon = 2 * d.random_eps(rng) - 1;
        state.theta_parcouru_refresh = 0.0;
        state.event_counter = 0;
        state.lift_counter = 0;
        state.initialized = true;
    }

    int site_current = state.site;
    int mu_current = state.mu;
    int epsilon_current = state.epsilon;
    size_t event_counter = 0;
    size_t lift_counter = 0;

    double theta_sample = params.theta_sample;
    double theta_refresh = params.theta_refresh;
    double theta_parcouru_sample = 0.0;
    double theta_parcouru_refresh = state.theta_parcouru_refresh;

    std::array<double, 2> reject_angles;
    std::array<double, 2> list_plaquettes;

    while (true) {
        compute_plaquettes(field, geo, site_current, mu_current, list_plaquettes);
        compute_reject_angles_fast(list_plaquettes, epsilon_current, beta, reject_angles, rng);

        int j = 0;
        double theta_reject = reject_angles[0];
        if (reject_angles[1] < theta_reject) {
            theta_reject = reject_angles[1];
            j = 1;
        }

        double dist_to_sample = theta_sample - theta_parcouru_sample;
        double dist_to_refresh = theta_refresh - theta_parcouru_refresh;
        double theta_step = std::min({theta_reject, dist_to_sample, dist_to_refresh});

        if (theta_step == dist_to_sample) {
            update(field, site_current, mu_current, dist_to_sample, epsilon_current);
            event_counter++;
            state.site = site_current;
            state.mu = mu_current;
            state.epsilon = epsilon_current;
            state.theta_parcouru_refresh = theta_parcouru_refresh + dist_to_sample;
            state.event_counter = event_counter;
            state.lift_counter = lift_counter;
            return;
        } else if (theta_step == dist_to_refresh) {
            update(field, site_current, mu_current, dist_to_refresh, epsilon_current);
            event_counter++;
            theta_parcouru_sample += dist_to_refresh;
            theta_parcouru_refresh = 0.0;
            site_current = random_site(geo, rng);
            mu_current = d.random_dir(rng);
            epsilon_current = 2 * d.random_eps(rng) - 1;
        } else {
            update(field, site_current, mu_current, theta_reject, epsilon_current);
            event_counter++;
            theta_parcouru_sample += theta_reject;
            theta_parcouru_refresh += theta_reject;

            std::pair<std::pair<int, int>, int> l;
            if (params.algo == 1) {
                l = lift_topological(field, geo, site_current, mu_current, j, params, rng);
            } else if (params.algo == 0) {
                l = lift_improved_fast_norev(geo, site_current, mu_current, j, rng);
            }
            lift_counter++;
            site_current = l.first.first;
            mu_current = l.first.second;
            epsilon_current = epsilon_current * l.second;
        }
    }
}

void solve_reject_fast(double A, double B, double& gamma, double& reject, int epsilon) {
    // Change A depending on epsilon
    A = (epsilon == -1) ? -A : A;

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

void algo2::compute_plaquettes(const GaugeField& field, const Geometry& geo, int site,
                               std::array<double, 4>& list_plaquettes) {
    int site_right = geo.get_neighbor(site, 0, pos);
    int site_top = geo.get_neighbor(site, 1, pos);
    int site_left = geo.get_neighbor(site, 0, neg);
    int site_bottom = geo.get_neighbor(site, 1, neg);
    list_plaquettes[0] = field.plaquette(geo, site_right);
    list_plaquettes[1] = field.plaquette(geo, site_top);
    list_plaquettes[2] = field.plaquette(geo, site_left);
    list_plaquettes[3] = field.plaquette(geo, site_bottom);
}

void algo2::compute_reject_angles_fast(const std::array<double, 4>& list_plaquettes,
                                       const double& beta, std::array<double, 4>& reject_angles,
                                       std::mt19937_64& rng, int epsilon) {
    static std::uniform_real_distribution<double> unif01_g(0.0, 1.0);
    std::array<double, 4> coeffs_plaquette = {-1.0, 1.0, 1.0, -1.0};
    for (int i = 0; i < 4; i++) {
        double gamma = -std::log(1.0 - unif01_g(rng));
        double A = -beta * std::cos(list_plaquettes[i]);
        double B = epsilon * coeffs_plaquette[i] * beta * std::sin(list_plaquettes[i]);
        solve_reject_fast(A, B, gamma, reject_angles[i], epsilon);
    }
}

void algo2::update(GaugeField& field, const Geometry& geo, int site, double theta_update) {
    int site_px = geo.get_neighbor(site, 0, pos);
    int site_pt = geo.get_neighbor(site, 1, pos);
    field.add_to_link(site, 0, theta_update);
    field.add_to_link(site_px, 1, theta_update);
    field.add_to_link(site, 1, theta_update);
    field.add_to_link(site_pt, 0, theta_update);
}

void algo2::ecmc_sample(LocalChainState& state, GaugeField& field, double beta, Distributions& d, const Geometry& geo,
                        const ECMCParams& params, std::mt19937_64& rng) {
    if (!state.initialized) {
        state.site = random_site(geo, rng);
        state.mu = -1;
        state.epsilon = 2 * d.random_eps(rng) - 1;
        state.theta_parcouru_refresh = 0.0;
        state.event_counter = 0;
        state.lift_counter = 0;
        state.initialized = true;
    }

    int site_current = state.site;
    int epsilon_current = state.epsilon;
    size_t event_counter = 0;
    size_t lift_counter = 0;

    double theta_sample = params.theta_sample;
    double theta_refresh = params.theta_refresh;
    double theta_parcouru_sample = 0.0;
    double theta_parcouru_refresh = state.theta_parcouru_refresh;

    std::array<double, 4> reject_angles;
    std::array<double, 4> list_plaquettes;

    while (true) {
        compute_plaquettes(field, geo, site_current, list_plaquettes);
        compute_reject_angles_fast(list_plaquettes, beta, reject_angles, rng, epsilon_current);

        int j = 0;
        double theta_reject = reject_angles[0];
        for (int i = 1; i < 4; i++) {
            if (reject_angles[i] < theta_reject) {
                theta_reject = reject_angles[i];
                j = i;
            }
        }

        double dist_to_sample = theta_sample - theta_parcouru_sample;
        double dist_to_refresh = theta_refresh - theta_parcouru_refresh;
        double theta_step = std::min({theta_reject, dist_to_sample, dist_to_refresh});

        if (theta_step == dist_to_sample) {
            //Sample
            update(field, geo, site_current, epsilon_current*dist_to_sample);
            event_counter++;
            state.site = site_current;
            state.epsilon = epsilon_current;
            state.theta_parcouru_refresh = theta_parcouru_refresh + dist_to_sample;
            state.event_counter = event_counter;
            state.lift_counter = lift_counter;
            return;
        } else if (theta_step == dist_to_refresh) {
            //Refresh
            update(field, geo, site_current, epsilon_current*dist_to_refresh);
            event_counter++;
            theta_parcouru_sample += dist_to_refresh;
            theta_parcouru_refresh = 0.0;
            site_current = random_site(geo, rng);
            epsilon_current = 2*d.random_eps(rng)-1;
        } else {
            //Update+lift
            update(field, geo, site_current, epsilon_current*theta_reject);
            event_counter++;
            theta_parcouru_sample += theta_reject;
            theta_parcouru_refresh += theta_reject;

            // Lift
            if (j == 0) {
                site_current = geo.get_neighbor(site_current, 0, pos);
                site_current = geo.get_neighbor(site_current, 0, pos);
            }
            if (j == 1) {
                site_current = geo.get_neighbor(site_current, 1, pos);
                site_current = geo.get_neighbor(site_current, 1, pos);
            }
            if (j == 2) {
                site_current = geo.get_neighbor(site_current, 0, neg);
                site_current = geo.get_neighbor(site_current, 0, neg);
            }
            if (j == 3) {
                site_current = geo.get_neighbor(site_current, 1, neg);
                site_current = geo.get_neighbor(site_current, 1, neg);
            }
            lift_counter++;
        }
    }
}
void algo2::solve_reject_fast(double A, double B, double& gamma, double& reject, int epsilon) {
    // Change A depending on epsilon
    A = (epsilon == -1) ? -A : A;

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
