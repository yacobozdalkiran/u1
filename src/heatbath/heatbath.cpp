#include "heatbath.h"

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Von Mises sampling using Best-Fisher algorithm.
 * Samples from P(theta) ~ exp(k * cos(theta)).
 */
static double sample_von_mises(double k, std::mt19937_64& rng) {
    if (k < 1e-8) {
        std::uniform_real_distribution<double> dist(-M_PI, M_PI);
        return dist(rng);
    }
    double a = 1.0 + std::sqrt(1.0 + 4.0 * k * k);
    double b = (a - std::sqrt(2.0 * a)) / (2.0 * k);
    double r = (1.0 + b * b) / (2.0 * b);

    std::uniform_real_distribution<double> dist_u(0.0, 1.0);

    while (true) {
        double u1 = dist_u(rng);
        double z = std::cos(M_PI * u1);
        double f = (1.0 + r * z) / (r + z);
        double c = k * (r - f);

        double u2 = dist_u(rng);
        if (u2 < c * (2.0 - c) || std::log(c / u2) + 1.0 - c >= 0.0) {
            double u3 = dist_u(rng);
            double theta = std::acos(f);
            if (u3 < 0.5) theta = -theta;
            return theta;
        }
    }
}

static void heatbath_hit(GaugeField& field, const Geometry& geo, int site, int mu, double beta,
                         std::mt19937_64& rng) {
    double X = 0, Y = 0;
    if (mu == 0) {
        // Staples for U_0(x)
        // Term 1 from P(x) = theta_0(x) + Delta1
        int x_p0 = geo.get_neighbor(site, 0, pos);
        int x_p1 = geo.get_neighbor(site, 1, pos);

        double delta1 = field.get_link(x_p0, 1) - field.get_link(x_p1, 0) - field.get_link(site, 1);

        // Term 2 from P(x-1) = Delta_old - theta_0(x)
        int x_m1 = geo.get_neighbor(site, 1, neg);
        int x_p0_m1 = geo.get_neighbor(x_p0, 1, neg);

        double delta_old =
            field.get_link(x_m1, 0) + field.get_link(x_p0_m1, 1) - field.get_link(x_m1, 1);

        X = std::cos(delta1) + std::cos(delta_old);
        Y = std::sin(delta_old) - std::sin(delta1);
    } else {
        // Staples for U_1(x)
        // Term 1 from P(x) = Delta_p1 - theta_1(x)
        int x_p0 = geo.get_neighbor(site, 0, pos);
        int x_p1 = geo.get_neighbor(site, 1, pos);

        double delta_p1 =
            field.get_link(site, 0) + field.get_link(x_p0, 1) - field.get_link(x_p1, 0);

        // Term 2 from P(x-0) = Delta_p2 + theta_1(x)
        int x_m0 = geo.get_neighbor(site, 0, neg);
        int x_m0_p1 = geo.get_neighbor(x_m0, 1, pos);

        double delta_p2 =
            field.get_link(x_m0, 0) - field.get_link(x_m0_p1, 0) - field.get_link(x_m0, 1);

        X = std::cos(delta_p1) + std::cos(delta_p2);
        Y = std::sin(delta_p1) - std::sin(delta_p2);
    }

    double rho = std::sqrt(X * X + Y * Y);
    double alpha = std::atan2(Y, X);

    double theta_new = alpha + sample_von_mises(beta * rho, rng);

    // Normalize to [-PI, PI]
    theta_new = std::atan2(std::sin(theta_new), std::cos(theta_new));

    field.set_link(site, mu, theta_new);
}

void heatbath_update(GaugeField& field, const Geometry& geo, double beta, std::mt19937_64& rng,
                     const HBParams& params) {
    int V = geo.V;
    for (int site = 0; site < V; ++site) {
        for (int mu = 0; mu < 2; ++mu) {
            for (int hits = 0; hits < params.n_hit; hits++)
                heatbath_hit(field, geo, site, mu, beta, rng);
        }
    }
}
