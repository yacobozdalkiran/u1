#include "field.h"

#include <cmath>

GaugeField::GaugeField(const Geometry& geo) : L(geo.L), V(L * L), links(2 * V) {}

void GaugeField::hot_start(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> dist_mpipi(-M_PI, M_PI);
    for (int site = 0; site < V; site++) {
        set_link(site, 0, dist_mpipi(rng));
        set_link(site, 1, dist_mpipi(rng));
    }
}

void GaugeField::cold_start() {
    for (int site = 0; site < V; site++) {
        set_link(site, 0, 0);
        set_link(site, 1, 0);
    }
}

double GaugeField::get_link(int site, int mu) const { return links[site * 2 + mu]; }

void GaugeField::set_link(int site, int mu, double new_value) { links[site * 2 + mu] = new_value; }

double GaugeField::plaquette(const Geometry& geo, int site) const {
    int site_px = geo.get_neighbor(site, 0, pos);
    int site_pt = geo.get_neighbor(site, 1, pos);

    double phi =
        get_link(site, 0) + get_link(site_px, 1) - get_link(site_pt, 0) - get_link(site, 1);
    return std::cos(phi);
}

void GaugeField::add_to_link(int site, int mu, double value) { links[site * 2 + mu] += value; }

double GaugeField::plaquette_phi(const Geometry& geo, int site, int mu, bool forward) const {
    int x = site;  // x
    int nu = 1 - mu;
    int xmu = geo.get_neighbor(x, mu, pos);  // x+mu
    if (forward) {
        // Plaquette forward
        int xnu = geo.get_neighbor(x, nu, pos);  // x+nu
        const auto& U0 = get_link(xmu, nu);
        const auto& U1 = get_link(xnu, mu);
        const auto& U2 = get_link(x, nu);
        return get_link(site, mu) + U0 - U2 - U1;
    } else {
        // Plaquette backward
        int xmunu = geo.get_neighbor(xmu, nu, neg);  // x+mu-nu
        int xmnu = geo.get_neighbor(x, nu, neg);     // x-nu
        auto V0 = get_link(xmunu, nu);
        auto V1 = get_link(xmnu, mu);
        auto V2 = get_link(xmnu, nu);
        return get_link(site, mu) - V1 - V0 + V2;
    }
}
