#include "observables.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double average_plaquette(const GaugeField& field, const Geometry& geo) {
    double sum = 0.0;
    for (int site = 0; site < geo.V; ++site) {
        sum += field.plaquette(geo, site);
    }
    return sum / static_cast<double>(geo.V);
}

double plaquette_phase(const GaugeField& field, const Geometry& geo, int site) {
    int site_px = geo.get_neighbor(site, 0, pos);
    int site_pt = geo.get_neighbor(site, 1, pos);
    
    // Angle total de la plaquette : theta_0(x) + theta_1(x+0) - theta_0(x+1) - theta_1(x)
    double phi = field.get_link(site, 0) + field.get_link(site_px, 1) - field.get_link(site_pt, 0) - field.get_link(site, 1);
    
    // Ramener dans [-PI, PI[
    while (phi >= M_PI) phi -= 2.0 * M_PI;
    while (phi < -M_PI) phi += 2.0 * M_PI;
    
    return phi;
}

double topological_charge(const GaugeField& field, const Geometry& geo) {
    double total_phase = 0.0;
    for (int site = 0; site < geo.V; ++site) {
        total_phase += plaquette_phase(field, geo, site);
    }
    return total_phase / (2.0 * M_PI);
}

double topological_susceptibility(const GaugeField& field, const Geometry& geo) {
    double Q = topological_charge(field, geo);
    return (Q * Q) / static_cast<double>(geo.V);
}
