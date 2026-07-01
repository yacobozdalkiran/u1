#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include "../field/field.h"
#include "../geometry/geometry.h"

/**
 * @brief Calculates the average plaquette over the entire lattice.
 */
double average_plaquette(const GaugeField& field, const Geometry& geo);

/**
 * @brief Calculates the phase of a specific plaquette, normalized to [-PI, PI[.
 */
double plaquette_phase(const GaugeField& field, const Geometry& geo, int site);

/**
 * @brief Calculates the total topological charge Q.
 * Q = 1/(2*PI) * sum(normalized_plaquette_phases)
 * Returns a double to monitor numerical precision.
 */
double topological_charge(const GaugeField& field, const Geometry& geo);

/**
 * @brief Calculates the topological susceptibility for the current configuration (Q^2 / V).
 */
double topological_susceptibility(const GaugeField& field, const Geometry& geo);

#endif
