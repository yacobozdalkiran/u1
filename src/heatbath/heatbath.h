#ifndef HEATBATH_H
#define HEATBATH_H

#include "../field/field.h"
#include "../geometry/geometry.h"
#include <random>


/**
 * @brief Structure to hot heatbath parameters.
 */
struct HBParams {
    int n_sweep;
    int n_hit;
};

/**
 * @brief Performs a Heatbath update on the entire gauge field.
 * 
 * @param field The gauge field to update.
 * @param geo The geometry of the lattice.
 * @param beta The inverse coupling constant.
 * @param rng Random number generator.
 */
void heatbath_update(GaugeField& field, const Geometry& geo, double beta, std::mt19937_64& rng, const HBParams& params);



#endif
