#ifndef IO_H
#define IO_H

#include <string>
#include <vector>
#include <random>
#include "../field/field.h"
#include "../geometry/geometry.h"
#include "../ecmc/ecmc.h"
#include "../heatbath/heatbath.h"

/**
 * @brief Structure to hold simulation parameters.
 */
struct SimulationParameters {
    //Core parameters
    int L;
    double beta;
    int n_updates;
    unsigned long long seed;
    //Measures params
    int plaq_each;
    int topo_each;
    //Save parameters
    int save_each;
    std::string output_dir;
    std::string name;
};

/**
 * @brief Structure to hold ECMC simulation parameters.
 */
struct ECMCSimulationParameters {
    SimulationParameters sim_params;
    // ECMC specific
    ECMCParams ecmc_params;
};

/**
 * @brief Structure to hold Heatbath simulation parameters.
 */
struct HBSimulationParameters {
    SimulationParameters sim_params;
    // Heatbath specific
    HBParams hb_params;
};

/**
 * @brief Reads simulation parameters from a text file.
 */
bool read_input_file_simparams(const std::string& filename, SimulationParameters& params);

/**
 * @brief Reads simulation parameters from a text file.
 */
bool read_input_file_ecmcparams(const std::string& filename, ECMCParams& params);

/**
 * @brief Reads simulation parameters from a text file.
 */
bool read_input_file_hbparams(const std::string& filename, HBParams& params);

/**
 * @brief Writes a vector of doubles to a file.
 * @param append If true, appends to the file instead of overwriting.
 */
bool write_vector_to_file(const std::string& filename, const std::vector<double>& data, bool append = false);

/**
 * @brief Loads a vector of doubles from a file.
 */
bool load_vector_from_file(const std::string& filename, std::vector<double>& data);

/**
 * @brief Saves the current gauge field configuration to a file.
 */
bool save_configuration(const std::string& filename, const GaugeField& field, const Geometry& geo);

/**
 * @brief Loads a gauge field configuration from a file.
 */
bool load_configuration(const std::string& filename, GaugeField& field, const Geometry& geo);

/**
 * @brief Saves the state of the random number generator to a file.
 */
bool save_rng_state(const std::string& filename, const std::mt19937_64& rng);

/**
 * @brief Loads the state of the random number generator from a file.
 */
bool load_rng_state(const std::string& filename, std::mt19937_64& rng);

/**
 * @brief Saves current iteration info.
 */
bool save_checkpoint_info(const std::string& filename, int iter);

/**
 * @brief Loads current iteration info.
 */
bool load_checkpoint_info(const std::string& filename, int& iter);

/**
 * @brief Saves the ECMC local chain state.
 */
bool save_ecmc_state(const std::string& filename, const LocalChainState& state);

/**
 * @brief Loads the ECMC local chain state.
 */
bool load_ecmc_state(const std::string& filename, LocalChainState& state);

/**
 * @brief Prints the Simulation paramaters.
 */
void print_parameters_sim(const SimulationParameters& sim_params);

/**
 * @brief Prints the ECMC paramaters.
 */
void print_parameters_ecmc(const ECMCParams& ecmc_params);

/**
 * @brief Prints the Heatbath paramaters.
 */
void print_parameters_hb(const HBParams& hb_params);
#endif
