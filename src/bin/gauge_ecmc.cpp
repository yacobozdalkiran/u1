#include <filesystem>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../ecmc/ecmc.h"
#include "../field/field.h"
#include "../geometry/geometry.h"
#include "../io/io.h"
#include "../obs/observables.h"

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.txt>" << std::endl;
        return 1;
    }

    ECMCSimulationParameters params;

    if (!read_input_file_simparams(argv[1], params.sim_params)) {
        std::cerr << "Error : could not read simulation params from inputs.\n";
        return 1;
    }

    if (!read_input_file_ecmcparams(argv[1], params.ecmc_params)) {
        std::cerr << "Error : could not read ecmc params from inputs.\n";
        return 1;
    }

    // Vérification explicite pour ECMC
    if (params.ecmc_params.theta_sample <= 0.0 || params.ecmc_params.theta_refresh <= 0.0) {
        std::cerr << "Error: ECMC requires 'theta_sample' and 'theta_refresh' to be defined in the "
                     "input file and positive.\n";
        return 1;
    }

    SimulationParameters sim_params = params.sim_params;
    ECMCParams ecmc_params = params.ecmc_params;

    print_parameters_sim(sim_params);
    print_parameters_ecmc(ecmc_params);
    Geometry geo(sim_params.L);
    GaugeField field(geo);
    std::mt19937_64 rng(sim_params.seed);
    field.cold_start();

    Distributions dists;
    LocalChainState state;

    std::vector<double> plaquettes;
    std::vector<double> charges;
    std::vector<double> susceptibilities;
    std::vector<double> event_counts;
    std::vector<double> lift_counts;
    std::vector<double> lambdas;

    int start_iter = 0;
    std::string base_path =
        sim_params.output_dir + "/" + sim_params.name + "/" + sim_params.name + "_";
    std::string cp_info_path = base_path + "checkpoint_info.txt";
    std::string cp_config_path = base_path + "checkpoint_config.txt";
    std::string cp_rng_path = base_path + "checkpoint_rng.txt";
    std::string cp_ecmc_path = base_path + "checkpoint_ecmc.txt";

    std::string res_plaq_path = base_path + "plaquettes.txt";
    std::string res_charge_path = base_path + "top_charge.txt";
    std::string res_susc_path = base_path + "top_susc.txt";
    std::string res_event_path = base_path + "event_counts.txt";
    std::string res_lift_path = base_path + "lift_counts.txt";
    std::string res_lambda_path = base_path + "lambdas.txt";

    bool resuming = false;

    if (fs::exists(cp_info_path) && fs::exists(cp_config_path) && fs::exists(cp_rng_path) &&
        fs::exists(cp_ecmc_path)) {
        std::cout << "Checkpoint found! Resuming simulation...\n";
        if (load_checkpoint_info(cp_info_path, start_iter) &&
            load_configuration(cp_config_path, field, geo) && load_rng_state(cp_rng_path, rng) &&
            load_ecmc_state(cp_ecmc_path, state)) {
            resuming = true;
            std::cout << "Resuming from iteration " << start_iter << "\n";
        } else {
            std::cerr << "Error loading checkpoint files.\n";
            return 1;
        }
    }

    if (!resuming) {
        write_vector_to_file(res_plaq_path, {}, false);
        write_vector_to_file(res_charge_path, {}, false);
        write_vector_to_file(res_susc_path, {}, false);
        write_vector_to_file(res_event_path, {}, false);
        write_vector_to_file(res_lift_path, {}, false);
        write_vector_to_file(res_lambda_path, {}, false);
    }

    int target_iter = start_iter + sim_params.n_updates;
    std::cout << "Target total configurations : " << target_iter << "\n\n";

    plaquettes.reserve(sim_params.save_each);
    charges.reserve(sim_params.save_each);
    susceptibilities.reserve(sim_params.save_each);
    event_counts.reserve(sim_params.save_each);
    lift_counts.reserve(sim_params.save_each);
    lambdas.reserve(sim_params.save_each);

    std::cout << std::fixed << std::setprecision(6);
    for (int i = start_iter + 1; i <= target_iter; ++i) {
        std::cout << "\n========== Configuration " << i << "/" << target_iter << " ==========\n";
        if (ecmc_params.algo==0 or ecmc_params.algo==1) ecmc_sample(state, field, sim_params.beta, dists, geo, ecmc_params, rng);
        if (ecmc_params.algo==2) algo2::ecmc_sample(state, field, sim_params.beta, dists, geo, ecmc_params, rng);

        double lambda = ecmc_params.theta_sample / static_cast<double>(state.event_counter);
        std::cout << "Lambda : " << lambda << "\n";
        event_counts.push_back(static_cast<double>(state.event_counter));
        lift_counts.push_back(static_cast<double>(state.lift_counter));
        lambdas.push_back(lambda);

        if (i % sim_params.plaq_each == 0) {
            double plaq = average_plaquette(field, geo);
            std::cout << "<P> = " << plaq << "\n";
            plaquettes.push_back(plaq);
        }
        if (i % sim_params.topo_each == 0) {
            double Q = topological_charge(field, geo);
            double chi = (Q * Q) / static_cast<double>(geo.V);
            std::cout << "Q = " << Q << "\n";
            std::cout << "Chi = " << chi << "\n";
            charges.push_back(Q);
            susceptibilities.push_back(chi);
        }


        if (i % sim_params.save_each == 0) {
            save_configuration(cp_config_path, field, geo);
            save_rng_state(cp_rng_path, rng);
            save_checkpoint_info(cp_info_path, i);
            save_ecmc_state(cp_ecmc_path, state);

            write_vector_to_file(res_plaq_path, plaquettes, true);
            write_vector_to_file(res_charge_path, charges, true);
            write_vector_to_file(res_susc_path, susceptibilities, true);
            write_vector_to_file(res_event_path, event_counts, true);
            write_vector_to_file(res_lift_path, lift_counts, true);
            write_vector_to_file(res_lambda_path, lambdas, true);

            plaquettes.clear();
            charges.clear();
            susceptibilities.clear();
            event_counts.clear();
            lift_counts.clear();
            lambdas.clear();
        }
    }

    if (!plaquettes.empty()) {
        std::cout << "\nSaving remaining results...\n";
        write_vector_to_file(res_plaq_path, plaquettes, true);
        write_vector_to_file(res_charge_path, charges, true);
        write_vector_to_file(res_susc_path, susceptibilities, true);
        write_vector_to_file(res_event_path, event_counts, true);
        write_vector_to_file(res_lift_path, lift_counts, true);
        write_vector_to_file(res_lambda_path, lambdas, true);

        save_configuration(cp_config_path, field, geo);
        save_rng_state(cp_rng_path, rng);
        save_checkpoint_info(cp_info_path, target_iter);
        save_ecmc_state(cp_ecmc_path, state);
    }

    std::cout << "Done!\n";
    return 0;
}
