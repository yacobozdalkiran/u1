#include <filesystem>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../field/field.h"
#include "../geometry/geometry.h"
#include "../heatbath/heatbath.h"
#include "../io/io.h"
#include "../obs/observables.h"

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.txt>" << std::endl;
        return 1;
    }

    HBSimulationParameters params;
    if (!read_input_file_simparams(argv[1], params.sim_params)) {
        return 1;
    }
    if (!read_input_file_hbparams(argv[1], params.hb_params)) {
        return 1;
    }

    SimulationParameters sim_params = params.sim_params;
    HBParams hb_params = params.hb_params;

    print_parameters_sim(params.sim_params);
    print_parameters_hb(hb_params);

    Geometry geo(sim_params.L);
    GaugeField field(geo);
    std::mt19937_64 rng(sim_params.seed);
    field.cold_start();

    std::vector<double> plaquettes;
    std::vector<double> charges;
    std::vector<double> susceptibilities;

    int start_iter = 0;
    std::string base_path =
        sim_params.output_dir + "/" + sim_params.name + "/" + sim_params.name + "_";
    std::string cp_info_path = base_path + "checkpoint_info.txt";
    std::string cp_config_path = base_path + "checkpoint_config.txt";
    std::string cp_rng_path = base_path + "checkpoint_rng.txt";

    std::string res_plaq_path = base_path + "plaquettes.txt";
    std::string res_charge_path = base_path + "top_charge.txt";
    std::string res_susc_path = base_path + "top_susc.txt";

    bool resuming = false;

    if (fs::exists(cp_info_path) && fs::exists(cp_config_path) && fs::exists(cp_rng_path)) {
        std::cout << "Checkpoint found! Resuming simulation...\n";
        if (load_checkpoint_info(cp_info_path, start_iter) &&
            load_configuration(cp_config_path, field, geo) && load_rng_state(cp_rng_path, rng)) {
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
    }

    int target_iter = start_iter + sim_params.n_updates;
    std::cout << "Target total configurations : " << target_iter << "\n\n";

    plaquettes.reserve(sim_params.save_each);
    charges.reserve(sim_params.save_each);
    susceptibilities.reserve(sim_params.save_each);

    std::cout << std::fixed << std::setprecision(6);
    for (int i = start_iter + 1; i <= target_iter; ++i) {
        std::cout << "\n========== Configuration " << i << "/" << target_iter << " ==========\n";
        heatbath_update(field, geo, sim_params.beta, rng, hb_params);

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

            write_vector_to_file(res_plaq_path, plaquettes, true);
            write_vector_to_file(res_charge_path, charges, true);
            write_vector_to_file(res_susc_path, susceptibilities, true);

            plaquettes.clear();
            charges.clear();
            susceptibilities.clear();
        }
    }

    if (!plaquettes.empty()) {
        std::cout << "\nSaving remaining results...\n";
        write_vector_to_file(res_plaq_path, plaquettes, true);
        write_vector_to_file(res_charge_path, charges, true);
        write_vector_to_file(res_susc_path, susceptibilities, true);

        // Final checkpoint to update info to target_iter
        save_configuration(cp_config_path, field, geo);
        save_rng_state(cp_rng_path, rng);
        save_checkpoint_info(cp_info_path, target_iter);
    }

    std::cout << "Done!\n";
    return 0;
}
