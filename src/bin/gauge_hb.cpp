#include <iostream>
#include <random>
#include <iomanip>
#include <vector>
#include <string>
#include <filesystem>

#include "../field/field.h"
#include "../geometry/geometry.h"
#include "../heatbath/heatbath.h"
#include "../obs/observables.h"
#include "../io/io.h"

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.txt>" << std::endl;
        return 1;
    }

    SimulationParameters params;
    if (!read_input_file(argv[1], params)) {
        return 1;
    }

    print_parameters(params);

    Geometry geo(params.L);
    GaugeField field(geo);
    std::mt19937_64 rng(params.seed);
    field.cold_start();

    std::vector<double> plaquettes;
    std::vector<double> charges;
    std::vector<double> susceptibilities;

    int start_iter = 0;
    std::string base_path = params.output_dir + "/" + params.name + "/" + params.name + "_";
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
            load_configuration(cp_config_path, field, geo) &&
            load_rng_state(cp_rng_path, rng)) {
            resuming = true;
            std::cout << "Resuming from iteration " << start_iter << "\n";
        } else {
            std::cerr << "Error loading checkpoint files. Starting from scratch.\n";
            start_iter = 0;
        }
    }

    if (!resuming) {
        write_vector_to_file(res_plaq_path, {}, false);
        write_vector_to_file(res_charge_path, {}, false);
        write_vector_to_file(res_susc_path, {}, false);
    }

    int target_iter = start_iter + params.n_updates;
    std::cout << "Target total configurations : " << target_iter << "\n\n";

    plaquettes.reserve(params.save_each);
    charges.reserve(params.save_each);
    susceptibilities.reserve(params.save_each);

    std::cout << std::fixed << std::setprecision(6);
    for (int i = start_iter + 1; i <= target_iter; ++i) {
        std::cout << "========== Configuration " << i << "/" << target_iter << " ==========\n";
        heatbath_update(field, geo, params.beta, rng);
        
        double plaq = average_plaquette(field, geo);
        double Q = topological_charge(field, geo);
        double chi = topological_susceptibility(field, geo);
        
        plaquettes.push_back(plaq);
        charges.push_back(Q);
        susceptibilities.push_back(chi);
        
        if (i % 10 == 0 || i == target_iter) {
            std::cout << "<P> = " << plaq << " | Q=" << Q << "\n";
        }

        if (i % params.save_each == 0) {
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
