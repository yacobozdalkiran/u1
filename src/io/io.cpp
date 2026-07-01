#include "io.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace fs = std::filesystem;

bool read_input_file_simparams(const std::string& filename, SimulationParameters& params) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open input file " << filename << std::endl;
        return false;
    }

    std::string line;

    bool has_L = false;
    bool has_beta = false;
    bool has_updates = false;
    bool has_out = false;
    bool has_seed = false;
    bool has_save = false;
    bool has_name = false;
    bool has_plaq_each = false;
    bool has_topo_each = false;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key)) continue;

        if (key == "L") {
            if (iss >> params.L) has_L = true;
        } else if (key == "beta") {
            if (iss >> params.beta) has_beta = true;
        } else if (key == "n_updates") {
            if (iss >> params.n_updates) has_updates = true;
        } else if (key == "seed") {
            if (iss >> params.seed) has_seed = true;

        } else if (key == "plaq_each") {
            if (iss >> params.plaq_each) has_plaq_each = true;
        } else if (key == "topo_each") {
            if (iss >> params.topo_each) has_topo_each = true;

        } else if (key == "name") {
            if (iss >> params.name) has_name = true;
        } else if (key == "output_dir") {
            if (iss >> params.output_dir) has_out = true;

        } else if (key == "save_each") {
            if (iss >> params.save_each) has_save = true;
        }
    }

    if (!has_L || !has_beta || !has_updates || !has_out || !has_seed || !has_save || !has_name) {
        std::cerr << "Error: Missing parameters in input file. Need L, beta, n_updates, "
                     "output_dir, name, seed, and save_each."
                  << std::endl;
        return false;
    }

    return true;
}

bool read_input_file_ecmcparams(const std::string& filename, ECMCParams& params) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open input file " << filename << std::endl;
        return false;
    }

    std::string line;
    bool has_param_theta_sample = false;
    bool has_param_theta_refresh = false;
    bool has_use_topological_lifting = false;
    bool has_eta = false;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key)) continue;

        if (key == "theta_sample") {
            if (iss >> params.theta_sample) has_param_theta_sample = true;
        } else if (key == "theta_refresh") {
            if (iss >> params.theta_refresh) has_param_theta_refresh = true;
        } else if (key == "use_topological_lifting") {
            if (iss >> params.use_topological_lifting) has_use_topological_lifting = true;
        } else if (key == "eta") {
            if (iss >> params.eta) has_eta = true;
        }
    }

    if (!has_param_theta_sample || !has_param_theta_refresh || !has_use_topological_lifting ||
        !has_eta) {
        std::cerr << "Error: Missing parameters in input file. Need theta_sample, "
                     "theta_refresh, use_topological_lifting, eta."
                  << std::endl;
        return false;
    }

    return true;
}

bool read_input_file_hbparams(const std::string& filename, HBParams& params) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open input file " << filename << std::endl;
        return false;
    }

    std::string line;
    bool has_n_hit = false;
    bool has_n_sweep = false;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key)) continue;

        if (key == "n_hit") {
            if (iss >> params.n_hit) has_n_hit = true;
        } else if (key == "n_sweep") {
            if (iss >> params.n_sweep) has_n_sweep = true;
        } else if (key == "use_topological_lifting") {
        }
    }

    if (!has_n_hit || !has_n_sweep) {
        std::cerr << "Error: Missing parameters in input file. Need n_sweep, n_hit." << std::endl;
        return false;
    }

    return true;
}

bool write_vector_to_file(const std::string& filename, const std::vector<double>& data,
                          bool append) {
    fs::path filepath(filename);
    if (filepath.has_parent_path()) {
        fs::create_directories(filepath.parent_path());
    }

    std::ofstream outfile;
    if (append) {
        outfile.open(filename, std::ios::app);
    } else {
        outfile.open(filename);
    }

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open output file " << filename << std::endl;
        return false;
    }

    if (data.empty()) return true;  // Just opening/creating the file if empty

    outfile << std::scientific << std::setprecision(10);
    for (const auto& val : data) {
        outfile << val << "\n";
    }
    return true;
}

bool load_vector_from_file(const std::string& filename, std::vector<double>& data) {
    std::ifstream infile(filename);
    if (!infile.is_open()) return false;

    data.clear();
    double val;
    while (infile >> val) {
        data.push_back(val);
    }
    return true;
}

bool save_configuration(const std::string& filename, const GaugeField& field, const Geometry& geo) {
    fs::path filepath(filename);
    if (filepath.has_parent_path()) {
        fs::create_directories(filepath.parent_path());
    }

    std::ofstream outfile(filename);
    if (!outfile.is_open()) return false;

    outfile << std::scientific << std::setprecision(15);
    for (int site = 0; site < geo.V; ++site) {
        for (int mu = 0; mu < 2; ++mu) {
            outfile << field.get_link(site, mu) << "\n";
        }
    }
    return true;
}

bool load_configuration(const std::string& filename, GaugeField& field, const Geometry& geo) {
    std::ifstream infile(filename);
    if (!infile.is_open()) return false;

    double val;
    for (int site = 0; site < geo.V; ++site) {
        for (int mu = 0; mu < 2; ++mu) {
            if (!(infile >> val)) return false;
            field.set_link(site, mu, val);
        }
    }
    return true;
}

bool save_rng_state(const std::string& filename, const std::mt19937_64& rng) {
    fs::path filepath(filename);
    if (filepath.has_parent_path()) {
        fs::create_directories(filepath.parent_path());
    }

    std::ofstream outfile(filename);
    if (!outfile.is_open()) return false;

    outfile << rng;
    return true;
}

bool load_rng_state(const std::string& filename, std::mt19937_64& rng) {
    std::ifstream infile(filename);
    if (!infile.is_open()) return false;

    infile >> rng;
    return true;
}

bool save_checkpoint_info(const std::string& filename, int iter) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) return false;
    outfile << iter;
    return true;
}

bool load_checkpoint_info(const std::string& filename, int& iter) {
    std::ifstream infile(filename);
    if (!infile.is_open()) return false;
    infile >> iter;
    return true;
}

bool save_ecmc_state(const std::string& filename, const LocalChainState& state) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) return false;

    outfile << state.site << " " << state.mu << " " << state.epsilon << " " << state.initialized
            << " " << std::scientific << std::setprecision(15) << state.theta_parcouru_refresh
            << " " << state.event_counter << " " << state.lift_counter << "\n";
    return true;
}

bool load_ecmc_state(const std::string& filename, LocalChainState& state) {
    std::ifstream infile(filename);
    if (!infile.is_open()) return false;

    if (!(infile >> state.site >> state.mu >> state.epsilon >> state.initialized >>
          state.theta_parcouru_refresh >> state.event_counter >> state.lift_counter)) {
        return false;
    }
    return true;
}
void print_parameters_sim(const SimulationParameters& sim_params) {
    std::cout << "--- Simulation Parameters ---\n";
    std::cout << "L                 : " << sim_params.L << "\n";
    std::cout << "Beta              : " << sim_params.beta << "\n";
    std::cout << "Updates (Sweeps)  : " << sim_params.n_updates << "\n";
    std::cout << "Seed              : " << sim_params.seed << "\n";
    std::cout << "Plaq. meas every  : " << sim_params.plaq_each << "\n";
    std::cout << "Topo. meas every  : " << sim_params.topo_each << "\n";
    std::cout << "Save Every        : " << sim_params.save_each << "\n";
    std::cout << "Name              : " << sim_params.name << "\n";
    std::cout << "Output Dir        : " << sim_params.output_dir << "/" << sim_params.name << "\n";
    std::cout << "----------------------------------\n\n";
}
void print_parameters_ecmc(const ECMCParams& ecmc_params) {
    std::cout << "--- ECMC Parameters ---\n";
    std::cout << "Theta sample  : " << ecmc_params.theta_sample << "\n";
    std::cout << "Theta refresh : " << ecmc_params.theta_refresh<< "\n";
    std::cout << "----------------------------------\n\n";
}
