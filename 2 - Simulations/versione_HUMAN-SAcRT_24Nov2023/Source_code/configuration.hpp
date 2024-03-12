#pragma once

#include <string>
#include <fstream>
#include <iostream>

class Configuration {
    const char* RAND_COND_PREFIX = "rand_cond";
    const char* SIM_PARAM_PREFIX = "sim_param";
    const char* ATRIAL_GEOMETRY_PREFIX = "atrial_geometry";
    const char* GAP_JUNC_PREFIX = "gap_junc";
    const char* SAVED_STATES_PREFIX = "saved_states";
    const char* VM_PREFIX = "Vm";
    const char* T_PREFIX = "t";

    std::string files_dir;
    std::string suffix;

    std::fstream get_configuration_file_stream(const std::string& file_prefix, bool create_if_none=false) const {
        std::string filepath = files_dir + file_prefix + suffix;
        std::fstream filestream = create_if_none
            ? std::fstream(filepath, std::fstream::out)
            : std::fstream(filepath);

        if (!filestream)
        {
            std::cout << "Cannot open param file " << filepath << "\n";
            abort();
        }

        return filestream;
    }

    /*
    std::ofstream get_configuration_file_ostream(const std::string& file_prefix, bool create_if_none=false) const {
        std::string filepath = files_dir + file_prefix + suffix;
        std::ofstream filestream = create_if_none
            ? std::ofstream(filepath, std::ofstream::out)
            : std::ofstream(filepath);

        if (!filestream)
        {
            std::cout << "Cannot open param file " << filepath << "\n";
            abort();
        }

        return filestream;
    }
*/

public:
    Configuration(const std::string& _files_dir, const std::string& _suffix)
        : files_dir(_files_dir), suffix(_suffix) {};

    std::fstream get_sim_param_file() const         { return get_configuration_file_stream(SIM_PARAM_PREFIX);       }
    std::fstream get_atrial_geometry_file() const   { return get_configuration_file_stream(ATRIAL_GEOMETRY_PREFIX); }
    std::fstream get_rand_cond_file() const         { return get_configuration_file_stream(RAND_COND_PREFIX);       }
    std::fstream get_gap_junc_file() const          { return get_configuration_file_stream(GAP_JUNC_PREFIX);        }
    std::fstream get_Vm_file() const                { return get_configuration_file_stream(VM_PREFIX, true);        }
    std::fstream get_t_file() const                 { return get_configuration_file_stream(T_PREFIX, true);         }
    std::fstream get_saved_states_file() const      { return get_configuration_file_stream(SAVED_STATES_PREFIX, true);    }
};

