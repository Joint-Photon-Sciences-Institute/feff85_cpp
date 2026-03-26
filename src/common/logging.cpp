// Logging — converted from src/COMMON/wlog.f, src/COMMON/chopen.f

#include "logging.hpp"
#include "../par/parallel.hpp"
#include <iostream>
#include <stdexcept>

namespace feff::common {

static Logger g_logger;

Logger& logger() {
    return g_logger;
}

Logger::~Logger() {
    close();
}

void Logger::open(const std::string& log_filename) {
    if (!log_file_.is_open()) {
        log_file_.open(log_filename, std::ios::out);
    }
}

void Logger::close() {
    if (log_file_.is_open()) {
        log_file_.close();
    }
}

void Logger::wlog(const std::string& msg) {
    auto& ps = feff::par::state();

    // Suppress output in sequential loops (par_type == 2)
    if (ps.par_type == 2) return;

    if (msg.empty()) {
        std::cout << "\n";
        if (ps.par_type != 3 && log_file_.is_open()) {
            log_file_ << "\n";
        }
    } else {
        std::cout << msg << "\n";
        if (ps.par_type != 3 && log_file_.is_open()) {
            log_file_ << msg << "\n";
        }
    }
}

void check_file_open(const std::ifstream& f, const std::string& filename,
                     const std::string& module_name) {
    if (!f.is_open()) {
        throw std::runtime_error("Error opening file '" + filename +
                                 "' in module " + module_name);
    }
}

void check_file_open(const std::ofstream& f, const std::string& filename,
                     const std::string& module_name) {
    if (!f.is_open()) {
        throw std::runtime_error("Error opening file '" + filename +
                                 "' in module " + module_name);
    }
}

} // namespace feff::common
