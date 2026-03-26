#pragma once
// Logging and file-open error checking.
// Converted from: src/COMMON/wlog.f, src/COMMON/chopen.f

#include <string>
#include <fstream>
#include <memory>

namespace feff::common {

/// Logger that writes to stdout and a log file (replaces Fortran wlog).
/// Respects par_type from ParallelState:
///   par_type == 2: suppress all output (sequential loop mode)
///   par_type == 3: write to stdout only (no log file)
class Logger {
public:
    Logger() = default;
    ~Logger();

    /// Open the log file. Call once at startup. Does nothing if already open.
    void open(const std::string& log_filename = "log.dat");

    /// Close the log file.
    void close();

    /// Write a message to stdout and log file (like Fortran wlog).
    void wlog(const std::string& msg);

    /// Check if log file is open.
    bool is_open() const { return log_file_.is_open(); }

private:
    std::ofstream log_file_;
};

/// Global logger accessor
Logger& logger();

/// Check that a file stream opened successfully; throw on failure.
/// Replaces Fortran chopen(ios, fname, mod).
void check_file_open(const std::ifstream& f, const std::string& filename,
                     const std::string& module_name);

/// Overload for ofstream
void check_file_open(const std::ofstream& f, const std::string& filename,
                     const std::string& module_name);

} // namespace feff::common
