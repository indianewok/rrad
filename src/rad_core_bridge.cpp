#include <Rcpp.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef RRAD_EMBEDDED
#include "vendor/rad_core/include/rad/embedded_console.hpp"
#endif

// RAD uses Striped Smith-Waterman for its adapter fallback.  R only compiles
// translation units directly under src/, so compile the vendored C and C++
// implementations into this one bridge translation unit.
extern "C" {
#include "vendor/rad_core/include/ssw/ssw.c"
}
#include "vendor/rad_core/include/ssw/ssw_cpp.cpp"

#include "vendor/rad_core/include/rad/rad_headers.h"

namespace fs = std::filesystem;
namespace core = rrad_rad_core;

namespace {

constexpr const char* kCoreVersion = "1.0.2";
constexpr int kCoreApiVersion = 1;
constexpr int kResourceSchema = 1;
constexpr const char* kSourceCommit =
    "24d5ce47e172222c96cd7dd19e94758d32aa5009";
constexpr const char* kEmbeddedSourceDigest =
    "ab7272977588d6c33db6d27ebae43010a75507e1ecc18be69faa1dce23f5c5a8";
constexpr const char* kResourceDigest =
    "b3081d5fd40d2d9309108fa79dd591d286e3ec1422582b0c520712659ffec4fb";

bool has_value(const Rcpp::List& options, const char* name) {
    return options.containsElementNamed(name) &&
           !Rf_isNull(options[name]);
}

std::string string_option(const Rcpp::List& options, const char* name,
                          const std::string& fallback = "") {
    if (!has_value(options, name)) return fallback;
    return Rcpp::as<std::string>(options[name]);
}

bool bool_option(const Rcpp::List& options, const char* name, bool fallback) {
    if (!has_value(options, name)) return fallback;
    return Rcpp::as<bool>(options[name]);
}

int int_option(const Rcpp::List& options, const char* name, int fallback) {
    if (!has_value(options, name)) return fallback;
    const double value = Rcpp::as<double>(options[name]);
    if (!std::isfinite(value) || std::floor(value) != value ||
        value < static_cast<double>(std::numeric_limits<int>::min()) ||
        value > static_cast<double>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument(std::string("`") + name +
                                    "` must be a finite whole number");
    }
    return static_cast<int>(value);
}

double double_option(const Rcpp::List& options, const char* name,
                     double fallback) {
    if (!has_value(options, name)) return fallback;
    const double value = Rcpp::as<double>(options[name]);
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string("`") + name +
                                    "` must be finite");
    }
    return value;
}

std::optional<int> optional_int(const Rcpp::List& options,
                                const char* name) {
    if (!has_value(options, name)) return std::nullopt;
    return int_option(options, name, 0);
}

std::int64_t int64_option(const Rcpp::List& options, const char* name,
                          std::int64_t fallback) {
    if (!has_value(options, name)) return fallback;
    const double value = Rcpp::as<double>(options[name]);
    // INT64_MAX rounds up to 2^63 when converted to double on platforms where
    // long double has no extra precision. Compare against the exact exclusive
    // upper bound so the subsequent integer cast is always defined.
    const double int64_exclusive_upper = std::ldexp(1.0, 63);
    if (!std::isfinite(value) || std::floor(value) != value ||
        value < -1.0 || value >= int64_exclusive_upper) {
        throw std::invalid_argument(std::string("`") + name +
                                    "` must be -1 or a non-negative whole number");
    }
    return static_cast<std::int64_t>(value);
}

class EnvironmentGuard {
public:
    EnvironmentGuard(std::string name, const std::string& value)
        : name_(std::move(name)) {
        if (const char* current = std::getenv(name_.c_str())) {
            previous_ = std::string(current);
        }
#ifdef _WIN32
        if (_putenv_s(name_.c_str(), value.c_str()) != 0) {
#else
        if (setenv(name_.c_str(), value.c_str(), 1) != 0) {
#endif
            throw std::runtime_error("could not set environment variable " +
                                     name_);
        }
    }

    ~EnvironmentGuard() {
#ifdef _WIN32
        (void)_putenv_s(name_.c_str(), previous_ ? previous_->c_str() : "");
#else
        if (previous_) {
            (void)setenv(name_.c_str(), previous_->c_str(), 1);
        } else {
            (void)unsetenv(name_.c_str());
        }
#endif
    }

private:
    std::string name_;
    std::optional<std::string> previous_;
};

class TemporaryFileGuard {
public:
    TemporaryFileGuard() = default;
    TemporaryFileGuard(const TemporaryFileGuard&) = delete;
    TemporaryFileGuard& operator=(const TemporaryFileGuard&) = delete;

    ~TemporaryFileGuard() {
        if (!path_) return;
        std::error_code error;
        fs::remove(*path_, error);
    }

    void reset(const std::string& path) { path_ = fs::path(path); }

private:
    std::optional<fs::path> path_;
};

class OpenMpThreadGuard {
public:
    explicit OpenMpThreadGuard(int requested)
        : previous_(omp_get_max_threads()) {
        omp_set_num_threads(requested);
    }

    ~OpenMpThreadGuard() noexcept { omp_set_num_threads(previous_); }

private:
    int previous_;
};

class StreamCapture {
public:
    std::string str() const { return capture_.str(); }

private:
    rrad_embedded_console::ScopedCapture capture_;
};

std::string installed_package_resource_dir() {
    Rcpp::Environment base = Rcpp::Environment::base_namespace();
    Rcpp::Function system_file = base["system.file"];

    Rcpp::CharacterVector found =
        system_file("rad/resources", Rcpp::_["package"] = "rrad",
                    Rcpp::_["mustWork"] = false);
    if (found.size() == 1 && found[0] != NA_STRING) {
        const std::string path = Rcpp::as<std::string>(found[0]);
        if (!path.empty() && fs::is_directory(path)) {
            return fs::canonical(path).string();
        }
    }
    return "";
}

std::string canonical_bundled_resource_dir() {
    // Only an installed package tree has the identity covered by the fixed
    // resource digest.  A cwd-relative source tree could be modified without
    // changing that handshake, so it must never be used as a fallback.
    return installed_package_resource_dir();
}

std::string normalize_resource_dir(const std::string& requested) {
    const std::string bundled = canonical_bundled_resource_dir();
    if (bundled.empty()) {
        throw std::runtime_error(
            "rrad's bundled RAD resources are unavailable; reinstall rrad");
    }

    fs::path bundled_path(bundled);
    if (!fs::is_directory(bundled_path / "read_layout") ||
        !fs::is_directory(bundled_path / "wl")) {
        throw std::runtime_error(
            "rrad's bundled RAD resource directory is incomplete: " +
            bundled_path.string());
    }

    if (!requested.empty()) {
        fs::path path(requested);
        if (path.filename() != "resources" &&
            fs::is_directory(path / "resources")) {
            path /= "resources";
        }
        std::error_code error;
        const fs::path requested_path = fs::canonical(path, error);
        if (error || requested_path != bundled_path) {
            throw std::invalid_argument(
                "`resource_dir` must identify rrad's canonical bundled RAD "
                "resource directory");
        }
    }
    return bundled;
}

using NamedPath = std::pair<std::string, fs::path>;

void append_whitelist_inputs(std::vector<NamedPath>& inputs,
                             const std::string& spec,
                             const std::string& description) {
    if (spec.empty()) return;
    const auto parts = core::whitelist_utils::parse_whitelist_specs(spec);
    for (std::size_t index = 0; index < parts.size(); ++index) {
        const fs::path resolved =
            core::whitelist_utils::kit_to_path(parts[index]);
        std::error_code error;
        if (!fs::is_regular_file(resolved, error) || error) {
            throw std::runtime_error(
                "whitelist input is not a regular file: " +
                resolved.string());
        }
        inputs.emplace_back(
            description + (parts.size() == 1
                               ? ""
                               : " component " + std::to_string(index + 1)),
            resolved);
    }
}

void reject_input_output_collisions(
    const std::vector<NamedPath>& inputs,
    const std::vector<NamedPath>& outputs) {
    for (const auto& input : inputs) {
        const fs::path input_path =
            fs::absolute(input.second).lexically_normal();
        for (const auto& output : outputs) {
            const fs::path output_path =
                fs::absolute(output.second).lexically_normal();
            bool collision = input_path == output_path;

            if (!collision) {
                std::error_code input_error;
                std::error_code output_error;
                const bool input_exists = fs::exists(input_path, input_error);
                const bool output_exists = fs::exists(output_path, output_error);
                if (input_error || output_error) {
                    throw std::runtime_error(
                        "could not verify RAD input/output path separation");
                }
                if (input_exists && output_exists) {
                    std::error_code equivalent_error;
                    collision = fs::equivalent(
                        input_path, output_path, equivalent_error);
                    if (equivalent_error) {
                        throw std::runtime_error(
                            "could not verify RAD input/output file identity");
                    }
                }
            }

            if (collision) {
                throw std::invalid_argument(
                    "RAD output " + output.first + " would overwrite " +
                    input.first + ": " + input_path.string());
            }
        }
    }
}

void require_nonempty_regular_file(const fs::path& path,
                                   const char* description) {
    std::error_code error;
    if (!fs::is_regular_file(path, error) || error ||
        fs::file_size(path, error) == 0 || error) {
        throw std::runtime_error(std::string(description) +
                                 " is missing or empty: " + path.string());
    }
}

struct TransactionOutput {
    std::string label;
    fs::path staged;
    fs::path target;
    bool publish = true;
};

void preflight_transaction_targets(const std::vector<NamedPath>& outputs);

// The lock pathname is intentionally persistent. Removing it after unlocking
// permits two waiters to lock different inodes at the same pathname.
class DemuxCommitLock {
public:
    explicit DemuxCommitLock(const fs::path& path) : path_(path) {
        int flags = O_CREAT | O_RDWR;
#ifdef O_CLOEXEC
        flags |= O_CLOEXEC;
#endif
#ifdef O_NOFOLLOW
        flags |= O_NOFOLLOW;
#endif
        descriptor_ = ::open(path_.string().c_str(), flags,
                             S_IRUSR | S_IWUSR);
        if (descriptor_ < 0) {
            throw std::runtime_error(
                "could not open RAD demux commit lock " + path_.string() +
                ": " + std::strerror(errno));
        }

        struct stat information {};
        if (::fstat(descriptor_, &information) != 0) {
            const std::string reason = std::strerror(errno);
            close_descriptor();
            throw std::runtime_error(
                "invalid RAD demux commit lock " + path_.string() + ": " +
                reason);
        }
        if (!S_ISREG(information.st_mode)) {
            close_descriptor();
            throw std::runtime_error(
                "invalid RAD demux commit lock " + path_.string() +
                ": lock path is not a regular file");
        }
        if (::fchmod(descriptor_, S_IRUSR | S_IWUSR) != 0) {
            const std::string reason = std::strerror(errno);
            close_descriptor();
            throw std::runtime_error(
                "could not secure RAD demux commit lock " + path_.string() +
                ": " + reason);
        }

        struct flock lock {};
        lock.l_type = F_WRLCK;
        lock.l_whence = SEEK_SET;
        lock.l_start = 0;
        lock.l_len = 0;
        while (::fcntl(descriptor_, F_SETLKW, &lock) != 0) {
            if (errno == EINTR) {
                try {
                    Rcpp::checkUserInterrupt();
                } catch (...) {
                    close_descriptor();
                    throw;
                }
                continue;
            }
            const std::string reason = std::strerror(errno);
            close_descriptor();
            throw std::runtime_error(
                "could not lock RAD demux outputs with " + path_.string() +
                ": " + reason);
        }
        locked_ = true;
    }

    DemuxCommitLock(const DemuxCommitLock&) = delete;
    DemuxCommitLock& operator=(const DemuxCommitLock&) = delete;

    ~DemuxCommitLock() { close_descriptor(); }

private:
    void close_descriptor() noexcept {
        if (descriptor_ < 0) return;
        if (locked_) {
            struct flock lock {};
            lock.l_type = F_UNLCK;
            lock.l_whence = SEEK_SET;
            lock.l_start = 0;
            lock.l_len = 0;
            (void)::fcntl(descriptor_, F_SETLK, &lock);
        }
        (void)::close(descriptor_);
        descriptor_ = -1;
        locked_ = false;
    }

    fs::path path_;
    int descriptor_ = -1;
    bool locked_ = false;
};

class SameFilesystemOutputTransaction {
public:
    explicit SameFilesystemOutputTransaction(const fs::path& output_base) {
        const fs::path normalized = fs::absolute(output_base).lexically_normal();
        if (normalized.filename().empty() || normalized.filename() == "." ||
            normalized.filename() == "..") {
            throw std::invalid_argument(
                "RAD output prefix must have a non-empty filename");
        }

        const fs::path parent = normalized.parent_path();
        std::error_code error;
        fs::create_directories(parent, error);
        if (error || !fs::is_directory(parent)) {
            throw std::runtime_error(
                "could not create RAD output directory: " + parent.string() +
                (error ? " (" + error.message() + ")" : ""));
        }

        static std::atomic<unsigned long long> counter{0};
        const auto nonce = static_cast<unsigned long long>(
            std::chrono::steady_clock::now().time_since_epoch().count());
        for (unsigned int attempt = 0; attempt < 100; ++attempt) {
            const auto serial = counter.fetch_add(1, std::memory_order_relaxed);
            stage_dir_ =
                parent /
                (".rrad-demux-stage-" + std::to_string(nonce) + "-" +
                 std::to_string(serial));
            errno = 0;
            if (::mkdir(stage_dir_.string().c_str(), S_IRWXU) == 0) {
                if (::chmod(stage_dir_.string().c_str(), S_IRWXU) != 0) {
                    const std::string reason = std::strerror(errno);
                    std::error_code cleanup_error;
                    fs::remove(stage_dir_, cleanup_error);
                    throw std::runtime_error(
                        "could not secure RAD output staging directory: " +
                        reason);
                }
                stage_base_ = stage_dir_ / normalized.filename();
                lock_path_ =
                    parent /
                    ("." + normalized.filename().string() +
                     ".rrad-demux-commit.lock");
                return;
            }
            const int mkdir_error = errno;
            if (mkdir_error != EEXIST) {
                throw std::runtime_error(
                    "could not create RAD output staging directory: " +
                    std::string(std::strerror(mkdir_error)));
            }
        }
        throw std::runtime_error(
            "could not allocate a unique RAD output staging directory");
    }

    SameFilesystemOutputTransaction(
        const SameFilesystemOutputTransaction&) = delete;
    SameFilesystemOutputTransaction& operator=(
        const SameFilesystemOutputTransaction&) = delete;

    ~SameFilesystemOutputTransaction() {
        if (preserve_stage_) return;
        std::error_code error;
        fs::remove_all(stage_dir_, error);
    }

    const fs::path& stage_base() const { return stage_base_; }
    const fs::path& lock_path() const { return lock_path_; }

    void commit(const std::vector<TransactionOutput>& outputs) {
        struct Backup {
            fs::path target;
            fs::path path;
            bool active = false;
        };
        std::vector<Backup> backups;
        backups.reserve(outputs.size());
        for (std::size_t index = 0; index < outputs.size(); ++index) {
            backups.push_back({
                outputs[index].target,
                stage_dir_ / ("previous-" + std::to_string(index)), false});
        }
        // Fixed-size flags make the post-rename bookkeeping non-allocating.
        // This matters most under memory pressure: every original output must
        // already have a known recovery path before the first mutation.
        std::vector<unsigned char> installed(outputs.size(), 0);

        // Everything must be complete before the first published pathname is
        // touched. A publish=false item explicitly removes a stale artifact as
        // part of the same commit, preventing debug/scan leftovers from being
        // mixed with a newer run.
        for (const auto& output : outputs) {
            std::error_code error;
            const bool staged_exists = fs::exists(output.staged, error);
            if (error) {
                throw std::runtime_error(
                    "could not inspect staged RAD output " + output.label +
                    ": " + error.message());
            }
            if (output.publish) {
                if (!staged_exists ||
                    !fs::is_regular_file(output.staged, error) || error) {
                    throw std::runtime_error(
                        "staged RAD output is missing or invalid (" +
                        output.label + "): " + output.staged.string());
                }
            } else if (staged_exists) {
                throw std::runtime_error(
                    "unexpected staged RAD output was produced (" +
                    output.label + "): " + output.staged.string());
            }
        }

        DemuxCommitLock commit_lock(lock_path_);
        std::vector<NamedPath> targets;
        targets.reserve(outputs.size());
        for (const auto& output : outputs) {
            targets.emplace_back(output.label, output.target);
        }
        // Generation may take hours, so repeat target validation while holding
        // the publication lock before moving any existing pathname.
        preflight_transaction_targets(targets);

        std::exception_ptr original_error;
        // Until either commit or rollback completes, the staging directory is
        // the recovery record. Keep it even if rollback reporting itself fails.
        preserve_stage_ = true;
        try {
            for (std::size_t index = 0; index < outputs.size(); ++index) {
                std::error_code error;
                if (!fs::exists(outputs[index].target, error)) {
                    if (error) {
                        throw std::runtime_error(
                            "could not inspect existing RAD output " +
                            outputs[index].label + ": " + error.message());
                    }
                    continue;
                }
                if (std::rename(outputs[index].target.string().c_str(),
                                backups[index].path.string().c_str()) != 0) {
                    throw std::runtime_error(
                        "could not preserve existing RAD output " +
                        outputs[index].target.string() + ": " +
                        std::strerror(errno));
                }
                backups[index].active = true;
            }

            for (std::size_t index = 0; index < outputs.size(); ++index) {
                const auto& output = outputs[index];
                if (!output.publish) continue;
                if (std::rename(output.staged.string().c_str(),
                                output.target.string().c_str()) != 0) {
                    throw std::runtime_error(
                        "could not atomically publish RAD output " +
                        output.target.string() + ": " +
                        std::strerror(errno));
                }
                installed[index] = 1;
            }
        } catch (...) {
            original_error = std::current_exception();
        }

        if (original_error) {
            std::vector<std::string> rollback_failures;
            for (std::size_t offset = outputs.size(); offset > 0; --offset) {
                const std::size_t index = offset - 1;
                if (!installed[index]) continue;
                std::error_code error;
                const bool removed = fs::remove(outputs[index].target, error);
                std::error_code inspect_error;
                const bool still_exists =
                    fs::exists(outputs[index].target, inspect_error);
                if (error || inspect_error || (!removed && still_exists)) {
                    rollback_failures.push_back(
                        "could not remove newly installed output " +
                        outputs[index].target.string() +
                        (error ? ": " + error.message()
                               : inspect_error
                                     ? ": " + inspect_error.message()
                                     : ""));
                }
            }
            for (std::size_t offset = backups.size(); offset > 0; --offset) {
                Backup& backup = backups[offset - 1];
                if (!backup.active) continue;
                if (std::rename(backup.path.string().c_str(),
                                backup.target.string().c_str()) != 0) {
                    rollback_failures.push_back(
                        "could not restore original output " +
                        backup.target.string() + " from preserved backup " +
                        backup.path.string() + ": " +
                        std::strerror(errno));
                } else {
                    backup.active = false;
                }
            }
            if (!rollback_failures.empty()) {
                std::string original_message = "unknown publication failure";
                try {
                    std::rethrow_exception(original_error);
                } catch (const std::exception& error) {
                    original_message = error.what();
                } catch (...) {
                }
                std::ostringstream message;
                message << original_message
                        << "; RAD demux rollback was incomplete";
                for (const auto& failure : rollback_failures) {
                    message << "; " << failure;
                }
                message << "; preserved staging directory "
                        << stage_dir_.string()
                        << " contains any unrestored original outputs and "
                           "requires manual recovery";
                throw std::runtime_error(message.str());
            }
            preserve_stage_ = false;
            std::rethrow_exception(original_error);
        }
        preserve_stage_ = false;
    }

private:
    fs::path stage_dir_;
    fs::path stage_base_;
    fs::path lock_path_;
    bool preserve_stage_ = false;
};

void preflight_transaction_targets(const std::vector<NamedPath>& outputs) {
    for (const auto& output : outputs) {
        const fs::path target = fs::absolute(output.second).lexically_normal();
        std::error_code error;
        const bool exists = fs::exists(target, error);
        if (error) {
            throw std::runtime_error("could not inspect RAD output " +
                                     output.first + ": " + error.message());
        }
        if (!exists) continue;
        if (!fs::is_regular_file(target, error) || error) {
            if (output.first == "layout cache") {
                throw std::runtime_error(
                    "Failed to install layout output: target is not a regular "
                    "file: " + target.string());
            }
            throw std::runtime_error(
                "output write failed: primary demultiplexed gzip is missing or "
                "empty because RAD output " + output.first +
                " is not a regular file: " + target.string());
        }
    }
}

std::string replace_all_copy(std::string value, const std::string& from,
                             const std::string& to) {
    if (from.empty() || from == to) return value;
    std::size_t offset = 0;
    while ((offset = value.find(from, offset)) != std::string::npos) {
        value.replace(offset, from.size(), to);
        offset += to.size();
    }
    return value;
}

class GzipReadGuard {
public:
    explicit GzipReadGuard(gzFile file) : file_(file) {}
    GzipReadGuard(const GzipReadGuard&) = delete;
    GzipReadGuard& operator=(const GzipReadGuard&) = delete;

    ~GzipReadGuard() {
        if (file_) (void)gzclose(file_);
    }

    gzFile get() const { return file_; }

    int close() {
        if (!file_) return Z_OK;
        gzFile file = file_;
        file_ = nullptr;
        return gzclose(file);
    }

private:
    gzFile file_;
};

void verify_primary_gzip(const fs::path& path, bool require_payload) {
    std::error_code error;
    if (!fs::is_regular_file(path, error) || error ||
        fs::file_size(path, error) == 0 || error) {
        throw std::runtime_error(
            "primary demultiplexed gzip is missing or empty: " +
            path.string());
    }

    GzipReadGuard input(gzopen(path.string().c_str(), "rb"));
    if (!input.get()) {
        throw std::runtime_error(
            "could not open primary demultiplexed gzip: " + path.string());
    }

    std::array<unsigned char, 64 * 1024> buffer{};
    bool has_payload = false;
    int count = 0;
    std::size_t bytes_since_interrupt = 0;
    constexpr std::size_t kInterruptInterval = 32U * 1024U * 1024U;
    while ((count = gzread(input.get(), buffer.data(),
                           static_cast<unsigned int>(buffer.size()))) > 0) {
        has_payload = true;
        bytes_since_interrupt += static_cast<std::size_t>(count);
        if (bytes_since_interrupt >= kInterruptInterval) {
            Rcpp::checkUserInterrupt();
            bytes_since_interrupt = 0;
        }
    }
    if (count < 0) {
        int zlib_code = Z_OK;
        const char* detail = gzerror(input.get(), &zlib_code);
        const std::string message = detail ? detail : "";
        (void)input.close();
        throw std::runtime_error(
            "primary demultiplexed gzip failed validation: " + path.string() +
            (message.empty() ? "" : " (" + message + ")"));
    }
    const int close_code = input.close();
    if (close_code != Z_OK) {
        throw std::runtime_error(
            "primary demultiplexed gzip failed to close after validation: " +
            path.string() + " (zlib code " + std::to_string(close_code) + ")");
    }
    if (require_payload && !has_payload) {
        throw std::runtime_error(
            "primary demultiplexed gzip contains no records: " +
            path.string());
    }
}

std::string companion_position_map(const std::string& layout_path) {
    static const std::string suffix = "_layout.csv";
    if (layout_path.size() <= suffix.size() ||
        layout_path.compare(layout_path.size() - suffix.size(), suffix.size(),
                            suffix) != 0) {
        return "";
    }
    return layout_path.substr(0, layout_path.size() - suffix.size()) +
           "_position_map.csv";
}

#include "rad_scan_wl_bridge.hpp"

struct AutoWhitelistResult {
    std::string detected_path;
    std::string reference_spec;
    NativeScanRun scan;
    NativeScanConfig config;
    std::string captured_log;

    std::string demux_spec() const {
        return reference_spec.empty() ? detected_path
                                      : reference_spec + ":" + detected_path;
    }
};

AutoWhitelistResult run_auto_whitelist(
    const core::ReadLayout& layout, const std::string& fastq_path,
    const std::string& outbase, const std::string& reference_whitelist,
    const std::string& adapter_override,
    const std::optional<int>& barcode_length_override,
    double max_error, size_t max_reads, size_t chunk_size,
    int requested_threads, int effective_threads,
    core::scan_whitelist_selection selection, bool verbose,
    std::function<void()> boundary_poll) {
    const core::ReadElement* barcode = nullptr;
    int barcode_count = 0;
    for (const auto& element : layout.by_order()) {
        if (element.global_class == "barcode" &&
            !core::seq_utils::is_rc(element.class_id)) {
            ++barcode_count;
            if (!barcode) barcode = &element;
        }
    }
    if (barcode_count == 0) {
        throw std::runtime_error(
            "auto-whitelist requires a layout with a barcode element");
    }
    if (barcode_count > 1) {
        throw std::runtime_error(
            "auto-whitelist currently supports single-barcode layouts only");
    }

    std::string adapter = adapter_override;
    if (adapter.empty()) {
        bool adapter_found = false;
        int best_order = 0;
        for (const auto& element : layout.by_order()) {
            if (element.type == "static" &&
                !core::seq_utils::is_rc(element.class_id) &&
                !element.seq.empty() && element.order < barcode->order &&
                (!adapter_found || element.order > best_order)) {
                adapter_found = true;
                best_order = element.order;
                adapter = element.seq;
            }
        }
    }
    if (adapter.empty()) {
        throw std::runtime_error(
            "auto-whitelist could not derive an upstream adapter from the layout");
    }

    const int barcode_length = barcode_length_override.value_or(
        barcode->expected_length.value_or(
            barcode->length_candidates.empty()
                ? 16
                : barcode->length_candidates.back()));
    const std::string reference = reference_whitelist.empty()
                                      ? barcode->whitelist_path
                                      : reference_whitelist;

    const int scan_limit =
        max_reads == static_cast<size_t>(-1) ||
                max_reads > static_cast<size_t>(std::numeric_limits<int>::max())
            ? 0
            : static_cast<int>(max_reads);

    NativeScanConfig config;
    config.input = fastq_path;
    config.output_prefix = outbase + "_scanwl";
    config.log_path = outbase + "_scan_wl.log";
    config.adapter_seq = adapter;
    config.whitelist = reference;
    config.barcode_length = barcode_length;
    config.max_reads = scan_limit;
    config.max_error = max_error;
    config.requested_threads = requested_threads;
    config.effective_threads = effective_threads;
    config.chunk_size = static_cast<int>(chunk_size);
    config.selection = selection;
    config.verbose = verbose;

    NativeScanRun scan;
    std::string captured_log;
    {
        OpenMpThreadGuard scan_thread_guard(effective_threads);
        StreamCapture scan_capture;
        scan = run_scan_whitelist(config, std::move(boundary_poll));
        captured_log = public_scan_log(scan_capture.str(), scan);
    }
    if (scan.stats.selected_barcodes == 0 ||
        !fs::is_regular_file(scan.whitelist) ||
        fs::file_size(scan.whitelist) == 0) {
        throw std::runtime_error("auto-whitelist selected no barcodes");
    }
    return {scan.whitelist, reference, std::move(scan), std::move(config),
            std::move(captured_log)};
}

Rcpp::List reformat_stats_as_list(const core::reformat_result& run,
                                  int requested_threads,
                                  int effective_threads) {
    const auto& stats = run.stats;
    return Rcpp::List::create(
        Rcpp::_["threads_requested"] = requested_threads,
        Rcpp::_["threads_effective"] = effective_threads,
        Rcpp::_["records_read"] =
            static_cast<double>(stats.records_read),
        Rcpp::_["records_written"] =
            static_cast<double>(stats.records_written),
        Rcpp::_["records_skipped_missing_cb"] =
            static_cast<double>(stats.records_skipped_missing_cb),
        Rcpp::_["records_reformatted"] =
            static_cast<double>(stats.records_reformatted),
        Rcpp::_["coordinate_mapped"] =
            static_cast<double>(stats.coordinate_mapped),
        Rcpp::_["coordinate_unmapped"] =
            static_cast<double>(stats.coordinate_unmapped),
        Rcpp::_["empty_sequences_skipped"] =
            static_cast<double>(stats.empty_sequences_skipped),
        Rcpp::_["split_barcodes"] =
            static_cast<double>(stats.split_barcodes),
        Rcpp::_["gzip_members_written"] =
            static_cast<double>(stats.gzip_members_written));
}

Rcpp::CharacterVector reformat_split_paths(
    const core::reformat_result& run) {
    Rcpp::CharacterVector paths(run.split_paths.size());
    Rcpp::CharacterVector names(run.split_paths.size());
    for (R_xlen_t index = 0;
         index < static_cast<R_xlen_t>(run.split_paths.size()); ++index) {
        const auto& item = run.split_paths[static_cast<std::size_t>(index)];
        names[index] = item.first;
        paths[index] = item.second;
    }
    paths.attr("names") = names;
    return paths;
}

bool fail_reformat_result_construction() noexcept {
    const char* value = std::getenv(
        "RRAD_TEST_REFORMAT_FAIL_RESULT_CONSTRUCTION");
    return value && std::strcmp(value, "1") == 0;
}

Rcpp::List core_info_impl() {
    const std::string resources = canonical_bundled_resource_dir();
#if defined(__APPLE__)
    const char* platform = "darwin";
#elif defined(__linux__)
    const char* platform = "linux";
#elif defined(_WIN32)
    const char* platform = "windows";
#else
    const char* platform = "unknown";
#endif

    return Rcpp::List::create(
        Rcpp::_["core_version"] = kCoreVersion,
        Rcpp::_["core_api_version"] = kCoreApiVersion,
        Rcpp::_["source_commit"] = kSourceCommit,
        Rcpp::_["embedded_source_digest"] = kEmbeddedSourceDigest,
        Rcpp::_["resource_schema"] = kResourceSchema,
        Rcpp::_["resource_digest"] = kResourceDigest,
        Rcpp::_["features"] = Rcpp::CharacterVector::create(
            "demux", "file-fastq", "gzip", "builtin-layouts",
            "custom-layout", "custom-whitelist", "auto-whitelist",
            "debug-output", "seqspec", "scan-wl", "reformat",
            "remote-resource-cache"),
        Rcpp::_["build"] = Rcpp::List::create(
            Rcpp::_["platform"] = platform,
            Rcpp::_["compiler"] = __VERSION__,
            Rcpp::_["cxx_standard"] = 17,
#ifdef _OPENMP
            Rcpp::_["openmp"] = true,
            Rcpp::_["openmp_version"] = _OPENMP,
#else
            Rcpp::_["openmp"] = false,
            Rcpp::_["openmp_version"] = R_NilValue,
#endif
            Rcpp::_["resources_available"] = !resources.empty(),
            Rcpp::_["resource_dir"] = resources),
        Rcpp::_["vendored"] = true);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List rad_core_info_cpp() {
    return core_info_impl();
}

// [[Rcpp::export(name = ".rad_any_file_equivalent_cpp")]]
bool rad_any_file_equivalent_cpp(Rcpp::CharacterVector candidates,
                                 Rcpp::CharacterVector protected_paths) {
    try {
        for (R_xlen_t candidate_index = 0;
             candidate_index < candidates.size(); ++candidate_index) {
            if (Rcpp::CharacterVector::is_na(candidates[candidate_index])) {
                throw std::invalid_argument("candidate paths cannot contain NA");
            }
            const fs::path candidate = fs::absolute(
                Rcpp::as<std::string>(candidates[candidate_index]))
                                           .lexically_normal();
            for (R_xlen_t protected_index = 0;
                 protected_index < protected_paths.size(); ++protected_index) {
                if (Rcpp::CharacterVector::is_na(
                        protected_paths[protected_index])) {
                    throw std::invalid_argument(
                        "protected paths cannot contain NA");
                }
                const fs::path protected_path = fs::absolute(
                    Rcpp::as<std::string>(protected_paths[protected_index]))
                                                    .lexically_normal();
                if (candidate == protected_path) return true;

                std::error_code candidate_error;
                std::error_code protected_error;
                const bool candidate_exists =
                    fs::exists(candidate, candidate_error);
                const bool protected_exists =
                    fs::exists(protected_path, protected_error);
                if (candidate_error || protected_error) {
                    throw std::runtime_error(
                        "could not inspect batch scan paths");
                }
                if (!candidate_exists || !protected_exists) continue;

                std::error_code equivalent_error;
                const bool equivalent = fs::equivalent(
                    candidate, protected_path, equivalent_error);
                if (equivalent_error) {
                    throw std::runtime_error(
                        "could not verify batch scan file identity");
                }
                if (equivalent) return true;
            }
        }
        return false;
    } catch (const std::exception& error) {
        Rcpp::stop("embedded RAD batch path preflight failed: %s",
                   error.what());
    }
    return false;
}

// [[Rcpp::export(name = ".rad_scan_wl_cpp")]]
Rcpp::List rad_scan_wl_cpp(Rcpp::List options) {
    const std::string input = string_option(options, "input");
    const std::string output_prefix = string_option(options, "output_prefix");
    if (input.empty()) Rcpp::stop("`input` is required");
    if (output_prefix.empty()) Rcpp::stop("`output_prefix` is required");

    try {
        NativeScanConfig config;
        config.input = input;
        config.output_prefix = output_prefix;
        config.adapter_seq = string_option(options, "adapter_seq");
        config.whitelist = string_option(options, "whitelist");
        config.bc1_whitelist = string_option(options, "bc1_whitelist");
        config.bc2_whitelist = string_option(options, "bc2_whitelist");
        config.barcode_length = int_option(options, "barcode_length", 16);
        config.left_margin = int_option(options, "left_margin", 0);
        config.right_margin = int_option(options, "right_margin", 0);
        config.max_reads = int_option(options, "max_reads", 0);
        config.max_error = double_option(
            options, "max_error", core::kDefaultScanWlMaxErrorRatio);
        config.requested_threads = int_option(options, "threads", 1);
        config.chunk_size = int_option(options, "chunk_size", 10000);
        config.rescan = bool_option(options, "rescan", false);
        config.bc1_length = int_option(options, "bc1_length", 0);
        config.bc2_length = int_option(options, "bc2_length", 0);
        config.umi_length = int_option(options, "umi_length", 9);
        config.offset_min = int_option(options, "offset_min", 0);
        config.offset_max = int_option(options, "offset_max", 3);
        config.verbose = bool_option(options, "verbose", false);

        const std::string selection =
            string_option(options, "selection", "high_specificity");
        if (selection == "high_specificity") {
            config.selection =
                core::scan_whitelist_selection::high_specificity;
        } else if (selection == "above_floor") {
            config.selection = core::scan_whitelist_selection::above_floor;
        } else {
            throw std::invalid_argument(
                "`selection` must be 'high_specificity' or 'above_floor'");
        }

        if (config.barcode_length < 1 || config.barcode_length > 32) {
            throw std::invalid_argument(
                "`barcode_length` must be between 1 and 32");
        }
        if (config.left_margin < 0 || config.right_margin < 0) {
            throw std::invalid_argument(
                "scan margins must be non-negative");
        }
        if (config.max_reads < 0) {
            throw std::invalid_argument("`max_reads` must be positive or zero");
        }
        if (!std::isfinite(config.max_error) || config.max_error < 0.0 ||
            config.max_error > 1.0) {
            throw std::invalid_argument(
                "`max_error` must be between 0 and 1");
        }
        if (config.requested_threads < 1) {
            throw std::invalid_argument("`threads` must be >= 1");
        }
        if (config.chunk_size < 1) {
            throw std::invalid_argument("`chunk_size` must be >= 1");
        }
        if (config.bc1_length < 0 || config.bc1_length > 32 ||
            config.bc2_length < 0 || config.bc2_length > 32) {
            throw std::invalid_argument(
                "split-barcode packed lengths must be between 1 and 32");
        }
        if (config.umi_length < 0 || config.offset_min < 0 ||
            config.offset_max < config.offset_min) {
            throw std::invalid_argument(
                "split-barcode UMI/offset values are invalid");
        }
        const bool partial_split = config.bc1_whitelist.empty() !=
                                   config.bc2_whitelist.empty();
        if (partial_split) {
            throw std::invalid_argument(
                "`bc1_whitelist` and `bc2_whitelist` must be supplied together");
        }
        if (config.rescan && config.split_mode()) {
            throw std::invalid_argument(
                "rescan is not compatible with split-barcode options");
        }
        if (config.rescan && !config.whitelist.empty()) {
            throw std::invalid_argument(
                "rescan is not compatible with a reference whitelist");
        }
        if (!config.rescan && config.adapter_seq.empty()) {
            throw std::invalid_argument(
                "`adapter_seq` is required for a FASTQ/FASTA scan");
        }
        if (config.split_mode() && !config.whitelist.empty()) {
            throw std::invalid_argument(
                "`whitelist` is not used in split-barcode mode");
        }

#ifdef _OPENMP
        config.effective_threads = config.requested_threads;
#else
        config.effective_threads = 1;
#endif
        OpenMpThreadGuard openmp_thread_guard(config.effective_threads);
        const std::string resources = normalize_resource_dir(
            string_option(options, "resource_dir"));
        EnvironmentGuard resource_guard("RAD_RESOURCES", resources);
        EnvironmentGuard embedded_only_guard("RAD_EMBEDDED_ONLY", "1");
        EnvironmentGuard pigz_guard("RAD_NO_PIGZ", "1");
        StreamCapture capture;
        const std::function<void()> interrupt_poll = []() {
            Rcpp::checkUserInterrupt();
        };
        const NativeScanRun run =
            run_scan_whitelist(config, interrupt_poll);

        Rcpp::List value = scan_run_as_list(run, config);
        value["success"] = true;
        value["backend"] = "embedded";
        value["core"] = core_info_impl();
        value["log"] = public_scan_log(capture.str(), run);
        return value;
    } catch (const std::exception& error) {
        Rcpp::stop("embedded RAD whitelist scan failed: %s", error.what());
    }
    return Rcpp::List::create();
}

// [[Rcpp::export(name = ".rad_reformat_cpp")]]
SEXP rad_reformat_cpp(Rcpp::List options) {
    const std::string input = string_option(options, "input_path");
    if (input.empty()) Rcpp::stop("`input_path` is required");

    try {
        core::reformat_options config;
        config.input_path = input;
        config.split_output_dir =
            string_option(options, "split_output_dir");
        config.split_by_barcode =
            bool_option(options, "split_by_barcode", false);
        config.reformat_header =
            bool_option(options, "reformat_header", false);
        config.parse_collapsed_id =
            bool_option(options, "parse_collapsed_id", false);
        const std::string delimiter =
            string_option(options, "delimiter", "_");
        if (delimiter.size() != 1) {
            throw std::invalid_argument(
                "`delimiter` must be exactly one single-byte character");
        }
        config.delimiter = delimiter.front();
        if (has_value(options, "coordinate_mode")) {
            config.coordinate_mode =
                string_option(options, "coordinate_mode");
        }
        config.bin_size_um = int_option(options, "bin_size_um", 2);
        config.output_path = string_option(options, "output_path");
        const int requested_threads = int_option(options, "threads", 1);
        if (requested_threads < 1) {
            throw std::invalid_argument("`threads` must be >= 1");
        }
#ifdef _OPENMP
        const int effective_threads = requested_threads;
#else
        const int effective_threads = 1;
#endif
        config.threads = effective_threads;
        const std::int64_t chunk_size =
            int64_option(options, "chunk_size", 5000);
        if (chunk_size < 1 ||
            static_cast<std::uint64_t>(chunk_size) >
                static_cast<std::uint64_t>(
                    std::numeric_limits<std::size_t>::max())) {
            throw std::invalid_argument(
                "`chunk_size` must be a positive whole number supported by this build");
        }
        config.chunk_size = static_cast<std::size_t>(chunk_size);
        config.verbose = bool_option(options, "verbose", false);
        config.overwrite = false;

        OpenMpThreadGuard openmp_thread_guard(effective_threads);
        const std::string resources = normalize_resource_dir(
            string_option(options, "resource_dir"));
        EnvironmentGuard resource_guard("RAD_RESOURCES", resources);
        EnvironmentGuard embedded_only_guard("RAD_EMBEDDED_ONLY", "1");
        EnvironmentGuard pigz_guard("RAD_NO_PIGZ", "1");
        StreamCapture capture;
        const std::function<void()> interrupt_poll = []() {
            Rcpp::checkUserInterrupt();
        };
        SEXP prepared_result = R_NilValue;
        const std::function<void(const core::reformat_result&)>
            prepare_result = [&](const core::reformat_result& run) {
                const std::string output = config.split_by_barcode
                                               ? run.split_output_dir
                                               : run.output_path;
                if (output.empty()) {
                    throw std::runtime_error(
                        "reformat core returned no output destination");
                }
                Rcpp::List result = Rcpp::List::create(
                    Rcpp::_["success"] = true,
                    Rcpp::_["backend"] = "embedded",
                    Rcpp::_["core"] = core_info_impl(),
                    Rcpp::_["input"] = run.input_path,
                    Rcpp::_["output"] = output,
                    Rcpp::_["mode"] = config.split_by_barcode
                                          ? "split"
                                          : "aggregate",
                    Rcpp::_["files"] = Rcpp::List::create(
                        Rcpp::_["fastq"] = config.split_by_barcode
                                                  ? std::string()
                                                  : run.output_path,
                        Rcpp::_["split_fastq"] = reformat_split_paths(run)),
                    Rcpp::_["stats"] = reformat_stats_as_list(
                        run, requested_threads, effective_threads),
                    Rcpp::_["log"] = capture.str());
                if (fail_reformat_result_construction()) {
                    throw std::bad_alloc();
                }
                // `result` protects the complete SEXP until this callback
                // returns. Publication performs no R allocation, and the
                // direct .Call wrapper returns this SEXP without wrapping it.
                prepared_result = result;
            };
        core::reformat_fastx(config, interrupt_poll, prepare_result);
        return prepared_result;
    } catch (const std::exception& error) {
        Rcpp::stop("embedded RAD reformat failed: %s", error.what());
    }
    return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List rad_demux_cpp(Rcpp::List options) {
    const std::string layout_requested = string_option(options, "layout");
    const std::string fastq_requested = string_option(options, "fastq");
    if (layout_requested.empty()) {
        Rcpp::stop("`layout` is required");
    }
    if (fastq_requested.empty()) {
        Rcpp::stop("`fastq` is required");
    }

    try {
        const fs::path input_path =
            fs::absolute(fastq_requested).lexically_normal();
        if (!fs::is_regular_file(input_path)) {
            throw std::invalid_argument("FASTQ file does not exist: " +
                                        input_path.string());
        }

        const int requested_threads = int_option(options, "threads", 1);
        const std::int64_t chunk_value =
            int64_option(options, "chunk_size", 5000);
        const std::int64_t max_reads_value =
            int64_option(options, "max_reads", -1);
        if (requested_threads < 1) {
            throw std::invalid_argument("`threads` must be >= 1");
        }
#ifdef _OPENMP
        const int threads = requested_threads;
#else
        const int threads = 1;
#endif
        OpenMpThreadGuard openmp_thread_guard(threads);
        if (chunk_value < 1) {
            throw std::invalid_argument("`chunk_size` must be >= 1");
        }

        const std::string correction_mode =
            string_option(options, "bc_correction_mode", "offensive");
        if (correction_mode != "offensive" && correction_mode != "defensive") {
            throw std::invalid_argument(
                "`bc_correction_mode` must be 'offensive' or 'defensive'");
        }
        const std::string joint_mode =
            string_option(options, "joint_bc_mode", "default");
        if (joint_mode != "default" && joint_mode != "strict") {
            throw std::invalid_argument(
                "`joint_bc_mode` must be 'default' or 'strict'");
        }

        const int min_read_length =
            int_option(options, "min_read_length", -1);
        if (min_read_length < -1) {
            throw std::invalid_argument("`min_read_length` must be >= 0 or NULL");
        }

        const bool verbose = bool_option(options, "verbose", false);
        const bool max_verbose = bool_option(options, "max_verbose", false);
        const bool write_debug = bool_option(options, "write_debug", false);
        const bool reuse_cache = bool_option(options, "reuse_cache", false);
        const bool auto_whitelist =
            bool_option(options, "auto_whitelist", false);
        const bool rc_umi = bool_option(options, "rc_umi", true);
        const std::string kit = string_option(options, "kit");
        const std::string global_whitelist =
            string_option(options, "global_whitelist");
        const std::string custom_whitelist =
            string_option(options, "custom_whitelist");
        if (!kit.empty() && !global_whitelist.empty()) {
            throw std::invalid_argument(
                "`kit` and `global_whitelist` are alternatives; supply only one");
        }
        if (auto_whitelist && !custom_whitelist.empty()) {
            throw std::invalid_argument(
                "`auto_whitelist` and `custom_whitelist` are alternatives; supply only one");
        }
        const std::string reference = kit.empty() ? global_whitelist : kit;
        const std::optional<int> whitelist_mutation =
            optional_int(options, "whitelist_mutation");
        const std::optional<int> generated_mutation =
            optional_int(options, "generated_mutation");

        const std::string scan_adapter =
            string_option(options, "scan_adapter");
        const std::optional<int> scan_barcode_length =
            optional_int(options, "scan_barcode_length");
        if (scan_barcode_length &&
            (*scan_barcode_length < 1 || *scan_barcode_length > 32)) {
            throw std::invalid_argument(
                "`scan_barcode_length` must be between 1 and 32");
        }
        const double scan_max_error = double_option(
            options, "scan_max_error", core::kDefaultScanWlMaxErrorRatio);
        if (!std::isfinite(scan_max_error) || scan_max_error < 0.0 ||
            scan_max_error > 1.0) {
            throw std::invalid_argument(
                "`scan_max_error` must be between 0 and 1");
        }
        const std::int64_t scan_max_reads_value = has_value(
            options, "scan_max_reads")
            ? int64_option(options, "scan_max_reads", -1)
            : max_reads_value;
        const std::int64_t scan_chunk_value = has_value(
            options, "scan_chunk_size")
            ? int64_option(options, "scan_chunk_size", chunk_value)
            : chunk_value;
        const int scan_requested_threads = has_value(options, "scan_threads")
                                               ? int_option(
                                                     options, "scan_threads",
                                                     requested_threads)
                                               : requested_threads;
        if (auto_whitelist &&
            (scan_max_reads_value < -1 || scan_max_reads_value == 0 ||
             scan_max_reads_value > std::numeric_limits<int>::max())) {
            throw std::invalid_argument(
                "`scan_max_reads` must be -1 or between 1 and INT_MAX");
        }
        if (auto_whitelist &&
            (scan_chunk_value < 1 ||
             scan_chunk_value > std::numeric_limits<int>::max())) {
            throw std::invalid_argument(
                "`scan_chunk_size` must be between 1 and INT_MAX");
        }
        if (scan_requested_threads < 1) {
            throw std::invalid_argument("`scan_threads` must be >= 1");
        }
#ifdef _OPENMP
        const int scan_effective_threads = scan_requested_threads;
#else
        const int scan_effective_threads = 1;
#endif
        const std::string scan_selection_name =
            string_option(options, "scan_selection", "high_specificity");
        core::scan_whitelist_selection scan_selection;
        if (scan_selection_name == "high_specificity") {
            scan_selection =
                core::scan_whitelist_selection::high_specificity;
        } else if (scan_selection_name == "above_floor") {
            scan_selection = core::scan_whitelist_selection::above_floor;
        } else {
            throw std::invalid_argument(
                "`scan_selection` must be 'high_specificity' or 'above_floor'");
        }

        const std::string resources = normalize_resource_dir(
            string_option(options, "resource_dir"));
        EnvironmentGuard resource_guard("RAD_RESOURCES", resources);
        // Ignore ~/.rad layout/whitelist overrides while the embedded core is
        // running.  Those files belong to the standalone CLI and otherwise
        // make an rrad run depend on unreported, per-user state.
        EnvironmentGuard embedded_only_guard("RAD_EMBEDDED_ONLY", "1");
        // The embedded backend is deliberately self-contained.  This also
        // avoids shelling out from an R process and guarantees zlib output.
        EnvironmentGuard pigz_guard("RAD_NO_PIGZ", "1");

        const std::string layout_csv =
            core::config_utils::check_if_custom_rl(layout_requested)
                ? layout_requested
                : core::config_utils::get_read_layout(layout_requested);

        fs::path output_dir = fs::absolute(
            string_option(options, "output_dir", "."));
        fs::path output_name =
            string_option(options, "output_prefix", "output");
        fs::path outbase = output_name.is_absolute()
                               ? output_name
                               : output_dir / output_name;
        outbase = outbase.lexically_normal();

        const std::string fastq_suffix =
            core::path_utils::get_fastqa_type(input_path.string());
        const fs::path fastq_output =
            outbase.string() + fastq_suffix + ".gz";
        const fs::path layout_cache = outbase.string() + "_layout.csv";
        const fs::path position_cache =
            outbase.string() + "_position_map.csv";
        const fs::path true_whitelist =
            outbase.string() + "_whitelist_true.csv";
        const fs::path global_whitelist_output =
            outbase.string() + "_whitelist_global.csv";
        const fs::path demux_log = outbase.string() + "_demux.log";
        const fs::path scan_csv = outbase.string() + "_scanwl.csv";
        const fs::path scan_txt = outbase.string() + "_scanwl.txt";
        const fs::path scan_log = outbase.string() + "_scan_wl.log";
        const fs::path debug_sig = outbase.string() + "_dbg.sig.gz";
        const fs::path debug_csv = outbase.string() + "_dbg.csv.gz";
        const fs::path debug_fastq =
            outbase.string() + "_dbg" + fastq_suffix + ".gz";
        const fs::path debug_metrics = outbase.string() + ".metrics.tsv";

        std::vector<NamedPath> input_paths{
            {"FASTQ input", input_path},
            {"layout input", fs::path(layout_csv)}};
        if (!reference.empty()) {
            append_whitelist_inputs(input_paths, reference,
                                    "global whitelist input");
        }
        if (!custom_whitelist.empty()) {
            append_whitelist_inputs(input_paths, custom_whitelist,
                                    "custom whitelist input");
        }

        std::vector<NamedPath> output_paths{
            {"primary FASTQ", fastq_output},
            {"layout cache", layout_cache},
            {"position-map cache", position_cache},
            {"true-whitelist summary", true_whitelist},
            {"global-whitelist summary", global_whitelist_output},
            {"demultiplexing log", demux_log},
            {"auto-whitelist table", scan_csv},
            {"auto-whitelist list", scan_txt},
            {"auto-whitelist log", scan_log},
            {"debug signature", debug_sig},
            {"debug table", debug_csv},
            {"debug FASTQ", debug_fastq},
            {"debug metrics", debug_metrics}};
        reject_input_output_collisions(input_paths, output_paths);
        preflight_transaction_targets(output_paths);

        SameFilesystemOutputTransaction publication(outbase);
        const std::vector<NamedPath> publication_control_paths{
            {"demux publication lock", publication.lock_path()}};
        reject_input_output_collisions(input_paths,
                                       publication_control_paths);
        const fs::path staged_outbase = publication.stage_base();
        const fs::path staged_fastq =
            staged_outbase.string() + fastq_suffix + ".gz";
        const fs::path staged_layout =
            staged_outbase.string() + "_layout.csv";
        const fs::path staged_position =
            staged_outbase.string() + "_position_map.csv";
        const fs::path staged_true_whitelist =
            staged_outbase.string() + "_whitelist_true.csv";
        const fs::path staged_global_whitelist =
            staged_outbase.string() + "_whitelist_global.csv";
        const fs::path staged_demux_log =
            staged_outbase.string() + "_demux.log";
        const fs::path staged_scan_csv =
            staged_outbase.string() + "_scanwl.csv";
        const fs::path staged_scan_txt =
            staged_outbase.string() + "_scanwl.txt";
        const fs::path staged_scan_log =
            staged_outbase.string() + "_scan_wl.log";
        const fs::path staged_debug_sig =
            staged_outbase.string() + "_dbg.sig.gz";
        const fs::path staged_debug_csv =
            staged_outbase.string() + "_dbg.csv.gz";
        const fs::path staged_debug_fastq =
            staged_outbase.string() + "_dbg" + fastq_suffix + ".gz";
        const fs::path staged_debug_metrics =
            staged_outbase.string() + ".metrics.tsv";

        StreamCapture capture;
        const auto run_started = std::chrono::steady_clock::now();
        const std::function<void()> interrupt_poll = []() {
            // This callback is invoked only by chunk_streaming's serial
            // boundary path, never by an OpenMP or writer worker.
            Rcpp::checkUserInterrupt();
        };

        core::ReadLayout read_layout;
        // A calibration cache is one logical artifact. Reusing either half
        // independently can combine a newly requested layout with an orphaned
        // position map (or vice versa), silently changing read interpretation.
        // Writers protect publication with this same per-prefix lock. Snapshot
        // both files into this run's private staging directory while holding
        // it, then parse and republish only those immutable copies.
        bool have_cache_pair = false;
        if (reuse_cache) {
            DemuxCommitLock cache_snapshot_lock(publication.lock_path());
            preflight_transaction_targets({
                {"layout cache", layout_cache},
                {"position-map cache", position_cache}});

            auto regular_file_if_present = [](const fs::path& path,
                                              const char* description) {
                std::error_code error;
                const bool exists = fs::exists(path, error);
                if (error) {
                    throw std::runtime_error(
                        std::string("could not inspect reusable ") +
                        description + ": " + error.message());
                }
                if (!exists) return false;
                const bool regular = fs::is_regular_file(path, error);
                if (error) {
                    throw std::runtime_error(
                        std::string("could not inspect reusable ") +
                        description + ": " + error.message());
                }
                return regular;
            };
            const bool layout_available = regular_file_if_present(
                layout_cache, "layout cache");
            const bool position_available = regular_file_if_present(
                position_cache, "position-map cache");

            have_cache_pair = layout_available && position_available;
            if (have_cache_pair) {
                std::error_code copy_error;
                fs::copy_file(layout_cache, staged_layout,
                              fs::copy_options::overwrite_existing,
                              copy_error);
                if (copy_error) {
                    throw std::runtime_error(
                        "could not snapshot reused layout cache: " +
                        copy_error.message());
                }
                copy_error.clear();
                fs::copy_file(position_cache, staged_position,
                              fs::copy_options::overwrite_existing,
                              copy_error);
                if (copy_error) {
                    throw std::runtime_error(
                        "could not snapshot reused position-map cache: " +
                        copy_error.message());
                }
                require_nonempty_regular_file(
                    staged_layout, "snapshotted layout cache");
                require_nonempty_regular_file(
                    staged_position, "snapshotted position-map cache");
            }
        }
        const bool have_layout = have_cache_pair;
        const bool have_position = have_cache_pair;
        if (have_cache_pair) {
            auto reject_cache_alias = [&](const std::string& input_label,
                                          const fs::path& input_path,
                                          const std::string& own_output_label) {
                std::vector<NamedPath> other_outputs;
                other_outputs.reserve(output_paths.size());
                for (const auto& output_path : output_paths) {
                    if (output_path.first != own_output_label) {
                        other_outputs.push_back(output_path);
                    }
                }
                reject_input_output_collisions(
                    {{input_label, input_path}}, other_outputs);
            };
            reject_cache_alias("reused layout cache input", layout_cache,
                               "layout cache");
            reject_cache_alias("reused position-map cache input",
                               position_cache, "position-map cache");
        }
        bool imported_companion = false;
        std::string imported_companion_path;
        TemporaryFileGuard seqspec_temp;

        if (have_layout) {
            read_layout.import_read_layout(staged_layout.string(), max_verbose,
                                           layout_csv);
        } else {
            // Match RAD's demux cache-miss path: translate seqspec YAML into a
            // temporary native layout, then feed the unchanged layout engine.
            std::string layout_source = layout_csv;
            if (core::seqspec::is_seqspec_file(layout_source)) {
                if (verbose || max_verbose) {
                    RAD_COUT <<
                        "[demux] Detected seqspec; translating to RAD layout...\n";
                }
                layout_source = core::seqspec::convert_to_temp_layout_csv(
                    layout_source, verbose || max_verbose,
                    auto_whitelist || !reference.empty() ||
                        !custom_whitelist.empty());
                seqspec_temp.reset(layout_source);
            }

            if (core::ReadLayout::is_generated_cache(layout_source)) {
                read_layout.import_read_layout(layout_source, max_verbose,
                                               layout_source);
                const std::string companion =
                    companion_position_map(layout_source);
                if (!companion.empty() && fs::is_regular_file(companion)) {
                    read_layout.import_position_map(companion, max_verbose);
                    imported_companion = true;
                    imported_companion_path = companion;
                }
            } else {
                read_layout.prep_new_layout(layout_source,
                                            verbose || max_verbose);
            }
        }

        std::vector<NamedPath> deferred_inputs;
        for (const auto& element : read_layout.by_order()) {
            const bool layout_whitelist_is_consumed =
                (reference.empty() && custom_whitelist.empty()) ||
                element.flags.find("joint_barcode") != std::string::npos;
            if (layout_whitelist_is_consumed &&
                element.global_class == "barcode" &&
                element.type == "variable" &&
                !element.whitelist_path.empty()) {
                append_whitelist_inputs(
                    deferred_inputs, element.whitelist_path,
                    "layout whitelist input for " + element.class_id);
            }
        }
        if (!imported_companion_path.empty()) {
            deferred_inputs.emplace_back("companion position-map input",
                                         imported_companion_path);
        }
        if (!deferred_inputs.empty()) {
            reject_input_output_collisions(deferred_inputs, output_paths);
            reject_input_output_collisions(deferred_inputs,
                                           publication_control_paths);
        }

        if (have_position) {
            read_layout.import_position_map(staged_position.string(),
                                            max_verbose);
        } else if (!imported_companion) {
            read_layout.generate_position_mapping();
        }

        if ((!have_layout || !have_position) && !imported_companion) {
            read_layout.write_to_csv(staged_outbase.string(), "layout");
            core::Misalignment_Setup misalignment(read_layout);
            misalignment.generate_misalignment_data(
                input_path.string(), read_layout, threads, 50000,
                interrupt_poll);
            read_layout.write_to_csv(staged_outbase.string(), "both");
        } else if (imported_companion) {
            read_layout.write_to_csv(staged_outbase.string(), "both");
        }
        require_nonempty_regular_file(staged_layout,
                                      "generated layout cache");
        require_nonempty_regular_file(staged_position,
                                      "generated position-map cache");

        std::optional<AutoWhitelistResult> auto_result;
        if (auto_whitelist) {
            auto_result = run_auto_whitelist(
                read_layout, input_path.string(), staged_outbase.string(),
                reference,
                scan_adapter, scan_barcode_length, scan_max_error,
                scan_max_reads_value < 0
                    ? static_cast<size_t>(-1)
                    : static_cast<size_t>(scan_max_reads_value),
                static_cast<size_t>(scan_chunk_value),
                scan_requested_threads, scan_effective_threads,
                scan_selection, verbose, interrupt_poll);
            read_layout.load_wl(std::nullopt, whitelist_mutation, verbose,
                                threads,
                                auto_result->reference_spec.empty()
                                    ? std::optional<std::string>{}
                                    : std::optional<std::string>{
                                          auto_result->reference_spec},
                                auto_result->detected_path);
        } else if (!custom_whitelist.empty()) {
            read_layout.load_wl(
                std::nullopt, whitelist_mutation, verbose, threads,
                reference.empty() ? std::optional<std::string>{}
                                  : std::optional<std::string>{reference},
                custom_whitelist);
        } else if (!reference.empty()) {
            read_layout.load_wl(std::nullopt, whitelist_mutation, verbose,
                                threads, reference, std::nullopt);
        } else {
            read_layout.load_wl(std::nullopt, whitelist_mutation, verbose,
                                threads);
        }

        std::size_t loaded_true_barcodes = 0;
        std::size_t loaded_global_barcodes = 0;
        for (const auto& item : read_layout.wl_map.lists) {
            loaded_true_barcodes += item.second.true_bcs.unique_val_size();
            loaded_global_barcodes += item.second.global_bcs.size();
        }

        const auto demux_started = std::chrono::steady_clock::now();
        core::sigalign_run_stats stats = core::SigString::sigalign(
            input_path.string(), read_layout, staged_outbase.string(),
            generated_mutation, max_verbose, threads,
            static_cast<size_t>(chunk_value),
            max_reads_value < 0 ? static_cast<size_t>(-1)
                                : static_cast<size_t>(max_reads_value),
            write_debug, correction_mode, joint_mode, rc_umi,
            min_read_length, interrupt_poll);
        stats.wall_time_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          demux_started)
                .count();
        stats.overhead_time_seconds = std::max(
            0.0, stats.wall_time_seconds - stats.process_time_seconds -
                     stats.output_staging_time_seconds);

        verify_primary_gzip(staged_fastq, stats.records_serialized > 0);

        const std::string whitelist_base =
            staged_outbase.string() + "_whitelist.csv";
        read_layout.save_wl(whitelist_base, verbose, true);

        const double total_wall_time =
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          run_started)
                .count();
        const double pass_rate =
            stats.total_reads == 0
                ? 0.0
                : 100.0 * static_cast<double>(stats.reads_demultiplexed) /
                      static_cast<double>(stats.total_reads);

        {
            std::ofstream log(staged_demux_log);
            if (!log) {
                throw std::runtime_error("could not write demux summary: " +
                                         staged_demux_log.string());
            }
            log << "RAD embedded demultiplexing summary\n"
                << "status=success\n"
                << "rad_version=" << kCoreVersion << "\n"
                << "rad_source_commit=" << kSourceCommit << "\n"
                << "rad_core_api=" << kCoreApiVersion << "\n"
                << "layout=" << layout_requested << "\n"
                << "input=" << input_path.string() << "\n"
                << "output_prefix=" << outbase.string() << "\n"
                << "threads_requested=" << requested_threads << "\n"
                << "threads_effective=" << threads << "\n"
                << "reuse_cache=" << (reuse_cache ? "true" : "false") << "\n"
                << "auto_whitelist=" << (auto_whitelist ? "true" : "false")
                << "\n";
            if (auto_result) {
                log << "scan_adapter=" << auto_result->config.adapter_seq
                    << "\n"
                    << "scan_barcode_length="
                    << auto_result->config.barcode_length << "\n"
                    << "scan_reference="
                    << auto_result->config.whitelist << "\n"
                    << "scan_max_error=" << auto_result->config.max_error
                    << "\n"
                    << "scan_max_reads=" << auto_result->config.max_reads
                    << "\n"
                    << "scan_chunk_size=" << auto_result->config.chunk_size
                    << "\n"
                    << "scan_threads_requested="
                    << auto_result->config.requested_threads << "\n"
                    << "scan_threads_effective="
                    << auto_result->config.effective_threads << "\n"
                    << "scan_selection="
                    << core::scan_whitelist_selection_name(
                           auto_result->config.selection)
                    << "\n"
                    << "scan_reads_processed="
                    << auto_result->scan.stats.reads_processed << "\n"
                    << "scan_total_extractions="
                    << auto_result->scan.stats.total_extractions << "\n"
                    << "scan_selected_barcodes="
                    << auto_result->scan.stats.selected_barcodes << "\n";
            }
            log << "total_reads=" << stats.total_reads << "\n"
                << "reads_passing_filter=" << stats.reads_passing_filter
                << "\n"
                << "reads_demultiplexed=" << stats.reads_demultiplexed
                << "\n"
                << "records_serialized=" << stats.records_serialized << "\n"
                << "demultiplex_rate_percent=" << pass_rate << "\n"
                << "chunks_processed=" << stats.chunks_processed << "\n"
                << "sigalign_wall_time_seconds=" << stats.wall_time_seconds
                << "\n"
                << "process_time_seconds=" << stats.process_time_seconds
                << "\n"
                << "output_staging_time_seconds="
                << stats.output_staging_time_seconds << "\n"
                << "overhead_time_seconds=" << stats.overhead_time_seconds
                << "\n"
                << "total_wall_time_seconds=" << total_wall_time << "\n";
            log.flush();
            if (!log) {
                throw std::runtime_error("could not write demux summary: " +
                                         staged_demux_log.string());
            }
            log.close();
            if (log.fail()) {
                throw std::runtime_error("could not close demux summary: " +
                                         staged_demux_log.string());
            }
        }

        auto staged_exists = [](const fs::path& path) {
            std::error_code error;
            const bool exists = fs::exists(path, error);
            if (error) {
                throw std::runtime_error("could not inspect staged RAD output: " +
                                         path.string());
            }
            if (!exists) return false;
            const bool regular = fs::is_regular_file(path, error);
            if (error) {
                throw std::runtime_error("could not inspect staged RAD output: " +
                                         path.string());
            }
            return regular;
        };
        auto published_path = [&](const fs::path& staged,
                                  const fs::path& target) {
            return staged_exists(staged) ? fs::absolute(target).string()
                                         : std::string();
        };

        const bool publish_true_whitelist =
            staged_exists(staged_true_whitelist);
        const bool publish_global_whitelist =
            staged_exists(staged_global_whitelist);

        Rcpp::List files = Rcpp::List::create(
            Rcpp::_["fastq"] = published_path(staged_fastq, fastq_output),
            Rcpp::_["layout"] = published_path(staged_layout, layout_cache),
            Rcpp::_["position_map"] =
                published_path(staged_position, position_cache),
            Rcpp::_["whitelist_true"] = published_path(
                staged_true_whitelist, true_whitelist),
            Rcpp::_["whitelist_global"] = published_path(
                staged_global_whitelist, global_whitelist_output),
            Rcpp::_["demux_log"] =
                published_path(staged_demux_log, demux_log),
            Rcpp::_["scan_whitelist_csv"] =
                auto_result ? published_path(staged_scan_csv, scan_csv) : "",
            Rcpp::_["scan_whitelist_txt"] =
                auto_result ? published_path(staged_scan_txt, scan_txt) : "",
            Rcpp::_["scan_log"] =
                auto_result ? published_path(staged_scan_log, scan_log) : "",
            Rcpp::_["debug_sig"] =
                write_debug ? published_path(staged_debug_sig, debug_sig) : "",
            Rcpp::_["debug_csv"] =
                write_debug ? published_path(staged_debug_csv, debug_csv) : "",
            Rcpp::_["debug_fastq"] =
                write_debug
                    ? published_path(staged_debug_fastq, debug_fastq)
                    : "",
            Rcpp::_["metrics"] =
                write_debug
                    ? published_path(staged_debug_metrics, debug_metrics)
                    : "");

        Rcpp::List stats_output = Rcpp::List::create(
            Rcpp::_["threads_requested"] = requested_threads,
            Rcpp::_["threads_effective"] = threads,
            Rcpp::_["total_reads"] = static_cast<double>(stats.total_reads),
            Rcpp::_["reads_passing_filter"] =
                static_cast<double>(stats.reads_passing_filter),
            Rcpp::_["reads_demultiplexed"] =
                static_cast<double>(stats.reads_demultiplexed),
            Rcpp::_["records_serialized"] =
                static_cast<double>(stats.records_serialized),
            Rcpp::_["chunks_processed"] =
                static_cast<double>(stats.chunks_processed),
            Rcpp::_["loaded_true_barcodes"] =
                static_cast<double>(loaded_true_barcodes),
            Rcpp::_["loaded_global_barcodes"] =
                static_cast<double>(loaded_global_barcodes),
            Rcpp::_["demultiplex_rate_percent"] = pass_rate,
            Rcpp::_["wall_time_seconds"] = stats.wall_time_seconds,
            Rcpp::_["process_time_seconds"] = stats.process_time_seconds,
            Rcpp::_["output_staging_time_seconds"] =
                stats.output_staging_time_seconds,
            Rcpp::_["overhead_time_seconds"] = stats.overhead_time_seconds,
            Rcpp::_["total_wall_time_seconds"] = total_wall_time);

        Rcpp::RObject scan_output = R_NilValue;
        if (auto_result) {
            auto_result->scan.output_prefix = outbase.string() + "_scanwl";
            auto_result->scan.csv = scan_csv.string();
            auto_result->scan.whitelist = scan_txt.string();
            auto_result->scan.log = scan_log.string();
            auto_result->config.output_prefix = outbase.string() + "_scanwl";
            auto_result->config.log_path = scan_log.string();
            auto_result->captured_log = replace_all_copy(
                auto_result->captured_log, staged_outbase.string(),
                outbase.string());
            Rcpp::List nested_scan = scan_run_as_list(
                auto_result->scan, auto_result->config);
            nested_scan["success"] = true;
            nested_scan["backend"] = "embedded";
            nested_scan["core"] = core_info_impl();
            nested_scan["log"] = auto_result->captured_log;
            scan_output = nested_scan;
        }

        // Build the complete R result before publication so allocation or
        // interrupt failures cannot occur after the old artifact set has been
        // replaced.
        Rcpp::checkUserInterrupt();
        const std::string public_log = replace_all_copy(
            capture.str(), staged_outbase.string(), outbase.string());
        Rcpp::List result = Rcpp::List::create(
            Rcpp::_["success"] = true,
            Rcpp::_["backend"] = "embedded",
            Rcpp::_["core"] = core_info_impl(),
            Rcpp::_["input"] = input_path.string(),
            Rcpp::_["layout"] = layout_requested,
            Rcpp::_["output_prefix"] = outbase.string(),
            Rcpp::_["files"] = files,
            Rcpp::_["stats"] = stats_output,
            Rcpp::_["scan"] = scan_output,
            Rcpp::_["log"] = public_log);

        const std::vector<TransactionOutput> transaction_outputs{
            {"primary FASTQ", staged_fastq, fastq_output, true},
            {"layout cache", staged_layout, layout_cache, true},
            {"position-map cache", staged_position, position_cache, true},
            {"true-whitelist summary", staged_true_whitelist, true_whitelist,
             publish_true_whitelist},
            {"global-whitelist summary", staged_global_whitelist,
             global_whitelist_output, publish_global_whitelist},
            {"demultiplexing log", staged_demux_log, demux_log, true},
            {"auto-whitelist table", staged_scan_csv, scan_csv,
             auto_whitelist},
            {"auto-whitelist list", staged_scan_txt, scan_txt,
             auto_whitelist},
            {"auto-whitelist log", staged_scan_log, scan_log,
             auto_whitelist},
            {"debug signature", staged_debug_sig, debug_sig, write_debug},
            {"debug table", staged_debug_csv, debug_csv, write_debug},
            {"debug FASTQ", staged_debug_fastq, debug_fastq, write_debug},
            {"debug metrics", staged_debug_metrics, debug_metrics,
             write_debug}};
        publication.commit(transaction_outputs);
        return result;
    } catch (const std::exception& error) {
        Rcpp::stop("embedded RAD demultiplexing failed: %s", error.what());
    }
    return Rcpp::List::create();  // unreachable; placates some compilers
}
