// Included by rad_core_bridge.cpp inside its anonymous namespace.  Keeping the
// bridge in the same translation unit is important because RAD's vendored core
// is header-defined and must only be instantiated once.

struct NativeScanConfig {
    std::string input;
    std::string output_prefix;
    std::string log_path;
    std::string adapter_seq;
    std::string whitelist;
    std::string bc1_whitelist;
    std::string bc2_whitelist;
    int barcode_length = 16;
    int left_margin = 0;
    int right_margin = 0;
    int max_reads = 0;
    double max_error = core::kDefaultScanWlMaxErrorRatio;
    int requested_threads = 1;
    int effective_threads = 1;
    int chunk_size = 10000;
    core::scan_whitelist_selection selection =
        core::scan_whitelist_selection::high_specificity;
    bool rescan = false;
    int bc1_length = 0;
    int bc2_length = 0;
    int umi_length = 9;
    int offset_min = 0;
    int offset_max = 3;
    bool verbose = false;

    bool split_mode() const {
        return !bc1_whitelist.empty() && !bc2_whitelist.empty();
    }
};

struct NativeScanRun {
    std::string input;
    std::string output_prefix;
    std::string csv;
    std::string whitelist;
    std::string log;
    std::string valid_pairs;
    std::string spatial_mask;
    std::string mode;
    core::whitelist_scan_stats stats;
    std::uint64_t input_rows = 0;
    std::size_t bc1_whitelist_size = 0;
    std::size_t bc2_whitelist_size = 0;
    int effective_bc1_length = 0;
    int effective_bc2_length = 0;
    std::string staging_base;
};

std::string public_scan_log(std::string log, const NativeScanRun& run) {
    auto replace_all = [&](const std::string& from, const std::string& to) {
        if (from.empty()) return;
        size_t offset = 0;
        while ((offset = log.find(from, offset)) != std::string::npos) {
            log.replace(offset, from.size(), to);
            offset += to.size();
        }
    };
    replace_all(run.staging_base + ".csv", run.csv);
    replace_all(run.staging_base + ".txt", run.whitelist);
    replace_all(run.staging_base + "_scan_wl.log", run.log);
    replace_all(run.staging_base + "_valid_pairs.csv", run.valid_pairs);
    replace_all(run.staging_base + "_spat_mask.csv", run.spatial_mask);
    return log;
}

class ScanStagingDirectory {
public:
    explicit ScanStagingDirectory(const fs::path& output_prefix) {
        const fs::path parent = output_prefix.parent_path();
        std::string pattern =
            (parent /
             (".rrad-scan-stage-" + output_prefix.filename().string() +
              "-XXXXXX"))
                .string();
        std::vector<char> writable_pattern(pattern.begin(), pattern.end());
        writable_pattern.push_back('\0');

        char* created = ::mkdtemp(writable_pattern.data());
        if (created == nullptr) {
            throw std::runtime_error(
                "could not allocate private whitelist scan staging "
                "directory beside the output: " +
                std::string(std::strerror(errno)));
        }
        path_ = fs::path(created);

        // mkdtemp creates mode 0700. Enforce the exact mode as a defensive
        // postcondition before any RAD writer opens a staged artifact.
        if (::chmod(path_.string().c_str(), S_IRWXU) != 0) {
            const std::string reason = std::strerror(errno);
            std::error_code cleanup_error;
            fs::remove_all(path_, cleanup_error);
            path_.clear();
            throw std::runtime_error(
                "could not secure whitelist scan staging directory: " +
                reason);
        }
    }

    ScanStagingDirectory(const ScanStagingDirectory&) = delete;
    ScanStagingDirectory& operator=(const ScanStagingDirectory&) = delete;

    ~ScanStagingDirectory() {
        if (path_.empty()) return;
        std::error_code error;
        fs::remove_all(path_, error);
    }

    const fs::path& path() const { return path_; }

private:
    fs::path path_;
};

void require_scan_output(const fs::path& path, const char* description) {
    std::error_code error;
    if (!fs::is_regular_file(path, error) || error) {
        throw std::runtime_error(std::string(description) +
                                 " was not created: " + path.string());
    }
}

struct ScanCommitItem {
    fs::path staging;
    fs::path target;
    bool publish = true;
};

// The lock file is intentionally persistent.  Removing an advisory-lock file
// after unlocking creates an inode race in which a waiter can hold the old
// inode while a third process locks a newly-created one at the same path.
class ScanCommitLock {
public:
    explicit ScanCommitLock(const fs::path& path) : path_(path) {
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
                "could not open whitelist scan commit lock " +
                path_.string() + ": " + std::strerror(errno));
        }

        struct stat information {};
        if (::fstat(descriptor_, &information) != 0) {
            const std::string reason = std::strerror(errno);
            close_descriptor();
            throw std::runtime_error(
                "invalid whitelist scan commit lock " + path_.string() +
                ": " + reason);
        }
        if (!S_ISREG(information.st_mode)) {
            close_descriptor();
            throw std::runtime_error(
                "invalid whitelist scan commit lock " + path_.string() +
                ": lock path is not a regular file");
        }
        if (::fchmod(descriptor_, S_IRUSR | S_IWUSR) != 0) {
            const std::string reason = std::strerror(errno);
            close_descriptor();
            throw std::runtime_error(
                "could not secure whitelist scan commit lock " +
                path_.string() + ": " + reason);
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
                "could not lock whitelist scan outputs with " +
                path_.string() + ": " + reason);
        }
        locked_ = true;
    }

    ScanCommitLock(const ScanCommitLock&) = delete;
    ScanCommitLock& operator=(const ScanCommitLock&) = delete;

    ~ScanCommitLock() { close_descriptor(); }

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

void preflight_scan_outputs(const std::vector<fs::path>& targets) {
    for (std::size_t index = 0; index < targets.size(); ++index) {
        const fs::path normalized =
            fs::absolute(targets[index]).lexically_normal();
        const fs::path parent = normalized.parent_path();
        if (!parent.empty()) {
            std::error_code error;
            fs::create_directories(parent, error);
            if (error) {
                throw std::runtime_error(
                    "could not create whitelist scan output directory: " +
                    error.message());
            }
        }

        std::error_code error;
        const bool exists = fs::exists(normalized, error);
        if (error) {
            throw std::runtime_error(
                "could not inspect whitelist scan output: " +
                normalized.string());
        }
        if (exists && fs::is_directory(normalized, error)) {
            throw std::invalid_argument(
                "whitelist scan output identifies a directory: " +
                normalized.string());
        }
        if (error) {
            throw std::runtime_error(
                "could not inspect whitelist scan output: " +
                normalized.string());
        }

        for (std::size_t prior = 0; prior < index; ++prior) {
            const fs::path prior_normalized =
                fs::absolute(targets[prior]).lexically_normal();
            bool same = normalized == prior_normalized;
            if (!same && exists) {
                std::error_code prior_error;
                if (fs::exists(prior_normalized, prior_error) &&
                    !prior_error) {
                    std::error_code equivalent_error;
                    same = fs::equivalent(normalized, prior_normalized,
                                          equivalent_error);
                    if (equivalent_error) {
                        throw std::runtime_error(
                            "could not compare whitelist scan outputs");
                    }
                }
            }
            if (same) {
                throw std::invalid_argument(
                    "whitelist scan outputs resolve to the same path: " +
                    normalized.string());
            }
        }
    }
}

void commit_scan_outputs(const std::vector<ScanCommitItem>& items,
                         const std::string& nonce,
                         const fs::path& lock_path) {
    struct Backup {
        fs::path target;
        fs::path path;
        bool active = false;
    };
    std::vector<Backup> backups;
    backups.reserve(items.size());
    for (std::size_t index = 0; index < items.size(); ++index) {
        backups.push_back({
            items[index].target,
            items[index].target.string() + ".rrad-backup-" + nonce + "-" +
                std::to_string(index),
            false});
    }
    // Keep all recovery paths and install flags allocated before publication.
    std::vector<unsigned char> installed(items.size(), 0);

    ScanCommitLock commit_lock(lock_path);
    std::vector<fs::path> targets;
    targets.reserve(items.size());
    for (const auto& item : items) targets.push_back(item.target);
    // The original preflight happens before the potentially expensive scan.
    // Repeating it under the commit lock closes the writer/writer TOCTOU gap.
    preflight_scan_outputs(targets);

    // Validate every recovery pathname before moving any current output.
    for (const auto& backup : backups) {
        std::error_code error;
        if (fs::exists(backup.path, error) || error) {
            throw std::runtime_error(
                "could not allocate whitelist scan output backup: " +
                backup.path.string());
        }
    }

    try {
        for (std::size_t index = 0; index < items.size(); ++index) {
            std::error_code error;
            if (!fs::exists(items[index].target, error)) {
                if (error) {
                    throw std::runtime_error(
                        "could not inspect existing whitelist scan output: " +
                        items[index].target.string());
                }
                continue;
            }
            if (std::rename(items[index].target.string().c_str(),
                            backups[index].path.string().c_str()) != 0) {
                throw std::runtime_error(
                    "could not stage existing whitelist scan output " +
                    items[index].target.string() + ": " +
                    std::strerror(errno));
            }
            backups[index].active = true;
        }

        for (std::size_t index = 0; index < items.size(); ++index) {
            const auto& item = items[index];
            if (!item.publish) continue;
            if (std::rename(item.staging.string().c_str(),
                            item.target.string().c_str()) != 0) {
                throw std::runtime_error(
                    "could not commit whitelist scan output " +
                    item.target.string() + ": " + std::strerror(errno));
            }
            installed[index] = 1;
        }
    } catch (...) {
        const std::exception_ptr original = std::current_exception();
        std::vector<std::string> rollback_failures;

        for (std::size_t offset = items.size(); offset > 0; --offset) {
            const std::size_t index = offset - 1;
            if (!installed[index]) continue;
            std::error_code error;
            const bool removed = fs::remove(items[index].target, error);
            std::error_code inspect_error;
            const bool still_exists =
                fs::exists(items[index].target, inspect_error);
            if (error || inspect_error || (!removed && still_exists)) {
                rollback_failures.push_back(
                    "could not remove newly installed output " +
                    items[index].target.string() +
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
            std::string original_message = "unknown commit failure";
            try {
                std::rethrow_exception(original);
            } catch (const std::exception& error) {
                original_message = error.what();
            } catch (...) {
            }
            std::ostringstream message;
            message << original_message
                    << "; whitelist scan rollback was incomplete";
            for (const auto& failure : rollback_failures) {
                message << "; " << failure;
            }
            message << "; preserved backup paths above contain the original "
                       "outputs and require manual recovery";
            throw std::runtime_error(message.str());
        }
        std::rethrow_exception(original);
    }

    for (const auto& backup : backups) {
        if (!backup.active) continue;
        std::error_code error;
        fs::remove(backup.path, error);
        if (error) {
            throw std::runtime_error(
                "whitelist scan outputs were committed, but an old-output "
                "backup could not be removed: " + backup.path.string());
        }
    }
}

Rcpp::List scan_config_as_list(const NativeScanConfig& config) {
    auto optional_string = [](const std::string& input) -> Rcpp::RObject {
        return input.empty() ? Rcpp::RObject(R_NilValue)
                             : Rcpp::RObject(Rcpp::wrap(input));
    };
    Rcpp::List value;
    value["input"] = config.input;
    value["output_prefix"] = config.output_prefix;
    value["adapter_seq"] = config.adapter_seq;
    value["barcode_length"] = config.barcode_length;
    value["left_margin"] = config.left_margin;
    value["right_margin"] = config.right_margin;
    value["max_reads"] = config.max_reads > 0
                             ? Rcpp::RObject(Rcpp::wrap(config.max_reads))
                             : Rcpp::RObject(R_NilValue);
    value["max_error"] = config.max_error;
    value["whitelist"] = optional_string(config.whitelist);
    value["threads"] = config.effective_threads;
    value["threads_requested"] = config.requested_threads;
    value["threads_effective"] = config.effective_threads;
    value["chunk_size"] = config.chunk_size;
    value["selection"] =
        core::scan_whitelist_selection_name(config.selection);
    value["rescan"] = config.rescan;
    value["bc1_whitelist"] = optional_string(config.bc1_whitelist);
    value["bc2_whitelist"] = optional_string(config.bc2_whitelist);
    value["bc1_length"] = config.bc1_length > 0
                              ? Rcpp::RObject(Rcpp::wrap(config.bc1_length))
                              : Rcpp::RObject(R_NilValue);
    value["bc2_length"] = config.bc2_length > 0
                              ? Rcpp::RObject(Rcpp::wrap(config.bc2_length))
                              : Rcpp::RObject(R_NilValue);
    value["umi_length"] = config.umi_length;
    value["offset_min"] = config.offset_min;
    value["offset_max"] = config.offset_max;
    value["verbose"] = config.verbose;
    return value;
}

Rcpp::List scan_stats_as_list(const NativeScanRun& run,
                              int requested_threads,
                              int effective_threads) {
    const auto& stats = run.stats;
    return Rcpp::List::create(
        Rcpp::_["threads_requested"] = requested_threads,
        Rcpp::_["threads_effective"] = effective_threads,
        Rcpp::_["reads_processed"] =
            static_cast<double>(stats.reads_processed),
        Rcpp::_["input_rows"] = static_cast<double>(run.input_rows),
        Rcpp::_["total_extractions"] =
            static_cast<double>(stats.total_extractions),
        Rcpp::_["unique_sequences"] =
            static_cast<double>(stats.unique_sequences),
        Rcpp::_["total_perfect_matches"] =
            static_cast<double>(stats.total_perfect_matches),
        Rcpp::_["unique_perfect_matches"] =
            static_cast<double>(stats.unique_perfect_matches),
        Rcpp::_["match_rate_percent"] = stats.match_rate_percent,
        Rcpp::_["floor"] = stats.floor,
        Rcpp::_["threshold"] = stats.threshold,
        Rcpp::_["threshold_rule"] = stats.threshold_rule,
        Rcpp::_["above_floor_barcodes"] =
            static_cast<double>(stats.above_floor_barcodes),
        Rcpp::_["high_specificity_barcodes"] =
            static_cast<double>(stats.final_barcodes),
        Rcpp::_["selected_barcodes"] =
            static_cast<double>(stats.selected_barcodes),
        Rcpp::_["selection"] =
            core::scan_whitelist_selection_name(stats.selection),
        Rcpp::_["bc1_whitelist_size"] =
            static_cast<double>(run.bc1_whitelist_size),
        Rcpp::_["bc2_whitelist_size"] =
            static_cast<double>(run.bc2_whitelist_size));
}

Rcpp::List scan_run_as_list(const NativeScanRun& run,
                            const NativeScanConfig& config) {
    NativeScanConfig effective = config;
    effective.input = run.input;
    effective.output_prefix = run.output_prefix;
    if (run.effective_bc1_length > 0) {
        effective.bc1_length = run.effective_bc1_length;
    }
    if (run.effective_bc2_length > 0) {
        effective.bc2_length = run.effective_bc2_length;
    }
    return Rcpp::List::create(
        Rcpp::_["input"] = run.input,
        Rcpp::_["output_prefix"] = run.output_prefix,
        Rcpp::_["mode"] = run.mode,
        Rcpp::_["config"] = scan_config_as_list(effective),
        Rcpp::_["files"] = Rcpp::List::create(
            Rcpp::_["csv"] = run.csv,
            Rcpp::_["whitelist"] = run.whitelist,
            Rcpp::_["scan_log"] = run.log,
            Rcpp::_["valid_pairs"] = run.valid_pairs,
            Rcpp::_["spatial_mask"] = run.spatial_mask),
        Rcpp::_["stats"] = scan_stats_as_list(
            run, effective.requested_threads, effective.effective_threads));
}

NativeScanRun run_scan_whitelist(
    const NativeScanConfig& config,
    std::function<void()> boundary_poll = {}) {
    const fs::path input = fs::absolute(config.input).lexically_normal();
    if (!fs::is_regular_file(input)) {
        throw std::invalid_argument("scan input does not exist: " +
                                    input.string());
    }
    if (config.output_prefix.empty()) {
        throw std::invalid_argument("`output_prefix` is required");
    }

    fs::path output_prefix =
        fs::absolute(config.output_prefix).lexically_normal();
    if (!output_prefix.parent_path().empty()) {
        std::error_code error;
        fs::create_directories(output_prefix.parent_path(), error);
        if (error) {
            throw std::runtime_error("could not create scan output directory: " +
                                     error.message());
        }
    }

    const fs::path final_csv = output_prefix.string() + ".csv";
    const fs::path final_txt = output_prefix.string() + ".txt";
    const fs::path final_log = config.log_path.empty()
                                   ? fs::path(output_prefix.string() +
                                              "_scan_wl.log")
                                   : fs::absolute(config.log_path)
                                         .lexically_normal();
    const fs::path final_pairs = output_prefix.string() + "_valid_pairs.csv";
    const fs::path final_mask = output_prefix.string() + "_spat_mask.csv";
    const fs::path commit_lock =
        output_prefix.parent_path() /
        ("." + output_prefix.filename().string() + ".rrad-commit.lock");
    if (final_log.parent_path() != output_prefix.parent_path()) {
        throw std::invalid_argument(
            "whitelist scan log must be in the output-prefix directory");
    }
    if (final_log == commit_lock) {
        throw std::invalid_argument(
            "whitelist scan log cannot replace its publication lock");
    }

    std::vector<NamedPath> inputs{{"scan input", input}};
    if (!config.whitelist.empty()) {
        inputs.emplace_back(
            "reference whitelist",
            fs::path(core::whitelist_utils::kit_to_path(config.whitelist)));
    }
    if (!config.bc1_whitelist.empty()) {
        inputs.emplace_back(
            "BC1 whitelist",
            fs::path(core::whitelist_utils::kit_to_path(config.bc1_whitelist)));
    }
    if (!config.bc2_whitelist.empty()) {
        inputs.emplace_back(
            "BC2 whitelist",
            fs::path(core::whitelist_utils::kit_to_path(config.bc2_whitelist)));
    }
    std::vector<NamedPath> outputs{
        {"scan statistics", final_csv},
        {"selected whitelist", final_txt},
        {"scan log", final_log},
        {"valid barcode pairs", final_pairs},
        {"spatial mask", final_mask}};
    reject_input_output_collisions(inputs, outputs);
    reject_input_output_collisions(
        inputs, {{"scan publication lock", commit_lock}});
    std::vector<fs::path> final_targets{
        final_csv, final_txt, final_log, final_pairs, final_mask};
    preflight_scan_outputs(final_targets);

    ScanStagingDirectory staging(output_prefix);
    const std::string stage_base =
        (staging.path() / output_prefix.filename()).string();
    const fs::path stage_csv = stage_base + ".csv";
    const fs::path stage_txt = stage_base + ".txt";
    const fs::path stage_log = stage_base + "_scan_wl.log";
    const fs::path stage_pairs = stage_base + "_valid_pairs.csv";
    const fs::path stage_mask = stage_base + "_spat_mask.csv";

    NativeScanRun run;
    run.input = input.string();
    run.output_prefix = output_prefix.string();
    run.staging_base = stage_base;
    const auto started = std::chrono::steady_clock::now();

    if (config.rescan) {
        run.mode = "rescan";
        const core::scan_wl_rescan_input retained =
            core::read_scan_wl_rescan_counts(input.string(), boundary_poll);
        run.input_rows = retained.input_rows;
        const std::unordered_set<core::int64_seq> no_filter;
        if (!core::count_perfect_matches_with_stats(
                retained.counts, no_filter, stage_csv.string(),
                stage_txt.string(), config.verbose, config.selection,
                &run.stats, boundary_poll)) {
            throw std::runtime_error("whitelist rescan failed");
        }
    } else if (config.split_mode()) {
        run.mode = "split";
        const std::string bc1_path =
            core::whitelist_utils::kit_to_path(config.bc1_whitelist);
        const std::string bc2_path =
            core::whitelist_utils::kit_to_path(config.bc2_whitelist);
        auto packed_length = [](const std::string& spec,
                                const std::string& resolved, int supplied)
            -> std::optional<std::uint16_t> {
            if (supplied > 0) {
                return static_cast<std::uint16_t>(supplied);
            }
            const std::string filename = fs::path(resolved).filename().string();
            if (spec == "splitseq_bc1" || spec == "splitseq_bc2" ||
                filename == "splitseq_bc1_bitlist.csv.gz" ||
                filename == "splitseq_bc2_bitlist.csv.gz") {
                return static_cast<std::uint16_t>(8);
            }
            return std::nullopt;
        };
        const std::optional<std::uint16_t> bc1_length = packed_length(
            config.bc1_whitelist, bc1_path, config.bc1_length);
        const std::optional<std::uint16_t> bc2_length = packed_length(
            config.bc2_whitelist, bc2_path, config.bc2_length);
        const core::split_whitelist bc1 =
            core::load_split_whitelist_csv(bc1_path, bc1_length);
        const core::split_whitelist bc2 =
            core::load_split_whitelist_csv(bc2_path, bc2_length);
        if (bc1.seqs.empty() || bc2.seqs.empty()) {
            throw std::runtime_error(
                "two-part barcode whitelists contain no valid sequences");
        }
        auto validate_declared_length = [](
            const core::split_whitelist& whitelist,
            const std::optional<std::uint16_t>& declared,
            const char* option_name) {
            if (declared &&
                (whitelist.lengths.size() != 1 ||
                 whitelist.lengths.front() != static_cast<int>(*declared))) {
                throw std::invalid_argument(
                    std::string("`") + option_name +
                    "` does not match the barcode lengths in its whitelist");
            }
        };
        validate_declared_length(bc1, bc1_length, "bc1_length");
        validate_declared_length(bc2, bc2_length, "bc2_length");
        run.bc1_whitelist_size = bc1.seqs.size();
        run.bc2_whitelist_size = bc2.seqs.size();
        run.effective_bc1_length = bc1.lengths.size() == 1
                                       ? bc1.lengths.front()
                                       : 0;
        run.effective_bc2_length = bc2.lengths.size() == 1
                                       ? bc2.lengths.front()
                                       : 0;
        const core::split_barcode_count_result barcodes =
            core::process_fastq_split_barcode(
                input.string(), config.adapter_seq, config.umi_length,
                bc1.seqs, bc1.lengths, bc2.seqs, bc2.lengths,
                config.offset_min, config.offset_max, config.max_reads,
                config.max_error, config.chunk_size,
                config.effective_threads, boundary_poll);
        if (!core::write_split_pair_counts(barcodes,
                                           stage_pairs.string())) {
            throw std::runtime_error("could not write valid barcode pairs");
        }
        if (!core::write_spat_mask_csv(barcodes, bc1, bc2,
                                       stage_mask.string())) {
            throw std::runtime_error("could not write spatial mask");
        }
        const std::unordered_set<core::int64_seq> no_filter;
        if (!core::count_perfect_matches_with_stats(
                barcodes, no_filter, stage_csv.string(), stage_txt.string(),
                config.verbose, config.selection, &run.stats,
                boundary_poll)) {
            throw std::runtime_error("two-part whitelist selection failed");
        }
    } else {
        run.mode = "single";
        core::scan_whitelist_filter filter;
        if (!config.whitelist.empty()) {
            filter = core::load_scan_whitelist(
                config.whitelist,
                static_cast<std::uint16_t>(config.barcode_length),
                config.verbose);
            if (filter.empty()) {
                throw std::runtime_error(
                    "reference whitelist contains no valid barcodes");
            }
        }
        core::barcode_count_result barcodes = core::process_fastq(
            input.string(), config.adapter_seq, config.barcode_length,
            config.left_margin, config.right_margin, config.max_reads,
            config.max_error, config.chunk_size, config.effective_threads,
            filter.empty() ? nullptr : &filter, boundary_poll);
        if (!barcodes.succeeded) {
            throw std::runtime_error("whitelist scan could not read the input");
        }
        if (!core::count_perfect_matches_with_stats(
                barcodes, filter.long_barcodes, stage_csv.string(),
                stage_txt.string(), config.verbose, config.selection,
                &run.stats, boundary_poll)) {
            throw std::runtime_error("whitelist selection failed");
        }
    }

    const double wall_seconds = std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() - started)
                                    .count();
    if (!core::write_whitelist_scan_log(stage_log.string(), input.string(),
                                        run.stats, wall_seconds)) {
        throw std::runtime_error("could not write whitelist scan summary");
    }

    require_scan_output(stage_csv, "whitelist statistics");
    require_scan_output(stage_txt, "selected whitelist");
    require_scan_output(stage_log, "whitelist scan log");
    if (config.split_mode()) {
        require_scan_output(stage_pairs, "valid barcode pair table");
        require_scan_output(stage_mask, "spatial mask");
    }

    const bool split = config.split_mode();
    std::vector<ScanCommitItem> commits{
        {stage_csv, final_csv, true},
        {stage_txt, final_txt, true},
        {stage_log, final_log, true},
        {stage_pairs, final_pairs, split},
        {stage_mask, final_mask, split}};
    commit_scan_outputs(commits, staging.path().filename().string(),
                        commit_lock);

    run.csv = final_csv.string();
    run.whitelist = final_txt.string();
    run.log = final_log.string();
    if (config.split_mode()) {
        run.valid_pairs = final_pairs.string();
        run.spatial_mask = final_mask.string();
    }
    return run;
}
