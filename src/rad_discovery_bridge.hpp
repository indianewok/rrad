// Included after the shared bridge helpers. No R API is called by worker threads.
SEXP rad_discover_cpp(Rcpp::List values) {
    namespace discovery = core::rad_discovery;
    namespace detail = discovery::run_detail;
    try {
        discovery::Options options;
        options.input = string_option(values, "input");
        options.output = string_option(values, "output_prefix");
        options.pattern = string_option(values, "pattern");
        options.parts_csv = string_option(values, "parts");
        options.barcode_length = int_option(values, "barcode_length", 0);
        options.sample_reads = int_option(values, "sample_reads", 10000);
        options.sample_bases = int_option(values, "sample_bases", 50000000);
        options.error_ratio = double_option(values, "max_error", 0.15);
        const int requested_threads = int_option(values, "threads", 1);
        options.threads = requested_threads;
#ifndef _OPENMP
        options.threads = 1;
#endif
        options.layout = bool_option(values, "layout", true);
        options.wl = bool_option(values, "whitelist", true);
        const bool overwrite = bool_option(values, "overwrite", false);
        options.boundary_poll = []() { Rcpp::checkUserInterrupt(); };
        if (options.input.empty() || options.output.empty() ||
            requested_threads < 1 || options.sample_reads < 1 ||
            options.sample_reads > static_cast<size_t>(INT_MAX) ||
            options.sample_bases < 1 || options.sample_bases > static_cast<size_t>(INT_MAX) ||
            options.barcode_length < 0 || options.barcode_length > 64 ||
            options.error_ratio < 0 || options.error_ratio > 1 ||
            (!options.layout && !options.wl)) {
            throw std::invalid_argument("invalid discovery input, mode, limits or alignment settings");
        }
        if (options.layout && !options.pattern.empty()) {
            throw std::invalid_argument("an explicit pattern is supported only in whitelist mode");
        }
        if (!options.layout && options.pattern.empty()) {
            throw std::invalid_argument("whitelist-only discovery requires a pattern");
        }
        require_nonempty_regular_file(options.input, "discovery input");
        std::vector<NamedPath> inputs{{"discovery input", options.input}};
        if (!options.parts_csv.empty()) {
            require_nonempty_regular_file(options.parts_csv, "discovery parts");
            inputs.emplace_back("discovery parts", options.parts_csv);
        }
        if (has_value(values, "catalogs")) {
            const Rcpp::CharacterVector catalogs(values["catalogs"]);
            const Rcpp::CharacterVector labels(values["catalog_names"]);
            const Rcpp::IntegerVector lengths(values["catalog_lengths"]);
            if (catalogs.size() != labels.size() || catalogs.size() != lengths.size())
                throw std::invalid_argument("catalog paths, names and lengths must agree");
            for (R_xlen_t i = 0; i < catalogs.size(); ++i) {
                if (Rcpp::CharacterVector::is_na(catalogs[i]) ||
                    Rcpp::CharacterVector::is_na(labels[i]) || lengths[i] < 0 || lengths[i] > 64)
                    throw std::invalid_argument("invalid discovery catalog");
                discovery::CatalogInput catalog{Rcpp::as<std::string>(labels[i]),
                    Rcpp::as<std::string>(catalogs[i]), lengths[i]};
                require_nonempty_regular_file(catalog.path, "discovery catalog");
                inputs.emplace_back("discovery catalog", catalog.path);
                options.catalogs.push_back(std::move(catalog));
            }
        }
        const std::string resources = normalize_resource_dir("");
        EnvironmentGuard resource_guard("RAD_RESOURCES", resources);
        EnvironmentGuard embedded_guard("RAD_EMBEDDED_ONLY", "1");
        EnvironmentGuard pigz_guard("RAD_NO_PIGZ", "1");
        OpenMpThreadGuard thread_guard(options.threads);
        StreamCapture capture;
        const fs::path final_prefix = fs::absolute(options.output).lexically_normal();
        const std::vector<std::string> suffixes{
            "_layout.csv", "_layout.partial.csv", "_parts.csv",
            "_whitelists.csv", "_whitelist_match.csv", "_discover.log"};
        std::vector<NamedPath> targets;
        for (const auto& suffix : suffixes)
            targets.emplace_back("discovery" + suffix, final_prefix.string() + suffix);
        reject_input_output_collisions(inputs, targets);
        preflight_transaction_targets(targets);
        if (!overwrite) for (const auto& target : targets) {
            if (fs::exists(target.second) || fs::is_symlink(target.second))
                throw std::invalid_argument("output already exists: " + target.second.string());
        }
        SameFilesystemOutputTransaction transaction(final_prefix);
        reject_input_output_collisions(inputs,
            {{"discovery publication lock", transaction.lock_path()}});
        options.output = transaction.stage_base().string();
        std::vector<std::string> unavailable;
        options.unavailable = &unavailable;
        size_t bases = 0;
        const auto reads = detail::sample_reads(options, bases);
        const auto sample_path = detail::write_sample(options.output, reads);
        discovery::LayoutInference inference;
        if (options.layout) {
            inference = discovery::infer_layout(reads, options);
            detail::write_parts(options.output + "_parts.csv", inference);
        } else {
            inference.recipes.push_back(detail::explicit_recipe(
                options, has_value(values, "barcode_length")));
        }
        std::vector<discovery::CatalogMatch> matches;
        std::vector<discovery::WhitelistChoice> choices;
        bool whitelist_complete = true;
        if (options.wl) {
            matches = discovery::identify_whitelists(sample_path, inference.recipes, options);
            detail::write_matches(options.output + "_whitelists.csv", matches);
            choices = discovery::select_whitelists(matches, inference.recipes);
            detail::write_choices(options.output + "_whitelist_match.csv", choices, matches);
            whitelist_complete = !choices.empty() && std::all_of(choices.begin(), choices.end(),
                [](const discovery::WhitelistChoice& choice) { return choice.status == "matched"; });
        }
        const std::string layout_suffix = inference.complete ? "_layout.csv" : "_layout.partial.csv";
        if (options.layout)
            discovery::write_layout(options.output + layout_suffix, inference, matches);
        const bool complete = (!options.layout || inference.complete) &&
                              (!options.wl || whitelist_complete);
        auto log = detail::output_file(options.output + "_discover.log");
        log << "command=discover\ninput=" << options.input
            << "\nsampling=prefix_whole_reads\nsample_reads=" << reads.size()
            << "\nsample_bases=" << bases << "\nthreads=" << options.threads
            << "\nlayout_complete=" << inference.complete
            << "\nwhitelist_scope=nonspatial_catalogs_only"
            << "\nwhitelist_ranking=highest_independent_exact_read_hits"
            << "\nstatus=" << (complete ? "complete" : "incomplete") << '\n';
        for (const auto& note : inference.notes) log << "note=" << note << '\n';
        for (const auto& name : unavailable) log << "unavailable_catalog=" << name << '\n';
        for (const auto& choice : choices)
            log << "whitelist_status." << choice.recipe_id << '=' << choice.status << '\n';
        log.close();
        std::vector<TransactionOutput> outputs;
        Rcpp::List files;
        for (const auto& suffix : suffixes) {
            const fs::path staged = options.output + suffix;
            const fs::path target = final_prefix.string() + suffix;
            const bool present = fs::exists(staged);
            outputs.push_back({"discovery" + suffix, staged, target, present});
            if (present) files[suffix.substr(1)] = target.string();
        }
        Rcpp::List result = Rcpp::List::create(
            Rcpp::_["success"] = true, Rcpp::_["complete"] = complete,
            Rcpp::_["layout_complete"] = options.layout && inference.complete,
            Rcpp::_["whitelist_complete"] = options.wl && whitelist_complete,
            Rcpp::_["sample_reads"] = static_cast<double>(reads.size()),
            Rcpp::_["sample_bases"] = static_cast<double>(bases),
            Rcpp::_["threads_requested"] = requested_threads,
            Rcpp::_["threads_effective"] = options.threads,
            Rcpp::_["unavailable_catalogs"] = unavailable,
            Rcpp::_["notes"] = inference.notes,
            Rcpp::_["files"] = files, Rcpp::_["core"] = core_info_impl(),
            Rcpp::_["log"] = capture.str());
        options.boundary_poll();
        transaction.commit(outputs, overwrite);
        return result;
    } catch (const std::exception& error) {
        Rcpp::stop("embedded RAD discovery failed: %s", error.what());
    }
    return R_NilValue;
}

SEXP rad_spatial_reformat_cpp(Rcpp::List values) {
    try {
        core::spatial_reformat_options options;
        options.input = string_option(values, "input");
        options.output = string_option(values, "output_dir");
        options.sample = string_option(values, "sample");
        options.map_path = string_option(values, "map");
        options.bc1_path = string_option(values, "bc1_axis");
        options.bc2_path = string_option(values, "bc2_axis");
        options.metadata_path = string_option(values, "metadata");
        options.unit = string_option(values, "unit", "bin");
        options.bin_size = int_option(values, "bin_size", 2);
        options.tissue_only = bool_option(values, "tissue_only", false);
        options.strict = bool_option(values, "strict", false);
        options.boundary_poll = []() { Rcpp::checkUserInterrupt(); };
        if (options.output.empty()) throw std::invalid_argument("output_dir is required");
        const fs::path destination = fs::absolute(options.output).lexically_normal();
        if (fs::exists(destination) || fs::is_symlink(destination))
            throw std::invalid_argument("output directory already exists: " + destination.string());
        std::vector<NamedPath> inputs;
        for (const auto& path : {options.input, options.map_path, options.bc1_path,
                                 options.bc2_path, options.metadata_path}) {
            if (path.empty()) continue;
            require_nonempty_regular_file(path, "spatial input");
            inputs.emplace_back("spatial input", path);
        }
        StreamCapture capture;
        SameFilesystemOutputTransaction transaction(destination);
        reject_input_output_collisions(inputs,
            {{"spatial publication lock", transaction.lock_path()}});
        options.output = transaction.stage_base().string();
        const auto run = core::reformat_spatial(options);
        std::vector<TransactionOutput> outputs;
        Rcpp::List files;
        for (const auto& name : {"reads.fastq.gz", "assignments.tsv.gz", "summary.tsv"}) {
            fs::path staged = fs::path(options.output) / name;
            if (std::string(name).find(".gz") != std::string::npos)
                verify_primary_gzip(staged, false);
            outputs.push_back({name, staged, destination / name, true});
            files[name] = (destination / name).string();
        }
        Rcpp::List counts;
        for (const auto& entry : run.counts) counts[entry.first] = static_cast<double>(entry.second);
        Rcpp::List result = Rcpp::List::create(
            Rcpp::_["success"] = true, Rcpp::_["files"] = files,
            Rcpp::_["records_read"] = static_cast<double>(run.records_read),
            Rcpp::_["records_written"] = static_cast<double>(run.records_written),
            Rcpp::_["counts"] = counts, Rcpp::_["core"] = core_info_impl(),
            Rcpp::_["log"] = capture.str());
        options.boundary_poll();
        if (!fs::create_directory(destination))
            throw std::invalid_argument("output directory already exists: " + destination.string());
        try {
            transaction.commit(outputs, false);
        } catch (...) {
            // Remove only our empty directory. Never remove unrelated files.
            std::error_code ignored;
            fs::remove(destination, ignored);
            throw;
        }
        return result;
    } catch (const std::exception& error) {
        Rcpp::stop("embedded RAD spatial export failed: %s", error.what());
    }
    return R_NilValue;
}
