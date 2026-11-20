#pragma once
namespace rad_discovery {
namespace run_detail {
std::string csv_string(const std::string& value) {
    std::string result = "\"";
    for (char c : value) {
        if (c == '"') result += '"';
        result += c;
    }
    return result + '"';
}

std::string join_aliases(const std::vector<std::string>& aliases) {
    std::string result;
    for (const auto& alias : aliases) {
        if (!result.empty()) result += ';';
        result += alias;
    }
    return result;
}

std::ofstream output_file(const std::string& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write output: " + path);
    out.exceptions(std::ios::badbit | std::ios::failbit);
    return out;
}

std::vector<Read> sample_reads(const Options& options, std::size_t& bases) {
    // zlib keeps prefix sampling bounded without a decompressor helper process.
    file_streaming files(options.input, options.threads, false);
    read_streaming reader(files);
    std::vector<Read> reads;
    reads.reserve(std::min<std::size_t>(options.sample_reads, 10000));
    bases = 0;
    while (reads.size() < options.sample_reads && bases < options.sample_bases) {
        if (options.boundary_poll) options.boundary_poll();
        auto record = reader.next_sequence();
        if (!record) break;
        if (record->seq.size() > options.sample_bases - bases) break;
        if (record->seq.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("sampled read exceeds the supported alignment length");
        }
        std::transform(record->seq.begin(), record->seq.end(), record->seq.begin(), [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
        bases += record->seq.size();
        reads.push_back({std::move(record->id), std::move(record->seq), std::move(record->qual)});
    }
    if (reader.read_failed()) {
        throw std::runtime_error("malformed sequence input in " + reader.read_error_path() +
                                 " (reader error " + std::to_string(reader.read_error_code()) + ")");
    }
    if (reads.empty()) {
        throw std::runtime_error("sample contains no nonempty reads within the base limit; check input or increase --sample-bases");
    }
    return reads;
}

std::string write_sample(const std::string& prefix, const std::vector<Read>& reads) {
    const bool fastq = std::all_of(reads.begin(), reads.end(), [](const Read& read) {
        return read.quality.size() == read.sequence.size();
    });
    const auto path = prefix + (fastq ? "_sample.fastq" : "_sample.fasta");
    auto out = output_file(path);
    for (std::size_t i = 0; i < reads.size(); ++i) {
        const auto& read = reads[i];
        // Sample-local identifiers avoid header ambiguities; the cached Read retains the original ID.
        if (fastq) out << "@sample_" << i << '\n' << read.sequence << "\n+\n" << read.quality << '\n';
        else out << ">sample_" << i << '\n' << read.sequence << '\n';
    }
    out.close();
    return path;
}

void write_parts(const std::string& path, const LayoutInference& inference) {
    auto out = output_file(path);
    out << "id,sequence,reads,forward,reverse\n";
    for (const auto& part : inference.parts) {
        out << csv_string(part.id) << ',' << csv_string(part.sequence) << ','
            << part.reads << ',' << part.forward << ',' << part.reverse << '\n';
    }
    out.close();
}

void write_matches(const std::string& path, const std::vector<CatalogMatch>& matches) {
    auto out = output_file(path);
    out << "recipe_id,catalog,aliases,length,catalog_size,matched_reads,matched_distinct,extracted_reads,observed_distinct,sampled_reads\n";
    for (const auto& match : matches) {
        out << csv_string(match.recipe_id) << ',' << csv_string(match.path) << ','
            << csv_string(join_aliases(match.aliases)) << ',' << match.length << ','
            << match.catalog_size << ',' << match.reads << ',' << match.distinct << ','
            << match.extracted << ',' << match.observed_distinct << ',' << match.sampled_reads << '\n';
    }
    out.close();
}

void write_choices(const std::string& path, const std::vector<WhitelistChoice>& choices,
                   const std::vector<CatalogMatch>& matches) {
    auto out = output_file(path);
    out << "status,recipe_id,catalog,aliases,matched_reads,extracted_reads,matched_distinct\n";
    if (choices.empty()) out << "no_hits,\"\",\"\",\"\",0,,0\n";
    for (const auto& choice : choices) {
        if (choice.winners.empty()) {
            // If multiple geometries for this recipe have different extraction
            // totals, no single denominator exists. An absent eligible catalog
            // also leaves this denominator unknown rather than establishing 0.
            std::optional<std::uint64_t> extracted;
            bool ambiguous = false;
            for (const auto& match : matches) {
                if (match.recipe_id != choice.recipe_id) continue;
                if (!extracted) extracted = match.extracted;
                else if (*extracted != match.extracted) ambiguous = true;
            }
            out << "no_hits," << csv_string(choice.recipe_id) << ",\"\",\"\",0,";
            if (extracted && !ambiguous) out << *extracted;
            out << ",0\n";
            continue;
        }
        for (const auto& match : choice.winners) {
            out << choice.status << ',' << csv_string(choice.recipe_id) << ','
                << csv_string(match.path) << ',' << csv_string(join_aliases(match.aliases)) << ','
                << match.reads << ',' << match.extracted << ',' << match.distinct << '\n';
        }
    }
    out.close();
}

Recipe explicit_recipe(Options& options, bool length_given) {
    std::size_t capture = 0, literal = 0;
    for (char& c : options.pattern) {
        if (c == 'N') { ++capture; continue; }
        if (c == 'n') continue;
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            throw std::invalid_argument("--pattern accepts A, C, G, T, N (capture), and n (match-only)");
        }
        ++literal;
    }
    if (literal == 0) throw std::invalid_argument("--pattern needs at least one literal anchor base");
    if (capture > 64) throw std::invalid_argument("--pattern captures more than 64 barcode bases");
    if (capture > 0) {
        if (length_given && options.barcode_length != static_cast<int>(capture)) {
            throw std::invalid_argument("--barcode-length must match the uppercase N capture count");
        }
        options.barcode_length = static_cast<int>(capture);
    } else if (!length_given) {
        throw std::invalid_argument("an unmasked --pattern requires --barcode-length");
    }
    return {"explicit_pattern", options.pattern, options.barcode_length, 0, 0};
}


} // namespace run_detail
} // namespace rad_discovery
