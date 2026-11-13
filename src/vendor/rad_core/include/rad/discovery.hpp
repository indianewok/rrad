#pragma once


namespace rad_discovery {

struct CatalogInput { std::string name, path; int length = 0; };

struct Options {
    std::string input;
    std::string output;
    std::string parts_csv;
    std::string pattern;
    int barcode_length = 0;
    int threads = 1;
    std::size_t sample_reads = 10000;
    std::size_t sample_bases = 50000000;
    double error_ratio = 0.15;
    bool layout = true;
    bool wl = true;
    std::vector<CatalogInput> catalogs;
    std::function<void()> boundary_poll;
    std::vector<std::string>* unavailable = nullptr;
};

struct Read {
    std::string id, sequence, quality;
};

struct Row {
    std::string id, sequence;
    int length = 0;
    std::string type, klass, whitelist;
    std::string length_spec, flags;
};

// A scan-wl extraction recipe, independent of any reference catalog.
struct Recipe {
    std::string id, pattern;
    int length = 0;
    int left_margin = 0, right_margin = 0;
};

struct PartSupport {
    std::string id, sequence;
    std::size_t reads = 0, forward = 0, reverse = 0;
};

struct LayoutInference {
    std::vector<Row> rows;
    std::vector<Recipe> recipes;
    std::vector<PartSupport> parts;
    std::vector<std::string> notes;
    bool complete = false;
};

struct CatalogMatch {
    std::string recipe_id, path;
    std::vector<std::string> aliases;
    int length = 0;
    std::uint64_t reads = 0, distinct = 0;
    std::uint64_t extracted = 0, observed_distinct = 0, sampled_reads = 0;
    std::uint64_t catalog_size = 0;
};

struct WhitelistChoice {
    std::string recipe_id;
    std::string status; // matched, tie, no_hits
    std::uint64_t matched_reads = 0;
    std::vector<CatalogMatch> winners;
};

// Selection is deliberately just the largest exact matching read count.
// Distinct barcode counts and catalog sizes are descriptive, never tie breakers.
inline std::vector<WhitelistChoice> select_whitelists(
        const std::vector<CatalogMatch>& matches, const std::vector<Recipe>& recipes) {
    std::map<std::string, std::map<std::string, CatalogMatch>> groups;
    for (const auto& recipe : recipes) groups[recipe.id];
    for (const auto& match : matches) {
        auto& catalogs = groups[match.recipe_id];
        auto found = catalogs.find(match.path);
        if (found == catalogs.end() || match.reads > found->second.reads)
            catalogs[match.path] = match;
    }
    std::vector<WhitelistChoice> selected;
    for (const auto& group : groups) {
        WhitelistChoice choice;
        choice.recipe_id = group.first;
        for (const auto& catalog : group.second)
            choice.matched_reads = std::max(choice.matched_reads, catalog.second.reads);
        if (choice.matched_reads) {
            for (const auto& catalog : group.second)
                if (catalog.second.reads == choice.matched_reads)
                    choice.winners.push_back(catalog.second);
        }
        choice.status = choice.winners.empty() ? "no_hits" : choice.winners.size() == 1 ? "matched" : "tie";
        selected.push_back(std::move(choice));
    }
    return selected;
}

LayoutInference infer_layout(const std::vector<Read>& reads, const Options& options);
std::vector<CatalogMatch> identify_whitelists(const std::string& sample_path,
                                            const std::vector<Recipe>& recipes,
                                            const Options& options);
void write_layout(const std::string& path, const LayoutInference& inference,
                  const std::vector<CatalogMatch>& matches);


} // namespace rad_discovery
